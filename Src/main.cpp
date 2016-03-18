/** @file man.cpp

  Convert OpenEXR image file to PNG image file.
*/
#include <ImfRgbaFile.h>
#include <ImfArray.h>
#include <png.h>
#include <iostream>
#include <algorithm>
#include <memory>
#include <vector>

/** The result of the conversion.
*/
enum Result {
  RESULT_SUCCESS = 0, ///< conversion success.
  RESULT_ERROR = 1, ///< conversion failure.
};

/** Read OpenEXR image from file.

  @param filename  OpenEXR file path.
  @param pPixels   The pointer to storing OpenEXR image if it is not nullptr.
  @param pWidth    The pointer to storing OpenEXR image width if it is not nullptr.
  @param pHeight   The pointer to storing OpenEXR image height if it is not nullptr.

  @retval RESULT_SUCCESS  Succeeded to read image file.
  @retval RESULT_ERROR    Failed  to read image file.
                          pPixels, pWidth, pHeight has invalid data, so shouldn't refer.
                          The file is broken or isn't OpenEXR format, probably.
*/
Result ReadRgbaExrFile(const char* filename, Imf::Array2D<Imf::Rgba>* pPixels, int* pWidth, int* pHeight) {
  Imf::RgbaInputFile file (filename);
  if (!file.isComplete()) {
	return RESULT_ERROR;
  }
  Imath::Box2i dw = file.dataWindow();

  if (pWidth) {
    *pWidth  = dw.max.x - dw.min.x + 1;
  }
  if (pHeight) {
    *pHeight = dw.max.y - dw.min.y + 1;
  }
  if (pPixels) {
    pPixels->resizeErase (*pHeight, *pWidth);
    file.setFrameBuffer(&(*pPixels)[0][0] - dw.min.x - dw.min.y * *pWidth, 1, *pWidth);
    file.readPixels(dw.min.y, dw.max.y);
  }
  return RESULT_SUCCESS;
}

/*** Print program usage.
*/
void PrintUsage() {
  std::cout << "exr2png.exe [-s scale] [infile] [outfile]" << std::endl;
  std::cout << std::endl;
  std::cout << "  Convert to PNG image from OpenEXR image." << std::endl;
  std::cout << "  The alpha component in PNG has the reciprocal of the color strength." << std::endl;
  std::cout << "  To restore the color, divide RGB by alpha." << std::endl;
  std::cout << std::endl;
  std::cout << "  -s scale: The low dynamic range scale." << std::endl;
  std::cout << "            0.5 maps the range of 0-0.5 to 0-255, 2.0 maps 0-2.0 to 0-255." << std::endl;
  std::cout << "            Note that higher scale reduce the maximum brightness." << std::endl;
  std::cout << "            If not passed this option, the default is 1.0." << std::endl;
  std::cout << "  infile  : OpenEXR image file." << std::endl;
  std::cout << "  outfile : PNG image file that converted from infile." << std::endl;
  std::cout << "            If not passed this option, use the infile that has replaced" << std::endl;
  std::cout << "            extension to '.png.'" << std::endl;
}

/** Convert half precision color format to uint8_t format.

  @param  n         The input of the half precision color.
  @param  strength  The strength of the color that is biggest of the half color of RGB.

  @return The uint8_t color that is converted from the half precision color.
*/
uint8_t HalfToUint8(half n, float strength) {
  return std::max<uint8_t>(0, std::min<uint8_t>(255, static_cast<uint8_t>((255.0f * n) / strength + 0.5f)));
}

/** Convert to the png_byte color from the half precision RGBA color.

  @param buf       The pointer to the png_byte color buffer.
  @param rgba      The source of the half precision RGBA color. The element of A is ignored.
  @param strength  The scale of the color strength. The default value is 1.0.
                   If it is 0.5, the color is doubled. If it is 2.0, the color is halved.

  This function stores 4 byte(RGBA) to buf.
  buf should have large enough to store it.
  The element of A in PNG has the reciprocal of strength of the color.
*/
void SetPixel(png_bytep buf, const Imf::Rgba& rgba, float strength) {
  uint8_t* p = static_cast<uint8_t*>(buf);
  half biggest = rgba.r;
  if (biggest < rgba.g) biggest = rgba.g;
  if (biggest < rgba.b) biggest = rgba.b;
  if (biggest >= strength) {
	strength *= biggest;
  }
  p[0] = HalfToUint8(rgba.r, strength);
  p[1] = HalfToUint8(rgba.g, strength);
  p[2] = HalfToUint8(rgba.b, strength);
  p[3] = std::max<uint8_t>(1, std::min<uint8_t>(255, static_cast<uint8_t>(255.0f / strength + 0.5f)));
}

/** The entry point.

  @param argc  argument count.
  @param argv  argument vector.

  @retval RESULT_SUCCESS  the conversion is success.
  @retval RESULT_ERROR    the conversion is failure.
*/
int main(int argc, char** argv) {
  const char* infilename = nullptr;
  std::string outfilename;
  float strengthScale = 1.0f;
  for (int i = 1; i < argc; ++i) {
	if (argv[i][0] == '-') {
	  if (argv[i][1] == 's' || argv[i][1] == 'S' && (argc >= i + 1)) {
		strengthScale = atof(argv[i + 1]);
		++i;
	  }
	  continue;
	}
    if (!infilename) {
      infilename = argv[i];
    } else if (outfilename.empty()) {
      outfilename = argv[i];  
      break;
    }
  }
  if (!infilename) {
    PrintUsage();
    return RESULT_SUCCESS;
  }
  if (outfilename.empty()) {
    outfilename = infilename;
    auto dotPos = outfilename.find_last_of('.');
    if (dotPos != std::string::npos) {
      outfilename.erase(dotPos, std::string::npos);
    }
    outfilename += ".png";
  }

  Imf::Array2D<Imf::Rgba> pixels;
  int w, h;
  if (ReadRgbaExrFile(infilename, &pixels, &w, &h) == RESULT_ERROR) {
	std::cerr << "Error: can't read '" << infilename << "' file." << std::endl;
	return RESULT_ERROR;
  }

  auto del = [](png_structp p) { png_destroy_write_struct(&p, nullptr); };
  std::unique_ptr<png_struct, decltype(del)> pngWriter(png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr), del);
  if (!pngWriter) {
	std::cerr << "Error: png_create_write_struct failed." << std::endl;
	return RESULT_ERROR;
  }
  png_infop pngInfo = png_create_info_struct(pngWriter.get());
  if (!pngInfo) {
	std::cerr << "Error: png_create_info_struct failed." << std::endl;
	return RESULT_ERROR;
  }

  std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(outfilename.c_str(), "wb"), &fclose);
  if (!fp) {
	std::cerr << "Error: can't open '" << outfilename << "' file." << std::endl;
	return RESULT_ERROR;
  }

  if (setjmp(png_jmpbuf(pngWriter.get()))) {
	std::cerr << "Error: png output failed." << std::endl;
    return RESULT_ERROR;
  }
  png_init_io(pngWriter.get(), fp.get());
  png_set_IHDR(pngWriter.get(), pngInfo, w, h, 8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  png_write_info(pngWriter.get(), pngInfo);
  std::vector<png_byte> row;
  row.resize(4/*8bit RGBA*/ * w * sizeof(png_byte));
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      SetPixel(&row[x * 4], pixels[y][x], strengthScale);
    }
    png_write_row(pngWriter.get(), row.data());
  }

  return RESULT_SUCCESS;
}
