/** @file man.cpp

  Convert OpenEXR image file to PNG image file.
*/
#include "tga.h"
#include "cubemapgen.h"
#include <ImfRgbaFile.h>
#include <ImfArray.h>
#include <png.h>
#include <iostream>
#include <algorithm>
#include <memory>
#include <vector>
#include <array>
#include <cassert>

#define EXR2PNG_VERSION 1.4
#define EXR2PNG_TO_STR_I(x) #x
#define EXR2PNG_TO_STR(x) EXR2PNG_TO_STR_I(x)

static const int maxMipLevels = 16;

/** The result of the conversion.
*/
enum Result {
  RESULT_SUCCESS = 0, ///< conversion success.
  RESULT_ERROR = 1, ///< conversion failure.
};

/** The output file format.
*/
enum Format {
  FORMAT_PNG, ///< Output PNG format.
  FORMAT_TGA, ///< Output TGA format.
};

/** The conversion mode.
*/
enum ConversionMode {
  CONVERSIONMODE_SINGLE, ///< Output a file(default).
  CONVERSIONMODE_CUBEMAP_FILTERED, ///< Output the filtered cubemap.
  CONVERSIONMODE_CUBEMAP_IRRADIANCE, ///< Output the irradiance cubemap.
  CONVERSIONMODE_CUBEMAP_NONFILTERD, ///< Output the non filtered cubemap.
};

enum ImageDir {
  ImageDir_Top,
  ImageDir_Left,
  ImageDir_Bottom,
  ImageDir_Right,
};

struct ImageRect {
  int left, top;
  int right, bottom;
  ImageDir dir;
};

typedef std::array<ImageRect, 6> ImageRectArray;
ImageRectArray cubemapRectArray = { {
  { 256, 256, 512, 512, ImageDir_Left }, // +X
  {   0,   0, 256, 256, ImageDir_Top }, // -X
  {   0, 256, 256, 512, ImageDir_Left }, // +Y
  { 512,   0, 768, 256, ImageDir_Bottom }, // -Y
  { 256,   0, 512, 256, ImageDir_Top }, // +Z
  { 512, 256, 768, 512, ImageDir_Left }, // -Z
} };

/** Reduce the image rects.

  @param array     The array of the image rect.
  @param miplevel  A mipmap level of the reducion target.

  @return The array of the image rect that was reduced.
*/
ImageRectArray Reduce(const ImageRectArray& array, int miplevel) {
  ImageRectArray  ret = array;
  for (auto& e : ret) {
	e.left >>= (miplevel - 1);
	e.top >>= (miplevel - 1);
	e.right = std::max(1, e.right >> (miplevel - 1));
	e.bottom = std::max(1, e.bottom >> (miplevel - 1));
  }
  return ret;
}

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
  try {
	Imf::RgbaInputFile file(filename);
	if (!file.isComplete()) {
	  return RESULT_ERROR;
	}
	Imath::Box2i dw = file.dataWindow();

	if (pWidth) {
	  *pWidth = dw.max.x - dw.min.x + 1;
	}
	if (pHeight) {
	  *pHeight = dw.max.y - dw.min.y + 1;
	}
	if (pPixels) {
	  pPixels->resizeErase(*pHeight, *pWidth);
	  file.setFrameBuffer(&(*pPixels)[0][0] - dw.min.x - dw.min.y * *pWidth, 1, *pWidth);
	  file.readPixels(dw.min.y, dw.max.y);
	}
	return RESULT_SUCCESS;
  }
  catch (std::exception& e) {
	std::cout << "Error: " << e.what() << std::endl;
  }
  return RESULT_ERROR;
}

typedef Imf::Array2D<Imf::Rgba> ExrImage;

Imf::Rgba Add(const Imf::Rgba& lhs, const Imf::Rgba& rhs) {
  return Imf::Rgba(lhs.r + rhs.r, lhs.g + rhs.g, lhs.b + rhs.b, lhs.a + rhs.a);
}

Imf::Rgba Mul(const Imf::Rgba& lhs, float rhs) {
  return Imf::Rgba(lhs.r * rhs, lhs.g * rhs, lhs.b * rhs, lhs.a * rhs);
}

/** Reduce EXR image to half.

  @param source    The source image buffer.
  @param dest      The destination image buffer.
*/
void ReduceExrImage2(const ExrImage& source, ExrImage& dest) {
  const int sw = source.width();
  const int sh = source.height();
  const int dw = sw / 2;
  const int dh = sh / 2;
  dest.resizeEraseUnsafe(dh, dw);
  for (int y = 0; y < dh; ++y) {
	for (int x = 0; x < dw; ++x) {
	  const int x2 = x * 2;
	  const int y2 = y * 2;
	  Imf::Rgba color = source[y2 + 0][x2 + 0];
	  color = Add(color, source[y2 + 0][x2 + 1]);
	  color = Add(color, source[y2 + 1][x2 + 0]);
	  color = Add(color, source[y2 + 1][x2 + 1]);
	  color = Mul(color, 1.0f / 4.0f);
	  dest[y][x] = color;
	}
  }
}

/** Get an image of the miplevel.

  @param source    The source image buffer.
  @param miplevel  The miplevel.
  @param dest      The destination image buffer.
*/
void ReduceExrImage(const ExrImage& source, int miplevel, ExrImage& dest) {
  const int sw = source.width();
  const int sh = source.height();
  if (miplevel < 2 || sw == 1 || sh == 1) {
	dest.resizeEraseUnsafe(sh, sw);
	for (int y = 0; y < sh; ++y) {
	  std::copy(source[y], source[y] + sw, dest[y]);
	}
  } else if (miplevel == 2) {
	ReduceExrImage2(source, dest);
  } else {
	ExrImage tmp;
	ReduceExrImage2(source, tmp);
	ReduceExrImage(tmp, miplevel - 1, dest);
  }
}

/*** Print program usage.
*/
void PrintUsage() {
  std::cout << "exr2png Ver. " << EXR2PNG_TO_STR(EXR2PNG_VERSION) << std::endl;
  std::cout << "Convert to PNG/TGA image from OpenEXR image." << std::endl;
  std::cout << "The alpha component in PNG has the reciprocal of the color strength." << std::endl;
  std::cout << "To restore the color, divide RGB by alpha." << std::endl;
  std::cout << std::endl;
  std::cout << "usage: exr2png.exe [-s scale] [infile] [outfile]" << std::endl;
  std::cout << std::endl;
  std::cout << "  -s scale: The low dynamic range scale." << std::endl;
  std::cout << "            0.5 maps the range of 0-0.5 to 0-255, 2.0 maps 0-2.0 to 0-255." << std::endl;
  std::cout << "            Note that higher scale reduce the maximum brightness." << std::endl;
  std::cout << "            If not passed this option, the default is 1.0." << std::endl;
  std::cout << "  -f type : The output format type." << std::endl;
  std::cout << "            'png' or 'tga' is accepted." << std::endl;
  std::cout << "  -c type : Enable cubemap mode." << std::endl;
  std::cout << "            You can set 6 region with the order of +x -x +y -y +z -z." << std::endl;
  std::cout << "            And output 6 files that have the postfix of '_[pn][xyz]'" << std::endl;
  std::cout << "            'type' is following type:" << std::endl;
  std::cout << "            filtered  : generate the filtered cubemap." << std::endl;
  std::cout << "            irradiance: generate the irradiance cubemap." << std::endl;
  std::cout << "            none      : generate the non filtered cubemap." << std::endl;
  std::cout << " -a angle : A half of the filter corn angle. the default is 1.0." << std::endl;
  std::cout << "            In general, when the mip level goes up one, it will double." << std::endl;
  std::cout << " -m level : A maximum mipmap level. level number of texture is  generated." << std::endl;
  std::cout << "            An image is reduced to '1 / 2^(level - 1)'. The default is 1." << std::endl;
  std::cout << " [+-][xyz] left top right bottom:" << std::endl;
  std::cout << "            Set the source cubemap region." << std::endl;
  std::cout << "            The reagion is the rectangle of pixels." << std::endl;
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
  const float value = std::max<float>(0.0f, n);
  return std::max<uint8_t>(0, std::min<uint8_t>(255, static_cast<uint8_t>((255.0f * value) / strength + 0.5f)));
}

/** Convert to the png_byte color from the half precision RGBA color.

  @param buf            The pointer to the png_byte color buffer.
  @param rgba           The source of the half precision RGBA color. The element of A is ignored.
  @param strengthScale  The scale of the color strength. The default value is 1.0.
                        If it is 0.5, the color is doubled. If it is 2.0, the color is halved.

  This function stores 4 byte(RGBA) to buf.
  buf should have large enough to store it.
  The element of A in PNG has the reciprocal of strength of the color.
*/
template<typename T>
void SetPixel(png_bytep buf, const T& rgba, float strengthScale) {
  uint8_t* p = static_cast<uint8_t*>(buf);
  float strength = 1.0f;
  half biggest = rgba.r;
  if (biggest < rgba.g) biggest = rgba.g;
  if (biggest < rgba.b) biggest = rgba.b;
  if (biggest > strengthScale) {
	strength *= biggest / strengthScale;
  }
  p[0] = HalfToUint8(rgba.r, strength * strengthScale);
  p[1] = HalfToUint8(rgba.g, strength * strengthScale);
  p[2] = HalfToUint8(rgba.b, strength * strengthScale);
  p[3] = std::max<uint8_t>(1, std::min<uint8_t>(255, static_cast<uint8_t>(255.0f / strength + 0.5f)));
}

/** Get Imf::Rgba object.
*/
Imf::Rgba Get(const Imf::Array2D<Imf::Rgba>& pixels, int x, int y) {
  return pixels[y][x];
}

struct FloatRgb {
  float r, g, b;
};
FloatRgb Get(const CubeMapGen::ImageSurface& surface, int x, int y) {
  const float* p = surface.GetSurfaceTexelPtr(x, y);
  return{ p[0], p[1], p[2] };
}

/** Write the pixels to file.
*/
template<typename T>
Result WriteFile(const T& pixels, int w, int h, float strengthScale, const std::string& outfilename, Format outputFormat) {
  switch (outputFormat) {
  case FORMAT_PNG: {
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
		SetPixel(&row[x * 4], Get(pixels, x, y), strengthScale);
	  }
	  png_write_row(pngWriter.get(), row.data());
	}
	break;
  }

  case FORMAT_TGA: {
	EXR2PNG::Tga32Image image(w, h);
	for (int y = 0; y < h; ++y) {
	  for (int x = 0; x < w; ++x) {
		uint8_t buf[4];
		SetPixel(buf, Get(pixels, x, h - y - 1), strengthScale);
		image.SetPixel(x, y, buf[0], buf[1], buf[2], buf[3]);
	  }
	}
	std::unique_ptr<FILE, decltype(&fclose)> fp(fopen(outfilename.c_str(), "wb"), &fclose);
	if (!fp) {
	  std::cerr << "Error: can't open '" << outfilename << "' file." << std::endl;
	  return RESULT_ERROR;
	}
	const EXR2PNG::TgaHeader header = image.GetHeader();
	fwrite(&header, sizeof(EXR2PNG::TgaHeader), 1, fp.get());
	fwrite(image.buf.data(), 4, image.buf.size(), fp.get());
	break;
  }
  }
  return RESULT_SUCCESS;
}

/// Result of IsInsideImage().
struct IsInsideImageResult {
  /**
    true : all rect is inside. false: one or more rect isn't inside.
    false: One or more region is outside the source image,
	       or the left(or top) of the region is greater than the right(or bottom).
  */
  bool isInside;

  /**
    If 'isInside' is 'false', the index of first rect that isn't inside,
    otherwize it shouldn't to reference.
  */
  uint8_t indexOfInvalidRect;

  operator bool() const { return isInside; }
};

/** Check that all regions is inside the source image.

  @param rectArray  The cubemap region array.
  @param pixels     The source EXR image.

  @return true  All regions is inside the source image.
*/
IsInsideImageResult IsInsideImage(const ImageRectArray& rectArray, const Imf::Array2D<Imf::Rgba>& pixels) {
  const int w = pixels.width();
  const int h = pixels.height();
  uint8_t i = 0;
  for (auto& e : rectArray) {
	if (e.left < 0 || e.left >= e.right || e.right > w) {
	  return{ false, i };
	}
	if (e.top < 0 || e.top >= e.bottom || e.bottom > h) {
	  return{ false, i };
	}
	++i;
  }
  return{ true, 0 };
}

/** Create the cubemap image array from an EXR image.

  @param pixels     The source EXR image.
  @param rectArray  The cubemap region array.

  @return The cubemap image array.
*/
CubeMapGen::CubemapImageSurface CreateCubemapFromEXR(const Imf::Array2D<Imf::Rgba>& pixels, const ImageRectArray& rectArray) {
  CubeMapGen::CubemapImageSurface cubemap;
  for (int i = 0; i < 6; ++i) {
	const auto& rect = rectArray[i];
	CubeMapGen::ImageSurface& surface = cubemap[i];
	const int w = rect.right - rect.left;
	const int h = rect.bottom - rect.top;
	//std::cout << "cubemap surface " << i << "(" << w << "," << h << ")" << std::endl;
	surface.Init(w, h);
	float* p = surface.GetSurfacePtr();
	switch (rect.dir) {
	default:
	case ImageDir_Top:
	  for (int y = rect.top; y < rect.bottom; ++y) {
		for (int x = rect.left; x < rect.right; ++x) {
		  const auto& e = pixels[y][x];
		  p[CubeMapGen::ImageChannel_Red] = e.r;
		  p[CubeMapGen::ImageChannel_Green] = e.g;
		  p[CubeMapGen::ImageChannel_Blue] = e.b;
		  p += CubeMapGen::numberOfImageChannels;
		}
	  }
	  break;
	case ImageDir_Left:
	  for (int x = rect.left; x < rect.right; ++x) {
		for (int y = rect.bottom - 1; y >= rect.top; --y) {
		  const auto& e = pixels[y][x];
		  p[CubeMapGen::ImageChannel_Red] = e.r;
		  p[CubeMapGen::ImageChannel_Green] = e.g;
		  p[CubeMapGen::ImageChannel_Blue] = e.b;
		  p += CubeMapGen::numberOfImageChannels;
		}
	  }
	  break;
	case ImageDir_Bottom:
	  for (int y = rect.bottom - 1; y >= rect.top; --y) {
		for (int x = rect.right - 1; x >= rect.left; --x) {
		  const auto& e = pixels[y][x];
		  p[CubeMapGen::ImageChannel_Red] = e.r;
		  p[CubeMapGen::ImageChannel_Green] = e.g;
		  p[CubeMapGen::ImageChannel_Blue] = e.b;
		  p += CubeMapGen::numberOfImageChannels;
		}
	  }
	  break;
	case ImageDir_Right:
	  for (int x = rect.right - 1; x >= rect.left; --x) {
		for (int y = rect.top; y < rect.bottom; ++y) {
		  const auto& e = pixels[y][x];
		  p[CubeMapGen::ImageChannel_Red] = e.r;
		  p[CubeMapGen::ImageChannel_Green] = e.g;
		  p[CubeMapGen::ImageChannel_Blue] = e.b;
		  p += CubeMapGen::numberOfImageChannels;
		}
	  }
	  break;
	}
	assert(static_cast<size_t>(p - surface.GetSurfacePtr()) <= surface.buf.size());
  }
  return cubemap;
}

/** The entry point.

  @param argc  argument count.
  @param argv  argument vector.

  @retval RESULT_SUCCESS  the conversion is success.
  @retval RESULT_ERROR    the conversion is failure.
*/
int main(int argc, char** argv) {
  bool verbose = false;
  const char* infilename = nullptr;
  std::string outfilename;
  float filterAngle = 1.0f;
  float strengthScale = 1.0f;
  int miplevel = 1;
  Format outputFormat = FORMAT_PNG;
  ConversionMode conversionMode = CONVERSIONMODE_SINGLE;
  for (int i = 1; i < argc; ++i) {
	if (argv[i][0] == '-' || argv[i][0] == '+') {
	  if ((argv[i][1] == 'x' || argv[i][1] == 'y' || argv[i][1] == 'z') && (argc >= i + 5)) {
		char* p;
		const int l = strtol(argv[i + 1], &p, 10);
		if (*p) {
		  std::cerr << "Error: Invalid number in " << argv[i] << ": " << argv[i + 1] << std::endl;
		  return RESULT_ERROR;
		}
		const int t = strtol(argv[i + 2], &p, 10);
		if (*p) {
		  std::cerr << "Error: Invalid number in " << argv[i] << ": " << argv[i + 2] << std::endl;
		  return RESULT_ERROR;
		}
		const int r = strtol(argv[i + 3], &p, 10);
		if (*p) {
		  std::cerr << "Error: Invalid number in " << argv[i] << ": " << argv[i + 3] << std::endl;
		  return RESULT_ERROR;
		}
		const int b = strtol(argv[i + 4], &p, 10);
		if (*p) {
		  std::cerr << "Error: Invalid number in " << argv[i] << ": " << argv[i + 4] << std::endl;
		  return RESULT_ERROR;
		}
		ImageDir dir;
		switch (argv[i + 5][0]) {
		case 'l': dir = ImageDir_Left; break;
		case 't': dir = ImageDir_Top; break;
		case 'r': dir = ImageDir_Right; break;
		case 'b': dir = ImageDir_Bottom; break;
		default:
		  std::cerr << "Error: Invalid direction in " << argv[i] << ": " << argv[i + 5] << std::endl;
		  return RESULT_ERROR;
		}
		const int offset = argv[i][0] == '+' ? 0 : 1;
		switch (argv[i][1]) {
		case 'x': cubemapRectArray[CubeMapGen::CP_FACE_X_POS + offset] = { l, t, r, b, dir }; break;
		case 'y': cubemapRectArray[CubeMapGen::CP_FACE_Y_POS + offset] = { l, t, r, b, dir }; break;
		case 'z': cubemapRectArray[CubeMapGen::CP_FACE_Z_POS + offset] = { l, t, r, b, dir }; break;
		}
		if (verbose) {
		  std::cout << "Set " << argv[i] << std::endl;
		}
		i += 5;
		continue;
	  }
	}
	if (argv[i][0] == '-') {
	  if ((argv[i][1] == 's' || argv[i][1] == 'S') && (argc >= i + 1)) {
		strengthScale = std::max(0.0f, static_cast<float>(atof(argv[i + 1])));
		++i;
	  } else if ((argv[i][1] == 'a'|| argv[i][1] == 'A') && (argc >= i + 1)) {
		filterAngle = std::max(0.0f, std::min(90.0f, static_cast<float>(atof(argv[i + 1]))));
		++i;
	  } else if ((argv[i][1] == 'm' || argv[i][1] == 'm') && (argc >= i + 1)) {
		miplevel = std::min(maxMipLevels, std::max(1, atoi(argv[i + 1])));
		++i;
	  } else if ((argv[i][1] == 'f' || argv[i][1] == 'F') && (argc >= i + 1)) {
		if (strcmp("png", argv[i + 1]) == 0) {
		  outputFormat = FORMAT_PNG;
		} else if (strcmp("tga", argv[i + 1]) == 0) {
		  outputFormat = FORMAT_TGA;
		} else {
		  std::cout << "Error: Unknown format: " << argv[i + 1] << std::endl;
		  std::cout << "       -f option only accepts 'png' or 'tga'" << std::endl;
		  return RESULT_ERROR;
		}
		++i;
	  } else if ((argv[i][1] == 'c' || argv[i][1] == 'C') && (argc >= i + 1)) {
		if (strcmp("filtered", argv[i + 1]) == 0) {
		  conversionMode = CONVERSIONMODE_CUBEMAP_FILTERED;
		} else if (strcmp("irradiance", argv[i + 1]) == 0) {
		  conversionMode = CONVERSIONMODE_CUBEMAP_IRRADIANCE;
		} else if (strcmp("none", argv[i + 1]) == 0) {
		  conversionMode = CONVERSIONMODE_CUBEMAP_NONFILTERD;
		} else {
		  std::cout << "Error: Unknown format in -c: " << argv[i + 1] << std::endl;
		  return RESULT_ERROR;
		}
		++i;
	  }
	  continue;
	}
	if (!infilename) {
	  if (verbose) {
		std::cout << "infile: " << argv[i] << std::endl;
	  }
      infilename = argv[i];
    } else if (outfilename.empty()) {
	  if (verbose) {
		std::cout << "outfile: " << argv[i] << std::endl;
	  }
	  outfilename = argv[i];
      break;
    }
  }

  if (!infilename) {
	PrintUsage();
	return RESULT_SUCCESS;
  }

  if (conversionMode == CONVERSIONMODE_SINGLE) {
	if (outfilename.empty()) {
	  outfilename = infilename;
	  auto dotPos = outfilename.find_last_of('.');
	  if (dotPos != std::string::npos) {
		outfilename.erase(dotPos, std::string::npos);
	  }
	  if (outputFormat == FORMAT_PNG) {
		outfilename += ".png";
	  } else {
		outfilename += ".tga";
	  }
	}

	Imf::Array2D<Imf::Rgba> pixels;
	int w, h;
	if (ReadRgbaExrFile(infilename, &pixels, &w, &h) == RESULT_ERROR) {
	  return RESULT_ERROR;
	}
	return WriteFile(pixels, w, h, strengthScale, outfilename, outputFormat);
  } else {
	if (outfilename.empty()) {
	  outfilename = infilename;
	  auto dotPos = outfilename.find_last_of('.');
	  if (dotPos != std::string::npos) {
		outfilename.erase(dotPos, std::string::npos);
	  }
	}

	try {
	  Imf::Array2D<Imf::Rgba> pixels;
	  int w, h;
	  if (ReadRgbaExrFile(infilename, &pixels, &w, &h) == RESULT_ERROR) {
		return RESULT_ERROR;
	  }
	  const auto result = IsInsideImage(cubemapRectArray, pixels);
	  if (!result) {
		static const char optionName[][3] = { "+x", "-x", "+y", "-y", "+z", "-z" };
		const int i = result.indexOfInvalidRect;
		std::cout << "Error: Invalid region in " << optionName[i] <<
		  ":(" << cubemapRectArray[i].left << ", " << cubemapRectArray[i].top <<
		  "-" << cubemapRectArray[i].right << ", " << cubemapRectArray[i].bottom << ")" << std::endl;
		return RESULT_ERROR;
	  }

	  CubeMapGen::CubemapImageSurface srcCubemap = CreateCubemapFromEXR(pixels, cubemapRectArray);
	  std::vector<CubeMapGen::CubemapImageSurface> destCubemapList;
	  destCubemapList.resize(miplevel);
	  int width = srcCubemap[0].m_Width;
	  int height = srcCubemap[0].m_Height;
	  for (auto& destCubemap : destCubemapList) {
		for (auto& e : destCubemap) {
		  e.Init(width, height);
		}
		width = std::max(1, width / 2);
		height = std::max(1, height / 2);
	  }
	  float specularPower = 2048.0f;
	  switch (conversionMode) {
	  default:
	  case CONVERSIONMODE_CUBEMAP_FILTERED: {
		CubeMapGen::CubemapImageSurface* pSrc = &srcCubemap;
		for (int i = 0; i < miplevel; ++i) {
		  CubeMapGen::CubemapImageSurface* pDest = &destCubemapList[i];
		  CubeMapGen::FilterCubeSurfaces(
			pSrc->data(), pDest->data(),
			filterAngle,
			CubeMapGen::CP_FILTER_TYPE_COSINE_POWER,
			true,
			0, CubeMapGen::numberOfCubemapFaces - 1,
			specularPower
		  );
		  filterAngle *= 2.0f;
		  specularPower *= 0.25f;
		  pSrc = pDest;
		}
		break;
	  }
	  case CONVERSIONMODE_CUBEMAP_IRRADIANCE:
		CubeMapGen::SHFilterCubeMap(srcCubemap.data(), destCubemapList[0].data());
		break;
	  case CONVERSIONMODE_CUBEMAP_NONFILTERD:
		destCubemapList[0] = srcCubemap;
		break;
	  }

	  static const char* const postfixArray[] = { "_px", "_nx", "_py", "_ny", "_pz", "_nz" };
	  for (int m = 0; m < miplevel; ++m) {
		for (int i = 0; i < CubeMapGen::numberOfCubemapFaces; ++i) {
		  std::string filename = outfilename + static_cast<char>('1' + m) + postfixArray[i];
		  if (outputFormat == FORMAT_PNG) {
			filename += ".png";
		  } else {
			filename += ".tga";
		  }
		  if (WriteFile(destCubemapList[m][i], destCubemapList[m][i].m_Width, destCubemapList[m][i].m_Height, strengthScale, filename, outputFormat) == RESULT_ERROR) {
			return RESULT_ERROR;
		  }
		}
	  }
	}
	catch (std::exception& e) {
	  std::cerr << "Error: " << e.what() << std::endl;
	  return RESULT_ERROR;
	}
  }
  return RESULT_SUCCESS;
}
