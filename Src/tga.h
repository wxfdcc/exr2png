/**
  @file tga.h

  This code creates the uncompressed TARGA32 image.
*/
#ifndef EXR2PNG_TGA_H_INCLUDED
#define EXR2PNG_TGA_H_INCLUDED
#include <cstdint>
#include <vector>

namespace EXR2PNG {

  /** The little endian 16 bit unsigned integer.

    TGA header should be packed format.
	This class provides the integer represented packed format.
  */
  struct ShortInt {
	uint8_t p[2];
	ShortInt& operator=(uint16_t value) {
	  p[0] = static_cast<uint8_t>(value & 0xff);
	  p[1] = static_cast<uint8_t>((value >> 8) & 0xff);
	  return *this;
	}
	operator uint16_t() const {
	  return static_cast<uint16_t>(p[0]) + static_cast<uint16_t>(p[1]) << 8;
	}
  };

  /** The header of TARGA image file.
  */
  struct TgaHeader {
	uint8_t  idlength; ///< Number of characters in identification field following this header.
	uint8_t  colourmaptype; ///< 0 means no color map is included. 1 means a color map is included,
	uint8_t  datatypecode; ///< Image type code.
	ShortInt colourmaporigin;
	ShortInt colourmaplength;
	uint8_t  colourmapdepth;
	ShortInt x_origin; ///< Integer ( lo-hi ) X coordinate of the lower left corner of the image.
	ShortInt y_origin; ///< Integer ( lo-hi ) Y coordinate of the lower left corner of the image.
	ShortInt width; ///< Integer ( lo-hi ) width of the image in pixels.
	ShortInt height; ///< Integer ( lo-hi ) height of the image in pixels.
	uint8_t  bitsperpixel; ///< Number of bits in a pixel.
	uint8_t  imagedescriptor; ///< bit 3-0: number of attribute bits associated with each pixel. bit 7-4: zero.
  };
  static_assert(sizeof(TgaHeader) == 18, "TGA Header has wrong size. it should be 18 bytes.");

  /** Set the information of TARGA32 to TARGA header.

    @param tga  TARGA header.
	@param w    The pixel width of the image.
	@param h    The pixel height of the image.
  */
  void SetUncompressedRgbInfo(TgaHeader& tga, int w, int h) {
	tga.idlength = 0;
	tga.colourmaptype = 0; // 0: no color map.
	tga.datatypecode = 2; // 2: Uncompressed, RGB images.

	// Ignored if color map type is 0.
	tga.colourmaporigin = 0;
	tga.colourmaplength = 0;
	tga.colourmapdepth = 0;

	tga.x_origin = 0;
	tga.y_origin = 0;
	tga.width = w;
	tga.height = h;
	tga.bitsperpixel = 32;
	tga.imagedescriptor = 8;
  }

  /** The TARGA32 image format data class.
  */
  struct Tga32Image {
	std::vector<uint32_t> buf;
	int width;
	int height;

	/** Constructor.

	  @param w  The pixel width of the image.
	  @param h  The pixel height of the image.
	*/
	Tga32Image(int w, int h) {
	  width = w;
	  height = h;
	  buf.resize(w * h);
	}

	/** Set the pixel data.

	  @param x  The X offset of target pixel.
	  @param y  The Y offset of target pixel.
	  @param r  The pixel color of red.
	  @param g  The pixel color of green.
	  @param b  The pixel color of blue.
	  @param a  The pixel attribute of alpha.

	  @retval true  Success.
	  @retval false Failure. The parameter x or y is in the out of image.
	*/
	bool SetPixel(int x, int y, uint32_t r, uint32_t g, uint32_t b, uint32_t a) {
#ifndef NDEBUG
	  if (x < 0 || x >= width) {
		return false;
	  }
	  if (y < 0 || y >= height) {
		return false;
	  }
#endif // NDEBUG
	  buf[y * width + x] = b + (g << 8) + (r << 16) + (a << 24);
	  return true;
	}

	/** Get TARGA32 Header information.
	*/
	TgaHeader GetHeader() const {
	  TgaHeader header;
	  SetUncompressedRgbInfo(header, width, height);
	  return header;
	}
  };

} // namespace EXR2PNG

#endif // EXR2PNG_TGA_H_INCLUDED