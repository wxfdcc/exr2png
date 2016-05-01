/**
  @file cubemapgen.h
*/
#ifndef EXR2PNG_CUBEMAPGEN_H_INCLUDED
#define EXR2PNG_CUBEMAPGEN_H_INCLUDED
#include <cstdint>
#include <vector>
#include <array>

namespace CubeMapGen {

/// used to index cube faces
enum Face {
  CP_FACE_X_POS,
  CP_FACE_X_NEG,
  CP_FACE_Y_POS,
  CP_FACE_Y_NEG,
  CP_FACE_Z_POS,
  CP_FACE_Z_NEG,
};
static const int32_t numberOfCubemapFaces = 6;

enum ImageChannel {
  ImageChannel_Red,
  ImageChannel_Green,
  ImageChannel_Blue,
};
static const int32_t numberOfImageChannels = 3;

/// Filter type.
enum FilterType {
  CP_FILTER_TYPE_DISC,
  CP_FILTER_TYPE_CONE,
  CP_FILTER_TYPE_COSINE,
  CP_FILTER_TYPE_ANGULAR_GAUSSIAN,
  CP_FILTER_TYPE_COSINE_POWER,
};

template<int C>
struct BaseImageSurface {
  int32_t m_Width;
  int32_t m_Height;
  static const int32_t m_NumChannels = C;
  std::vector<float> buf;

  void Init(int32_t w, int32_t h) {
	m_Width = w;
	m_Height = h;
	buf.resize(w * h * m_NumChannels);
  }

  float* GetSurfacePtr() { return buf.data(); }
  const float* GetSurfacePtr() const { return buf.data(); }

  float* GetSurfaceTexelPtr(int32_t u, int32_t v) {
	return const_cast<float*>(const_cast<const BaseImageSurface*>(this)->GetSurfaceTexelPtr(u, v));
  }
  const float* GetSurfaceTexelPtr(int32_t u, int32_t v) const {
	return buf.data() + ((m_Width * v) + u) * m_NumChannels;
  }
};
typedef BaseImageSurface<numberOfImageChannels> ImageSurface;
typedef std::array<ImageSurface, numberOfCubemapFaces> CubemapImageSurface;
typedef BaseImageSurface<4> NormalSurface;

void FilterCubeSurfaces(
  const ImageSurface *a_SrcCubeMap,
  ImageSurface *a_DstCubeMap,
  float a_FilterConeAngle = 1.0f,
  int32_t a_FilterType = CP_FILTER_TYPE_COSINE_POWER,
  bool a_bUseSolidAngle = true,
  int32_t a_FaceIdxStart = 0,
  int32_t a_FaceIdxEnd = 5,
  float a_SpecularPower = 2048.0f
);

/** Filter the cubemap using the SH filter.

  @param  inputSurface           The pointer to the array of six input color surface.
                                 @sa numberOfCubemapFaces

  @param  outputSurface          The pointer to the array of output color surface array.
                                 It must have been initialized.
								 @sa numberOfCubemapFaces

  @param  useSolidAngleWeighing  - true  Use the solid angle weight.
                                 - false Unuse.
*/
void SHFilterCubeMap(
  const ImageSurface* m_InputSurface,
  ImageSurface* m_OutputSurface,
  bool a_bUseSolidAngleWeighting = true
);

} // namespace CubeMapGen

#endif // EXR2PNG_CUBEMAPGEN_H_INCLUDED