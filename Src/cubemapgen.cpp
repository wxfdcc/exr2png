/**
  @file cubemapgen.cpp
*/
#include "cubemapgen.h"
#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>

namespace CubeMapGen {

#define VM_ADD3_UNTYPED(d, sa, sb) { d[0]=sa[0]+sb[0]; d[1]=sa[1]+sb[1]; d[2]=sa[2]+sb[2]; }
#define VM_ADD3(d, sa, sb) VM_ADD3_UNTYPED(((float*)(d)), ((float*)(sa)), ((float*)(sb)))

#define VM_SCALE3_UNTYPED(d, s, f) {d[0]=s[0]*f; d[1]=s[1]*f; d[2]=s[2]*f; }
#define VM_SCALE3(d, s, f) VM_SCALE3_UNTYPED(((float*)(d)), ((float*)(s)), ((float)(f)))

#define VM_DOTPROD3_UNTYPED(sa, sb) (sa[0]*sb[0]+ sa[1]*sb[1]+ sa[2]*sb[2])
#define VM_DOTPROD3(sa, sb) VM_DOTPROD3_UNTYPED(((float*)(sa)), ((float*)(sb)))

#define VM_NORM3_UNTYPED(d, s) {double __idsq; __idsq=1.0/sqrt(VM_DOTPROD3_UNTYPED(s,s)); d[0]=s[0]*__idsq; d[1]=s[1]*__idsq; d[2]=s[2]*__idsq; }
#define VM_NORM3_UNTYPED_F32(d, s) {float __idsq; __idsq=static_cast<float>(1.0/sqrt(VM_DOTPROD3_UNTYPED(s,s))); d[0]=s[0]*__idsq; d[1]=s[1]*__idsq; d[2]=s[2]*__idsq; }
#define VM_NORM3(d, s)  VM_NORM3_UNTYPED_F32(((float*)(d)), ((float*)(s)))

#define VM_ABS3_UNTYPED(d, s) { d[0] = fabs(s[0]); d[1] = fabs(s[1]); d[2] = fabs(s[2]); }
#define VM_ABS3(d, s) VM_ABS3_UNTYPED(((float*)(d)), ((float*)(s)))

const int32_t m_NumFilterLUTEntries = 4096;
const int32_t g_numChannels = 3;
const double PI_DOUBLE = 3.14159265358979323846;
const float PI_FLOAT = static_cast<float>(PI_DOUBLE);

//------------------------------------------------------------------------------
// D3D cube map face specification
//   mapping from 3D x,y,z cube map lookup coordinates 
//   to 2D within face u,v coordinates
//
//   --------------------> U direction 
//   |                   (within-face texture space)
//   |         _____
//   |        |     |
//   |        | +Y  |
//   |   _____|_____|_____ _____
//   |  |     |     |     |     |
//   |  | -X  | +Z  | +X  | -Z  |
//   |  |_____|_____|_____|_____|
//   |        |     |
//   |        | -Y  |
//   |        |_____|
//   |
//   v   V direction
//      (within-face texture space)
//------------------------------------------------------------------------------

//used to index image edges
// NOTE.. the actual number corresponding to the edge is important
//  do not change these, or the code will break
//
// CP_EDGE_LEFT   is u = 0
// CP_EDGE_RIGHT  is u = width-1
// CP_EDGE_TOP    is v = 0
// CP_EDGE_BOTTOM is v = height-1
#define CP_EDGE_LEFT   0
#define CP_EDGE_RIGHT  1
#define CP_EDGE_TOP    2
#define CP_EDGE_BOTTOM 3

//information about cube maps neighboring face after traversing
// across an edge
struct CPCubeMapNeighbor
{
  uint8_t m_Face;    //index of neighboring face
  uint8_t m_Edge;    //edge in neighboring face that abuts this face
};

//Information about neighbors and how texture coorrdinates change across faces 
//  in ORDER of left, right, top, bottom (e.g. edges corresponding to u=0, 
//  u=1, v=0, v=1 in the 2D coordinate system of the particular face.
//Note this currently assumes the D3D cube face ordering and orientation
CPCubeMapNeighbor sg_CubeNgh[6][4] =
{
  //XPOS face
  { { CP_FACE_Z_POS, CP_EDGE_RIGHT },
  { CP_FACE_Z_NEG, CP_EDGE_LEFT },
  { CP_FACE_Y_POS, CP_EDGE_RIGHT },
  { CP_FACE_Y_NEG, CP_EDGE_RIGHT } },
  //XNEG face
  { { CP_FACE_Z_NEG, CP_EDGE_RIGHT },
  { CP_FACE_Z_POS, CP_EDGE_LEFT },
  { CP_FACE_Y_POS, CP_EDGE_LEFT },
  { CP_FACE_Y_NEG, CP_EDGE_LEFT } },
  //YPOS face
  { { CP_FACE_X_NEG, CP_EDGE_TOP },
  { CP_FACE_X_POS, CP_EDGE_TOP },
  { CP_FACE_Z_NEG, CP_EDGE_TOP },
  { CP_FACE_Z_POS, CP_EDGE_TOP } },
  //YNEG face
  { { CP_FACE_X_NEG, CP_EDGE_BOTTOM },
  { CP_FACE_X_POS, CP_EDGE_BOTTOM },
  { CP_FACE_Z_POS, CP_EDGE_BOTTOM },
  { CP_FACE_Z_NEG, CP_EDGE_BOTTOM } },
  //ZPOS face
  { { CP_FACE_X_NEG, CP_EDGE_RIGHT },
  { CP_FACE_X_POS, CP_EDGE_LEFT },
  { CP_FACE_Y_POS, CP_EDGE_BOTTOM },
  { CP_FACE_Y_NEG, CP_EDGE_TOP } },
  //ZNEG face
  { { CP_FACE_X_POS, CP_EDGE_RIGHT },
  { CP_FACE_X_NEG, CP_EDGE_LEFT },
  { CP_FACE_Y_POS, CP_EDGE_TOP },
  { CP_FACE_Y_NEG, CP_EDGE_BOTTOM } }
};

//3x2 matrices that map cube map indexing vectors in 3d 
// (after face selection and divide through by the 
//  _ABSOLUTE VALUE_ of the max coord)
// into NVC space
//Note this currently assumes the D3D cube face ordering and orientation
#define CP_UDIR     0
#define CP_VDIR     1
#define CP_FACEAXIS 2

float sgFace2DMapping[6][3][3] = {
  //XPOS face
  { { 0,  0, -1 },   //u towards negative Z
  { 0, -1,  0 },   //v towards negative Y
  { 1,  0,  0 } },  //pos X axis  
					//XNEG face
  { { 0,  0,  1 },   //u towards positive Z
  { 0, -1,  0 },   //v towards negative Y
  { -1,  0,  0 } },  //neg X axis       
					 //YPOS face
  { { 1, 0, 0 },     //u towards positive X
  { 0, 0, 1 },     //v towards positive Z
  { 0, 1 , 0 } },   //pos Y axis  
					//YNEG face
  { { 1, 0, 0 },     //u towards positive X
  { 0, 0 , -1 },   //v towards negative Z
  { 0, -1 , 0 } },  //neg Y axis  
					//ZPOS face
  { { 1, 0, 0 },     //u towards positive X
  { 0, -1, 0 },    //v towards negative Y
  { 0, 0,  1 } },   //pos Z axis  
					//ZNEG face
  { { -1, 0, 0 },    //u towards negative X
  { 0, -1, 0 },    //v towards negative Y
  { 0, 0, -1 } },   //neg Z axis  
};

struct BoundingBox {
  int32_t m_minCoord[3];  //upper left back corner
  int32_t m_maxCoord[3];  //lower right front corner

  /// Clear bounding box extents
  void Clear() {
	m_minCoord[0] = std::numeric_limits<int32_t>::max();
	m_minCoord[1] = std::numeric_limits<int32_t>::max();
	m_minCoord[2] = std::numeric_limits<int32_t>::max();
	m_maxCoord[0] = std::numeric_limits<int32_t>::min();
	m_maxCoord[1] = std::numeric_limits<int32_t>::min();
	m_maxCoord[2] = std::numeric_limits<int32_t>::min();
  }

  /// Augment bounding box extents by specifying point to include in bounding box
  void Augment(int32_t aX, int32_t aY, int32_t aZ) {
	m_minCoord[0] = std::min(m_minCoord[0], aX);
	m_minCoord[1] = std::min(m_minCoord[1], aY);
	m_minCoord[2] = std::min(m_minCoord[2], aZ);
	m_maxCoord[0] = std::max(m_maxCoord[0], aX);
	m_maxCoord[1] = std::max(m_maxCoord[1], aY);
	m_maxCoord[2] = std::max(m_maxCoord[2], aZ);
  }

  /// Clamp minimum values in bbox to be no larger than aX, aY, aZ
  void ClampMin(int32_t aX, int32_t aY, int32_t aZ) {
	m_minCoord[0] = std::max(m_minCoord[0], aX);
	m_minCoord[1] = std::max(m_minCoord[1], aY);
	m_minCoord[2] = std::max(m_minCoord[2], aZ);
  }

  /// Clamp maximum values in bbox to be no larger than aX, aY, aZ
  void ClampMax(int32_t aX, int32_t aY, int32_t aZ) {
	m_maxCoord[0] = std::min(m_maxCoord[0], aX);
	m_maxCoord[1] = std::min(m_maxCoord[1], aY);
	m_maxCoord[2] = std::min(m_maxCoord[2], aZ);
  }
};

//--------------------------------------------------------------------------------------
// Convert cubemap face texel coordinates and face idx to 3D vector
// note the U and V coords are integer coords and range from 0 to size-1
//  this routine can be used to generate a normalizer cube map
//--------------------------------------------------------------------------------------
void TexelCoordToVect(int32_t a_FaceIdx, float a_U, float a_V, int32_t a_Size, float *a_XYZ)
{
  float nvcU, nvcV;
  float tempVec[3];

  // Change from original AMD code
  // transform from [0..res - 1] to [- (1 - 1 / res) .. (1 - 1 / res)]
  // + 0.5f is for texel center addressing
  nvcU = (2.0f * ((float)a_U + 0.5f) / (float)a_Size) - 1.0f;
  nvcV = (2.0f * ((float)a_V + 0.5f) / (float)a_Size) - 1.0f;

  // Code from Nvtt:
  // http://code.google.com/p/nvidia-texture-tools/source/browse/trunk/src/nvtt/CubeSurface.cpp
  float a = powf(float(a_Size), 2.0f) / powf(float(a_Size - 1), 3.0f);
  nvcU = a * powf(nvcU, 3) + nvcU;
  nvcV = a * powf(nvcV, 3) + nvcV;

  // Get current vector
  //generate x,y,z vector (xform 2d NVC coord to 3D vector)
  //U contribution
  VM_SCALE3(a_XYZ, sgFace2DMapping[a_FaceIdx][CP_UDIR], nvcU);
  //V contribution
  VM_SCALE3(tempVec, sgFace2DMapping[a_FaceIdx][CP_VDIR], nvcV);
  VM_ADD3(a_XYZ, tempVec, a_XYZ);
  //add face axis
  VM_ADD3(a_XYZ, sgFace2DMapping[a_FaceIdx][CP_FACEAXIS], a_XYZ);
  //normalize vector
  VM_NORM3(a_XYZ, a_XYZ);
}

//--------------------------------------------------------------------------------------
// Convert 3D vector to cubemap face texel coordinates and face idx 
// note the U and V coords are integer coords and range from 0 to size-1
//  this routine can be used to generate a normalizer cube map
//
// returns face IDX and texel coords
//--------------------------------------------------------------------------------------
/*
Mapping Texture Coordinates to Cube Map Faces
Because there are multiple faces, the mapping of texture coordinates to positions on cube map faces
is more complicated than the other texturing targets.  The EXT_texture_cube_map extension is purposefully
designed to be consistent with DirectX 7's cube map arrangement.  This is also consistent with the cube
map arrangement in Pixar's RenderMan package.
For cube map texturing, the (s,t,r) texture coordinates are treated as a direction vector (rx,ry,rz)
emanating from the center of a cube.  (The q coordinate can be ignored since it merely scales the vector
without affecting the direction.) At texture application time, the interpolated per-fragment (s,t,r)
selects one of the cube map face's 2D mipmap sets based on the largest magnitude coordinate direction
the major axis direction). The target column in the table below explains how the major axis direction
maps to the 2D image of a particular cube map target.

major axis
direction     target                              sc     tc    ma
----------    ---------------------------------   ---    ---   ---
+rx          GL_TEXTURE_CUBE_MAP_POSITIVE_X_EXT   -rz    -ry   rx
-rx          GL_TEXTURE_CUBE_MAP_NEGATIVE_X_EXT   +rz    -ry   rx
+ry          GL_TEXTURE_CUBE_MAP_POSITIVE_Y_EXT   +rx    +rz   ry
-ry          GL_TEXTURE_CUBE_MAP_NEGATIVE_Y_EXT   +rx    -rz   ry
+rz          GL_TEXTURE_CUBE_MAP_POSITIVE_Z_EXT   +rx    -ry   rz
-rz          GL_TEXTURE_CUBE_MAP_NEGATIVE_Z_EXT   -rx    -ry   rz

Using the sc, tc, and ma determined by the major axis direction as specified in the table above,
an updated (s,t) is calculated as follows
s   =   ( sc/|ma| + 1 ) / 2
t   =   ( tc/|ma| + 1 ) / 2
If |ma| is zero or very nearly zero, the results of the above two equations need not be defined
(though the result may not lead to GL interruption or termination).  Once the cube map face's 2D mipmap
set and (s,t) is determined, texture fetching and filtering proceeds like standard OpenGL 2D texturing.
*/
// Note this method return U and V in range from 0 to size-1
void VectToTexelCoord(const float *a_XYZ, int32_t a_Size, int32_t *a_FaceIdx, int32_t *a_U, int32_t *a_V)
{
  float nvcU, nvcV;
  float absXYZ[3];
  float maxCoord;
  float onFaceXYZ[3];
  int32_t   faceIdx;
  int32_t   u, v;

  //absolute value 3
  VM_ABS3(absXYZ, a_XYZ);

  if ((absXYZ[0] >= absXYZ[1]) && (absXYZ[0] >= absXYZ[2])) {
	maxCoord = absXYZ[0];

	if (a_XYZ[0] >= 0) //face = XPOS
	{
	  faceIdx = CP_FACE_X_POS;
	} else {
	  faceIdx = CP_FACE_X_NEG;
	}
  } else if ((absXYZ[1] >= absXYZ[0]) && (absXYZ[1] >= absXYZ[2])) {
	maxCoord = absXYZ[1];

	if (a_XYZ[1] >= 0) //face = XPOS
	{
	  faceIdx = CP_FACE_Y_POS;
	} else {
	  faceIdx = CP_FACE_Y_NEG;
	}
  } else  // if( (absXYZ[2] > absXYZ[0]) && (absXYZ[2] > absXYZ[1]) )
  {
	maxCoord = absXYZ[2];

	if (a_XYZ[2] >= 0) //face = XPOS
	{
	  faceIdx = CP_FACE_Z_POS;
	} else {
	  faceIdx = CP_FACE_Z_NEG;
	}
  }

  //divide through by max coord so face vector lies on cube face
  VM_SCALE3(onFaceXYZ, a_XYZ, 1.0f / maxCoord);
  nvcU = VM_DOTPROD3(sgFace2DMapping[faceIdx][CP_UDIR], onFaceXYZ);
  nvcV = VM_DOTPROD3(sgFace2DMapping[faceIdx][CP_VDIR], onFaceXYZ);

  // Modify original AMD code to return value from 0 to Size - 1
  u = (int32_t)floor((a_Size - 1) * 0.5f * (nvcU + 1.0f));
  v = (int32_t)floor((a_Size - 1) * 0.5f * (nvcV + 1.0f));

  *a_FaceIdx = faceIdx;
  *a_U = u;
  *a_V = v;
}

//--------------------------------------------------------------------------------------
// gets texel ptr in a cube map given a direction vector, and an array of 
//  CImageSurfaces that represent the cube faces.
//   
//--------------------------------------------------------------------------------------
const float* GetCubeMapTexelPtr(const float* a_XYZ, const ImageSurface *a_Surface)
{
  int32_t u, v, faceIdx;

  //get face idx and u, v texel coordinate in face
  VectToTexelCoord(a_XYZ, a_Surface[0].m_Width, &faceIdx, &u, &v);

  return(a_Surface[faceIdx].GetSurfaceTexelPtr(u, v));
}

static float AreaElement(float x, float y)
{
  return atan2(x * y, sqrt(x * x + y * y + 1));
}

float TexelCoordSolidAngle(int32_t a_FaceIdx, float a_U, float a_V, int32_t a_Size)
{
  // transform from [0..res - 1] to [- (1 - 1 / res) .. (1 - 1 / res)]
  // (+ 0.5f is for texel center addressing)
  float U = (2.0f * ((float)a_U + 0.5f) / (float)a_Size) - 1.0f;
  float V = (2.0f * ((float)a_V + 0.5f) / (float)a_Size) - 1.0f;

  // Shift from a demi texel, mean 1.0f / a_Size with U and V in [-1..1]
  float InvResolution = 1.0f / a_Size;

  // U and V are the -1..1 texture coordinate on the current face.
  // Get projected area for this texel
  float x0 = U - InvResolution;
  float y0 = V - InvResolution;
  float x1 = U + InvResolution;
  float y1 = V + InvResolution;
  float SolidAngle = AreaElement(x0, y0) - AreaElement(x0, y1) - AreaElement(x1, y0) + AreaElement(x1, y1);

  return SolidAngle;
}

//--------------------------------------------------------------------------------------
//Builds a normalizer cubemap, with the texels solid angle stored in the fourth component
//
//Takes in a cube face size, and an array of 6 surfaces to write the cube faces into
//
//Note that this normalizer cube map stores the vectors in unbiased -1 to 1 range.
// if _bx2 style scaled and biased vectors are needed, uncomment the SCALE and BIAS
// below
//--------------------------------------------------------------------------------------
void BuildNormalizerSolidAngleCubemap(int32_t a_Size, NormalSurface *a_Surface)
{
  int32_t iCubeFace, u, v;

  //iterate over cube faces
  for (iCubeFace = 0; iCubeFace < 6; iCubeFace++) {
	a_Surface[iCubeFace].Init(a_Size, a_Size);  //First three channels for norm cube, and last channel for solid angle

	//fast texture walk, build normalizer cube map
	float *texelPtr = a_Surface[iCubeFace].GetSurfacePtr();

	for (v = 0; v < a_Surface[iCubeFace].m_Height; v++) {
	  for (u = 0; u < a_Surface[iCubeFace].m_Width; u++) {
		// SL_BEGIN
		TexelCoordToVect(iCubeFace, (float)u, (float)v, a_Size, texelPtr);
		// SL END
		//VM_SCALE3(texelPtr, texelPtr, 0.5f);
		//VM_BIAS3(texelPtr, texelPtr, 0.5f);

		*(texelPtr + 3) = TexelCoordSolidAngle(iCubeFace, (float)u, (float)v, a_Size);

		texelPtr += a_Surface[iCubeFace].m_NumChannels;
	  }
	}
  }
}

//--------------------------------------------------------------------------------------
//build filter lookup table
//
//--------------------------------------------------------------------------------------
std::vector<float> BuildAngleWeightLUT(int32_t a_FilterType, float a_FilterAngle)
{
  std::vector<float> m_FilterLUT;
  m_FilterLUT.resize(m_NumFilterLUTEntries);

  // note that CP_FILTER_TYPE_DISC weights all taps equally and does not need a lookup table    
  if (a_FilterType == CP_FILTER_TYPE_CONE) {
	//CP_FILTER_TYPE_CONE is a cone centered around the center tap and falls off to zero 
	//  over the filtering radius
	float filtAngleRad = a_FilterAngle * PI_FLOAT / 180.0f;
	for (int32_t iLUTEntry = 0; iLUTEntry<m_NumFilterLUTEntries; iLUTEntry++) {
	  float angle = acos((float)iLUTEntry / (float)(m_NumFilterLUTEntries - 1));
	  float filterVal = (filtAngleRad - angle) / filtAngleRad;
	  if (filterVal < 0) {
		filterVal = 0;
	  }
	  //note that gaussian is not weighted by 1.0 / (sigma* sqrt(2 * PI)) seen as weights
	  // weighted tap accumulation in filters is divided by sum of weights
	  m_FilterLUT[iLUTEntry] = filterVal;
	}
  } else if (a_FilterType == CP_FILTER_TYPE_ANGULAR_GAUSSIAN) {
	//fit 3 standard deviations within angular extent of filter
	float stdDev = (a_FilterAngle * PI_FLOAT / 180.0f) / 3.0f;
	float inv2Variance = 1.0f / (2.0f * stdDev * stdDev);
	for (int32_t iLUTEntry = 0; iLUTEntry<m_NumFilterLUTEntries; iLUTEntry++) {
	  float angle = acos((float)iLUTEntry / (float)(m_NumFilterLUTEntries - 1));
	  float filterVal = exp(-(angle * angle) * inv2Variance);
	  //note that gaussian is not weighted by 1.0 / (sigma* sqrt(2 * PI)) seen as weights
	  // weighted tap accumulation in filters is divided by sum of weights
	  m_FilterLUT[iLUTEntry] = filterVal;
	}
  }
  return m_FilterLUT;
}

//--------------------------------------------------------------------------------------
//Clear filter extents for the 6 cube map faces
//--------------------------------------------------------------------------------------
void ClearFilterExtents(BoundingBox *aFilterExtents)
{
  for (int32_t iCubeFaces = 0; iCubeFaces < 6; iCubeFaces++) {
	aFilterExtents[iCubeFaces].Clear();
  }
}


//--------------------------------------------------------------------------------------
//Define per-face bounding box filter extents
//
// These define conservative texel regions in each of the faces the filter can possibly 
// process.  When the pixels in the regions are actually processed, the dot product  
// between the tap vector and the center tap vector is used to determine the weight of 
// the tap and whether or not the tap is within the cone.
//
//--------------------------------------------------------------------------------------
void DetermineFilterExtents(float *a_CenterTapDir, int32_t a_SrcSize, int32_t a_BBoxSize,
  BoundingBox *a_FilterExtents)
{
  int32_t u, v;
  int32_t faceIdx;
  int32_t minU, minV, maxU, maxV;
  int32_t i;

  //neighboring face and bleed over amount, and width of BBOX for
  // left, right, top, and bottom edges of this face
  int32_t bleedOverAmount[4];
  int32_t bleedOverBBoxMin[4];
  int32_t bleedOverBBoxMax[4];

  int32_t neighborFace;
  int32_t neighborEdge;

  //get face idx, and u, v info from center tap dir
  VectToTexelCoord(a_CenterTapDir, a_SrcSize, &faceIdx, &u, &v);

  //define bbox size within face
  a_FilterExtents[faceIdx].Augment(u - a_BBoxSize, v - a_BBoxSize, 0);
  a_FilterExtents[faceIdx].Augment(u + a_BBoxSize, v + a_BBoxSize, 0);

  a_FilterExtents[faceIdx].ClampMin(0, 0, 0);
  a_FilterExtents[faceIdx].ClampMax(a_SrcSize - 1, a_SrcSize - 1, 0);

  //u and v extent in face corresponding to center tap
  minU = a_FilterExtents[faceIdx].m_minCoord[0];
  minV = a_FilterExtents[faceIdx].m_minCoord[1];
  maxU = a_FilterExtents[faceIdx].m_maxCoord[0];
  maxV = a_FilterExtents[faceIdx].m_maxCoord[1];

  //bleed over amounts for face across u=0 edge (left)    
  bleedOverAmount[0] = (a_BBoxSize - u);
  bleedOverBBoxMin[0] = minV;
  bleedOverBBoxMax[0] = maxV;

  //bleed over amounts for face across u=1 edge (right)    
  bleedOverAmount[1] = (u + a_BBoxSize) - (a_SrcSize - 1);
  bleedOverBBoxMin[1] = minV;
  bleedOverBBoxMax[1] = maxV;

  //bleed over to face across v=0 edge (up)
  bleedOverAmount[2] = (a_BBoxSize - v);
  bleedOverBBoxMin[2] = minU;
  bleedOverBBoxMax[2] = maxU;

  //bleed over to face across v=1 edge (down)
  bleedOverAmount[3] = (v + a_BBoxSize) - (a_SrcSize - 1);
  bleedOverBBoxMin[3] = minU;
  bleedOverBBoxMax[3] = maxU;

  //compute bleed over regions in neighboring faces
  for (i = 0; i < 4; i++) {
	if (bleedOverAmount[i] > 0) {
	  neighborFace = sg_CubeNgh[faceIdx][i].m_Face;
	  neighborEdge = sg_CubeNgh[faceIdx][i].m_Edge;

	  //For certain types of edge abutments, the bleedOverBBoxMin, and bleedOverBBoxMax need to 
	  //  be flipped: the cases are 
	  // if a left   edge mates with a left or bottom  edge on the neighbor
	  // if a top    edge mates with a top or right edge on the neighbor
	  // if a right  edge mates with a right or top edge on the neighbor
	  // if a bottom edge mates with a bottom or left  edge on the neighbor
	  //Seeing as the edges are enumerated as follows 
	  // left   =0 
	  // right  =1 
	  // top    =2 
	  // bottom =3            
	  // 
	  // so if the edge enums are the same, or the sum of the enums == 3, 
	  //  the bbox needs to be flipped
	  if ((i == neighborEdge) || ((i + neighborEdge) == 3)) {
		bleedOverBBoxMin[i] = (a_SrcSize - 1) - bleedOverBBoxMin[i];
		bleedOverBBoxMax[i] = (a_SrcSize - 1) - bleedOverBBoxMax[i];
	  }


	  //The way the bounding box is extended onto the neighboring face
	  // depends on which edge of neighboring face abuts with this one
	  switch (sg_CubeNgh[faceIdx][i].m_Edge) {
	  case CP_EDGE_LEFT:
		a_FilterExtents[neighborFace].Augment(0, bleedOverBBoxMin[i], 0);
		a_FilterExtents[neighborFace].Augment(bleedOverAmount[i], bleedOverBBoxMax[i], 0);
		break;
	  case CP_EDGE_RIGHT:
		a_FilterExtents[neighborFace].Augment((a_SrcSize - 1), bleedOverBBoxMin[i], 0);
		a_FilterExtents[neighborFace].Augment((a_SrcSize - 1) - bleedOverAmount[i], bleedOverBBoxMax[i], 0);
		break;
	  case CP_EDGE_TOP:
		a_FilterExtents[neighborFace].Augment(bleedOverBBoxMin[i], 0, 0);
		a_FilterExtents[neighborFace].Augment(bleedOverBBoxMax[i], bleedOverAmount[i], 0);
		break;
	  case CP_EDGE_BOTTOM:
		a_FilterExtents[neighborFace].Augment(bleedOverBBoxMin[i], (a_SrcSize - 1), 0);
		a_FilterExtents[neighborFace].Augment(bleedOverBBoxMax[i], (a_SrcSize - 1) - bleedOverAmount[i], 0);
		break;
	  }

	  //clamp filter extents in non-center tap faces to remain within surface
	  a_FilterExtents[neighborFace].ClampMin(0, 0, 0);
	  a_FilterExtents[neighborFace].ClampMax(a_SrcSize - 1, a_SrcSize - 1, 0);
	}

	//If the bleed over amount bleeds past the adjacent face onto the opposite face 
	// from the center tap face, then process the opposite face entirely for now. 
	//Note that the cases in which this happens, what usually happens is that 
	// more than one edge bleeds onto the opposite face, and the bounding box 
	// encompasses the entire cube map face.
	if (bleedOverAmount[i] > a_SrcSize) {
	  uint32_t oppositeFaceIdx;

	  //determine opposite face 
	  switch (faceIdx) {
	  case CP_FACE_X_POS:
		oppositeFaceIdx = CP_FACE_X_NEG;
		break;
	  case CP_FACE_X_NEG:
		oppositeFaceIdx = CP_FACE_X_POS;
		break;
	  case CP_FACE_Y_POS:
		oppositeFaceIdx = CP_FACE_Y_NEG;
		break;
	  case CP_FACE_Y_NEG:
		oppositeFaceIdx = CP_FACE_Y_POS;
		break;
	  case CP_FACE_Z_POS:
		oppositeFaceIdx = CP_FACE_Z_NEG;
		break;
	  case CP_FACE_Z_NEG:
		oppositeFaceIdx = CP_FACE_Z_POS;
		break;
	  default:
		break;
	  }

	  //just encompass entire face for now
	  a_FilterExtents[oppositeFaceIdx].Augment(0, 0, 0);
	  a_FilterExtents[oppositeFaceIdx].Augment((a_SrcSize - 1), (a_SrcSize - 1), 0);
	}
  }
  minV = minV;
}

// This function return the BaseFilterAngle require by cubemapgen to its FilterExtends
// It allow to optimize the texel to access base on the specular power.
static float GetBaseFilterAngle(float cosinePower)
{
  // We want to find the alpha such that:
  // cos(alpha)^cosinePower = epsilon
  // That's: acos(epsilon^(1/cosinePower))
  const float threshold = 0.000001f;  // Empirical threshold (Work perfectly, didn't check for a more big number, may get some performance and still god approximation)
  float Angle = 180.0f;
  if (Angle != 0.0f) {
	Angle = acosf(powf(threshold, 1.0f / cosinePower));
	Angle *= 180.0f / (float)PI_DOUBLE; // Convert to degree
	Angle *= 2.0f; // * 2.0f because cubemapgen divide by 2 later
  }
  return Angle;
}

float GetFilterAngle(int32_t a_FilterType, float specularPower = 2048.0f)
{
  // If we use SpecularPower, automatically calculate the a_BaseFilterAngle required, this will speed the process
  if (a_FilterType == CP_FILTER_TYPE_COSINE_POWER) {
	return GetBaseFilterAngle(specularPower);
  }
  return 0.0f;
}

//--------------------------------------------------------------------------------------
//ProcessFilterExtents 
//  Process bounding box in each cube face 
//
//--------------------------------------------------------------------------------------
void ProcessFilterExtents(
  const float *a_CenterTapDir,
  float a_DotProdThresh,
  const BoundingBox* a_FilterExtents,
  const NormalSurface *a_NormCubeMap,
  const ImageSurface *a_SrcCubeMap,
  float *a_DstVal,
  uint32_t a_FilterType,
  bool a_bUseSolidAngleWeighting,
  float a_SpecularPower
)
{
  int32_t iFaceIdx, u, v;
  int32_t faceWidth;
  int32_t k;

  //pointers used to walk across the image surface to accumulate taps
  const float *normCubeRowStartPtr;
  const float *srcCubeRowStartPtr;
  const float *texelVect;


  //accumulators are 64-bit floats in order to have the precision needed 
  // over a summation of a large number of pixels 
  double dstAccum[4];
  double weightAccum;

  float tapDotProd;   //dot product between center tap and current tap

  int32_t normCubePitch;
  int32_t srcCubePitch;
  int32_t normCubeRowWalk;
  int32_t srcCubeRowWalk;

  int32_t uStart, uEnd;
  int32_t vStart, vEnd;

  std::vector<float> m_FilterLUT = BuildAngleWeightLUT(a_FilterType, GetFilterAngle(a_FilterType));

  int32_t nSrcChannels = a_SrcCubeMap[0].m_NumChannels;

  //norm cube map and srcCubeMap have same face width
  faceWidth = a_NormCubeMap[0].m_Width;

  //amount to add to pointer to move to next scanline in images
  normCubePitch = faceWidth * a_NormCubeMap[0].m_NumChannels;
  srcCubePitch = faceWidth * a_SrcCubeMap[0].m_NumChannels;

  //dest accum
  for (k = 0; k < g_numChannels; k++) {
	dstAccum[k] = 0.0f;
  }

  weightAccum = 0.0f;

  // Add a more efficient path (without test and switch) for cosine power,
  // Basically just a copy past.
  if (a_FilterType != CP_FILTER_TYPE_COSINE_POWER) {
	//iterate over cubefaces
	for (iFaceIdx = 0; iFaceIdx < 6; iFaceIdx++) {
	  uStart = a_FilterExtents[iFaceIdx].m_minCoord[0];
	  vStart = a_FilterExtents[iFaceIdx].m_minCoord[1];
	  uEnd = a_FilterExtents[iFaceIdx].m_maxCoord[0];
	  vEnd = a_FilterExtents[iFaceIdx].m_maxCoord[1];

	  normCubeRowStartPtr = a_NormCubeMap[iFaceIdx].GetSurfacePtr() + (a_NormCubeMap[iFaceIdx].m_NumChannels * ((vStart * faceWidth) + uStart));
	  srcCubeRowStartPtr = a_SrcCubeMap[iFaceIdx].GetSurfacePtr() + (a_SrcCubeMap[iFaceIdx].m_NumChannels * ((vStart * faceWidth) + uStart));

	  //note that <= is used to ensure filter extents always encompass at least one pixel if bbox is non empty
	  for (v = vStart; v <= vEnd; v++) {
		normCubeRowWalk = 0;
		srcCubeRowWalk = 0;

		for (u = uStart; u <= uEnd; u++) {
		  //pointer to direction in cube map associated with texel
		  texelVect = (normCubeRowStartPtr + normCubeRowWalk);

		  //check dot product to see if texel is within cone
		  tapDotProd = VM_DOTPROD3(texelVect, a_CenterTapDir);

		  if (tapDotProd >= a_DotProdThresh) {
			float weight;

			//for now just weight all taps equally, but ideally
			// weight should be proportional to the solid angle of the tap
			if (a_bUseSolidAngleWeighting) {   //solid angle stored in 4th channel of normalizer/solid angle cube map
			  weight = *(texelVect + 3);
			} else {   //all taps equally weighted
			  weight = 1.0f;
			}

			switch (a_FilterType) {
			case CP_FILTER_TYPE_CONE:
			case CP_FILTER_TYPE_ANGULAR_GAUSSIAN:
			{
			  //weights are in same lookup table for both of these filter types
			  weight *= m_FilterLUT[(int32_t)(tapDotProd * (m_NumFilterLUTEntries - 1))];
			}
			break;
			case CP_FILTER_TYPE_COSINE:
			{
			  if (tapDotProd > 0.0f) {
				weight *= tapDotProd;
			  } else {
				weight = 0.0f;
			  }
			}
			break;
			case CP_FILTER_TYPE_DISC:
			default:
			  break;
			}

			//iterate over channels
			for (k = 0; k < nSrcChannels; k++)   //(aSrcCubeMap[iFaceIdx].m_NumChannels) //up to 4 channels 
			{
			  dstAccum[k] += weight * *(srcCubeRowStartPtr + srcCubeRowWalk);
			  srcCubeRowWalk++;
			}

			weightAccum += weight; //accumulate weight
		  } else {
			//step across source pixel
			srcCubeRowWalk += nSrcChannels;
		  }

		  normCubeRowWalk += a_NormCubeMap[iFaceIdx].m_NumChannels;
		}

		normCubeRowStartPtr += normCubePitch;
		srcCubeRowStartPtr += srcCubePitch;
	  }
	}
  } else // if (a_FilterType != CP_FILTER_TYPE_COSINE_POWER)
  {
	float IsPhongBRDF = 0.0f; //1.0f; // This value will be added to the specular power

	for (iFaceIdx = 0; iFaceIdx < 6; iFaceIdx++) {
	  uStart = a_FilterExtents[iFaceIdx].m_minCoord[0];
	  vStart = a_FilterExtents[iFaceIdx].m_minCoord[1];
	  uEnd = a_FilterExtents[iFaceIdx].m_maxCoord[0];
	  vEnd = a_FilterExtents[iFaceIdx].m_maxCoord[1];

	  normCubeRowStartPtr = a_NormCubeMap[iFaceIdx].GetSurfacePtr() + (a_NormCubeMap[iFaceIdx].m_NumChannels * ((vStart * faceWidth) + uStart));
	  srcCubeRowStartPtr = a_SrcCubeMap[iFaceIdx].GetSurfacePtr() + (a_SrcCubeMap[iFaceIdx].m_NumChannels * ((vStart * faceWidth) + uStart));

	  //note that <= is used to ensure filter extents always encompass at least one pixel if bbox is non empty
	  for (v = vStart; v <= vEnd; v++) {
		normCubeRowWalk = 0;
		srcCubeRowWalk = 0;

		for (u = uStart; u <= uEnd; u++) {
		  //pointer to direction in cube map associated with texel
		  texelVect = (normCubeRowStartPtr + normCubeRowWalk);

		  //check dot product to see if texel is within cone
		  tapDotProd = VM_DOTPROD3(texelVect, a_CenterTapDir);

		  if (tapDotProd >= a_DotProdThresh && tapDotProd > 0.0f) {
			float weight;

			//solid angle stored in 4th channel of normalizer/solid angle cube map
			weight = *(texelVect + 3);

			// Here we decide if we use a Phong/Blinn or a Phong/Blinn BRDF.
			// Phong/Blinn BRDF is just the Phong/Blinn model multiply by the cosine of the lambert law
			// so just adding one to specularpower do the trick.					   
			weight *= pow(tapDotProd, (a_SpecularPower + IsPhongBRDF));

			//iterate over channels
			for (k = 0; k < nSrcChannels; k++)   //(aSrcCubeMap[iFaceIdx].m_NumChannels) //up to 4 channels 
			{
			  dstAccum[k] += weight * *(srcCubeRowStartPtr + srcCubeRowWalk);
			  srcCubeRowWalk++;
			}

			weightAccum += weight; //accumulate weight
		  } else {
			//step across source pixel
			srcCubeRowWalk += nSrcChannels;
		  }

		  normCubeRowWalk += a_NormCubeMap[iFaceIdx].m_NumChannels;
		}

		normCubeRowStartPtr += normCubePitch;
		srcCubeRowStartPtr += srcCubePitch;
	  }
	}
  } // else // (a_FilterType != CP_FILTER_TYPE_COSINE_POWER)

	//divide through by weights if weight is non zero
  if (weightAccum != 0.0f) {
	for (k = 0; k < g_numChannels; k++) {
	  a_DstVal[k] = (float)(dstAccum[k] / weightAccum);
	}
  } else {   //otherwise sample nearest
	const float *texelPtr = GetCubeMapTexelPtr(a_CenterTapDir, a_SrcCubeMap);
	for (k = 0; k < g_numChannels; k++) {
	  a_DstVal[k] = texelPtr[k];
	}
  }
}

//--------------------------------------------------------------------------------------
// Fixup cube edges
//
// average texels on cube map faces across the edges
//--------------------------------------------------------------------------------------
void FixupCubeEdges(ImageSurface *a_CubeMap, int32_t a_FixupWidth)
{
  // If there is no fixup, or fixup width = 0, do nothing.
  // In case of Bent/Warp/Stretch Fixup and width of 1, we take the average of the texel color.
  if ((a_FixupWidth == 0) || (a_CubeMap[0].m_Width != 1)) {
	return;
  }

  //special case 1x1 cubemap, average face colors
  if (a_CubeMap[0].m_Width == 1) {
	//iterate over channels
	const int nChannels = a_CubeMap->m_NumChannels;
	for (int k = 0; k<nChannels; k++) {
	  float accum = 0.0f;
	  for (int iFace = 0; iFace<6; iFace++) {
		accum += *(a_CubeMap[iFace].GetSurfacePtr() + k);
	  }
	  accum /= 6.0f;
	  for (int iFace = 0; iFace<6; iFace++) {
		*(a_CubeMap[iFace].GetSurfacePtr() + k) = accum;
	  }
	}
	return;
  }
}

//--------------------------------------------------------------------------------------
//The key to the speed of these filtering routines is to quickly define a per-face 
//  bounding box of pixels which enclose all the taps in the filter kernel efficiently.  
//  Later these pixels are selectively processed based on their dot products to see if 
//  they reside within the filtering cone.
//
//This is done by computing the smallest per-texel angle to get a conservative estimate 
// of the number of texels needed to be covered in width and height order to filter the
// region.  the bounding box for the center taps face is defined first, and if the 
// filtereing region bleeds onto the other faces, bounding boxes for the other faces are 
// defined next
//--------------------------------------------------------------------------------------
void FilterCubeSurfaces(
  const ImageSurface *a_SrcCubeMap,
  ImageSurface *a_DstCubeMap,
  float a_FilterConeAngle,
  int32_t a_FilterType,
  bool a_bUseSolidAngle,
  int32_t a_FaceIdxStart,
  int32_t a_FaceIdxEnd,
  float a_SpecularPower
)
{
  //ImageSurface normCubeMap[6];     //
  BoundingBox   filterExtents[6];   //bounding box per face to specify region to process
									// note that pixels within these regions may be rejected 
									// based on the     
  int32_t iCubeFace, u, v;
  int32_t srcSize = a_SrcCubeMap[0].m_Width;
  int32_t dstSize = a_DstCubeMap[0].m_Width;

  float srcTexelAngle;
  float dotProdThresh;
  int32_t   filterSize;

  //angle about center tap to define filter cone
  float filterAngle;

  //min angle a src texel can cover (in degrees)
  srcTexelAngle = (180.0f / PI_FLOAT) * atan2f(1.0f, (float)srcSize);

  //filter angle is 1/2 the cone angle
  filterAngle = a_FilterConeAngle / 2.0f;

  //ensure filter angle is larger than a texel
  if (filterAngle < srcTexelAngle) {
	filterAngle = srcTexelAngle;
  }

  //ensure filter cone is always smaller than the hemisphere
  if (filterAngle > 90.0f) {
	filterAngle = 90.0f;
  }

  //the maximum number of texels in 1D the filter cone angle will cover
  //  used to determine bounding box size for filter extents
  filterSize = (int32_t)ceil(filterAngle / srcTexelAngle);

  //ensure conservative region always covers at least one texel
  if (filterSize < 1) {
	filterSize = 1;
  }

  //dotProdThresh threshold based on cone angle to determine whether or not taps 
  // reside within the cone angle
  dotProdThresh = cosf((PI_FLOAT / 180.0f) * filterAngle);

  NormalSurface normCubeMap[6];
  BuildNormalizerSolidAngleCubemap(a_SrcCubeMap[0].m_Width, normCubeMap);

  //process required faces
  for (iCubeFace = a_FaceIdxStart; iCubeFace <= a_FaceIdxEnd; iCubeFace++) {
	float *texelPtr = a_DstCubeMap[iCubeFace].GetSurfacePtr();

	//iterate over dst cube map face texel
	for (v = 0; v < dstSize; v++) {
	  for (u = 0; u < dstSize; u++) {
		float centerTapDir[3];  //direction of center tap
		//get center tap direction
		TexelCoordToVect(iCubeFace, (float)u, (float)v, dstSize, centerTapDir);

		//clear old per-face filter extents
		ClearFilterExtents(filterExtents);

		//define per-face filter extents
		DetermineFilterExtents(centerTapDir, srcSize, filterSize, filterExtents);

		//perform filtering of src faces using filter extents 
		ProcessFilterExtents(
		  centerTapDir,
		  dotProdThresh,
		  filterExtents,
		  normCubeMap,
		  a_SrcCubeMap,
		  texelPtr,
		  a_FilterType,
		  a_bUseSolidAngle,
		  a_SpecularPower
		);

		texelPtr += a_DstCubeMap[iCubeFace].m_NumChannels;
	  }
	}
  }

  FixupCubeEdges(a_DstCubeMap, 1);
}

// SH order use for approximation of irradiance cubemap is 5, mean 5*5 equals 25 coefficients
#define MAX_SH_ORDER 5
#define NUM_SH_COEFFICIENT (MAX_SH_ORDER * MAX_SH_ORDER)

// See Peter-Pike Sloan paper for these coefficients
static double SHBandFactor[NUM_SH_COEFFICIENT] = { 1.0,
2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0,
1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, // The 4 band will be zeroed
-1.0 / 24.0, -1.0 / 24.0, -1.0 / 24.0, -1.0 / 24.0, -1.0 / 24.0, -1.0 / 24.0, -1.0 / 24.0, -1.0 / 24.0, -1.0 / 24.0 };

void EvalSHBasis(const float* dir, double* res)
{
  // Can be optimize by precomputing constant.
  static const double SqrtPi = sqrt(PI_DOUBLE);

  double xx = dir[0];
  double yy = dir[1];
  double zz = dir[2];

  // x[i] == pow(x, i), etc.
  double x[MAX_SH_ORDER + 1], y[MAX_SH_ORDER + 1], z[MAX_SH_ORDER + 1];
  x[0] = y[0] = z[0] = 1.;
  for (int32_t i = 1; i < MAX_SH_ORDER + 1; ++i) {
	x[i] = xx * x[i - 1];
	y[i] = yy * y[i - 1];
	z[i] = zz * z[i - 1];
  }

  res[0] = (1 / (2.*SqrtPi));

  res[1] = -(sqrt(3 / PI_DOUBLE)*yy) / 2.;
  res[2] = (sqrt(3 / PI_DOUBLE)*zz) / 2.;
  res[3] = -(sqrt(3 / PI_DOUBLE)*xx) / 2.;

  res[4] = (sqrt(15 / PI_DOUBLE)*xx*yy) / 2.;
  res[5] = -(sqrt(15 / PI_DOUBLE)*yy*zz) / 2.;
  res[6] = (sqrt(5 / PI_DOUBLE)*(-1 + 3 * z[2])) / 4.;
  res[7] = -(sqrt(15 / PI_DOUBLE)*xx*zz) / 2.;
  res[8] = sqrt(15 / PI_DOUBLE)*(x[2] - y[2]) / 4.;

  res[9] = (sqrt(35 / (2.*PI_DOUBLE))*(-3 * x[2] * yy + y[3])) / 4.;
  res[10] = (sqrt(105 / PI_DOUBLE)*xx*yy*zz) / 2.;
  res[11] = -(sqrt(21 / (2.*PI_DOUBLE))*yy*(-1 + 5 * z[2])) / 4.;
  res[12] = (sqrt(7 / PI_DOUBLE)*zz*(-3 + 5 * z[2])) / 4.;
  res[13] = -(sqrt(21 / (2.*PI_DOUBLE))*xx*(-1 + 5 * z[2])) / 4.;
  res[14] = (sqrt(105 / PI_DOUBLE)*(x[2] - y[2])*zz) / 4.;
  res[15] = -(sqrt(35 / (2.*PI_DOUBLE))*(x[3] - 3 * xx*y[2])) / 4.;

  res[16] = (3 * sqrt(35 / PI_DOUBLE)*xx*yy*(x[2] - y[2])) / 4.;
  res[17] = (-3 * sqrt(35 / (2.*PI_DOUBLE))*(3 * x[2] * yy - y[3])*zz) / 4.;
  res[18] = (3 * sqrt(5 / PI_DOUBLE)*xx*yy*(-1 + 7 * z[2])) / 4.;
  res[19] = (-3 * sqrt(5 / (2.*PI_DOUBLE))*yy*zz*(-3 + 7 * z[2])) / 4.;
  res[20] = (3 * (3 - 30 * z[2] + 35 * z[4])) / (16.*SqrtPi);
  res[21] = (-3 * sqrt(5 / (2.*PI_DOUBLE))*xx*zz*(-3 + 7 * z[2])) / 4.;
  res[22] = (3 * sqrt(5 / PI_DOUBLE)*(x[2] - y[2])*(-1 + 7 * z[2])) / 8.;
  res[23] = (-3 * sqrt(35 / (2.*PI_DOUBLE))*(x[3] - 3 * xx*y[2])*zz) / 4.;
  res[24] = (3 * sqrt(35 / PI_DOUBLE)*(x[4] - 6 * x[2] * y[2] + y[4])) / 16.;
}

void SHFilterCubeMap(
  const ImageSurface* m_InputSurface,
  ImageSurface* m_OutputSurface,
  bool a_bUseSolidAngleWeighting
)
{
  const ImageSurface* SrcCubeImage = m_InputSurface;
  ImageSurface* DstCubeImage = m_OutputSurface;

  int32_t SrcSize = SrcCubeImage->m_Width;
  int32_t DstSize = DstCubeImage->m_Width;

  //pointers used to walk across the image surface
  const float *normCubeRowStartPtr;
  const float *srcCubeRowStartPtr;
  float *dstCubeRowStartPtr;

  const int32_t SrcCubeMapNumChannels = SrcCubeImage[0].m_NumChannels;
  const int32_t DstCubeMapNumChannels = DstCubeImage[0].m_NumChannels;

  //First step - Generate SH coefficient for the diffuse convolution

  //Regenerate normalization cubemap for source cubemap
  //Normalized vectors per cubeface and per-texel solid angle 
  NormalSurface normCubeMap[6];
  BuildNormalizerSolidAngleCubemap(m_InputSurface->m_Width, normCubeMap);

  const int32_t NormCubeMapNumChannels = normCubeMap[0].m_NumChannels; // This need to be init here after the generation of normCubeMap

  //This is a custom implementation of D3DXSHProjectCubeMap to avoid to deal with LPDIRECT3DSURFACE9 pointer
  //Use Sh order 2 for a total of 9 coefficient as describe in http://www.cs.berkeley.edu/~ravir/papers/envmap/
  //accumulators are 64-bit floats in order to have the precision needed 
  //over a summation of a large number of pixels 
  double SHr[NUM_SH_COEFFICIENT];
  double SHg[NUM_SH_COEFFICIENT];
  double SHb[NUM_SH_COEFFICIENT];
  double SHdir[NUM_SH_COEFFICIENT];

  memset(SHr, 0, NUM_SH_COEFFICIENT * sizeof(double));
  memset(SHg, 0, NUM_SH_COEFFICIENT * sizeof(double));
  memset(SHb, 0, NUM_SH_COEFFICIENT * sizeof(double));
  memset(SHdir, 0, NUM_SH_COEFFICIENT * sizeof(double));

  double weightAccum = 0.0;
  double weight = 0.0;

  for (int32_t iFaceIdx = 0; iFaceIdx < 6; iFaceIdx++) {
	for (int32_t y = 0; y < SrcSize; y++) {
	  normCubeRowStartPtr = normCubeMap[iFaceIdx].GetSurfaceTexelPtr(0, y);
	  srcCubeRowStartPtr = SrcCubeImage[iFaceIdx].GetSurfaceTexelPtr(0, y);

	  for (int32_t x = 0; x < SrcSize; x++) {
		//pointer to direction and solid angle in cube map associated with texel
		const float* texelVect = &normCubeRowStartPtr[NormCubeMapNumChannels * x];

		if (a_bUseSolidAngleWeighting) {   //solid angle stored in 4th channel of normalizer/solid angle cube map
		  weight = *(texelVect + 3);
		} else {   //all taps equally weighted
		  weight = 1.0;
		}

		EvalSHBasis(texelVect, SHdir);

		// Convert to float64
		double R = srcCubeRowStartPtr[(SrcCubeMapNumChannels * x) + 0];
		double G = srcCubeRowStartPtr[(SrcCubeMapNumChannels * x) + 1];
		double B = srcCubeRowStartPtr[(SrcCubeMapNumChannels * x) + 2];

		for (int32_t i = 0; i < NUM_SH_COEFFICIENT; i++) {
		  SHr[i] += R * SHdir[i] * weight;
		  SHg[i] += G * SHdir[i] * weight;
		  SHb[i] += B * SHdir[i] * weight;
		}

		weightAccum += weight;
	  }
	}
  }

  //Normalization - The sum of solid angle should be equal to the solid angle of the sphere (4 PI), so
  // normalize in order our weightAccum exactly match 4 PI.
  for (int32_t i = 0; i < NUM_SH_COEFFICIENT; ++i) {
	SHr[i] *= 4.0 * PI_DOUBLE / weightAccum;
	SHg[i] *= 4.0 * PI_DOUBLE / weightAccum;
	SHb[i] *= 4.0 * PI_DOUBLE / weightAccum;
  }

  //Second step - Generate cubemap from SH coefficient

  //Regenerate normalization cubemap for the destination cubemap
  //Normalized vectors per cubeface and per-texel solid angle 
  BuildNormalizerSolidAngleCubemap(DstCubeImage->m_Width, normCubeMap);

  for (int32_t iFaceIdx = 0; iFaceIdx < 6; iFaceIdx++) {
	for (int32_t y = 0; y < DstSize; y++) {
	  normCubeRowStartPtr = normCubeMap[iFaceIdx].GetSurfaceTexelPtr(0, y);
	  dstCubeRowStartPtr = DstCubeImage[iFaceIdx].GetSurfaceTexelPtr(0, y);

	  for (int32_t x = 0; x < DstSize; x++) {
		//pointer to direction and solid angle in cube map associated with texel
		const float* texelVect = &normCubeRowStartPtr[NormCubeMapNumChannels * x];

		EvalSHBasis(texelVect, SHdir);

		// get color value
		float R = 0.0f, G = 0.0f, B = 0.0f;

		for (int32_t i = 0; i < NUM_SH_COEFFICIENT; ++i) {
		  R += (float)(SHr[i] * SHdir[i] * SHBandFactor[i]);
		  G += (float)(SHg[i] * SHdir[i] * SHBandFactor[i]);
		  B += (float)(SHb[i] * SHdir[i] * SHBandFactor[i]);
		}

		dstCubeRowStartPtr[(DstCubeMapNumChannels * x) + 0] = R;
		dstCubeRowStartPtr[(DstCubeMapNumChannels * x) + 1] = G;
		dstCubeRowStartPtr[(DstCubeMapNumChannels * x) + 2] = B;
		if (DstCubeMapNumChannels > 3) {
		  dstCubeRowStartPtr[(DstCubeMapNumChannels * x) + 3] = 1.0f;
		}
	  }
	}
  }
}

}// namespace CubeMapGen
