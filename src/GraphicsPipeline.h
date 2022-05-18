#pragma   once  
#include "la.h"
#include "Camera.h"
#include "CImg.h"
#include "objparser.h"
#include <filesystem>
#include "platforms/framebuffer.h"
#include <mutex>
#include "isocontour.h"




enum MeshBuildConvention {
	LIKE_GL_TRIANGLE,
	LIKE_GL_TRIANGLE_STRIP
};

struct Texture {
	mVec3f* TextureArr;
	int  TW;
	int  TH;
};

typedef struct vertex {
	mVec3f Coordinate;
	//attributes
	mVec3f normal;
	mVec3f shading;
	mVec2<float> st;
}vtx;


//__declspec(align(1))
struct fragment {
	float depth;
	mVec2<int>XY;
	mVec3f RGBv;
	fragment() = default;
	__forceinline fragment(float depth,mVec2<int>XY, mVec3f RGBv):XY(XY),RGBv(RGBv),depth(depth) {

	}
	[[clang::always_inline]]fragment(fragment&& other)
		noexcept : XY(std::move(other.XY)), RGBv(std::move(other.RGBv)), depth(other.depth)
	{
	}

	fragment& operator=(const fragment& other) = default;
};



struct ShadingMaterial {
	mVec3f Ka;
	mVec3f Kd;
	mVec3f Ks;
	float f;
};



class TriangleWithAttributes {
public:
	vtx vtxa, vtxb, vtxc;
	inline TriangleWithAttributes(vtx a, vtx b, vtx c) :vtxa(a), vtxb(b), vtxc(c) {
	}
	inline TriangleWithAttributes(mVec3f ai, mVec3f bi, mVec3f ci) {
		vtxa.Coordinate = ai;
		vtxb.Coordinate = bi;
		vtxc.Coordinate = ci;
	}
	TriangleWithAttributes() = default;
	inline mVec3f Normal() {
		mVec3f L1(vtxa.Coordinate, vtxc.Coordinate);
		mVec3f L2(vtxb.Coordinate, vtxc.Coordinate);
		mVec3f out = L2.cross_product(L1);
		out.normalize();
		/*	mVec3f out2 = L1.cross_product(L2);
			if (dir*out > 0)return out2;*/
		return  out;

	}
	inline Triangle GetTriangleVertexes() const{
		return  Triangle(this->vtxa.Coordinate, this->vtxb.Coordinate, this->vtxc.Coordinate);

	}
	
};


class RenderableObject {
private:
	ShadingMaterial shma;	// material
	bool fixObject = false;
	Matrix4 ModleMatrix;
public:
	float texture_scaling = 1;
	unsigned int texturechannelID = 0;//{0-no texture, 1 :procedual strip, 2,3,...:real textures}
	
	std::vector<mVec3i> TrianglesIdx;
	std::vector<mVec3f> Vertices;
	std::vector<mVec3f> VtxNormals;
	std::vector<mVec2<float>>VtxTexUV;
	
	Matrix4 NormalMatrix ;
	inline void fixit() {
		fixObject = true;
	}
	inline void setModelMatrix(const Matrix4& mMat) {
		if (fixObject)return;
		this->ModleMatrix= mMat;
	}
	inline Matrix4 getModelMatrix() const{
		return this->ModleMatrix ;
	}
	inline void updateMaterial(const ShadingMaterial& newmat) {
		this->shma = newmat;
	}
	inline void updateMaterial(mVec3f a, mVec3f d, mVec3f s,float f) {
	if (!fixObject)
		this->shma = ShadingMaterial({ a, d, s,f });
	}
	inline ShadingMaterial Material() const{
		return this->shma;
	}
	inline void movePos(const mVec3f& STEP) {
		if (fixObject)return;
		this->ModleMatrix = translateMatrix(STEP)* this->ModleMatrix;
	}

	RenderableObject() {
		ModleMatrix = eye(1);
		this->shma = { mVec3f(1,1,1) * 0.002, mVec3f(1, 1, 1) * 600, mVec3f(1, 1, 1) * 20,8 };//  0.0005      8
	};

	inline RenderableObject(std::vector<mVec3i> TrianglesIdx, std::vector<mVec3f>Vertexes, std::vector<mVec3f>VertexesNormal, std::vector<mVec2<float>>VertexesT) {
		this->TrianglesIdx = TrianglesIdx;
		this->Vertices =Vertexes;
		this->VtxNormals = VertexesNormal;
		this->VtxTexUV = VertexesT;
		ModleMatrix = eye(1);
		this->shma = { mVec3f(1,1,1) * 0.002, mVec3f(1, 1, 1)* 600, mVec3f(1, 1, 1)*20,8};//  0.0005      8
	}

	inline void InplaceRotation(int axis, float degree) {
		if (fixObject)return;
		//auto MatInv = inv(ModleMatrix);
		auto rotate = rotateMatrix(axis,degree);
		this->ModleMatrix = this->ModleMatrix * rotate;
	}
	
	void BuildLikGLBegin(const std::vector<mVec3f>& vtxs, const std::vector<mVec3f>& normals, enum MeshBuildConvention Conv);


};


class Cube : public RenderableObject{
public:
	 inline Cube():RenderableObject() {
		float side = 1.0f;
		float side2 = side / 2.0f;

		float v[24 * 3] = {
			// Front
		   -side2, -side2, side2,
			side2, -side2, side2,
			side2,  side2, side2,
		   -side2,  side2, side2,
		   // Right
			side2, -side2, side2,
			side2, -side2, -side2,
			side2,  side2, -side2,
			side2,  side2, side2,
			// Back
			-side2, -side2, -side2,
			-side2,  side2, -side2,
			 side2,  side2, -side2,
			 side2, -side2, -side2,
			 // Left
			 -side2, -side2, side2,
			 -side2,  side2, side2,
			 -side2,  side2, -side2,
			 -side2, -side2, -side2,
			 // Bottom
			 -side2, -side2, side2,
			 -side2, -side2, -side2,
			  side2, -side2, -side2,
			  side2, -side2, side2,
			  // Top
			  -side2,  side2, side2,
			   side2,  side2, side2,
			   side2,  side2, -side2,
			  -side2,  side2, -side2
		};

		float n[24 * 3] = {
			// Front
			0.0f, 0.0f, 1.0f,
			0.0f, 0.0f, 1.0f,
			0.0f, 0.0f, 1.0f,
			0.0f, 0.0f, 1.0f,
			// Right
			1.0f, 0.0f, 0.0f,
			1.0f, 0.0f, 0.0f,
			1.0f, 0.0f, 0.0f,
			1.0f, 0.0f, 0.0f,
			// Back
			0.0f, 0.0f, -1.0f,
			0.0f, 0.0f, -1.0f,
			0.0f, 0.0f, -1.0f,
			0.0f, 0.0f, -1.0f,
			// Left
			-1.0f, 0.0f, 0.0f,
			-1.0f, 0.0f, 0.0f,
			-1.0f, 0.0f, 0.0f,
			-1.0f, 0.0f, 0.0f,
			// Bottom
			0.0f, -1.0f, 0.0f,
			0.0f, -1.0f, 0.0f,
			0.0f, -1.0f, 0.0f,
			0.0f, -1.0f, 0.0f,
			// Top
			0.0f, 1.0f, 0.0f,
			0.0f, 1.0f, 0.0f,
			0.0f, 1.0f, 0.0f,
			0.0f, 1.0f, 0.0f
		};

		float tex[24 * 2] = {
			// Front
			0.0f, 0.0f,
			1.0f, 0.0f,
			1.0f, 1.0f,
			0.0f, 1.0f,
			// Right
			0.0f, 0.0f,
			1.0f, 0.0f,
			1.0f, 1.0f,
			0.0f, 1.0f,
			// Back
			0.0f, 0.0f,
			1.0f, 0.0f,
			1.0f, 1.0f,
			0.0f, 1.0f,
			// Left
			0.0f, 0.0f,
			1.0f, 0.0f,
			1.0f, 1.0f,
			0.0f, 1.0f,
			// Bottom
			0.0f, 0.0f,
			1.0f, 0.0f,
			1.0f, 1.0f,
			0.0f, 1.0f,
			// Top
			0.0f, 0.0f,
			1.0f, 0.0f,
			1.0f, 1.0f,
			0.0f, 1.0f
		};

		int el[6*2*3] = {
			0,1,2,0,2,3,
			4,5,6,4,6,7,
			8,9,10,8,10,11,
			12,13,14,12,14,15,
			16,17,18,16,18,19,
			20,21,22,20,22,23
		};
		this->Vertices.reserve(36*3);
		this->VtxTexUV.reserve(36 * 3);
		this->VtxNormals.reserve(36 * 3);
		this->TrianglesIdx.reserve(36 * 3);
		for (int i=0;i<12;i++)
		{
			int faceid [3]={el[3 * i], el[3 * i + 1], el[3 * i + 2] };
			this->TrianglesIdx.emplace_back(faceid[0], faceid[1], faceid[2]);
		}
		for (int i = 0; i < 24; i++)
		{
			this->Vertices.emplace_back(v[i*3+0], v[i * 3 + 1], v[i * 3 + 2]);
			this->VtxNormals.emplace_back(n[i * 3 + 0], n[i * 3 + 1], n[i * 3 + 2]);
			this->VtxTexUV.emplace_back(tex[i*2], tex[i*2+1]);
		}

	}
};




class GraphicsPipeline {

public:
	bool Clip=true;
	bool UsePhongShading;//false for GouraudShading  true for phongshading 
	Camera  _Camera;
	int window_height, window_width;

	GraphicsPipeline() = default;
	GraphicsPipeline(int w, int h);
	void GraphicsPipeline::clearPipeline();

	void LoadTexture(std::string path);
	void LoadTexture(std::filesystem::path tex_path);
	STRONG_INLINE  mVec3f LookupTexel(mVec2<float> st, const int textureChannel) const  {
			auto linear = [](float alpha, mVec3f a, mVec3f b) {
				return a * alpha + b * (1 - alpha);
			};
			const auto& tex = mTextures[textureChannel];
			auto TW = tex.TW;
			auto TH = tex.TH;
			auto TextureArr = tex.TextureArr;

			float s = st.x;
			float t = st.y;
			s = s > 1 ? 1 : (s < 0 ? 0 : s);
			t = t > 1 ? 1 : (t < 0 ? 0 : t);
			s *= (TW - 1);
			t *= (TH - 1);
			int s1 = int(floorf(s));
			int s2 = int(ceilf(s));
			int t1 = int(floorf(t));
			int t2 = int(ceilf(t));
			mVec3f C1 = TextureArr[t1 * TW + s1];
			mVec3f C2 = TextureArr[t1 * TW + s2];
			mVec3f C3 = TextureArr[t2 * TW + s1];
			mVec3f C4 = TextureArr[t2 * TW + s2];
			mVec3f C12 = linear(s - s1, C1, C2);
			mVec3f C34 = linear(s - s1, C3, C4);
			mVec3f C = linear(t - t1, C12, C34);
			return C;
	}



	void Render(const std::vector<RenderableObject>& robjs1, const std::vector<RenderableObject>& robjs2, framebuffer_t* Fb);
	inline int getActiveTextureNumber() {
		return mTextures.size() ;
	}

private:
	std::vector<Texture>mTextures;

    STRONG_INLINE  void VertexesProcess(const RenderableObject& yu);
	STRONG_INLINE  void Rasterization(const RenderableObject& robj, framebuffer_t* Fb) const;
	STRONG_INLINE mVec4f FragmentShading(const Triangle& TriInEyeSpace, const mVec3i& TriIdx, const mVec2<float>& uv, mVec3f& recip_w, const RenderableObject& robj)const;

	float* depthbuffer = NULL;


	

	std::vector<mVec3f> vec3Buffers[MAXIMUM_GPIPELINE_BUFFERS];
	std::vector<mVec4f> vec4Buffers[MAXIMUM_GPIPELINE_BUFFERS];
	std::vector<mVec3i> vec3iBuffer;
	std::vector<mVec2f> vec2Buffer;

	std::vector<TriangleWithAttributes>TargetRenderTriangles;

	void clip_with_plane(const Plane<float>& c_plane);

};
















