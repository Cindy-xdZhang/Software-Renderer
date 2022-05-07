#pragma   once  
#include "la.h"
#include "Camera.h"
#include "CImg.h"
#include <filesystem>
#include "platforms/framebuffer.h"
#include <mutex>

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
	mVec3f EyeSpaceCoordinate;
	mVec2<float> st;
	//vertex operator+(vertex right) {
	//	return { Coordinate + right.Coordinate,normal + right.normal, shading + right.shading, EyeSpaceCoordinate + right.EyeSpaceCoordinate };
	//}
	//vertex operator*(float right) {
	//	return { Coordinate * right,normal * right, shading * right, EyeSpaceCoordinate * right };
	//}
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
	vtx a, b, c;
	inline TriangleWithAttributes(vtx a, vtx b, vtx c) :a(a), b(b), c(c) {
	}
	inline TriangleWithAttributes(mVec3f ai, mVec3f bi, mVec3f ci) {
		a.Coordinate = ai;
		b.Coordinate = bi;
		c.Coordinate = ci;
	}
	TriangleWithAttributes() = default;
	inline mVec3f Normal() {
		mVec3f L1(a.Coordinate, c.Coordinate);
		mVec3f L2(b.Coordinate, c.Coordinate);
		mVec3f out = L2.cross_product(L1);
		out.normalize();
		/*	mVec3f out2 = L1.cross_product(L2);
			if (dir*out > 0)return out2;*/
		return  out;

	}
	inline Triangle GetTriangleVertexes() const{
		return  Triangle(this->a.Coordinate, this->b.Coordinate, this->c.Coordinate);

	}
	
};


class RenderableObject {
private:
	ShadingMaterial shma;	// material

public:
	std::vector<mVec3i> TrianglesIdx;
	std::vector<mVec3f> Vertexes;
	std::vector<mVec3f> VertexesNormal;
	std::vector<mVec2<float>>VtxTexUV;
	Matrix4 ModleMatrix ;
	Matrix4 NormalMatrix ;
	inline void updateMaterial(mVec3f a, mVec3f d, mVec3f s,float f) {
		this->shma = ShadingMaterial({ a, d, s,f });
	}
	inline ShadingMaterial Material() const{
		return this->shma;
	}
	RenderableObject() = default;
	inline RenderableObject(std::vector<mVec3i> TrianglesIdx, std::vector<mVec3f>Vertexes, std::vector<mVec3f>VertexesNormal, std::vector<mVec2<float>>VertexesT) {
		this->TrianglesIdx = TrianglesIdx;
		this->Vertexes =Vertexes;
		this->VertexesNormal = VertexesNormal;
		this->VtxTexUV = VertexesT;
		ModleMatrix = translateMatrix(mVec3f(0, -12, -10)) ;
		this->shma = { mVec3f(1,1,1) * 0.002, mVec3f(1, 1, 1)* 600, mVec3f(1, 1, 1)*20,8};//  0.0005      8
	}

	inline void InplaceRotation(int axis, float degree) {
		//auto MatInv = inv(ModleMatrix);
		auto rotate = rotateMatrix(axis,degree);
		this->ModleMatrix = this->ModleMatrix * rotate;
	}


};




class GraphicsPipeline {

public:
	bool Clip;
	bool UsePhongShading;//false for GouraudShading  true for phongshading 
	Camera  _Camera;
	int window_height, window_width;
	
	int TextureMode=2;

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

	void Render(const RenderableObject& yu, framebuffer_t* Fb) ;

private:

	float* depthbuffer = NULL;
	STRONG_INLINE  void VertexesProcess(const RenderableObject& yu);
	STRONG_INLINE mVec3f ADSshading(mVec3f point, mVec3f normal, ShadingMaterial Mp) const;	
	STRONG_INLINE void Rasterization(const ShadingMaterial& myu, framebuffer_t* Fb) const;
	//inline void FragmentProcess(framebuffer_t* Fb);

	mVec3f LightSource;

	std::vector<TriangleWithAttributes> TargetRenderTriangles;
	//std::vector<std::vector<fragment>> FragmentsAfterRasterization;
	std::vector<Texture>mTextures;
};














