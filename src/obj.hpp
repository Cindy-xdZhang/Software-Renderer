#pragma   once  
#include "assert.h"
#include "la.h"
#include "matrix.h"
#include "Base.h"
#include <mutex>
#include <thread>
#include "timer.hpp"
#include "CImg.h"
#include "platforms/framebuffer.h"

#define LI 100
using namespace cimg_library;

static mVec3 InitEyePos = mVec3(0, 0, 200);
static mVec3 InitGazeDirection = mVec3(0, 0, -1);
static mVec3 InitTopDirection = mVec3(0, 1, 0);

struct Texture {
	mVec3* TextureArr;
	int TW;
	int TH;
};






typedef struct vertex {
	mVec3 Coordinate;

	//attributes
	mVec3 normal;
	mVec3 shading;
	mVec3 EyeSpaceCoordinate;
	mVec2<float> st;
	//vertex operator+(vertex right) {
	//	return { Coordinate + right.Coordinate,normal + right.normal, shading + right.shading, EyeSpaceCoordinate + right.EyeSpaceCoordinate };
	//}
	//vertex operator*(float right) {
	//	return { Coordinate * right,normal * right, shading * right, EyeSpaceCoordinate * right };
	//}
}vtx;

typedef struct fragment {
	mVec2<int>XY;
	mVec3 RGB;
	double depth;
}fgt;

struct ShadingMaterial {
	mVec3 Ka;
	mVec3 Kd;
	mVec3 Ks;
	float f;
};



class TriangleWithAttributes {
public:
	vtx a, b, c;
	inline TriangleWithAttributes(vtx a, vtx b, vtx c) :a(a), b(b), c(c) {
	}
	inline TriangleWithAttributes(mVec3 ai, mVec3 bi, mVec3 ci) {
		a.Coordinate = ai;
		b.Coordinate = bi;
		c.Coordinate = ci;
	}
	TriangleWithAttributes() = default;
	inline mVec3 Normal() {
		mVec3 L1(a.Coordinate, c.Coordinate);
		mVec3 L2(b.Coordinate, c.Coordinate);
		mVec3 out = L2.cross_product(L1);
		out.normalize();
		/*	mVec3 out2 = L1.cross_product(L2);
			if (dir*out > 0)return out2;*/
		return  out;

	}
	inline Triangle GetTriangleVertexes() {
		return  Triangle(this->a.Coordinate, this->b.Coordinate, this->c.Coordinate);

	}
	
};


#ifdef CLIP
typedef struct vertex4 {
	mVec4 Coordinate;

	//attributes
	mVec3 normal;
	mVec3 shading;
	mVec3 EyeSpaceCoordinate;
	//mVec2<float> st;
	vertex4 operator+(vertex4 right) {
		return { Coordinate + right.Coordinate,normal + right.normal, shading + right.shading, EyeSpaceCoordinate + right.EyeSpaceCoordinate };
	}
	vertex4 operator*(float right) {
		return { Coordinate * right,normal * right, shading * right, EyeSpaceCoordinate * right };
	}
};
class Vertex4TriangleWithAttributes {
public:
	vertex4  a, b, c;
	/*inline Triangle GetTriangleVertexes() {
		return  Triangle(this->a.Coordinate.tomVec3(), this->b.Coordinate.tomVec3(), this->c.Coordinate.tomVec3());

	}*/
	Vertex4TriangleWithAttributes(vertex4 a, vertex4 b, vertex4 c) :a(a), b(b), c(c) {

	}
	TriangleWithAttributes toTriangleWithAttributes() {
		vertex  a1 = { a.Coordinate.tomVec3(),a.normal,a.shading,a.EyeSpaceCoordinate };
		vertex	b1 = { b.Coordinate.tomVec3(),b.normal,b.shading,b.EyeSpaceCoordinate };
		vertex	c1 = { c.Coordinate.tomVec3(),a.normal,c.shading,c.EyeSpaceCoordinate };
		return TriangleWithAttributes(a1, b1, c1);
	}
};
//clip_plane为裁剪平面的自定义结构体，vert_list存储了待裁剪凸多边形的所有顶点
//num_vert为顶点个数，in_list为需要保留下来的裁剪平面内侧顶点的列表
void clip_with_plane(Plane c_plane, std::vector<Vertex4TriangleWithAttributes> in_list, std::vector<Vertex4TriangleWithAttributes>& out_list) {
	auto cal_project_distance = [](Plane ic_plane, mVec3 ptx) ->float {
		return (ptx - ic_plane.px)*ic_plane.n;
	};
	auto cal_insertRatio = [](Plane ic_plane, mVec3 ptx1, mVec3 ptx2) ->float {
		float a= ic_plane.n*ptx1 - ic_plane.n*ic_plane.px;
		float b = ic_plane.n*(ptx1 - ptx2);

		return a / b;
	};
	for (auto triangle : in_list) {
		int i;
		vertex4 vert_list[3] = { triangle.a,  triangle.b, triangle.c };

		int previous_index, current_index;
	    vertex4 remainPoints[4];
		int remaincount = 0;
		for (i = 0; i < 3; i++)
		{
			//从最后一个点开始，遍历所有边
			current_index = i;
			previous_index = (i - 1 + 3) % 3;
			mVec3 pre_vertex = vert_list[previous_index].Coordinate.tomVec3(); //边的起始点
			mVec3 cur_vertex = vert_list[current_index].Coordinate.tomVec3();  //边的终止点

			float d1 = cal_project_distance(c_plane, pre_vertex);
			float d2 = cal_project_distance(c_plane, cur_vertex);

			//如果该边与裁剪平面有交点，则计算交点并存入in_list
			if (d1 * d2 < 0)
			{
				float t = cal_insertRatio(c_plane, pre_vertex, cur_vertex);//求出t值
				vertex4 intersectionPostion =  vert_list[current_index]+  vert_list[previous_index] * (1-t);
				remainPoints[remaincount++] = intersectionPostion;

			}
			//如果终止点在内侧，直接存入in_list
			if (d2 < 0)
			{
				remainPoints[remaincount++] = vert_list[current_index];
			}

			

		}
		if (remaincount == 3) {
			out_list.emplace_back(Vertex4TriangleWithAttributes(remainPoints[0], remainPoints[1], remainPoints[2]));

		}
		else	if (remaincount == 4) {
		out_list.emplace_back(Vertex4TriangleWithAttributes(remainPoints[0], remainPoints[1], remainPoints[2]));
		out_list.emplace_back(Vertex4TriangleWithAttributes(remainPoints[1],remainPoints[2], remainPoints[3]));
		}
		else if (remaincount != 0) {
			std::cout << " Unknown Error!";
		}
	
	
	
	
	
	
	
	}


	return;
}


#endif


class RenderObject {
protected:
	ShadingMaterial shma;	// material

public:
	std::vector<mVec3i> TrianglesIdx;
	std::vector<mVec3> Vertexes;
	std::vector<mVec3> VertexesNormal;
	std::vector<mVec2<float>>VertexesTex;
	Matrix4 ModleMatrix ;
	Matrix4 NormalMatrix ;
	inline void updateMaterial(mVec3 a, mVec3 d, mVec3 s,float f) {
		this->shma = ShadingMaterial({ a, d, s,f });
	}
	inline ShadingMaterial Material() const{
		return this->shma;
	}
	RenderObject() = default;
	RenderObject(std::vector<mVec3i> TrianglesIdx, std::vector<mVec3>Vertexes, std::vector<mVec3>VertexesNormal, std::vector<mVec2<float>>VertexesT) {
		this->TrianglesIdx = TrianglesIdx;
		this->Vertexes =Vertexes;
		this->VertexesNormal = VertexesNormal;
		this->VertexesTex = VertexesT;
		ModleMatrix =  translateMatrix(mVec3(-125, -125, -125));
		this->shma = { mVec3(1,1,1) * 0.002, mVec3(1, 1, 1)* 600, mVec3(1, 1, 1)*20,8};//  0.0005      8
	}
};

class GraphicsPipeline {
private:
	
	MyTimer timer;

	float* depthbuffer = NULL;
	inline void VertexesProcess(const RenderObject& yu );
	inline void Rasterization(const ShadingMaterial& myu );
	inline void FragmentProcess(framebuffer_t* Fb);
	mVec3 LightSource;
	std::vector<TriangleWithAttributes> TargetRenderTriangles;
	std::vector<fragment> FragmentsAfterRasterization;

	std::vector<Texture>mTextures;
public:
	bool Clip;
	int TextureChanel;
	bool UsePhongShading;//false for GouraudShading  ture for phongshading 
	mVec3 InitLightPos = mVec3(0, 20, 0);
	Camera  _Camera;
	int h, w;
	GraphicsPipeline() = default;
	GraphicsPipeline(int h,int w): h(h),w(w){
		//mTextures.emplace_back(LoadTexture());
		TextureChanel = mTextures.size()-1;

		Clip = false;
		UsePhongShading = false;
		_Camera = Camera(InitEyePos, InitGazeDirection, InitTopDirection);
		_Camera.SetFrustm(1,-1,-1,1,90.0f, InitEyePos.z+200);
		
	};
	
	inline void Render(const RenderObject& yu, framebuffer_t* Fb);
	inline void GraphicsPipeline::clearPipeline();

	inline  Texture  LoadTexture() {
		CImg<unsigned  char> TextureImg("C:\\Users\\8\\source\\repos\\assignment3\\assignment3\\brick1.bmp");
		Texture  tX;

		tX.TW = TextureImg.width();
		tX.TH = TextureImg.height();
		tX.TextureArr = (mVec3*)malloc(tX.TW * tX.TH * sizeof(mVec3));
		for (int r = 0; r < tX.TH; r++)
			for (int c = 0; c < tX.TW; c++) {
				tX.TextureArr[r * tX.TW + c].x = float(TextureImg(c, r, 0)) / 255.0f;
				tX.TextureArr[r * tX.TW + c].y = float(TextureImg(c, r, 1)) / 255.0f;
				tX.TextureArr[r * tX.TW + c].z = float(TextureImg(c, r, 2)) / 255.0f;
			}
		return tX;
	}

	inline mVec3 LookupTexel(mVec2<float> st) {
		auto linear = [](float alpha, mVec3 a, mVec3 b) {
			return a * alpha + b * (1 - alpha);
		};
		const auto& tex = mTextures[TextureChanel];
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
		mVec3 C1 = TextureArr[t1 * TW + s1];
		mVec3 C2 = TextureArr[t1 * TW + s2];
		mVec3 C3 = TextureArr[t2 * TW + s1];
		mVec3 C4 = TextureArr[t2 * TW + s2];
		mVec3 C12 = linear(s - s1, C1, C2);
		mVec3 C34 = linear(s - s1, C3, C4);
		mVec3 C = linear(t - t1, C12, C34);
		return C;

	}
};

inline void GraphicsPipeline::Render(const RenderObject& yu,  framebuffer_t* Fb) {


	timer.begin();

	auto mt=yu.Material();

	VertexesProcess(yu);
	Rasterization(mt);
	FragmentProcess(Fb);
	timer.end();
	clearPipeline();

	//free(depthbuffer);

	if(UsePhongShading)
	std::cout << "PhongShading FPS=" << 1 / timer.interval << "." << std::endl;
	else 
	std::cout << "GauroudShading FPS=" << 1 / timer.interval << "." << std::endl;

}

inline void GraphicsPipeline::clearPipeline() {
	TargetRenderTriangles.clear();
	//TargetRenderTriangles.swap(vector<Triangle>());
	FragmentsAfterRasterization.clear();
	//	FragmentsAfterRasterization.swap(vector<fragment>());
}

inline void GraphicsPipeline::VertexesProcess(const RenderObject& yu ) {
	auto Vertexshading = [this](mVec3 point, mVec3 normal, ShadingMaterial Mp) {
		mVec3 Ia = Mp.Ka*LI;
		mVec3  Sourse = LightSource;
		mVec3 eye = mVec3(0,0,0);
		//Reflection
		mVec3 ReflectionD = mVec3(Sourse, point);
		float distance2 = pow((Sourse.x - point.x), 2) + pow((Sourse.y - point.y), 2) + pow((Sourse.z - point.z), 2);
		float Intensity = (LI / distance2);
		mVec3 light_dir(Sourse, point);
		mVec3 eye_dir(eye, point);
		light_dir.normalize();
		eye_dir.normalize();
		mVec3 h = (eye_dir + light_dir);
		h.normalize();
		float dotproduct = normal * light_dir;
		mVec3 Id = Mp.Kd*Intensity*(dotproduct > 0 ? dotproduct : 0);
		float dotproduct2 = normal * h;
		float tmp = dotproduct2 > 0.0f ? dotproduct2 : 0.0f;
		tmp = pow(tmp, Mp.f);
		mVec3 Is = Mp.Ks*Intensity*tmp;
		mVec3  If = (Is + Id + Ia);
		If.ColorClamp(); 
		return If;
	};

	//-------------------------------
	//**** Per-vertex Operation  ****
	//Note:
	//Matrix4 P = squashMatrix(_Camera.getNearPlane(), _Camera.getFarPlane());
	//Matrix4 Ortho = OrthoMatrix(_Camera.r, _Camera.l, _Camera.t, _Camera.b, _Camera.getNearPlane(), _Camera.getFarPlane());
	//Per=Ortho*P
	//-------------------------------
	//GET M V,P VP matrix
	Matrix4 M = yu.ModleMatrix;
	Matrix4 V = _Camera.genViewMat();
	Matrix4 Per = _Camera.genPerspectiveMat();
	Matrix4 Vp = ViewPortMatrix(w,h);

	Matrix4 ViewModel = ((V *M));
	Matrix4 PersViewPort =  (Vp *Per );
	Matrix4 NoramlViewModel = inv(ViewModel).Transpose();


	//-------------------------------
	// apply  MVP projction +Pervertex lighting
	//-------------------------------
	mVec4 LightSourceHomo = V * mVec4(InitLightPos, 1);
	LightSource = LightSourceHomo.tomVec3();

	int Vid = 0;
	std::vector<vertex>TargetRenderVertexes;
	TargetRenderVertexes.reserve(yu.Vertexes.size());
	for (const mVec3& Vertex : yu.Vertexes) {
		//-------------------------------
		//  Model + View  transformation
		//-------------------------------
		mVec4 EyeHomoCoordinates = (ViewModel*  mVec4(Vertex, 1));
		mVec3 vertex_normal = (NoramlViewModel * mVec4(yu.VertexesNormal[Vid], 0)).tomVec3();
		//-------------------------------
		//Pervertex  LightSource
		//-------------------------------
		mVec3 shading = Vertexshading(EyeHomoCoordinates.tomVec3(), yu.VertexesNormal[Vid], yu.Material());
		//-------------------------------
		//Perpective + Viewport transformation
	   //-------------------------------
		mVec4 HomoCoordinates = (PersViewPort*  EyeHomoCoordinates);
		HomoCoordinates = HomoCoordinates / HomoCoordinates.w;
		vertex tmp = { HomoCoordinates.tomVec3(),vertex_normal,shading, EyeHomoCoordinates.tomVec3(),yu.VertexesTex[Vid]};
		TargetRenderVertexes.emplace_back(tmp);
		Vid++;
	}
	//-------------------------------
	//primitives assembly
	//-------------------------------
	TargetRenderTriangles.reserve(yu.TrianglesIdx.size());
	for (auto idxs : yu.TrianglesIdx) {
		TargetRenderTriangles.emplace_back(TriangleWithAttributes(TargetRenderVertexes[idxs.x], TargetRenderVertexes[idxs.y], TargetRenderVertexes[idxs.z]));
	}


	

}

inline void  GraphicsPipeline::Rasterization(const ShadingMaterial& myu) {
	auto clamp = [](float in)->float {
		return in > 1.0f ? 1.0f : (in < 0.0f ? 0.0f : in);
	};
	auto fract = [](float x)->float {
		return x - floor(x);
	};
	//after perspective ortho viewport z=n+f/n-f - 2*f*n/z(n-f)
	float n = _Camera.getNearPlane();
	float f = _Camera.getFarPlane();
	float depthA = (n + f) / (n - f);
	float zcoord = f * n * 2;
	auto ADSshading = [this](mVec3 point, mVec3 normal, ShadingMaterial Mp) {
		mVec3 Ia = Mp.Ka*LI;
		mVec3  Sourse = LightSource;
		mVec3 eye = _Camera.position;
		//Reflection
		mVec3 ReflectionD = mVec3(Sourse, point);
		float distance2 = pow((Sourse.x - point.x), 2) + pow((Sourse.y - point.y), 2) + pow((Sourse.z - point.z), 2);
		float Intensity = (LI / distance2);
		mVec3 light_dir(Sourse, point);
		mVec3 eye_dir(eye, point);
		light_dir.normalize();
		eye_dir.normalize();
		//normal.normalize();
		mVec3 h = (eye_dir + light_dir);
		h.normalize();
		float dotproduct = normal * light_dir;
		mVec3 Id = Mp.Kd*Intensity*(dotproduct > 0 ? dotproduct : 0);
		float dotproduct2 = normal * h;
		float tmp = dotproduct2 > 0.0f ? dotproduct2 : 0.0f;
		tmp = pow(tmp, Mp.f);
		mVec3 Is = Mp.Ks*Intensity*tmp;
		mVec3  If = (Is + Id + Ia);
		If.ColorClamp();
		return If;
	};

	
	auto FramentShader = [&](TriangleWithAttributes& Tri, mVec2<float>uv)->mVec4 {
		float iaz = 1.0f / (Tri.a.Coordinate.z - depthA);
		float ibz = 1.0f / (Tri.b.Coordinate.z - depthA);
		float icz = 1.0f / (Tri.c.Coordinate.z - depthA);
		float interpolateIZ = iaz + uv.x  * (icz - iaz) + uv.y* (ibz - iaz);//P = A + u * (C - A) + v * (B - A)  
		float interpolateZ = (1.0f / interpolateIZ);
		mVec3 tmp;
		if (UsePhongShading) {
			//interpolate normal
			mVec3 normal = Tri.a.normal + ((Tri.c.normal - Tri.a.normal)*  uv.x + (Tri.b.normal - Tri.a.normal)* uv.y);
			normal.normalize();//normal in eyespace,light source in eyespace
			mVec3 point = Tri.a.EyeSpaceCoordinate  *(1 - uv.x - uv.y) + (Tri.c.EyeSpaceCoordinate*  uv.x) + (Tri.b.EyeSpaceCoordinate  *uv.y);
			tmp = ADSshading(point, normal, myu);

			//stripe texture
			if(TextureChanel ==1){
				mVec3 backcolor = { 0,0,0 };
				float scale = 100;
				float fuzz = 2;
				float width = 10;
				mVec2<float> Texcoords = Tri.a.st + ((Tri.c.st - Tri.a.st)*  uv.x + (Tri.b.st - Tri.a.st)* uv.y);
			//	Texcoords = Texcoords * interpolateZ;
				float scaleT = fract(Texcoords.y*scale);

				float frac1 = clamp(scaleT / fuzz);
				float frac2 = clamp((scaleT - width) / fuzz);

				frac1 = frac1 * (1.0f - frac2);
				frac1 = frac1 * frac1*(3.0f - (2.0f*frac1));
				tmp = backcolor * (frac1)+tmp * (1 - frac1);
			
			}
			else if (TextureChanel == 2) {
				float scale = 20;
				mVec2<float> Texcoords = Tri.a.st + ((Tri.c.st - Tri.a.st)*  uv.x + (Tri.b.st - Tri.a.st)* uv.y);
				//Texcoords = Texcoords* interpolateZ;
				Texcoords.x= fract(Texcoords.x*scale);
				Texcoords.y = fract(Texcoords.y*scale);
				//mVec3 backcolor =LookupTexel(Texcoords);
				mVec3 backcolor = { 0.5,0.5,0.5};
				tmp = { backcolor.x*tmp.x,  backcolor.y*tmp.y, backcolor.z*tmp.z };
				tmp.ColorClamp();
			}
			


		}
		else {
			//interpolate lighnting
			tmp = Tri.a.shading  *(1 - uv.x - uv.y) + (Tri.c.shading*  uv.x) + (Tri.b.shading  *uv.y);
		

		}
		return { tmp,interpolateZ };
	};





	
	auto RasterAtriangle = [&](TriangleWithAttributes& TriWithAtrib)->std::vector<fragment> {
		auto Tri = TriWithAtrib.GetTriangleVertexes();
		std::vector<fragment> RasterResult;
		int Xmin, Xmax, Ymin, Ymax;
		Xmin = int( floorf(MIN(MIN(Tri.a.x, Tri.b.x), Tri.c.x)));
		Xmax = int(ceilf(MAX(MAX(Tri.a.x, Tri.b.x), Tri.c.x)));
		Ymin = int(floorf(MIN(MIN(Tri.a.y, Tri.b.y), Tri.c.y)));
		Ymax = int(ceilf(MAX(MAX(Tri.a.y, Tri.b.y), Tri.c.y)));
		if(Xmin<0|| Xmin>=w || Xmax<0 || Xmax>=w  )return RasterResult;
		if (Ymin<0 || Ymin>=h|| Ymax<0 || Ymax>=h )return RasterResult;
		if (Tri.a.z > 0 || Tri.b.z > 0 || Tri.c.z > 0 )return RasterResult;
		RasterResult.reserve((Xmax - Xmin)*(Ymax - Ymin));


		mVec2<float>FirstTime;
		mVec2<float>Incremental_uv;
		mVec2<float>Lastuv;
	
		for (int x = Xmin; x <= Xmax; x++) {
			
			for (int y = Ymin; y <= Ymax; y++) {
				if (y == Ymin) {
					mVec2<float>  uv = Tri.PointIsInTriangle(mVec3(x, y, 0));
					if (uv.x >= 0 && uv.y >= 0 && uv.x + uv.y <= 1) {
						mVec4 attibs = FramentShader(TriWithAtrib, uv);
						RasterResult.emplace_back(fragment{ mVec2<int>{x, y } ,attibs.tomVec3(),attibs.w });
					}
					FirstTime = uv;


				}
				else if (y == Ymin + 1) {
					mVec2<float>  uv = Tri.PointIsInTriangle(mVec3(x, y, 0));
					if (uv.x >= 0 && uv.y >= 0 && uv.x + uv.y <= 1) {
					
						mVec4 attibs = FramentShader(TriWithAtrib, uv);
						RasterResult.emplace_back(fragment{ mVec2<int>{x, y } ,attibs.tomVec3(),attibs.w });
					}
					Incremental_uv = uv - FirstTime;
					Lastuv = uv;
				}
				else {
					mVec2<float> uv = Lastuv + Incremental_uv;
					Lastuv = uv;
					if (uv.x >= 0 && uv.y >= 0 && uv.x + uv.y <= 1) {
						mVec4 attibs = FramentShader(TriWithAtrib, uv);
						RasterResult.emplace_back(fragment{ mVec2<int>{x, y } ,attibs.tomVec3(),attibs.w });
					}
				}

			}
		}

		return RasterResult;

	};





	int TrianglesNumber = TargetRenderTriangles.size();
	FragmentsAfterRasterization.reserve(TrianglesNumber *50); 
#ifdef SingleThread
	for (auto& tri : TargetRenderTriangles) {
		auto Frags = RasterAtriangle(tri);

		FragmentsAfterRasterization.insert(FragmentsAfterRasterization.end(), Frags.begin(), Frags.end());
	}
#else
	//Rasterization
	std::mutex g_num_mutex;
	auto RastertrianglesPerThread = [&](int startId, int endID) {
		std::vector<fragment> ThisTheadFragments;
		ThisTheadFragments.reserve((endID - startId) * 50);
		for (int j = startId; j < endID; j++) {
			auto Frags = RasterAtriangle(TargetRenderTriangles[j]);
			ThisTheadFragments.insert(ThisTheadFragments.end(), Frags.begin(), Frags.end());
		}
		g_num_mutex.lock();
		FragmentsAfterRasterization.insert(FragmentsAfterRasterization.end(), ThisTheadFragments.begin(), ThisTheadFragments.end());
		g_num_mutex.unlock();
		return ;
	};
	int threadcount = 16;
	int WorkloadPerThread = floorf(TrianglesNumber / threadcount);
	std::vector<  std::thread >Threads;
	for (int i = 0; i < threadcount; i++) {
		std::thread t(RastertrianglesPerThread, WorkloadPerThread*i, (i + 1) * WorkloadPerThread);
		Threads.push_back(std::move(t));
	}
	for (auto& th : Threads) {
		th.join();
	}
#endif


	
		

	//interpolate
}

inline void GraphicsPipeline::FragmentProcess(framebuffer_t* Fb) {
	float* depthbuffer = Fb->depth_buffer;
	assert(depthbuffer != NULL);
	for (int i = 0; i < w * h; i++)depthbuffer[i] = (float)INT32_MAX;
	for (fragment TriFrag : FragmentsAfterRasterization) {
		mVec2<int> ScreenXY = { TriFrag.XY.x, h-(TriFrag.XY.y +1)};
		double distance = abs(TriFrag.depth -_Camera.position.z);
		if (depthbuffer[ScreenXY.y*w + ScreenXY.x] > distance ) {
			unsigned char fragColor[4] = { TriFrag.RGB.x*255,TriFrag.RGB.y * 255,TriFrag.RGB.z * 255,0 };

			Fb->setvalue(ScreenXY.x, ScreenXY.y,fragColor);
			//Fb->color_buffer[ScreenXY.y*w + ScreenXY.x] = TriFrag.RGB;

			depthbuffer[ScreenXY.y*w + ScreenXY.x] = distance;
		}
	} 

}


#ifdef CLIP
if (Clip) {
	int Vid = 0;
	for (mVec3& Vertex : yu.Vertexes) {
		//-------------------------------
		//  Model + View  transformation
		//-------------------------------
		mVec4 EyeHomoCoordinates = (ViewModel*  mVec4(Vertex, 1));
		yu.VertexesNormal[Vid] = (NoramlViewModel * mVec4(yu.VertexesNormal[Vid], 0)).tomVec3();
		//-------------------------------
		//Pervertex  LightSource
		//-------------------------------
		mVec3 shading = Vertexshading(EyeHomoCoordinates.tomVec3(), yu.VertexesNormal[Vid], yu.Material());

		vertex4 tmp = { EyeHomoCoordinates, yu.VertexesNormal[Vid],shading, EyeHomoCoordinates.tomVec3() };
		TargetRenderVertexes.emplace_back(tmp);
		Vid++;
	}
	//-------------------------------
	//primitives assembly
	//-------------------------------
	std::vector<Vertex4TriangleWithAttributes> tmpTargetRenderTriangles;
	tmpTargetRenderTriangles.reserve(yu.TrianglesIdx.size());
	TargetRenderTriangles.reserve(yu.TrianglesIdx.size());
	for (auto idxs : yu.TrianglesIdx) {
		tmpTargetRenderTriangles.emplace_back(Vertex4TriangleWithAttributes(TargetRenderVertexes[idxs.x], TargetRenderVertexes[idxs.y], TargetRenderVertexes[idxs.z]));
	}
	//-------------------------------
	//todo:clipping cullling
	//-------------------------------
	float thta = Radians(_Camera.fov / 2);
	Plane topPlane = { { 0,0,0}, {0,cos(thta),sin(thta) } };
	Plane bottomPlane = { { 0,0,0} ,{0,-cos(thta),sin(thta) } };
	Plane leftPlane = { { 0,0,0} ,{-cos(thta),0,sin(thta) } };
	Plane rightPlane = { { 0,0,0} ,{cos(thta),0,sin(thta) } };
	Plane nearPlane = { { 0,0,_Camera.getNearPlane()}, {0,0,1 } };
	Plane farPlane = { { 0,0,_Camera.getFarPlane()}, {0,0,-1} };
	std::vector<Vertex4TriangleWithAttributes> tmpTargetRenderTriangles1;
	std::vector<Vertex4TriangleWithAttributes> tmpTargetRenderTriangles2;
	std::vector<Vertex4TriangleWithAttributes> tmpTargetRenderTriangles3;
	std::vector<Vertex4TriangleWithAttributes> tmpTargetRenderTriangles4;
	std::vector<Vertex4TriangleWithAttributes> tmpTargetRenderTriangles5;
	std::vector<Vertex4TriangleWithAttributes> tmpTargetRenderTriangles6;
	tmpTargetRenderTriangles1.reserve(yu.TrianglesIdx.size() * 1.2);//1->2
	tmpTargetRenderTriangles2.reserve(yu.TrianglesIdx.size() * 1.2);//1->2
	tmpTargetRenderTriangles3.reserve(yu.TrianglesIdx.size() *  1.2);//1->2
	tmpTargetRenderTriangles4.reserve(yu.TrianglesIdx.size() *  1.2);//1->2
	tmpTargetRenderTriangles5.reserve(yu.TrianglesIdx.size() *  1.2);//1->2
	tmpTargetRenderTriangles6.reserve(yu.TrianglesIdx.size() *  1.2);//1->2
	TargetRenderTriangles.reserve(yu.TrianglesIdx.size() * 2);//1->2

	clip_with_plane(topPlane, tmpTargetRenderTriangles, tmpTargetRenderTriangles1);
	clip_with_plane(bottomPlane, tmpTargetRenderTriangles1, tmpTargetRenderTriangles2);
	clip_with_plane(leftPlane, tmpTargetRenderTriangles2, tmpTargetRenderTriangles3);
	clip_with_plane(rightPlane, tmpTargetRenderTriangles3, tmpTargetRenderTriangles4);
	clip_with_plane(nearPlane, tmpTargetRenderTriangles4, tmpTargetRenderTriangles5);
	clip_with_plane(farPlane, tmpTargetRenderTriangles, tmpTargetRenderTriangles6);
	//-------------------------------
	//Perpective + Viewport transformation
	//-------------------------------
	for (auto& tri : tmpTargetRenderTriangles6) {
		tri.a.Coordinate = (PersViewPort * tri.a.Coordinate);
		tri.a.Coordinate = (tri.a.Coordinate / tri.a.Coordinate.w);
		tri.b.Coordinate = (PersViewPort * tri.b.Coordinate);
		tri.b.Coordinate = (tri.b.Coordinate / tri.b.Coordinate.w);
		tri.c.Coordinate = (PersViewPort * tri.c.Coordinate);
		tri.c.Coordinate = (tri.c.Coordinate / tri.c.Coordinate.w);
		TargetRenderTriangles.emplace_back(tri.toTriangleWithAttributes());

	}
}
#else
#endif


