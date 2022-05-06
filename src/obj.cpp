#include <execution>
#include <mutex>
#include "matrix.h"
#include "assert.h"
#include "obj.h"
#include "timer.hpp"

#define LI 100
using namespace cimg_library;

static mVec3f InitEyePos = mVec3f(0, 0, 200);
static mVec3f InitGazeDirection = mVec3f(0, 0, -1);
static mVec3f InitTopDirection = mVec3f(0, 1, 0);


void GraphicsPipeline::LoadTexture(std::filesystem::path tex_path) {
	CImg<unsigned  char> TextureImg(tex_path.string().c_str());
	Texture  tX;
	tX.TW = TextureImg.width();
	tX.TH = TextureImg.height();
	tX.TextureArr = (mVec3f*)malloc(tX.TW * tX.TH * sizeof(mVec3f));
	for (int r = 0; r < tX.TH; r++)
		for (int c = 0; c < tX.TW; c++) {
			tX.TextureArr[r * tX.TW + c].x = float(TextureImg(c, r, 0)) / 255.0f;
			tX.TextureArr[r * tX.TW + c].y = float(TextureImg(c, r, 1)) / 255.0f;
			tX.TextureArr[r * tX.TW + c].z = float(TextureImg(c, r, 2)) / 255.0f;
		}
	mTextures.push_back(tX);
	return ;
}

void  GraphicsPipeline::LoadTexture(std::string path) {
	CImg<unsigned  char> TextureImg(path.c_str());
	Texture  tX;

	tX.TW = TextureImg.width();
	tX.TH = TextureImg.height();
	tX.TextureArr = (mVec3f*)malloc(tX.TW * tX.TH * sizeof(mVec3f));
	for (int r = 0; r < tX.TH; r++)
		for (int c = 0; c < tX.TW; c++) {
			tX.TextureArr[r * tX.TW + c].x = float(TextureImg(c, r, 0)) / 255.0f;
			tX.TextureArr[r * tX.TW + c].y = float(TextureImg(c, r, 1)) / 255.0f;
			tX.TextureArr[r * tX.TW + c].z = float(TextureImg(c, r, 2)) / 255.0f;
		}
	mTextures.push_back(tX);
	return ;
}




void GraphicsPipeline::VertexesProcess(const RenderableObject& yu) {

	auto Vertexshading = [this](mVec3f point, mVec3f normal, ShadingMaterial Mp) {
		mVec3f Ia = Mp.Ka * LI;
		mVec3f  Sourse = LightSource;
		mVec3f eye = mVec3f(0, 0, 0);
		//Reflection
		mVec3f ReflectionD = mVec3f(Sourse, point);
		float distance2 = pow((Sourse.x - point.x), 2) + pow((Sourse.y - point.y), 2) + pow((Sourse.z - point.z), 2);
		float Intensity = (LI / distance2);
		mVec3f light_dir(Sourse, point);
		mVec3f eye_dir(eye, point);
		light_dir.normalize();
		eye_dir.normalize();
		mVec3f h = (eye_dir + light_dir);
		h.normalize();
		float dotproduct = normal * light_dir;
		mVec3f Id = Mp.Kd * Intensity * (dotproduct > 0 ? dotproduct : 0);
		float dotproduct2 = normal * h;
		float tmp = dotproduct2 > 0.0f ? dotproduct2 : 0.0f;
		tmp = pow(tmp, Mp.f);
		mVec3f Is = Mp.Ks * Intensity * tmp;
		mVec3f  If = (Is + Id + Ia);
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
	Matrix4 Vp = ViewPortMatrix(window_width, window_height);

	Matrix4 ViewModel = ((V * M));
	Matrix4 PersViewPort = (Vp * Per);
	Matrix4 NoramlViewModel = inv(ViewModel).Transpose();


	//-------------------------------
	// apply  MVP projction +Pervertex lighting
	//-------------------------------
	mVec4f LightSourceHomo = V * mVec4f(InitLightPos, 1);
	LightSource = LightSourceHomo.tomVec3();

	int Vid = 0;
	std::vector<vertex>TargetRenderVertexes;
	TargetRenderVertexes.reserve(yu.Vertexes.size());
	for (const mVec3f& Vertex : yu.Vertexes) {
		//-------------------------------
		//  Model + View  transformation
		//-------------------------------
		mVec4f EyeHomoCoordinates = (ViewModel * mVec4f(Vertex, 1));
		mVec3f vertex_normal = (NoramlViewModel * mVec4f(yu.VertexesNormal[Vid], 0)).tomVec3();
		//-------------------------------
		//Pervertex  LightSource
		//-------------------------------
		mVec3f shading = Vertexshading(EyeHomoCoordinates.tomVec3(), yu.VertexesNormal[Vid], yu.Material());
		//-------------------------------
		//Perpective + Viewport transformation
	   //-------------------------------
		mVec4f HomoCoordinates = (PersViewPort * EyeHomoCoordinates);
		HomoCoordinates = HomoCoordinates / HomoCoordinates.w;
		vertex tmp = { HomoCoordinates.tomVec3(),vertex_normal,shading, EyeHomoCoordinates.tomVec3(),yu.VertexesTex[Vid] };
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


template<typename T>
static inline T PerspetiveCorrectDepth(T iaz, T ibz, T icz, T uv_u, T uv_v) {

	float interpolateIZ = iaz + uv_u * (icz - iaz) + uv_v * (ibz - iaz);//P = A + u * (C - A) + v * (B - A)  
	float  interpolateZ = (1.0f / interpolateIZ);
	return interpolateZ;
}



mVec3f GraphicsPipeline::ADSshading(mVec3f point, mVec3f normal, ShadingMaterial Mp) const {
	mVec3f Ia = Mp.Ka * LI;
	mVec3f  Sourse = LightSource;
	mVec3f eye = _Camera.position;
	//Reflection
	mVec3f ReflectionD = mVec3f(Sourse, point);
	float distance2 = pow((Sourse.x - point.x), 2) + pow((Sourse.y - point.y), 2) + pow((Sourse.z - point.z), 2);
	float Intensity = (LI / distance2);
	mVec3f light_dir(Sourse, point);
	mVec3f eye_dir(eye, point);
	light_dir.normalize();
	eye_dir.normalize();
	//normal.normalize();
	mVec3f h = (eye_dir + light_dir);
	h.normalize();
	float dotproduct = normal * light_dir;
	mVec3f Id = Mp.Kd * Intensity * (dotproduct > 0 ? dotproduct : 0);
	float dotproduct2 = normal * h;
	float tmp = dotproduct2 > 0.0f ? dotproduct2 : 0.0f;
	tmp = pow(tmp, Mp.f);
	mVec3f Is = Mp.Ks * Intensity * tmp;
	mVec3f  If = (Is + Id + Ia);
	If.ColorClamp();
	return If;
};


void GraphicsPipeline::Rasterization(const ShadingMaterial& myu, framebuffer_t* Fb) const
{
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



	auto FramentShader = [&](const TriangleWithAttributes& Tri, mVec2<float>uv)->mVec4f {


		mVec3f tmp;
		if (UsePhongShading) {
			//interpolate normal
			mVec3f normal = Tri.a.normal + ((Tri.c.normal - Tri.a.normal) * uv.x + (Tri.b.normal - Tri.a.normal) * uv.y);
			normal.normalize();//normal in eyespace,light source in eyespace
			mVec3f point = Tri.a.EyeSpaceCoordinate * (1 - uv.x - uv.y) + (Tri.c.EyeSpaceCoordinate * uv.x) + (Tri.b.EyeSpaceCoordinate * uv.y);
			tmp = ADSshading(point, normal, myu);

			//stripe texture
			if (TextureMode== 1) {
				mVec3f backcolor = { 0,0,0 };
				float scale = 100;
				float fuzz = 2;
				float width = 10;
				mVec2<float> Texcoords = Tri.a.st + ((Tri.c.st - Tri.a.st) * uv.x + (Tri.b.st - Tri.a.st) * uv.y);
				//	Texcoords = Texcoords * interpolateZ;
				float scaleT = fract(Texcoords.y * scale);

				float frac1 = clamp(scaleT / fuzz);
				float frac2 = clamp((scaleT - width) / fuzz);

				frac1 = frac1 * (1.0f - frac2);
				frac1 = frac1 * frac1 * (3.0f - (2.0f * frac1));
				tmp = backcolor * (frac1)+tmp * (1 - frac1);

			}
			else if (TextureMode == 2) {
				float scale = 20;
				mVec2<float> Texcoords = Tri.a.st + ((Tri.c.st - Tri.a.st) * uv.x + (Tri.b.st - Tri.a.st) * uv.y);
				//Texcoords = Texcoords* interpolateZ;
				Texcoords.x = fract(Texcoords.x * scale);
				Texcoords.y = fract(Texcoords.y * scale);
				//mVec3f backcolor =LookupTexel(Texcoords);
				mVec3f backcolor = { 0.5,0.5,0.5 };
				tmp = { backcolor.x * tmp.x,  backcolor.y * tmp.y, backcolor.z * tmp.z };
				tmp.ColorClamp();
			}



		}
		else {
			//interpolate lighnting
			tmp = Tri.a.shading * (1 - uv.x - uv.y) + (Tri.c.shading * uv.x) + (Tri.b.shading * uv.y);


		}
		return { tmp,1.0 };//opacity not used for now.
	};



	int TrianglesNumber = TargetRenderTriangles.size();


#ifdef SingleThread
	for (auto& tri : TargetRenderTriangles) {
		auto Frags = RasterAtriangle(tri);

		FragmentsAfterRasterization.insert(FragmentsAfterRasterization.end(), Frags.begin(), Frags.end());
	}
#else
	//Rasterization
	/*
		//std::mutex g_num_mutex;
	//auto RastertrianglesPerThread = [&](int startId, int endID) {
	//	std::vector<fragment> ThisTheadFragments;
	//	ThisTheadFragments.reserve((endID - startId) * 50);
	//	for (int j = startId; j < endID; j++) {
	//		 RasterAtriangle(TargetRenderTriangles[j], ThisTheadFragments);
	//		//hisTheadFragments.insert(ThisTheadFragments.end(), Frags.begin(), Frags.end());
	//	}
	//	g_num_mutex.lock();
	//	FragmentsAfterRasterization.insert(FragmentsAfterRasterization.end(), ThisTheadFragments.begin(), ThisTheadFragments.end());
	//	g_num_mutex.unlock();
	//	return ;
	//};
	int WorkloadPerThread = floorf(TrianglesNumber / threadcount);
	std::vector<  std::thread >Threads;
	for (int i = 0; i < threadcount; i++) {
		std::thread t(RastertrianglesPerThread, WorkloadPerThread * i, (i + 1) * WorkloadPerThread);
		Threads.push_back(std::move(t));
	}
	for (auto& th : Threads) {
		th.join();
	}

	*/

	/*const int nT = TargetRenderTriangles.size();
	FragmentsAfterRasterization.resize(nT);*/

	float* depthbuffer = Fb->depth_buffer;
	assert(depthbuffer != NULL);
	memset(depthbuffer, 0x7f, DEFAULT_WINDOW_HEIGHT * DEFAULT_WINDOW_WIDTH * sizeof(float));

	auto EarlyZtest = [&](int x, int y, float z)->bool {
		//mVec2<int> ScreenXY = { x, window_height - (y + 1) };
		mVec2<int> ScreenXY = { y ,x};
		double distance = abs(z - _Camera.position.z);
		if (depthbuffer[ScreenXY.x * window_width + ScreenXY.y] > distance)
		{
			depthbuffer[ScreenXY.x * window_width + ScreenXY.y] = distance;
			return  true;
		}
		return false;

	};


	auto outputshading = [&](const fragment& TriFrag) {
		//mVec2<int> ScreenXY = {  window_height - (TriFrag.XY.y + 1) ,TriFrag.XY.x };
		mVec2<int> ScreenXY = { TriFrag.XY.y ,TriFrag.XY.x };
		double distance = abs(TriFrag.depth - _Camera.position.z);

		{
			//std::lock_guard<std::mutex> lock(buffer_mutex);
			if (depthbuffer[ScreenXY.x * window_width + ScreenXY.y] > distance) {
				unsigned char fragColor[4] = { static_cast<unsigned char>(TriFrag.RGBv.x * 255),static_cast<unsigned char>(TriFrag.RGBv.y * 255),static_cast<unsigned char>(TriFrag.RGBv.z * 255),0 };
				Fb->setvalue(ScreenXY.x, ScreenXY.y, fragColor);
				//Fb->color_buffer[ScreenXY.y*w + ScreenXY.x] = TriFrag.RGB;
				depthbuffer[ScreenXY.x * window_width + ScreenXY.y] = distance;
			}

		}


	};


	auto RasterAtriangle = [&](const TriangleWithAttributes& TriWithAtrib) {
		constexpr bool UseEarlyZ = true;


		constexpr size_t MaxFragCount = DEFAULT_WINDOW_WIDTH * DEFAULT_WINDOW_WIDTH * 0.5;
		auto Tri = TriWithAtrib.GetTriangleVertexes();
		int Xmin, Xmax, Ymin, Ymax;
		Xmin = int(floorf(MIN(MIN(Tri.a.x, Tri.b.x), Tri.c.x)));
		Xmax = int(ceilf(MAX(MAX(Tri.a.x, Tri.b.x), Tri.c.x)));
		Ymin = int(floorf(MIN(MIN(Tri.a.y, Tri.b.y), Tri.c.y)));
		Ymax = int(ceilf(MAX(MAX(Tri.a.y, Tri.b.y), Tri.c.y)));
		if (Xmin < 0 || Xmin >= window_width || Xmax < 0 || Xmax >= window_width)return;
		if (Ymin < 0 || Ymin >= window_height || Ymax < 0 || Ymax >= window_height)return;
		if (Tri.a.z > 0 || Tri.b.z > 0 || Tri.c.z > 0)return;

		//std::vector<fragment>RasterResult;
		//RasterResult.reserve((Xmax - Xmin)* (Ymax - Ymin) * 8);

		//early-z
		const float iaz = 1.0f / (TriWithAtrib.a.Coordinate.z - depthA);
		const float ibz = 1.0f / (TriWithAtrib.b.Coordinate.z - depthA);
		const float icz = 1.0f / (TriWithAtrib.c.Coordinate.z - depthA);



		mVec2<float>FirstTime;
		mVec2<float>Incremental_uv;
		mVec2<float>Lastuv;

		for (int x = Xmin; x <= Xmax; x++) {
			for (int y = Ymin; y <= Ymax; y++) {
				if (y == Ymin) {
					mVec2<float>  uv = Tri.PointIsInTriangle(mVec3f(x, y, 0));
					if (uv.x >= 0 && uv.y >= 0 && uv.x + uv.y <= 1) {
						if constexpr (UseEarlyZ)
						{
							float interpolatez = PerspetiveCorrectDepth(iaz, ibz, icz, uv.x, uv.y);
							if (!EarlyZtest(x, y, interpolatez))continue;

						}
						mVec4f attibs = FramentShader(TriWithAtrib, uv);
						//RasterResult.emplace_back(attibs.w, mVec2<int>{x, y }, attibs.tomVec3f());
						outputshading({ attibs.w, mVec2<int>{x, y }, attibs.tomVec3() });

					}
					FirstTime = uv;


				}
				else if (y == Ymin + 1) {
					mVec2<float>  uv = Tri.PointIsInTriangle(mVec3f(x, y, 0));
					if (uv.x >= 0 && uv.y >= 0 && uv.x + uv.y <= 1) {

						if constexpr (UseEarlyZ)
						{
							float interpolatez = PerspetiveCorrectDepth(iaz, ibz, icz, uv.x, uv.y);
							if (!EarlyZtest(x, y, interpolatez))continue;

						}

						mVec4f attibs = FramentShader(TriWithAtrib, uv);
						//RasterResult.emplace_back(attibs.w ,mVec2<int>{x, y } ,attibs.tomVec3f());
						outputshading({ attibs.w, mVec2<int>{x, y }, attibs.tomVec3() });
					}
					Incremental_uv = uv - FirstTime;
					Lastuv = uv;
				}
				else {
					mVec2<float> uv = Lastuv + Incremental_uv;
					Lastuv = uv;
					if (uv.x >= 0 && uv.y >= 0 && uv.x + uv.y <= 1) {

						if constexpr (UseEarlyZ)
						{
							float interpolatez = PerspetiveCorrectDepth(iaz, ibz, icz, uv.x, uv.y);
							if (!EarlyZtest(x, y, interpolatez))continue;

						}

						mVec4f attibs = FramentShader(TriWithAtrib, uv);
						//RasterResult.emplace_back(attibs.w, mVec2<int>{x, y }, attibs.tomVec3f());
						outputshading({ attibs.w, mVec2<int>{x, y }, attibs.tomVec3() });
					}
				}

			}
		}

		//return 0;




	};

	std::vector<int>count(TargetRenderTriangles.size());
	auto policy = std::execution::par_unseq;
	//std::transform(policy, TargetRenderTriangles.begin(), TargetRenderTriangles.end(), FragmentsAfterRasterization.begin(), count.begin(),RasterAtriangle);//3fps
	std::for_each(policy, TargetRenderTriangles.begin(), TargetRenderTriangles.end(), RasterAtriangle);



#endif





	//interpolate
}

void GraphicsPipeline::Render(const RenderableObject& robj, framebuffer_t* Fb) {
	timer.begin();
	auto material = robj.Material();
	VertexesProcess(robj);
	Rasterization(material, Fb);
	//FragmentProcess(Fb);
	timer.end();
	clearPipeline();
	if (UsePhongShading)
		std::cout << "PhongShading FPS=" << 1 / timer.interval << "." << std::endl;
	else
		std::cout << "GauroudShading FPS=" << 1 / timer.interval << "." << std::endl;
}



void GraphicsPipeline::clearPipeline() {
	TargetRenderTriangles.clear();
	//TargetRenderTriangles.swap(vector<Triangle>());
	//FragmentsAfterRasterization.clear();
	//	FragmentsAfterRasterization.swap(vector<fragment>());
}


GraphicsPipeline::GraphicsPipeline(int h, int w) : window_height(h), window_width(w) {
	//mTextures.emplace_back(LoadTexture());

	Clip = false;
	UsePhongShading = false;
	_Camera = Camera(InitEyePos, InitGazeDirection, InitTopDirection);
	_Camera.SetFrustm(1, -1, -1, 1, 90.0f, InitEyePos.z + 200);

	//pre-allocate memory:
	//FragmentsAfterRasterization.resize(4096);


};






#ifdef CLIP
typedef struct vertex4 {
	mVec4 Coordinate;

	//attributes
	mVec3f normal;
	mVec3f shading;
	mVec3f EyeSpaceCoordinate;
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
		return  Triangle(this->a.Coordinate.tomVec3f(), this->b.Coordinate.tomVec3f(), this->c.Coordinate.tomVec3f());

	}*/
	Vertex4TriangleWithAttributes(vertex4 a, vertex4 b, vertex4 c) :a(a), b(b), c(c) {

	}
	TriangleWithAttributes toTriangleWithAttributes() {
		vertex  a1 = { a.Coordinate.tomVec3f(),a.normal,a.shading,a.EyeSpaceCoordinate };
		vertex	b1 = { b.Coordinate.tomVec3f(),b.normal,b.shading,b.EyeSpaceCoordinate };
		vertex	c1 = { c.Coordinate.tomVec3f(),a.normal,c.shading,c.EyeSpaceCoordinate };
		return TriangleWithAttributes(a1, b1, c1);
	}
};
//clip_plane为裁剪平面的自定义结构体，vert_list存储了待裁剪凸多边形的所有顶点
//num_vert为顶点个数，in_list为需要保留下来的裁剪平面内侧顶点的列表
void clip_with_plane(Plane c_plane, std::vector<Vertex4TriangleWithAttributes> in_list, std::vector<Vertex4TriangleWithAttributes>& out_list) {
	auto cal_project_distance = [](Plane ic_plane, mVec3f ptx) ->float {
		return (ptx - ic_plane.px) * ic_plane.n;
	};
	auto cal_insertRatio = [](Plane ic_plane, mVec3f ptx1, mVec3f ptx2) ->float {
		float a = ic_plane.n * ptx1 - ic_plane.n * ic_plane.px;
		float b = ic_plane.n * (ptx1 - ptx2);

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
			mVec3f pre_vertex = vert_list[previous_index].Coordinate.tomVec3f(); //边的起始点
			mVec3f cur_vertex = vert_list[current_index].Coordinate.tomVec3f();  //边的终止点

			float d1 = cal_project_distance(c_plane, pre_vertex);
			float d2 = cal_project_distance(c_plane, cur_vertex);

			//如果该边与裁剪平面有交点，则计算交点并存入in_list
			if (d1 * d2 < 0)
			{
				float t = cal_insertRatio(c_plane, pre_vertex, cur_vertex);//求出t值
				vertex4 intersectionPostion = vert_list[current_index] + vert_list[previous_index] * (1 - t);
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
			out_list.emplace_back(Vertex4TriangleWithAttributes(remainPoints[1], remainPoints[2], remainPoints[3]));
		}
		else if (remaincount != 0) {
			std::cout << " Unknown Error!";
		}







	}


	return;
}


#endif
