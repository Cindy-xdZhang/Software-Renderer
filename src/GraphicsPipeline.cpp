#include <execution>
#include <mutex>
#include "matrix.h"
#include "assert.h"
#include "GraphicsPipeline.h"
#include "timer.hpp"
#include <numeric>


#define LI 200
using namespace cimg_library;

//#TODO: 
//1.resize windows call_back; done
//2.show fps on windows title? done
// 3. marching cube for procedural land generation
// 3.2 debug why scoll mouse change texture?
// 3.3 alpha blend
//4.args to start ray-tracing or graphics pipeline
// 4.1.Clip
//5.sky box
//6.pbr ggx 
//7.water surface


mVec3f InitLightPos = mVec3f(0, 10, 0);
float texture_scaling = 100;

static mVec3f LightSourceEyeSpace ;
static mVec3f InitEyePos = mVec3f(0, 0, 20);
static mVec3f InitGazeDirection = mVec3f(0, 0, -1);
static mVec3f InitTopDirection = mVec3f(0, 1, 0);
 //static std::map<std::string, unsigned int >bufferConvention={
//"position": 0,
//"normal" : 1,
//"textureuv" : 2,
//"per-vertex lighting " : 3,
//"eyespace coordinate": 4,
//"screen space coordinate": 5,
//};

 STRONG_INLINE  mVec3f ADSshading(const mVec3f &eye, const mVec3f& aLightSource, const mVec3f& point, const mVec3f& normal, const ShadingMaterial& Mp) {
	 mVec3f Ia = Mp.Ka * LI;
	 mVec3f  Sourse = aLightSource;
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




void RenderableObject::BuildLikGLBegin(const std::vector<mVec3f>&vtxs,const std::vector<mVec3f>&normals, enum MeshBuildConvention Conv) {
	 assert(Conv == MeshBuildConvention::LIKE_GL_TRIANGLE);
	 Vertexes.clear();
	 Vertexes.resize(vtxs.size());
	 VertexesNormal.clear();
	 VertexesNormal.resize(vtxs.size());
	 std::copy(vtxs.begin(), vtxs.end(), Vertexes.begin());
	 std::copy(normals.begin(), normals.end(), VertexesNormal.begin());

	  TrianglesIdx.clear();
	  int nTriangles = vtxs.size() / 3;
	  TrianglesIdx.reserve(nTriangles);
	 for (int i=0;i< nTriangles; i+=3)
	 {
		 TrianglesIdx.emplace_back( i,i + 1,i + 2 );
	 }
	 VtxTexUV=OBjReader::AssignTextureCoordinates(vtxs);

 }














void GraphicsPipeline::LoadTexture(std::filesystem::path tex_path) {
	LoadTexture(tex_path.string());

}
	

void  GraphicsPipeline::LoadTexture(std::string stringPath) {
	
	auto iter = stringPath.find(".jpg");
	auto iter2 = stringPath.find(".bmp");
	if (iter != std::string::npos)
	{
		CImg<unsigned char > TextureImg(stringPath.c_str()),visu();
		Texture  tX;
		tX.TW = TextureImg.width();
		tX.TH = TextureImg.height();
		tX.TextureArr = (mVec3f*)malloc(static_cast<size_t>(tX.TW) * static_cast<size_t>(tX.TH) * static_cast<size_t>(sizeof(mVec3f)));
		for (int r = 0; r < tX.TH; r++)
			for (int c = 0; c < tX.TW; c++) {
				tX.TextureArr[r * tX.TW + c].x = float(TextureImg(c, r, 0));
				tX.TextureArr[r * tX.TW + c].y = float(TextureImg(c, r, 1));
				tX.TextureArr[r * tX.TW + c].z = float(TextureImg(c, r, 2));
			}
		mTextures.push_back(tX);
	}
	else if (iter2 != std::string::npos) {

		CImg<unsigned char > TextureImg(stringPath.c_str());
		Texture  tX;
		tX.TW = TextureImg.width();
		tX.TH = TextureImg.height();
		tX.TextureArr = (mVec3f*)malloc(static_cast<size_t>(tX.TW) * static_cast<size_t>(tX.TH) * static_cast<size_t>(sizeof(mVec3f)));
		for (int r = 0; r < tX.TH; r++)
			for (int c = 0; c < tX.TW; c++) {
				tX.TextureArr[r * tX.TW + c].x = float(TextureImg(c, r, 0)) / 255.0f;
				tX.TextureArr[r * tX.TW + c].y = float(TextureImg(c, r, 1)) / 255.0f;
				tX.TextureArr[r * tX.TW + c].z = float(TextureImg(c, r, 2)) / 255.0f;
			}

		mTextures.push_back(tX);
	}
	else
	{
		throw("invalid texture file format.");

	}
	return;
}





void GraphicsPipeline::VertexesProcess(const RenderableObject& yu) {

	

	//-------------------------------
	//**** Per-vertex Operation  ****
	//Matrix4 P = squashMatrix(_Camera.getNearPlane(), _Camera.getFarPlane());
	//Matrix4 Ortho = OrthoMatrix(_Camera.r, _Camera.l, _Camera.t, _Camera.b, _Camera.getNearPlane(), _Camera.getFarPlane());
	//Per=Ortho*P
	//-------------------------------

	Matrix4 M = yu.ModleMatrix;
	Matrix4 V = _Camera.genViewMat();
	Matrix4 Per = _Camera.genPerspectiveMat();
	Matrix4 Vp = ViewPortMatrix(window_width, window_height);

	Matrix4 ViewModel = ((V * M));
	Matrix4 PersViewPort = (Vp * Per);
	Matrix4 NoramlViewModel = inv(ViewModel).Transpose();

	// apply  MVP projction +Pervertex lighting
	mVec4f LightSourceHomo = V * mVec4f(InitLightPos, 1);
	LightSourceEyeSpace = (LightSourceHomo/ LightSourceHomo.w).tomVec3();


	int Vid = 0;
	std::vector<vertex>TargetRenderVertexes;
	TargetRenderVertexes.reserve(yu.Vertexes.size());
	for (const mVec3f& Vertex : yu.Vertexes) {
		//-------------------------------
		//  Model + View  transformation
		//-------------------------------
		mVec4f EyeSpaceHomoCoordinates = (ViewModel * mVec4f(Vertex, 1));
		mVec3f EyeSpaceCoordinates = EyeSpaceHomoCoordinates.HomoCordinates2InHomoVec3();

		mVec3f vertex_normal = (NoramlViewModel * mVec4f(yu.VertexesNormal[Vid], 0)).tomVec3();
		//-------------------------------
		//Pervertex  LightSource
		//-------------------------------
		mVec3f shading = ADSshading({0,0,0},LightSourceEyeSpace, EyeSpaceCoordinates, yu.VertexesNormal[Vid], yu.Material());
		//-------------------------------
		//Perpective + Viewport transformation
	   //-------------------------------
		mVec4f ScreenSpaceHomoCoordinates = (PersViewPort * EyeSpaceHomoCoordinates);
		ScreenSpaceHomoCoordinates = ScreenSpaceHomoCoordinates / ScreenSpaceHomoCoordinates.w;


		vertex tmp = { ScreenSpaceHomoCoordinates.tomVec3(),vertex_normal,shading, EyeSpaceHomoCoordinates.tomVec3(),yu.VtxTexUV[Vid] };
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
	//# TODO: Clip



}


template<typename T>
static inline T PerspetiveCorrectDepth(T iaz, T ibz, T icz, T uv_u, T uv_v) {

	float interpolateIZ = iaz + uv_u * (icz - iaz) + uv_v * (ibz - iaz);//P = A + u * (C - A) + v * (B - A)  
	float  interpolateZ = (1.0f / interpolateIZ);
	return interpolateZ;
}




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


	auto FramentShader = [&](const TriangleWithAttributes& Tri, mVec2<float>uv, float interpolateZ)->mVec4f {

		mVec3f tmp;
		if (UsePhongShading) {
			//interpolate normal
			mVec3f normal = Tri.vtxa.normal + ((Tri.vtxc.normal - Tri.vtxa.normal) * uv.x + (Tri.vtxb.normal - Tri.vtxa.normal) * uv.y);
			normal.normalize();//normal in eyespace,light source in eyespace
			mVec3f point = Tri.vtxa.EyeSpaceCoordinate * (1 - uv.x - uv.y) + (Tri.vtxc.EyeSpaceCoordinate * uv.x) + (Tri.vtxb.EyeSpaceCoordinate * uv.y);
			tmp = ADSshading({ 0,0,0 }, LightSourceEyeSpace,point, normal, myu);

			//stripe texture
			if (TextureMode== 1) {
				mVec3f backcolor = { 0,0,0 };
				float scale = texture_scaling;
				float fuzz = 2;
				float width = 10;
				mVec2<float> Texcoords = Tri.vtxa.st + ((Tri.vtxc.st - Tri.vtxa.st) * uv.x + (Tri.vtxb.st - Tri.vtxa.st) * uv.y);
				//	Texcoords = Texcoords * interpolateZ;
				float scaleT = fract(Texcoords.y * scale);

				float frac1 = clamp(scaleT / fuzz);
				float frac2 = clamp((scaleT - width) / fuzz);

				frac1 = frac1 * (1.0f - frac2);
				frac1 = frac1 * frac1 * (3.0f - (2.0f * frac1));
				tmp = backcolor * (frac1)+tmp * (1 - frac1);

			}
			else if (TextureMode == 2) {
				float scale = texture_scaling;
				mVec2<float> Texcoords = Tri.vtxa.st + ((Tri.vtxc.st - Tri.vtxa.st) * uv.x + (Tri.vtxb.st - Tri.vtxa.st) * uv.y);
				Texcoords = Texcoords* interpolateZ;
				Texcoords.x = fract(Texcoords.x * scale);
				Texcoords.y = fract(Texcoords.y * scale);
				mVec3f backcolor =LookupTexel(Texcoords, mTextures.size()-1);
				tmp = { backcolor.x * tmp.x,  backcolor.y * tmp.y, backcolor.z * tmp.z };
				tmp.ColorClamp();
			}



		}
		else {
			//interpolate lighnting
			tmp = Tri.vtxa.shading * (1 - uv.x - uv.y) + (Tri.vtxc.shading * uv.x) + (Tri.vtxb.shading * uv.y);


		}
		return { tmp,1.0 };//opacity not used for now.
	};



	int TrianglesNumber = TargetRenderTriangles.size();

	float* depthbuffer = Fb->depth_buffer;
	assert(depthbuffer != NULL);
	memset(depthbuffer, 0x7f, DEFAULT_WINDOW_HEIGHT * DEFAULT_WINDOW_WIDTH * sizeof(float));
	std::mutex buffer_mutex;
	//Rasterization
	
	auto EarlyZtest = [&](int x, int y, float z)->bool {
		//mVec2<int> ScreenXY = { x, window_height - (y + 1) };
		mVec2<int> ScreenXY = { y ,x };
		//double distance = abs(z - _Camera.position.z);
		//std::lock_guard<std::mutex> lock(buffer_mutex);
		if (depthbuffer[ScreenXY.x * window_width + ScreenXY.y] >= z)
		{
			depthbuffer[ScreenXY.x * window_width + ScreenXY.y] = z;
			return  true;
		}
		return false;

	};


	auto outputshading = [&](const fragment& TriFrag) {
		//mVec2<int> ScreenXY = {  window_height - (TriFrag.XY.y + 1) ,TriFrag.XY.x };
		mVec2<int> ScreenXY = { TriFrag.XY.y ,TriFrag.XY.x };

		{
			//std::lock_guard<std::mutex> lock(buffer_mutex);
			//look at -z direction, this depth is z in screen space 
			if ( depthbuffer[ScreenXY.x * window_width + ScreenXY.y] >= TriFrag.depth) {
				unsigned char fragColor[4] = { static_cast<unsigned char>(TriFrag.RGBv.x * 255),static_cast<unsigned char>(TriFrag.RGBv.y * 255),static_cast<unsigned char>(TriFrag.RGBv.z * 255),0 };
				Fb->setvalue(ScreenXY.x, ScreenXY.y, fragColor);
				//Fb->color_buffer[ScreenXY.y*w + ScreenXY.x] = TriFrag.RGB;
				depthbuffer[ScreenXY.x * window_width + ScreenXY.y] = TriFrag.depth;
			}

		}
	};


	auto RasterAtriangle = [&](const TriangleWithAttributes& TriWithAtrib) {
		constexpr bool UseEarlyZ = true;


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
		const float iaz = 1.0f / (TriWithAtrib.vtxa.Coordinate.z - depthA);
		const float ibz = 1.0f / (TriWithAtrib.vtxb.Coordinate.z - depthA);
		const float icz = 1.0f / (TriWithAtrib.vtxc.Coordinate.z - depthA);



		mVec2<float>FirstTime;
		mVec2<float>Incremental_uv;
		mVec2<float>Lastuv{};

		for (int x = Xmin; x <= Xmax; x++) {
			for (int y = Ymin; y <= Ymax; y++) {
				if (y == Ymin) {
					mVec2<float>  uv = Tri.PointIsInTriangle(mVec3f(x, y, 0));
					if (uv.x >= 0 && uv.y >= 0 && uv.x + uv.y <= 1) {
						float interpolatez = PerspetiveCorrectDepth(iaz, ibz, icz, uv.x, uv.y);
						if constexpr (UseEarlyZ)
							if (!EarlyZtest(x, y, interpolatez))continue;

						mVec4f attibs = FramentShader(TriWithAtrib, uv, interpolatez);
						//RasterResult.emplace_back(attibs.w, mVec2<int>{x, y }, attibs.tomVec3f());
						outputshading({ interpolatez, mVec2<int>{x, y }, attibs.tomVec3() });

					}
					FirstTime = uv;


				}
				else if (y == Ymin + 1) {
					mVec2<float>  uv = Tri.PointIsInTriangle(mVec3f(x, y, 0));
					if (uv.x >= 0 && uv.y >= 0 && uv.x + uv.y <= 1) {

						float interpolatez = PerspetiveCorrectDepth(iaz, ibz, icz, uv.x, uv.y);
						if constexpr (UseEarlyZ)
							if (!EarlyZtest(x, y, interpolatez))continue;

						mVec4f attibs = FramentShader(TriWithAtrib, uv, interpolatez);
						//RasterResult.emplace_back(attibs.w ,mVec2<int>{x, y } ,attibs.tomVec3f());
						outputshading({ interpolatez, mVec2<int>{x, y }, attibs.tomVec3() });
					}
					Incremental_uv = uv - FirstTime;
					Lastuv = uv;
				}
				else {
					mVec2<float> uv = Lastuv + Incremental_uv;
					Lastuv = uv;
					if (uv.x >= 0 && uv.y >= 0 && uv.x + uv.y <= 1) {

						float interpolatez = PerspetiveCorrectDepth(iaz, ibz, icz, uv.x, uv.y);
						if constexpr (UseEarlyZ)
							if (!EarlyZtest(x, y, interpolatez))continue;

						mVec4f attibs = FramentShader(TriWithAtrib, uv, interpolatez);
						//RasterResult.emplace_back(attibs.w, mVec2<int>{x, y }, attibs.tomVec3f());
						outputshading({ interpolatez, mVec2<int>{x, y }, attibs.tomVec3() });
					}
				}

			}
		}

		//return 0;




	};

//#define SingleThread
#ifdef SingleThread
	auto policy = std::execution::seq;
#else
	auto policy = std::execution::par_unseq;
#endif
	std::for_each(policy, TargetRenderTriangles.begin(), TargetRenderTriangles.end(), RasterAtriangle);





	//interpolate
}

void GraphicsPipeline::Render(const RenderableObject& robj, framebuffer_t* Fb) {
	
	auto material = robj.Material();
	VertexesProcess(robj);
	Rasterization(material, Fb);
	//FragmentProcess(Fb);
	
	clearPipeline();
	
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
	
	_Camera = Camera(InitEyePos, InitGazeDirection, InitTopDirection);
	_Camera.SetFrustm(1, -1, -1, 1, 90.0f, InitEyePos.z + 200);

	//pre-allocate memory:
	//FragmentsAfterRasterization.resize(4096);


};




#if 0
void clip_with_plane(const Plane<float>& c_plane, std::vector<TriangleWithAttributes>& in_list,
	std::vector<TriangleWithAttributes>& outer_list) {

	for (auto& triangle : in_list) {
		int in_vert_num = 0;
		int previous_index, current_index;
		mVec3f thisTriangleVertices[3] = {
			triangle.vtxa.EyeSpaceCoordinate,
			triangle.vtxb.EyeSpaceCoordinate,
			triangle.vtxc.EyeSpaceCoordinate
		};
		mVec4f remainPoints[4];
		int remain_vert_num = 0;
		for (int i = 0; i < 3; i++)
		{
			current_index = i;
			previous_index = (i - 1 + 3) % 3;
			mVec3f pre_vertex = thisTriangleVertices[previous_index]; //
			mVec3f cur_vertex = thisTriangleVertices[current_index];  //

			float d1 = c_plane.cal_project_distance(pre_vertex);
			float d2 = c_plane.cal_project_distance(cur_vertex);

			if (d1 * d2 < 0)
			{
				float t = c_plane.cal_intersectRatio(pre_vertex, cur_vertex);
				mVec3f intersectionPostion = cur_vertex + pre_vertex * (1 - t);
				remainPoints[remain_vert_num++] = { intersectionPostion ,t };

			}
			if (d2 < 0)
			{
				remainPoints[remain_vert_num++] = { cur_vertex ,1.0f };
			}
		}

		//if (remain_vert_num == 3) {
		//	TriangleWithAttributes
		//	outer_list.emplace_back(Vertex4TriangleWithAttributes(remainPoints[0], remainPoints[1], remainPoints[2]));

		//}
		//else	if (remain_vert_num == 4) {
		//	outer_list.emplace_back(Vertex4TriangleWithAttributes(remainPoints[0], remainPoints[1], remainPoints[2]));
		//	outer_list.emplace_back(Vertex4TriangleWithAttributes(remainPoints[1], remainPoints[2], remainPoints[3]));
		//}
		//else if (remain_vert_num != 0) {
		//	std::cout << "E: Unknown Error!";
		//	return;
		//}



	}


	return;
}

#endif


