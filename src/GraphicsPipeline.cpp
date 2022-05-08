#include <execution>
#include <mutex>
#include "matrix.h"
#include "assert.h"
#include "GraphicsPipeline.h"
#include "timer.hpp"
#include <numeric>


#define LI 100
using namespace cimg_library;

//#TODO: 
// //0.0 depth buffer has problem when there are multiple objects!!!!!
//0. use buffer remove one redundant perspective divide in vertex process 
//1.resize windows call_back;       done
//2.show fps on windows title?          done
// 3. marching cube for procedural land generation.      done
// 3.2 debug why scoll mouse change texture when use perspective divide?
// 3.3 alpha blend
//4.parse args to start ray-tracing or graphics pipeline
// 4.1.Clip done
// 4.2 homospace clip
//5.sky box
//6.pbr ggx 
//7.water surface


mVec3f InitLightPos = mVec3f(0, 10, 0);
float texture_scaling = 100;

static mVec3f LightSourceEyeSpace ;
static mVec3f InitEyePos = mVec3f(0, 10, 20);
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
	 Vertices.clear();
	 Vertices.resize(vtxs.size());
	 VtxNormals.clear();
	 VtxNormals.resize(vtxs.size());
	 std::copy(vtxs.begin(), vtxs.end(), Vertices.begin());
	 std::copy(normals.begin(), normals.end(), VtxNormals.begin());

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

	const int N_VERTEX = yu.Vertices.size();
	const int N_VERTEXNORMAL = yu.VtxNormals.size();
	const int N_VERTEXUV= yu.VtxTexUV.size();
	assert(N_VERTEX==N_VERTEXNORMAL&& N_VERTEX == N_VERTEXUV);
	/*auto& vertex_coordinate_buffer = Buffers[0];
	auto& vertex_normal_buffer = Buffers[1];
	auto& vertex_shading_buffer = Buffers[3];
	auto& vertex_EyeSpacePos_buffer = Buffers[4];
	vertex_coordinate_buffer.reserve(N_VERTEX);
	vertex_normal_buffer.reserve(N_VERTEX);
	vertex_shading_buffer.reserve(N_VERTEX);
	vertex_EyeSpacePos_buffer.reserve(N_VERTEX);*/
	//-------------------------------
	//**** Per-vertex Operation  ****
	//Matrix4 P = squashMatrix(_Camera.getNearPlane(), _Camera.getFarPlane());
	//Matrix4 Ortho = OrthoMatrix(_Camera.r, _Camera.l, _Camera.t, _Camera.b, _Camera.getNearPlane(), _Camera.getFarPlane());
	//Per=Ortho*P
	//-------------------------------

	Matrix4 M = yu.getModelMatrix();
	Matrix4 V = _Camera.genViewMat();

	Matrix4 ViewModel = ((V * M));
	Matrix4 NoramlViewModel = inv(ViewModel).Transpose();

	mVec4f LightSourceHomo = V * mVec4f(InitLightPos, 1);
	LightSourceEyeSpace = (LightSourceHomo/ LightSourceHomo.w).tomVec3();




	std::vector<vertex>TargetRenderVertexes;
	TargetRenderVertexes.reserve(N_VERTEX);
	int Vid = 0;

	// apply  MV transformation +Pervertex lighting

	for (const mVec3f& Vertex : yu.Vertices) {
		//-------------------------------
		//  Model + View  transformation
		//-------------------------------
		mVec4f EyeSpaceHomoCoordinates = (ViewModel * mVec4f(Vertex, 1));
		mVec3f EyeSpaceCoordinates = EyeSpaceHomoCoordinates.HomoCordinates2InHomoVec3();

		mVec3f vertex_normal = (NoramlViewModel * mVec4f(yu.VtxNormals[Vid], 0)).tomVec3();
		//-------------------------------
		//Pervertex  LightSource
		//-------------------------------
		mVec3f shading = ADSshading({0,0,0},LightSourceEyeSpace, EyeSpaceCoordinates, vertex_normal, yu.Material());
	

		vertex tmp = { EyeSpaceHomoCoordinates.tomVec3(),vertex_normal,shading, yu.VtxTexUV[Vid] };
		TargetRenderVertexes.push_back(std::move(tmp));
		/*vertex_coordinate_buffer.emplace_back(ScreenSpaceInHomo);
		vertex_normal_buffer.emplace_back(vertex_normal);
		vertex_shading_buffer.emplace_back(shading);
		vertex_EyeSpacePos_buffer.emplace_back(EyeSpaceCoordinates);
		Vid++;*/

	}



	//-------------------------------
	//primitives assembly
	//-------------------------------
	
	TargetRenderTriangles.reserve(yu.TrianglesIdx.size());
	for (auto idxs : yu.TrianglesIdx) {
			TargetRenderTriangles.emplace_back(TargetRenderVertexes[idxs.x], TargetRenderVertexes[idxs.y], TargetRenderVertexes[idxs.z]);
	}
	//# TODO: Clip
	float alphaOver2 = Radians(_Camera.fov * 0.5);
	float n = _Camera.getNearPlane();
	float f = _Camera.getFarPlane();
	Plane<float>Top = { mVec3f{ 0,0,0 }, mVec3f{0,std::cos(alphaOver2),std::sin(alphaOver2)}};
	Plane<float>Bottom = { mVec3f{ 0,0,0 },mVec3f {0,-std::cos(alphaOver2),std::sin(alphaOver2)} };
	Plane<float>Left= { mVec3f { 0,0,0 }, mVec3f{-std::cos(alphaOver2),0,std::sin(alphaOver2)} };
	Plane<float>Right= { mVec3f{ 0,0,0 }, mVec3f{std::cos(alphaOver2),0,std::sin(alphaOver2)} };
	Plane<float>Near= { mVec3f { 0,0,n }, mVec3f{0,0,1} };
	Plane<float>Far= { mVec3f{ 0,0,f }, mVec3f{0,0,-1} };
	
	clip_with_plane(Near, TargetRenderTriangles, TargetRenderTriangles);
	clip_with_plane(Top, TargetRenderTriangles, TargetRenderTriangles);
	clip_with_plane(Bottom, TargetRenderTriangles, TargetRenderTriangles);
	clip_with_plane(Left, TargetRenderTriangles, TargetRenderTriangles);
	clip_with_plane(Right, TargetRenderTriangles, TargetRenderTriangles);
	clip_with_plane(Far, TargetRenderTriangles, TargetRenderTriangles);

}


template<typename T>
static inline T PerspetiveCorrectDepth(T iaz, T ibz, T icz, T uv_u, T uv_v) {

	float interpolateIZ = iaz + uv_u * (icz - iaz) + uv_v * (ibz - iaz);//P = A + u * (C - A) + v * (B - A)  
	float  interpolateZ = (1.0f / interpolateIZ);
	return interpolateZ;
}




void GraphicsPipeline::Rasterization( ShadingMaterial& myu, framebuffer_t* Fb, const unsigned int texChannel) const
{
	auto clamp = [](float in)->float {
		return in > 1.0f ? 1.0f : (in < 0.0f ? 0.0f : in);
	};
	auto fract = [](float x)->float {
		return x - floor(x);
	};
	const Matrix4 Per = _Camera.genPerspectiveMat();
	const Matrix4 Vp = ViewPortMatrix(window_width, window_height);
	const Matrix4 PersViewPort = (Vp * Per);


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
			mVec3f point = Tri.vtxa.Coordinate * (1 - uv.x - uv.y) + (Tri.vtxc.Coordinate * uv.x) + (Tri.vtxb.Coordinate * uv.y);
			tmp = ADSshading({ 0,0,0 }, LightSourceEyeSpace,point, normal, myu);

			//stripe texture
			if (texChannel == 1) {
				mVec3f backcolor = { 0,0,0 };
				float scale = texture_scaling;
				float fuzz = 2;
				float width = 10;
				mVec2<float> Texcoords = Tri.vtxa.st + ((Tri.vtxc.st - Tri.vtxa.st) * uv.x + (Tri.vtxb.st - Tri.vtxa.st) * uv.y);
				Texcoords = Texcoords * interpolateZ;
				float scaleT = fract(Texcoords.y * scale);

				float frac1 = clamp(scaleT / fuzz);
				float frac2 = clamp((scaleT - width) / fuzz);

				frac1 = frac1 * (1.0f - frac2);
				frac1 = frac1 * frac1 * (3.0f - (2.0f * frac1));
				tmp = backcolor * (frac1)+tmp * (1 - frac1);

			}
			else if (texChannel >=2) {
				float scale = texture_scaling;
				mVec2<float> Texcoords = Tri.vtxa.st + ((Tri.vtxc.st - Tri.vtxa.st) * uv.x + (Tri.vtxb.st - Tri.vtxa.st) * uv.y);
				Texcoords = Texcoords* interpolateZ;
				Texcoords.x = fract(Texcoords.x * scale);
				Texcoords.y = fract(Texcoords.y * scale);
				mVec3f backcolor =LookupTexel(Texcoords, texChannel-2);
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
			if ( depthbuffer[ScreenXY.x * window_width + ScreenXY.y] > TriFrag.depth) {
				unsigned char fragColor[4] = { static_cast<unsigned char>(TriFrag.RGBv.x * 255),static_cast<unsigned char>(TriFrag.RGBv.y * 255),static_cast<unsigned char>(TriFrag.RGBv.z * 255),0 };
				Fb->setvalue(ScreenXY.x, ScreenXY.y, fragColor);
				//Fb->color_buffer[ScreenXY.y*w + ScreenXY.x] = TriFrag.RGB;
				depthbuffer[ScreenXY.x * window_width + ScreenXY.y] = TriFrag.depth;
			}

		}
	};


	auto RasterAtriangle = [&](const TriangleWithAttributes& TriWithAtrib) {
		constexpr bool UseEarlyZ = true;


		//TriWithAtrib store eyespace  coordinates
		auto Tri = TriWithAtrib.GetTriangleVertexes();
		//Perpective + Viewport transformation
		//Tri will store screen space coordinates
		Tri.a = (PersViewPort * mVec4f{ Tri.a,1.0 }).HomoCordinates2InHomoVec3();
		Tri.b = (PersViewPort * mVec4f{ Tri.b,1.0 }).HomoCordinates2InHomoVec3();
		Tri.c = (PersViewPort * mVec4f{ Tri.c,1.0 }).HomoCordinates2InHomoVec3();


		int Xmin, Xmax, Ymin, Ymax;
		Xmin = int(floorf(MIN(MIN(Tri.a.x, Tri.b.x), Tri.c.x)));
		Xmax = int(ceilf(MAX(MAX(Tri.a.x, Tri.b.x), Tri.c.x)));
		Ymin = int(floorf(MIN(MIN(Tri.a.y, Tri.b.y), Tri.c.y)));
		Ymax = int(ceilf(MAX(MAX(Tri.a.y, Tri.b.y), Tri.c.y)));
		if (Xmin < 0 || Xmin >= window_width || Xmax < 0 || Xmax >= window_width)return;
		if (Ymin < 0 || Ymin >= window_height || Ymax < 0 || Ymax >= window_height)return;
		if (Tri.a.z > 0 || Tri.b.z > 0 || Tri.c.z > 0)return;


		//early-z
		const float iaz = 1.0f / (Tri.a.z - depthA);
		const float ibz = 1.0f / (Tri.b.z - depthA);
		const float icz = 1.0f / (Tri.c.z - depthA);

		mVec2<float>FirstTime;
		mVec2<float>Incremental_uv;
		mVec2<float>Lastuv{};

		for (int x = Xmin; x <= Xmax; x++) {
			for (int y = Ymin; y <= Ymax; y++) {
				if (y == Ymin) {
					mVec2<float>  uv = Tri.PointIsInTriangle(mVec3f(x, y, 0));
					if (uv.x >= 0 && uv.y >= 0 && uv.x + uv.y <= 1) {
						float interpolatez = PerspetiveCorrectDepth(iaz, ibz, icz, uv.x, uv.y);
						/*if constexpr (UseEarlyZ)
							if (!EarlyZtest(x, y, interpolatez))continue;*/

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
					/*	if constexpr (UseEarlyZ)
							if (!EarlyZtest(x, y, interpolatez))continue;*/

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
						/*if constexpr (UseEarlyZ)
							if (!EarlyZtest(x, y, interpolatez))continue;*/

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

void GraphicsPipeline::Render(const std::vector<RenderableObject>& robjs, framebuffer_t* Fb) {
	float* depthbuffer = Fb->depth_buffer;
	assert(depthbuffer != NULL);
	memset(depthbuffer, 0x7f, DEFAULT_WINDOW_HEIGHT * DEFAULT_WINDOW_WIDTH * sizeof(float));

	for (const auto& r_object : robjs)
	{
		auto material = r_object.Material();
		VertexesProcess(r_object);
		Rasterization(material, Fb, r_object.texturechannelID);
		//FragmentProcess(Fb);
		clearPipeline();

	}

	
	
	
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
	_Camera.SetFrustm(1, -1, -1, 1, 90.0f, InitEyePos.z + 400);

	//pre-allocate memory:
	//FragmentsAfterRasterization.resize(4096);


};


void clip_with_plane(const Plane<float>& c_plane, std::vector<TriangleWithAttributes>& in_list,
	std::vector<TriangleWithAttributes>& outer_list) {
	std::vector<TriangleWithAttributes> tempOuterList;
	tempOuterList.reserve(in_list.size());
	for (auto& triangle : in_list) {
		int in_vert_num = 0;
		int previous_index, current_index;
		vertex thisTriangleVertices[3] = {
			triangle.vtxa,
			triangle.vtxb,
			triangle.vtxc
		};
		vertex remainPoints[4];
		int remain_vert_num = 0;
		for (int i = 0; i < 3; i++)
		{
			current_index = i;
			previous_index = (i - 1 + 3) % 3;
			mVec3f pre_vertex = thisTriangleVertices[previous_index].Coordinate; //
			mVec3f cur_vertex = thisTriangleVertices[current_index].Coordinate;  //

			float d1 = c_plane.cal_project_distance(pre_vertex);
			float d2 = c_plane.cal_project_distance(cur_vertex);

			if (d1 * d2 < 0)
			{
				float t = c_plane.cal_intersectRatio(pre_vertex, cur_vertex);
				mVec3f interpolateCordinate= cur_vertex + pre_vertex * (1 - t);
				//attributes
				mVec3f interpolatenormal = thisTriangleVertices[current_index].normal + thisTriangleVertices[previous_index].normal * (1 - t);
				mVec3f interpolateshading= thisTriangleVertices[current_index].shading + thisTriangleVertices[previous_index].shading * (1 - t);
				mVec2f interpolateUV= thisTriangleVertices[current_index].st + thisTriangleVertices[previous_index].st * (1 - t);
				mVec2<float> st;
				remainPoints[remain_vert_num++] = { interpolateCordinate,interpolatenormal,interpolateshading,interpolateUV  };

			}
			if (d2 < 0)
			{
				remainPoints[remain_vert_num++] = { cur_vertex ,1.0f };
			}
		}

		if (remain_vert_num == 3) {
			tempOuterList.emplace_back(triangle);
		}
		else	if (remain_vert_num == 4) {
			TriangleWithAttributes newTriangle1 = { remainPoints[0],remainPoints[1],remainPoints[2] };
			TriangleWithAttributes newTriangle2 = { remainPoints[2],remainPoints[1],remainPoints[3] };
			tempOuterList.emplace_back(std::move(newTriangle1));
			tempOuterList.emplace_back(std::move(newTriangle2));

		}
		else if (remain_vert_num != 0) {
			std::cout << "E: Unknown Error!";
			return;
		}



	}

	std::swap(tempOuterList, outer_list);
	return;
}
