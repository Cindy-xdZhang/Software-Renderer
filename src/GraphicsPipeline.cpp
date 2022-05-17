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
//1.resize windows call_back;       done
//2.show fps on windows title?          done
// 3. marching cube for procedural land generation.      done
// 3.2 debug why scroll mouse change texture when use perspective divide?
// 3.3 alpha blend 
//4.parse args to start ray-tracing or graphics pipeline
// 4.1.Clip 
// 4.2 clip cause zigzag?
// 4.3 frustum culling while do frustum  clipping
//5.sky box
//6.pbr ggx 
//7.water surface
//10. use buffer remove one redundant perspective divide in vertex process 


mVec3f InitLightPos = mVec3f(0, 10, 0);
float texture_scaling = 100;

static mVec3f LightSourceEyeSpace;
static mVec3f InitEyePos = mVec3f(0, 2, 20);
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

STRONG_INLINE  mVec3f ADSshading(const mVec3f& eye, const mVec3f& aLightSource, const mVec3f& point, const mVec3f& normal, const ShadingMaterial& Mp) {
	mVec3f Ia = Mp.Ka * LI;
	mVec3f  Sourse = aLightSource;
	//Reflection
	mVec3f ReflectionD = mVec3f(Sourse, point);
	//float distance2 = pow((Sourse.x - point.x), 2) + pow((Sourse.y - point.y), 2) + pow((Sourse.z - point.z), 2);
	float distance2 =(Sourse-point).getEuclideannNorms();
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




void RenderableObject::BuildLikGLBegin(const std::vector<mVec3f>& vtxs, const std::vector<mVec3f>& normals, enum MeshBuildConvention Conv) {
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
	for (int i = 0; i < nTriangles; i += 3)
	{
		TrianglesIdx.emplace_back(i, i + 1, i + 2);
	}
	VtxTexUV = OBjReader::AssignTextureCoordinates(vtxs);

}





void GraphicsPipeline::LoadTexture(std::filesystem::path tex_path) {
	LoadTexture(tex_path.string());

}


void  GraphicsPipeline::LoadTexture(std::string stringPath) {

	auto iter = stringPath.find(".jpg");
	auto iter2 = stringPath.find(".bmp");
	if (iter != std::string::npos)
	{
		CImg<unsigned char > TextureImg(stringPath.c_str()), visu();
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
	const int N_VERTEXUV = yu.VtxTexUV.size();
	const int N_FACES = yu.TrianglesIdx.size();
	assert(N_VERTEX == N_VERTEXNORMAL && N_VERTEX == N_VERTEXUV);

	

	auto& face_idx_buffer = vec3iBuffer;
	auto& vertex_texture_buffer = vec2Buffer;
	auto& vertex_coordinate_buffer = vec4Buffers[0];
	auto& vertex_normal_buffer = vec3Buffers[0];
	auto& vertex_shading_buffer = vec3Buffers[1];

	vertex_normal_buffer.reserve(N_VERTEX);
	vertex_coordinate_buffer.reserve(N_VERTEX);
	vertex_shading_buffer.reserve(N_VERTEX);
	face_idx_buffer.resize(N_VERTEX);
	vertex_texture_buffer.resize(N_VERTEX);
	std::copy(yu.TrianglesIdx.begin(), yu.TrianglesIdx.end(), face_idx_buffer.begin());
	std::copy(yu.VtxTexUV.begin(), yu.VtxTexUV.end(), vertex_texture_buffer.begin());
	
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
	LightSourceEyeSpace = (LightSourceHomo / LightSourceHomo.w).tomVec3();


	

	int Vid = 0;

	// apply  MV transformation +Pervertex lighting

	for (const mVec3f& Vertex : yu.Vertices) {
		//-------------------------------
		//  Model + View  transformation
		//-------------------------------
		mVec4f EyeSpaceHomoCoordinates = (ViewModel * mVec4f(Vertex, 1));
		mVec3f EyeSpaceCoordinates = EyeSpaceHomoCoordinates.HomoCordinates2InHomoVec3();

		mVec3f vertex_normal = (NoramlViewModel * mVec4f(yu.VtxNormals[Vid], 0)).tomVec3();
		vertex_normal.normalize();
		//-------------------------------
		//Pervertex  LightSource
		//-------------------------------
		mVec3f PerVertexShading = ADSshading({ 0,0,0 }, LightSourceEyeSpace, EyeSpaceCoordinates, vertex_normal, yu.Material());


		vertex_coordinate_buffer.push_back(std::move(EyeSpaceHomoCoordinates));
		vertex_normal_buffer.push_back(std::move(vertex_normal));
		vertex_shading_buffer.push_back(std::move(PerVertexShading));
	}


	//-------------------------------
	//primitives assembly
	//-------------------------------

	//Clip
	if(Clip) {
		float alphaOver2 = Radians(_Camera.fov * 0.5);
		float n = _Camera.getNearPlane();
		float f = _Camera.getFarPlane();
		Plane<float>Top = { mVec3f{ 0,0,0 }, mVec3f{0,std::cos(alphaOver2),std::sin(alphaOver2)} };
		Plane<float>Bottom = { mVec3f{ 0,0,0 },mVec3f {0,-std::cos(alphaOver2),std::sin(alphaOver2)} };
		Plane<float>Left = { mVec3f { 0,0,0 }, mVec3f{-std::cos(alphaOver2),0,std::sin(alphaOver2)} };
		Plane<float>Right = { mVec3f{ 0,0,0 }, mVec3f{std::cos(alphaOver2),0,std::sin(alphaOver2)} };
		Plane<float>Near = { mVec3f { 0,0,n }, mVec3f{0,0,1} };
		Plane<float>Far = { mVec3f{ 0,0,f }, mVec3f{0,0,-1} };

		clip_with_plane(Near);
		clip_with_plane(Top);
		clip_with_plane(Bottom);
		clip_with_plane(Left);
		clip_with_plane(Right);
		clip_with_plane(Far);
	}


}


template<typename T>
static inline T PerspetiveCorrectDepth(T iaz, T ibz, T icz, T uv_u, T uv_v) {

	float interpolateIZ = iaz + uv_u * (icz - iaz) + uv_v * (ibz - iaz);//P = A + u * (C - A) + v * (B - A)  
	float  interpolateZ = (1.0f / interpolateIZ);
	return interpolateZ;
}




void GraphicsPipeline::Rasterization(ShadingMaterial& myu, framebuffer_t* Fb, const unsigned int texChannel) const
{
	auto clamp = [](float in)->float {
		return in > 1.0f ? 1.0f : (in < 0.0f ? 0.0f : in);
	};
	auto fract = [](float x)->float {
		return x - floor(x);
	};
	auto& face_idx_buffer = vec3iBuffer;
	auto& vertex_texture_buffer = vec2Buffer;
	auto& vertex_coordinate_buffer = vec4Buffers[0];
	auto& vertex_normal_buffer = vec3Buffers[0];
	auto& vertex_shading_buffer = vec3Buffers[1];


	const Matrix4 Per = _Camera.genPerspectiveMat();
	const Matrix4 Vp = ViewPortMatrix(window_width, window_height);
	const Matrix4 PersViewPort = (Vp * Per);


	//after perspective ortho viewport z=n+f/n-f - 2*f*n/z(n-f)
	float n = _Camera.getNearPlane();
	float f = _Camera.getFarPlane();
	float depthA = (n + f) / (n - f);
	float zcoord = f * n * 2;




	auto FramentShader = [&](const Triangle& TriInEyeSpace,const mVec3i TriIdx, const mVec2<float>&uv, float interpolateZ)->mVec4f {

		mVec3f tmp;
		auto normal_a = vertex_normal_buffer[TriIdx.x];
		auto normal_b = vertex_normal_buffer[TriIdx.y];
		auto normal_c = vertex_normal_buffer[TriIdx.z];

		auto vtex_a= vertex_texture_buffer[TriIdx.x];
		auto vtex_b = vertex_texture_buffer[TriIdx.y];
		auto vtex_c = vertex_texture_buffer[TriIdx.z];

		auto perVshading_a = vertex_shading_buffer[TriIdx.x];
		auto perVshading_b = vertex_shading_buffer[TriIdx.y];
		auto perVshading_c = vertex_shading_buffer[TriIdx.z];

		if (UsePhongShading) {
			//interpolate normal
			mVec3f normal = normal_a+ ((normal_c- normal_a) * uv.x + (normal_b- normal_a) * uv.y);

			//already normalize it in vertex shader, no need re-normalize
			//normal.normalize();//normal in eyespace,light source in eyespace

			mVec3f point = TriInEyeSpace.a* (1 - uv.x - uv.y) + (TriInEyeSpace.c* uv.x) + (TriInEyeSpace.b* uv.y);
			tmp = ADSshading({ 0,0,0 }, LightSourceEyeSpace, point, normal, myu);

			mVec2<float> Texcoords = vtex_a + ((vtex_c - vtex_a) * uv.x + (vtex_b - vtex_a) * uv.y);
			//stripe texture
			if (texChannel == 1) {
				mVec3f backcolor = { 0,0,0 };
				float scale = texture_scaling;
				float fuzz = 2;
				float width = 10;
				Texcoords = Texcoords * interpolateZ;
				float scaleT = fract(Texcoords.y * scale);

				float frac1 = clamp(scaleT / fuzz);
				float frac2 = clamp((scaleT - width) / fuzz);

				frac1 = frac1 * (1.0f - frac2);
				frac1 = frac1 * frac1 * (3.0f - (2.0f * frac1));
				tmp = backcolor * (frac1)+tmp * (1 - frac1);

			}
			else if (texChannel >= 2) {
				float scale = texture_scaling;
				Texcoords = Texcoords * interpolateZ;
				Texcoords.x = fract(Texcoords.x * scale);
				Texcoords.y = fract(Texcoords.y * scale);
				mVec3f backcolor = LookupTexel(Texcoords, texChannel - 2);
				tmp = { backcolor.x * tmp.x,  backcolor.y * tmp.y, backcolor.z * tmp.z };
				tmp.ColorClamp();
			}



		}
		else {
			//interpolate lightning
			tmp = perVshading_a* (1 - uv.x - uv.y) + (perVshading_c * uv.x) + (perVshading_b* uv.y);


		}
		return { tmp,1.0 };//opacity not used for now.
	};



	int TrianglesNumber = vec3iBuffer.size();

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
		//mVec2<int> ScreenXY = { TriFrag.XY.y ,TriFrag.XY.x };

		{
			//std::lock_guard<std::mutex> lock(buffer_mutex);
			//look at -z direction, this depth is z in screen space 
			if (depthbuffer[TriFrag.XY.y * window_width + TriFrag.XY.x] <TriFrag.depth) {
				unsigned char fragColor[4] = { static_cast<unsigned char>(TriFrag.RGBv.x * 255),
					static_cast<unsigned char>(TriFrag.RGBv.y * 255),static_cast<unsigned char>(TriFrag.RGBv.z * 255),0 };
				Fb->setvalue(TriFrag.XY.y, TriFrag.XY.x, fragColor);
				depthbuffer[TriFrag.XY.y * window_width + TriFrag.XY.x] = TriFrag.depth;
			}

		}
	};


	auto RasterAtriangle = [&](const mVec3i& TriIndex) {
		constexpr bool UseEarlyZ = false;
		
		//buffer store eyespace  coordinates
		//Perpective + Viewport transformation
		auto Eye_tri_a = vertex_coordinate_buffer[TriIndex.x];
		auto Eye_tri_b = vertex_coordinate_buffer[TriIndex.y];
		auto Eye_tri_c = vertex_coordinate_buffer[TriIndex.z];

		auto Tri_a = (PersViewPort * Eye_tri_a).HomoCordinates2InHomoVec3();
		auto Tri_b = (PersViewPort * Eye_tri_b).HomoCordinates2InHomoVec3();
		auto Tri_c = (PersViewPort * Eye_tri_c).HomoCordinates2InHomoVec3();
		
		Triangle TriInEyeSpace= { Eye_tri_a.tomVec3(),Eye_tri_b.tomVec3(),Eye_tri_c.tomVec3()};
		Triangle  Tri = {Tri_a,Tri_b,Tri_c};


	
		int Xmin, Xmax, Ymin, Ymax;
		Xmin = int(floorf(MIN(MIN(Tri.a.x, Tri.b.x), Tri.c.x)));
		Xmax = int(ceilf(MAX(MAX(Tri.a.x, Tri.b.x), Tri.c.x)));
		Ymin = int(floorf(MIN(MIN(Tri.a.y, Tri.b.y), Tri.c.y)));
		Ymax = int(ceilf(MAX(MAX(Tri.a.y, Tri.b.y), Tri.c.y)));
		if ( Xmin > window_width ||Xmax <0 )return;
		if (Ymin > window_height || Ymax < 0)return;
	/*	if (Xmin < 0 || Xmin >= window_width || Xmax < 0 || Xmax >= window_width)return;
		if (Ymin < 0 || Ymin >= window_height || Ymax < 0 || Ymax >= window_height)return;*/
		Xmin = Xmin < 0 ? 0 : Xmin;
		Ymin = Ymin < 0 ? 0 : Ymin;
		Xmax = Xmax >= window_width ? window_width-1 : Xmax;
		Ymax = Ymax >= window_height ? window_height-1 : Ymax;

		//early-z
		/*const float iaz = 1.0f / (Tri.a.z - depthA);
		const float ibz = 1.0f / (Tri.b.z - depthA);
		const float icz = 1.0f / (Tri.c.z - depthA);*/

		mVec2<float>FirstTime;
		mVec2<float>Incremental_uv;
		mVec2<float>Lastuv{};

		for (int x = Xmin; x <= Xmax; x++) {
			for (int y = Ymin; y <= Ymax; y++) {
				if (y == Ymin) {
					mVec2<float>  uv = Tri.PointIsInTriangle(mVec3f(x, y, 0));
					if (uv.x >= 0 && uv.y >= 0 && uv.x + uv.y <= 1) {
						//float interpolatez = PerspetiveCorrectDepth(iaz, ibz, icz, uv.x, uv.y);
						float interpolatez = Tri.a.z + uv.x * (Tri.c.z - Tri.b.z ) + uv.y* (Tri.b.z - Tri.a.z);//P = A + u * (C - A) + v * (B - A)  ;
						if constexpr (UseEarlyZ)
							if (!EarlyZtest(x, y, interpolatez))continue;

						mVec4f attibs = FramentShader(TriInEyeSpace, TriIndex, uv, interpolatez);
						//RasterResult.emplace_back(attibs.w, mVec2<int>{x, y }, attibs.tomVec3f());
						outputshading({ interpolatez, mVec2<int>{x, y }, attibs.tomVec3() });

					}
					FirstTime = uv;


				}
				else if (y == Ymin + 1) {
					mVec2<float>  uv = Tri.PointIsInTriangle(mVec3f(x, y, 0));
					if (uv.x >= 0 && uv.y >= 0 && uv.x + uv.y <= 1) {

						float interpolatez = Tri.a.z + uv.x * (Tri.c.z - Tri.b.z) + uv.y * (Tri.b.z - Tri.a.z);//P = A + u * (C - A) + v * (B - A)  ;
							if constexpr (UseEarlyZ)
								if (!EarlyZtest(x, y, interpolatez))continue;

						mVec4f attibs = FramentShader(TriInEyeSpace, TriIndex, uv, interpolatez);
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

						float interpolatez = Tri.a.z + uv.x * (Tri.c.z - Tri.b.z) + uv.y * (Tri.b.z - Tri.a.z);//P = A + u * (C - A) + v * (B - A)  ;
						if constexpr (UseEarlyZ)
							if (!EarlyZtest(x, y, interpolatez))continue;

						mVec4f attibs = FramentShader(TriInEyeSpace, TriIndex, uv, interpolatez);
						//RasterResult.emplace_back(attibs.w, mVec2<int>{x, y }, attibs.tomVec3f());
						outputshading({ interpolatez, mVec2<int>{x, y }, attibs.tomVec3() });
					}
				}

			}
		}

		//return 0;




	};



#ifdef SingleThread
	auto policy = std::execution::seq;
#else
	auto policy = std::execution::par_unseq;
#endif
	std::for_each(policy, vec3iBuffer.begin(), vec3iBuffer.end(), RasterAtriangle);


	//interpolate
}

void GraphicsPipeline::Render(const std::vector<RenderableObject>& robjs, framebuffer_t* Fb) {
	float* depthbuffer = Fb->depth_buffer;
	assert(depthbuffer != NULL);
	memset(depthbuffer, 0xfe, static_cast<size_t>(window_height) * static_cast<size_t>(window_width) * sizeof(float));

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
	
	for (auto& vbuffer : vec4Buffers) {
		vbuffer.clear();
	}
	for (auto& vbuffer : vec3Buffers) {
		vbuffer.clear();
	}
	vec3iBuffer.clear();
	vec2Buffer.clear();
}

GraphicsPipeline::GraphicsPipeline(int h, int w) : window_height(h), window_width(w) {
	//mTextures.emplace_back(LoadTexture());

	_Camera = Camera(InitEyePos, InitGazeDirection, InitTopDirection);
	_Camera.SetFrustm(1, -1, -1, 1, 45.0f, InitEyePos.z + 200);

	//pre-allocate memory:
	//FragmentsAfterRasterization.resize(4096);


};




void  GraphicsPipeline::clip_with_plane(const Plane<float>& c_plane) {

	auto& face_idx_buffer = vec3iBuffer;
	auto& vertex_texture_buffer = vec2Buffer;
	auto& vertex_coordinate_buffer = vec4Buffers[0];
	auto& vertex_normal_buffer = vec3Buffers[0];
	auto& vertex_shading_buffer = vec3Buffers[1];


	std::vector<mVec3i> tempOuterList;
	tempOuterList.reserve(face_idx_buffer.size());

	for (auto& triangle : face_idx_buffer) {
		int in_vert_num = 0;
		int previous_index, current_index;
		int thisTriangleVertices[3] = {
			triangle.x,
			triangle.y,
			triangle.z
		};

		
		int remain_vert_num = 0;
		int remain_vert_id[4] = { -1,-1,-1,-1 };

		for (int i = 0; i < 3; i++)
		{
			current_index = i;
			previous_index = (i - 1 + 3) % 3;
			const int pre_v_id = thisTriangleVertices[previous_index];
			const int cur_v_id = thisTriangleVertices[current_index];

			mVec3f pre_vertex = vertex_coordinate_buffer[pre_v_id].tomVec3(); //
			mVec3f cur_vertex = vertex_coordinate_buffer[cur_v_id].tomVec3();  //

			float d1 = c_plane.cal_project_distance(pre_vertex);
			float d2 = c_plane.cal_project_distance(cur_vertex);

			if (d1 * d2 < 0)
			{
				float t = d1 / (d1 - d2);
				mVec3f interpolateCordinate = cur_vertex * t + pre_vertex * (1 - t);
				//attribute
				mVec3f interpolatenormal = vertex_normal_buffer[cur_v_id] * t + vertex_normal_buffer[pre_v_id] * (1 - t);
				mVec3f interpolateshading = vertex_shading_buffer[cur_v_id] * t + vertex_shading_buffer[pre_v_id] * (1 - t);

				mVec2f interpolateUV = vertex_texture_buffer[cur_v_id]* t + vertex_texture_buffer[pre_v_id]* (1 - t);
		
				int lastVid = vertex_coordinate_buffer.size() ;
				vertex_coordinate_buffer.emplace_back(interpolateCordinate,1.0);
				vertex_normal_buffer.emplace_back(interpolatenormal);
				vertex_shading_buffer.emplace_back(interpolateshading);
				vertex_texture_buffer.emplace_back(interpolateUV);
				remain_vert_id[remain_vert_num++] = lastVid;
			}
			if (d2 <= 0)
			{
				remain_vert_id[remain_vert_num++] = cur_v_id;
			}


		}

		
		if (remain_vert_num == 3) {
			tempOuterList.emplace_back(remain_vert_id[0], remain_vert_id[1], remain_vert_id[2]);
		}
		else    if (remain_vert_num == 4) {
			tempOuterList.emplace_back(remain_vert_id[0], remain_vert_id[1], remain_vert_id[2]);
			tempOuterList.emplace_back(remain_vert_id[0], remain_vert_id[2], remain_vert_id[3]);

		}
		else if (remain_vert_num != 0) {
			std::cout << "E: Unknown Error!";
			return;
		}

	}

	std::swap(tempOuterList, face_idx_buffer);
	return;
}
