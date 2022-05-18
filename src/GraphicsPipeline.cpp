#include <assert.h>
#include <execution>
#include <mutex>
#include <numeric>
#include <immintrin.h>
#include "matrix.h"
#include "GraphicsPipeline.h"
#include "timer.hpp"

#define LI 100
using namespace cimg_library;

//#TODO: 
//1.resize windows call_back;       done
//2.show fps on windows title?          done
// 3. marching cube for procedural land generation.      done
// 3.2 debug why scroll mouse change texture when use perspective divide?
// 3.3 alpha blend 
//4.parse args to start ray-tracing or graphics pipeline
// 4.1.Clip  done
// 4.3 frustum culling while do frustum  clipping done
//5.sky box
//6.pbr ggx 
//7.water surface
//10. use buffer remove one redundant perspective divide in vertex process 


mVec3f InitLightPos = mVec3f(0, 10, 0);


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
	float distance2 = (Sourse - point).getEuclideannNorms();
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
	face_idx_buffer.resize(yu.TrianglesIdx.size());
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
	if (Clip) {
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


void GraphicsPipeline::Render(const std::vector<RenderableObject>& robjs1, const std::vector<RenderableObject>& robjs2, framebuffer_t* Fb) {
	float* depthbuffer = Fb->depth_buffer;
	assert(depthbuffer != NULL);
	memset(depthbuffer, 0xfe, static_cast<size_t>(window_height) * static_cast<size_t>(window_width) * sizeof(float));

	for (const auto& r_object : robjs1)
	{
		auto material = r_object.Material();
		VertexesProcess(r_object);
		Rasterization(r_object, Fb);
		//FragmentProcess(Fb);
		clearPipeline();

	}
	for (const auto& r_object : robjs2)
	{
		auto material = r_object.Material();
		VertexesProcess(r_object);
		Rasterization(r_object, Fb);
		//FragmentProcess(Fb);
		clearPipeline();

	}

}


//PerspetiveCorrect
template<typename T>
static T interpolate_varyings(
	const T& src0, const T& src1, const T& src2, mVec3f weights, mVec3f recip_w) {

	float weight0 = recip_w.x * weights.x;
	float weight1 = recip_w.y * weights.y;
	float weight2 = recip_w.z * weights.z;
	float normalizer = 1 / (weight0 + weight1 + weight2);
	
	auto sum = src0 * weight0 + src1 * weight1 + src2 * weight2;
	T dst = sum * normalizer;
	return dst;
}




void GraphicsPipeline::Rasterization(const RenderableObject& robj,framebuffer_t* Fb) const
{
	auto myu = robj.Material();

	auto& face_idx_buffer = vec3iBuffer;
	auto& vertex_texture_buffer = vec2Buffer;
	auto& vertex_coordinate_buffer = vec4Buffers[0];
	auto& vertex_normal_buffer = vec3Buffers[0];
	auto& vertex_shading_buffer = vec3Buffers[1];


	const Matrix4 Per = _Camera.genPerspectiveMat();
	const Matrix4 Vp = ViewPortMatrix(window_width, window_height);
	const Matrix4 PersViewPort = (Vp * Per);


	//after perspective ortho viewport z'=n+f/n-f +  2fn/z(f-n)
	float n = _Camera.getNearPlane();
	float f = _Camera.getFarPlane();
	float depthAffine = (n + f) / (n - f);
	float ZCoeffcient = (f - n) /( f * n * 2);





	int TrianglesNumber = vec3iBuffer.size();

	float* depthbuffer = Fb->depth_buffer;
	std::mutex buffer_mutex;
	//Rasterization

	auto ZBufferTest = [&](const fragment& TriFrag) {
		{
			//std::lock_guard<std::mutex> lock(buffer_mutex);
			//look at -z direction, this depth is z in screen space 
			if (depthbuffer[TriFrag.XY.y * window_width + TriFrag.XY.x] < TriFrag.depth) {
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

		//screen space homo coordinate
		auto Homo_tri_a = PersViewPort * Eye_tri_a;
		auto Homo_tri_b = PersViewPort * Eye_tri_b;
		auto Homo_tri_c = PersViewPort * Eye_tri_c;

		auto inv_w_a = 1 / Homo_tri_a.w;
		auto inv_w_b = 1 / Homo_tri_b.w;
		auto inv_w_c = 1 / Homo_tri_c.w;
		mVec3f recip_w = {inv_w_a,inv_w_b,inv_w_c};
		//screen space in homo coordinate
		auto Tri_a = (Homo_tri_a * inv_w_a).tomVec3();
		auto Tri_b = (Homo_tri_b * inv_w_b).tomVec3();
		auto Tri_c = (Homo_tri_c * inv_w_c).tomVec3();

		Triangle TriInEyeSpace = { Eye_tri_a.HomoCordinates2InHomoVec3(),Eye_tri_b.HomoCordinates2InHomoVec3(),Eye_tri_c.HomoCordinates2InHomoVec3() };

		//after perspective ortho viewport z'=n+f/n-f +  2fn/z(f-n)
		Triangle  Tri_ScreenSpace = { Tri_a,Tri_b,Tri_c };



		int Xmin, Xmax, Ymin, Ymax;
		Xmin = int(floorf(MIN(MIN(Tri_ScreenSpace.a.x, Tri_ScreenSpace.b.x), Tri_ScreenSpace.c.x)));
		Xmax = int(ceilf(MAX(MAX(Tri_ScreenSpace.a.x, Tri_ScreenSpace.b.x), Tri_ScreenSpace.c.x)));
		Ymin = int(floorf(MIN(MIN(Tri_ScreenSpace.a.y, Tri_ScreenSpace.b.y), Tri_ScreenSpace.c.y)));
		Ymax = int(ceilf(MAX(MAX(Tri_ScreenSpace.a.y, Tri_ScreenSpace.b.y), Tri_ScreenSpace.c.y)));
		if (Xmin > window_width || Xmax < 0)return;
		if (Ymin > window_height || Ymax < 0)return;
	
		Xmin = Xmin < 0 ? 0 : Xmin;
		Ymin = Ymin < 0 ? 0 : Ymin;
		Xmax = Xmax >= window_width ? window_width - 1 : Xmax;
		Ymax = Ymax >= window_height ? window_height - 1 : Ymax;

		//early-z
		const float iaz = Tri_ScreenSpace.a.z;// z'=n+f/n-f +  2fn/z(f-n)
		const float ibz = Tri_ScreenSpace.b.z ;
		const float icz = Tri_ScreenSpace.c.z ;

		mVec2<float>FirstTime;
		mVec2<float>Incremental_uv;
		mVec2<float>Lastuv{};

		for (int x = Xmin; x <= Xmax; x++) {
			for (int y = Ymin; y <= Ymax; y++) {
				if (y == Ymin) {
					mVec2<float>  uv = Tri_ScreenSpace.PointIsInTriangle(mVec3f(x + 0.5, y + 0.5, 0));
					if (uv.x >= 0 && uv.y >= 0 && uv.x + uv.y <= 1) {
						float InterpolatedZ = iaz + uv.x * (ibz - iaz) + uv.y * (icz - iaz);
						/* early depth testing */
						if (InterpolatedZ >= depthbuffer[y * window_width + x]) {
							mVec4f attibs = FragmentShading(TriInEyeSpace, TriIndex, uv, recip_w, robj);
							//RasterResult.emplace_back(attibs.w, mVec2<int>{x, y }, attibs.tomVec3f());
							ZBufferTest({ InterpolatedZ, mVec2<int>{x, y }, attibs.tomVec3() });
						}

					}
					FirstTime = uv;
				}
				else if (y == Ymin + 1) {
					mVec2<float>  uv = Tri_ScreenSpace.PointIsInTriangle(mVec3f(x + 0.5, y + 0.5, 0));
					if (uv.x >= 0 && uv.y >= 0 && uv.x + uv.y <= 1) {
						float InterpolatedZ = iaz + uv.x * (ibz - iaz) + uv.y * (icz - iaz);
						/* early depth testing */
						if (InterpolatedZ >= depthbuffer[y * window_width + x]) {
							mVec4f attibs = FragmentShading(TriInEyeSpace, TriIndex, uv, recip_w, robj);
							//RasterResult.emplace_back(attibs.w, mVec2<int>{x, y }, attibs.tomVec3f());
							ZBufferTest({ InterpolatedZ, mVec2<int>{x, y }, attibs.tomVec3() });
						}

					}
					Incremental_uv = uv - FirstTime;
					Lastuv = uv;
				}
				else {
					mVec2<float> uv = Lastuv + Incremental_uv;
					Lastuv = uv;
					if (uv.x >= 0 && uv.y >= 0 && uv.x + uv.y <= 1) {
						float InterpolatedZ = iaz + uv.x * (ibz - iaz) + uv.y * (icz - iaz);
						/* early depth testing */
						if (InterpolatedZ >= depthbuffer[y * window_width + x]) {
							mVec4f attibs = FragmentShading(TriInEyeSpace, TriIndex, uv, recip_w, robj);
							//RasterResult.emplace_back(attibs.w, mVec2<int>{x, y }, attibs.tomVec3f());
							ZBufferTest({ InterpolatedZ, mVec2<int>{x, y }, attibs.tomVec3() });
						}
					}
				}

			}
		}


	};



#ifdef SingleThread
	auto policy = std::execution::seq;
#else
	auto policy = std::execution::par_unseq;
#endif
	std::for_each(policy, vec3iBuffer.begin(), vec3iBuffer.end(), RasterAtriangle);


	//interpolate
}


STRONG_INLINE mVec4f GraphicsPipeline::FragmentShading(const Triangle& TriInEyeSpace, const mVec3i& TriIdx, const mVec2<float>& uv, mVec3f& recip_w, const RenderableObject& robj) const {
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
	mVec3f tmp;
	_mm_prefetch(reinterpret_cast<char const*> (&vertex_normal_buffer[TriIdx.x]), 0);
	_mm_prefetch(reinterpret_cast<char const*> (&vertex_normal_buffer[TriIdx.y]), 0);
	_mm_prefetch(reinterpret_cast<char const*> (&vertex_normal_buffer[TriIdx.z]), 0);
	_mm_prefetch(reinterpret_cast<char const*> (&vertex_texture_buffer[TriIdx.x]), 0);
	_mm_prefetch(reinterpret_cast<char const*> (&vertex_texture_buffer[TriIdx.y]), 0);
	_mm_prefetch(reinterpret_cast<char const*> (&vertex_texture_buffer[TriIdx.z]), 0);

	auto normal_a = vertex_normal_buffer[TriIdx.x];
	auto normal_b = vertex_normal_buffer[TriIdx.y];
	auto normal_c = vertex_normal_buffer[TriIdx.z];

	auto vtex_a = vertex_texture_buffer[TriIdx.x];
	auto vtex_b = vertex_texture_buffer[TriIdx.y];
	auto vtex_c = vertex_texture_buffer[TriIdx.z];

	auto texChannel=  robj.texturechannelID;
	auto texture_scaling = robj.texture_scaling;

	if (UsePhongShading) {
		//interpolate normal
		mVec3f normal = normal_a + ((normal_c - normal_a) * uv.x + (normal_b - normal_a) * uv.y);

		//already normalize it in vertex shader, no need re-normalize
		//normal.normalize();//normal in eyespace,light source in eyespace

		mVec3f point = TriInEyeSpace.a * (1 - uv.x - uv.y) + (TriInEyeSpace.c * uv.x) + (TriInEyeSpace.b * uv.y);
		tmp = ADSshading({ 0,0,0 }, LightSourceEyeSpace, point, normal, robj.Material());

		mVec2<float> Texcoords = interpolate_varyings<mVec2f>(vtex_a,vtex_b,vtex_c , { (1 - uv.x - uv.y),uv.x,uv.y }, recip_w);
			

		//stripe texture
		if (texChannel == 1) {
			mVec3f backcolor = { 0,0,0 };
			float scale = texture_scaling;
			float fuzz = 2;
			float width = 10;
			float scaleT = fract(Texcoords.y * scale);

			float frac1 = clamp(scaleT / fuzz);
			float frac2 = clamp((scaleT - width) / fuzz);

			frac1 = frac1 * (1.0f - frac2);
			frac1 = frac1 * frac1 * (3.0f - (2.0f * frac1));
			tmp = backcolor * (frac1)+tmp * (1 - frac1);

		}
		else if (texChannel >= 2) {
			float scale = texture_scaling;
			Texcoords.x = fract(Texcoords.x * scale);
			Texcoords.y = fract(Texcoords.y * scale);
			mVec3f backcolor = LookupTexel(Texcoords, texChannel - 2);
			tmp = { backcolor.x * tmp.x,  backcolor.y * tmp.y, backcolor.z * tmp.z };
			tmp.ColorClamp();
		}

	}
	else {
		_mm_prefetch(reinterpret_cast<char const*> (&vertex_shading_buffer[TriIdx.x]), 0);
		_mm_prefetch(reinterpret_cast<char const*> (&vertex_shading_buffer[TriIdx.y]), 0);
		_mm_prefetch(reinterpret_cast<char const*> (&vertex_shading_buffer[TriIdx.z]), 0);
		auto perVshading_a = vertex_shading_buffer[TriIdx.x];
		auto perVshading_b = vertex_shading_buffer[TriIdx.y];
		auto perVshading_c = vertex_shading_buffer[TriIdx.z];
		//interpolate lightning
		tmp = perVshading_a * (1 - uv.x - uv.y) + (perVshading_c * uv.x) + (perVshading_b * uv.y);


	}
	return { tmp,1.0 };//opacity not used for now.

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
	_Camera.SetFrustm(1, -1, -1.0,1.0, 45.0f, InitEyePos.z + 400);

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
		mVec3f  Vtxs[3] = {
			vertex_coordinate_buffer[triangle.x].tomVec3(),
			vertex_coordinate_buffer[triangle.y].tomVec3(),
			vertex_coordinate_buffer[triangle.z].tomVec3()
		};

		float TriVtxDistance[3] = {
			c_plane.cal_project_distance(Vtxs[0]),
			c_plane.cal_project_distance(Vtxs[1]),
			c_plane.cal_project_distance(Vtxs[2]),
		};

		int remain_vert_num = 0;
		int remain_vert_id[4] = { -1,-1,-1,-1 };

		for (int i = 0; i < 3; i++)
		{
			current_index = i;
			previous_index = (i - 1 + 3) % 3;
			const int pre_v_id = thisTriangleVertices[previous_index];
			const int cur_v_id = thisTriangleVertices[current_index];
			float d1 = TriVtxDistance[previous_index];
			float d2 = TriVtxDistance[current_index];

			if (d1 * d2 < 0)
			{
				float t = d1 / (d1 - d2);
				mVec3f interpolateCordinate = Vtxs[current_index] * t + Vtxs[previous_index] * (1 - t);
				//attribute
				mVec3f interpolatenormal = vertex_normal_buffer[cur_v_id] * t + vertex_normal_buffer[pre_v_id] * (1 - t);
				mVec3f interpolateshading = vertex_shading_buffer[cur_v_id] * t + vertex_shading_buffer[pre_v_id] * (1 - t);

				mVec2f interpolateUV = vertex_texture_buffer[cur_v_id] * t + vertex_texture_buffer[pre_v_id] * (1 - t);

				int lastVid = vertex_coordinate_buffer.size();
				vertex_coordinate_buffer.emplace_back(interpolateCordinate, 1.0);
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
