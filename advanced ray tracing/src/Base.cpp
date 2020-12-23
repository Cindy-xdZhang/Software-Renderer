#include "Base.h"
#include "objparser.h"
#include "timer.hpp"
#include <thread>

using namespace std;
mVec3 initlightPos = mVec3(-7, 0, -12);


//eyepos should deprecate in future
mVec3 eyepos = mVec3(0, 6, -12);
int MousePosition_x=-1;
int MousePosition_y=-1;

Object_List objlist;
Camera* MyCamera;
BvhTopLevelStructure* my;
std::vector<Surface>AllObjects;
//int* ObjsBuffer=new int[WINDOW_W*WINDOW_H] ;
bool UseDistributedRaySSP = false;
bool UseDistributedRaySoftShadow = false;
Camera::Camera(mVec3 _position, mVec3 _front, mVec3 _up, float fov=45.0f)
{
	n= std::normal_distribution<double> (0.5, 0.5);


	Reset_Pos = _position;
	Reset_up = _up;
	Reset_front = _front;

	position = _position;

	//保存相机前方向
	front = normalize(_front);
	//	//计算相机上方向
	up = normalize(_up);
	//计算相机右方向

	right = front.cross_product(up);
	right.normalize();
	
	//up = right.cross_product(front);
	//up.normalize();
	//这里是计算相机投影范围的，暂时考虑的是正方形区域
	scale = tan(fov * 0.5 *PI / 180) * 2; 
	InitViewPlane();
	yaw = 90.0f;
	pitch = 0;
}
void Camera::UpdateFov(float z) {
	fov += z;
	if (fov > 120.0f)fov = 120.0f;
	if (fov < 45.0f)fov = 45.0f;
	scale = tan(fov * 0.5 *PI / 180) * 2;
}
void Camera::UpdateTarget(mVec3 _position) {
	

	mVec3 direction(_position, position);

	this->front = normalize(direction);

	right = front.cross_product(up);
	right.normalize();
}
//竖直
void Camera::UpdatePitchAngle(float dpitch) {
	pitch += dpitch;
	if (pitch > 89.0f)pitch = 89.0f;
	if (pitch < -89.0f)pitch = -89.0f;

	mVec3 direction(1.0f);
	direction.y = cos(Radians(yaw))*cos(Radians(pitch));
	direction.x = sin(Radians(pitch));
	direction.z= sin(Radians(yaw))*cos(Radians(pitch));
	this->front =normalize(direction) ;
	

	//计算相机上方向
	up = right.cross_product(front);
	up.normalize();
	//right = front.cross_product(up);
	//right.normalize();
}
//水平
void Camera::UpdateYawAngle( float dyaw) {
	yaw += dyaw;

	mVec3 direction(1.0f);
	direction.y = cos(Radians(yaw))*cos(Radians(pitch));
	direction.x = sin(Radians(pitch));
	direction.z = sin(Radians(yaw))*cos(Radians(pitch));
	this->front = direction;
	//计算相机上方向
	/*up = right.cross_product(front);
	up.normalize();*/
	right = front.cross_product(up);
	right.normalize();
	
}
void Camera:: InitViewPlane(int nx , int ny , float z , float l , float r, float b , float t ) {

	pixel_spacing_x = (r - l) / nx;
	pixel_spacing_y = (t - b) / ny;
}
Ray Camera::generatePerpectiveRay(int j, int i)
{
	i = i - WINDOW_H / 2;
	j = j - WINDOW_W / 2;
	mVec3 u = right * (i + 0.5) * (pixel_spacing_x);
	mVec3 v = up * (j + 0.5) * (pixel_spacing_y) ;
	

	return Ray(position, normalize(front*scale + u + v) );

}
Ray Camera::generateDistributedPerpectiveRay(int j, int i,mVec3 pqn)
{

	float randomKese = n(e);

	i = i - WINDOW_H / 2;
	j = j - WINDOW_W / 2;

	mVec3 u = right * (i + (pqn.x+randomKese)/ pqn.z) * (pixel_spacing_x);
	mVec3 v = up * (j + (pqn.y + randomKese) / pqn.z) * (pixel_spacing_y);


	return Ray(position, normalize(front*scale + u + v));

}

Ray Camera::generateParallelRay(int j, int i)
{
	i = i - WINDOW_H / 2;
	j = j - WINDOW_W / 2;
	mVec3 u = right * (i + 0.5) * (pixel_spacing_x);
	mVec3 v = up * (j + 0.5) * (pixel_spacing_y);


	return Ray(position + u + v, normalize(front));

}
void Camera::Reset() {
	position = Reset_Pos;
	//保存相机前方向
	front = normalize(Reset_front);
	//计算相机右方向

	right = front.cross_product(Reset_up);
	right.normalize();
	//计算相机上方向
	up = right.cross_product(front);
	up.normalize();
	//这里是计算相机投影范围的，暂时考虑的是正方形区域
	fov = 45.0f;
	scale = tan(fov * 0.5 *PI / 180) * 2;
}
//p的x y z分别带表 相机上下 左右 前后
void Camera::UpdatePos(mVec3 p) {
	position += (up*p.x + right * p.y + front * p.z);
}


class ViewPlane {
public:
	int nx;
	int ny;
private:
	float z;
	float r;
	float l;
	float b;
	float t;
	float pixel_spacing_x;
	float pixel_spacing_y;
public:
	ViewPlane(int nx = WINDOW_W, int ny = WINDOW_H, float z = 0, float l = -2, float r = 2, float b = -1.5, float t = 1.5) :nx(nx), ny(ny), z(z), l(l), r(r), b(b), t(t) {
		pixel_spacing_x = (r - l) / nx;
		pixel_spacing_y = (t - b) / ny;
	}
	//uv到cols rows像素坐标
	mVec3 Cordinate_Translate_cr2uv(int i, int j) {
		return mVec3(l + pixel_spacing_y*  (i + 0.5), b + pixel_spacing_x * (j + 0.5), this->z);
	}
};

void setupSence() {

	mVec3  Up = mVec3(-1, 0, 0);
	mVec3  Down = mVec3(1, 0, 0);

	//Sphere* spoRed = new Sphere(0.4f, 4.0f, 1.0f, 1.0f);
	//Sphere* spoGreen = new Sphere(mVec3(-4, 4, 4), 0.3);//head
	//Sphere* spoblue = new Sphere(mVec3(-1, 2, 6), 1);//budy
	//Sphere* spoyellow = new Sphere(mVec3(-4, 2, 5), 0.3);

	Plane* Base_plane = new Plane(mVec3(1.4, 0, 0), Up);

	//Triangular_pyramid* pyramid = new Triangular_pyramid(mVec3(-1, -1, 1.4), mVec3(1.4, -2, 1), mVec3(1.4, -4.4, 2), mVec3(1.4, -0.2, 2));
	//Cubid*pcubid = new Cubid(mVec3(1, 1, 2), 1.0f, 1.0f, 1.0f);

	//VertexInMesh v1 = mVec3(-1, -1, 1.4);
	//VertexInMesh v2 = mVec3(1.4, -2, 1);
	//VertexInMesh v3 = mVec3(1.4, -2, 2);
	//VertexInMesh v4 = mVec3(1.4, -0.2, 2);
	//VertexInMesh v5 = mVec3(1.4, -0.24, 1);
	//std::vector<VertexInMesh > tmpv = { v1, v2 , v3, v4 };
	//std::vector<TriangleInMesh > tmpt = { { 0,1,2 },{ 0,2,3 },{ 0,1,3 },{ 3,1,2 } };
	//MeshObjects* TestMesh = new MeshObjects(tmpt,tmpv);
	//TestMesh->boudleMaterial(mVec3(0.3, 0.9, 1) * 0.001, mVec3(1, 0.9, 0.7) * 1, mVec3(1, 1, 1) *1.0);
	//spoGreen->boudleMaterial(mVec3(0, 1, 0) * 0.01, mVec3(0, 1, 0)* 1.5, mVec3(1, 1, 1)*1.0);//8

	Base_plane->boudleMaterial(mVec3(0.5, 0.5, 0.5) * 0.01, mVec3(0.5, 0.5, 0.5) * 1, mVec3(0.5, 0.5, 0.5) * 1,64);
	//spoyellow->boudleMaterial(mVec3(1, 0.9, 0) * 0.01, mVec3(0.9, 0.9, 0.3) * 1, mVec3(1, 1, 1) *1.0);
	//spoblue->boudleMaterial(mVec3(0, 0.2, 1) * 0.01, mVec3(0, 0.2, 1) * 1, mVec3(1, 1, 1) * 1.0);
	////spoRed2->boudleMaterial(mVec3(1, 0, 0) * 0.01, mVec3(1, 0, 0)* 1.5, mVec3(1, 1, 1) *1.0);//6.3
	
	

	MyCamera = new Camera(eyepos, mVec3(0, -0.2, 1), Down, 45.0f);

	std::string path = "C:\\Users\\8\\source\\repos\\newell_teaset\\teapot.obj";
	OBjReader ojp;
	ShareVertexMesh Target=ojp.readObj2ShareMesh(path);
	MeshObjects TestMesh2 =  MeshObjects(Target);
	TestMesh2.boudleMaterial(mVec3(1.1, 1.0, 0.01) * 0.01, mVec3(1.0, 0.9, 1.0) * 1.5, mVec3(1, 1, 0.1) *1.0, 64);

	
	//objlist.AddObj(*TestMesh2);
	//objlist.AddObj(*TestMesh3);

	std::string path2 = "C:\\Users\\8\\source\\repos\\newell_teaset\\teacup.obj";
	OBjReader ojp2;
	ShareVertexMesh Target2 = ojp2.readObj2ShareMesh(path2);
	std::string path3 = "C:\\Users\\8\\source\\repos\\newell_teaset\\spoon.obj";
	OBjReader ojp3;
	ShareVertexMesh Target3 = ojp3.readObj2ShareMesh(path3);
	MeshObjects TestMesh4 =  MeshObjects(Target2);
	MeshObjects TestMesh5 =  MeshObjects(Target3);
	MeshObjects TestMesh6 =  MeshObjects(Target3);
	TestMesh4.updateModelMatrix(translateMatrix(mVec3(-1, 1, 6)));
	TestMesh5.updateModelMatrix(translateMatrix(mVec3(1, 3, 1)));
	TestMesh6.updateModelMatrix(translateMatrix(mVec3(-3, 4, -3)));
	TestMesh4.boudleMaterial(mVec3(0.0, 1.0, 1) * 0.005, mVec3(0, 1.0, 1.0) * 1.5, mVec3(1, 1, 1) *1.0, 1);
	TestMesh5.boudleMaterial(mVec3(0.0, 0.0, 1) * 0.005, mVec3(0, 0.0, 1.0) * 1.5, mVec3(1, 1, 1) *1.0, 1);
	TestMesh6.boudleMaterial(mVec3(1.0, 0.0, 1.1) * 0.006, mVec3(1, 0.5, 0.1) * 1.8, mVec3(1, 1, 1) *1.0, 1);

	/*objlist.AddObj(*TestMesh4);
	objlist.AddObj(*TestMesh5);
	objlist.AddObj(*TestMesh6);
	objlist.AddObj(*Base_plane);*/
	std::vector<MeshObjects>Meshes;
	Meshes.emplace_back(TestMesh2);
	Meshes.emplace_back(TestMesh5);
	for (int i = 0; i < 7; i++) {

		Meshes.emplace_back(TestMesh6);
		Meshes.emplace_back(TestMesh4);
	}

    my=new BvhTopLevelStructure(Meshes);
	objlist.AddObj(*my);
	objlist.AddObj(*Base_plane);

}

std::vector<mVec3> RandomLightSource(int n) {
	std::default_random_engine e; //引擎
	std::normal_distribution<double>nor(0.5,1);
	std::vector<mVec3> liht;
	mVec3 facebase1 = {0,0,0.25};
	mVec3 facebase2= { 0,0.25,0 };
	for(int n1=0;n1<n;n1++)
		for (int n2 = 0; n2 < n; n2++) {
			float seta1 = nor(e);
			float seta2= nor(e);
			mVec3 lightpose = initlightPos + facebase1 * seta1 + facebase2 * seta2;
			liht.emplace_back(lightpose);
		}

	return liht;

}

std::vector<mVec3> RandomEyeBias(int n) {
	std::default_random_engine e; //引擎
	std::normal_distribution<double>nor(0.5, 0.25);
	std::vector<mVec3> liht;
	mVec3 facebase1 = { 1,0,0 };
	mVec3 facebase2 = { 0,1,0 };
	for (int n1 = 0; n1 < n; n1++)
		for (int n2 = 0; n2 < n; n2++) {
			float seta1 = nor(e);
			float seta2 = nor(e);
			seta1 = seta1 > 1.0 ? 1.0 : seta1;
			seta1 = seta1 < 0.0 ? 0.0 : seta1;
			seta2 = seta2 > 1.0 ? 1.0 : seta2;
			seta2 = seta2 < 0.0 ? 0.0 : seta2;
			mVec3 lightpose =  facebase1 * seta1 + facebase2 * seta2;
			liht.emplace_back(lightpose);
		}

	return liht;

}

mVec3 ShadingARay(const Object_List& uobjlist, const Ray&cRay,mVec3 LightSourse) {
	float minIntersectionTime = (float)INT_MAX;
	mVec3 FinalShadingColor;
	int objs_len = uobjlist.Len();
	
	auto ADSshading_shadow = [&](ShadingMaterial Mp, mVec3 point, mVec3 normal,int cuid) {
		mVec3 Ia = Mp.Ka*LigntIntensity;
		//FloatRGB Imax = FloatRGB(1, 1, 1);
		//Reflection
		mVec3 ReflectionD = mVec3(LightSourse, point);
		ReflectionD.normalize();
		Ray Reflection = Ray(point, ReflectionD);
		if (objlist.OccludeByOther(cuid,Reflection)) { Ia.ColorClamp(); return  Ia; }


		float distance2 = pow((LightSourse.x - point.x), 2) + pow((LightSourse.y - point.y), 2) + pow((LightSourse.z - point.z), 2);
		float Intensity = (LigntIntensity / distance2);
		mVec3 light_dir(LightSourse, point);
		mVec3 eye_dir(MyCamera->position, point);
		light_dir.normalize();
		eye_dir.normalize();
		//normal.normalize(); already normaled after interpolation
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

	for (int objid = 0; objid < objs_len; objid++) {
		auto Cobjp = uobjlist.GetPointer(objid);
		int Cobjptype = uobjlist.GetType(objid);

		if (Cobjptype == XD_PLANE) {
			Plane Cobj = *((Plane*)Cobjp);
			Matrix4 M = Cobj.pIModleMatrix;
			Matrix4 normalM = Cobj.pNormalMatrix;
			Ray tRay = Ray((M*mVec4(cRay.ori, 1.0f)).tomVec3(), (M*mVec4(cRay.dir, 0.0f)).tomVec3());
			float time = intersection(Cobj, tRay);
			mVec3  intersectPoint = cRay.point(time);
			if (time > 0 && time < minIntersectionTime) {
				 mVec3 Normal = (normalM * mVec4(Cobj.n, 0.0)).tomVec3();
				 FinalShadingColor= ADSshading_shadow((Cobj.Material()), intersectPoint, Normal, objid);
				 minIntersectionTime = time;
			}
		}
		else if (Cobjptype == XD_MESH) {
			MeshObjects* Cobj = ((MeshObjects*)Cobjp);
			Matrix4 M = (Cobj->pIModleMatrix);
			Matrix4 normalM = (Cobj->pNormalMatrix);
			Ray tRay = Ray((M*mVec4(cRay.ori, 1.0f)).tomVec3(), (M*mVec4(cRay.dir, 0.0f)).tomVec3());

			intersectInfo time_n = Cobj->intersection(Cobj->pBVHmesh, tRay);
			float time = time_n.time;
			mVec3  intersectPoint = cRay.point(time);
	
			if (time > 0 && time < minIntersectionTime) {
				mVec3 Normal = (normalM * mVec4(time_n.n, 0.0)).tomVec3();
				FinalShadingColor = ADSshading_shadow((Cobj->Material()), intersectPoint, Normal, objid);
				minIntersectionTime = time;
			
			}
		}
		else if (Cobjptype == XD_TOPLEVELMESH) {
			BvhTopLevelStructure* Cobj = ((BvhTopLevelStructure*)Cobjp);
		

			intersectInfo time_n = Cobj->intersection(Cobj->pBVHmesh, cRay);
			float time = time_n.time;
			mVec3  intersectPoint = cRay.point(time);

			if (time > 0 && time < minIntersectionTime) {
			
				FinalShadingColor = ADSshading_shadow(time_n.shm, intersectPoint, time_n.n, objid);
				minIntersectionTime = time;

			}
		}
	}//for obj
	return FinalShadingColor;
}

//===============multi thread simple_ray_tracing===============

void RenderCols(int startcol, int endcol, ViewPlane vp, Object_List objlist, mVec3* framebuffer)
{
	int distributedSample = 5;
	float sampleWeight = 1.0f/((float)distributedSample*distributedSample);
	auto lihnts=RandomLightSource(distributedSample);
	auto eyebias = RandomEyeBias(distributedSample);

	for (int i = startcol; i < endcol; i++)
		for (int j = 0; j < WINDOW_H; j++) {
			mVec3  uv_pixel = vp.Cordinate_Translate_cr2uv(j, i);
			mVec3 Color;
			if (UseDistributedRaySSP) {
				for (int samle1 = 0; samle1 < distributedSample; samle1++) {
					for (int samle2 = 0; samle2 < distributedSample; samle2++){
						mVec3 pqn = { (float)samle1 ,(float)samle2 ,(float)distributedSample };

					Ray  cRay = MyCamera->generateDistributedPerpectiveRay(j, i, pqn);
					//MyCamera->position += eyebias[distributedSample*distributedSample-(samle1*distributedSample + samle2)-1];
	
					if (UseDistributedRaySoftShadow)
						Color += ShadingARay(objlist, cRay, lihnts[samle1*distributedSample + samle2])*sampleWeight;
					else Color += ShadingARay(objlist, cRay, lihnts[0])*sampleWeight;
					}
				}
			}
			else {
				Ray  cRay = MyCamera->generatePerpectiveRay(j, i);
				Color += ShadingARay(objlist, cRay, initlightPos);
			}
			
			framebuffer[j*WINDOW_W + i] = Color;
	
		}//for pixel

}



void multithread_simple_ray_tracing(mVec3* framebuffer) {

	//x,y,z -》y为列数 left and right（0-640），x为行数up and down（0-480）,x越大越向下； z越大越远
	//mVec  Up = mVec (-1, 0, 0)mVec  right = mVec(0, 1, 0);
	mVec3 Eye = eyepos;
	//mVec3 LightSourse = Lightpos;
	ViewPlane vp = ViewPlane();



	int threadcount = 16;
	int WorkloadPerThread = floorf(WINDOW_W / threadcount);
	std::vector<  std::thread >Threads;
	for (int i = 0; i < threadcount; i++) {
		std::thread th0(RenderCols, WorkloadPerThread*i, (i + 1) * WorkloadPerThread, vp, objlist, framebuffer);
		Threads.push_back(std::move(th0));
	}
	for (auto& th : Threads) {
		th.join();
	}
}
//===============single thread simple_ray_tracing===============
void simple_ray_tracing(mVec3* framebuffer, bool perspective, mVec3 Lightpos) {
	/*MyTimer timer;
	timer.begin();
	float* depth_buffer = new float[WINDOW_W*WINDOW_H]();
	mVec3 Eye = eyepos;
	mVec3 LightSourse = Lightpos;
	ViewPlane vp = ViewPlane();

	int objs_len = objlist.Len();
	for (int i = 0; i < WINDOW_H; i++)
		for (int j = 0; j < WINDOW_W; j++) {
			mVec3  uv_pixel = vp.Cordinate_Translate_cr2uv(i, j);
			Ray  cRay;
			//for prespective projction
			if (perspective == true) {
				 //cRay = Ray(Eye, mVec3(uv_pixel, Eye));
				cRay = MyCamera->generatePerpectiveRay(i, j);
			}
			//for parallel projction
			else {
				mVec3 uv_ = mVec3(uv_pixel.x, uv_pixel.y, 0);
				cRay = Ray(Eye + uv_, mVec3(0, 0, 1));
			}
			for (int objid = 0; objid < objs_len; objid++) {
				parameter opt = { objlist, objid, i, j, cRay, depth_buffer, framebuffer };
				RenderObject(&opt);
			}//for obj
		}//for pixel
	/*timer.end();
	cout << "Ray tracing on intel i5-6300 spends " << timer.interval << " s." << endl;
	delete depth_buffer; */
}





