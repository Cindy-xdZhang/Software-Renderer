#include "Base.h"
#include "timer.hpp"
#include <thread>
#define I 60
using namespace std;
float rotateTime = 0.0f;
//eyepos should deprecate in future
mVec3 eyepos = mVec3(0, 0, -1);
int MousePosition_x=-1;
int MousePosition_y=-1;
mVec3 lightPos = mVec3(-12, 0, -12);
Object_List objlist;
int MouseSelectObjID = -1;
int ObjsBuffer[WINDOW_H*WINDOW_W] = {-1};
Camera* MyCamera;

Camera::Camera(mVec3 _position, mVec3 _front, mVec3 _up, float fov=45.0f)
{
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



FloatRGB ADSshading_shadow(mVec3  Sourse, mVec3 eye, mVec3 point, mVec3 normal, ShadingMaterial Mp, Object_List ol, int skip_id, float f = 8) {

	FloatRGB Ia = Mp.Ka*I;
	//FloatRGB Imax = FloatRGB(1, 1, 1);
	//Reflection
	mVec3 ReflectionD = mVec3(Sourse, point);
	Ray Reflection = Ray(point, ReflectionD);
	if (ol.OccludeByOther(skip_id, Reflection)) { Ia.Clamp(); return  Ia; }


	float distance2 = pow((Sourse.x - point.x), 2) + pow((Sourse.y - point.y), 2) + pow((Sourse.z - point.z), 2);
	float Intensity = (I / distance2);
	mVec3 light_dir(Sourse, point);
	mVec3 eye_dir(eye, point);
	light_dir.normalize();
	eye_dir.normalize();
	normal.normalize();
	mVec3 h = (eye_dir + light_dir);
	h.normalize();
	float dotproduct = normal * light_dir;
	FloatRGB Id = Mp.Kd*Intensity*(dotproduct > 0 ? dotproduct : 0);
	float dotproduct2 = normal * h;
	float tmp = dotproduct2 > 0.0f ? dotproduct2 : 0.0f;
	tmp = pow(tmp, f);
	FloatRGB Is = Mp.Ks*Intensity*tmp;
	FloatRGB  If = (Is + Id + Ia);
	If.Clamp();
	return If;

}

typedef struct parameter {
	Object_List objlist;
	int objid; int i; int j;
	Ray cRay;
	mVec3 LightSourse;
	mVec3 Eye;
	float* depth_buffer;
	FloatRGB* framebuffer;
}*pparameter;
void RenderObject(pparameter uu) {
	int i = uu->i;
	int j = uu->j;
	Ray cRay = uu->cRay;
	float* depth_buffer = uu->depth_buffer;
	mVec3 LightSourse = uu->LightSourse;
	Object_List uobjlist = uu->objlist;
	int objid = uu->objid;
	FloatRGB* framebuffer = uu->framebuffer;
	auto Cobjp = uobjlist.GetPointer(uu->objid);
	int Cobjptype = uobjlist.GetType(uu->objid);

	if (Cobjptype == XD_SPHERE) {
		Sphere Cobj = *((Sphere*)Cobjp);
		Matrix4 M = *Cobj.pIModleMatrix;
		Matrix4 normalM = *Cobj.pNormalMatrix;
		Ray tRay = Ray((M*mVec4(cRay.ori, 1.0f)).tomVec3(), (M*mVec4(cRay.dir,0.0f)).tomVec3());

		float time = intersection(Cobj, tRay);
		mVec3  intersectPoint = tRay.point(time);
		float distance = intersectPoint.z - MyCamera->position.z;
		distance = distance > 0 ? distance : -distance;
		if (time > 0 && (depth_buffer[i*WINDOW_W + j] == 0 || (depth_buffer[i*WINDOW_W + j] > 0 && distance < depth_buffer[i*WINDOW_W + j]))) {
			mVec3 Normal = Cobj.Normal(intersectPoint);
			Normal = (normalM * mVec4(Normal,0.0)).tomVec3();
			framebuffer[i*WINDOW_W + j] = ADSshading_shadow(LightSourse, MyCamera->position, intersectPoint, Normal, (Cobj.Material()), uobjlist, objid, 64);
			depth_buffer[i*WINDOW_W + j] = distance;
			ObjsBuffer[i*WINDOW_W + j] = objid;
		}

	}
	else if (Cobjptype == XD_CONE) {
		Cone Cobj = *((Cone*)Cobjp);
		Matrix4 M = *Cobj.pIModleMatrix;
		Matrix4 normalM = *Cobj.pNormalMatrix;
		Ray tRay = Ray((M*mVec4(cRay.ori, 1.0f)).tomVec3(), (M*mVec4(cRay.dir, 0.0f)).tomVec3());

		float time = intersection(Cobj, tRay);
		mVec3  intersectPoint = tRay.point(time);
		float distance = intersectPoint.z - MyCamera->position.z;
		distance = distance > 0 ? distance : -distance;
		if (time > 0 && (depth_buffer[i*WINDOW_W + j] == 0 || (depth_buffer[i*WINDOW_W + j] > 0 && distance <= depth_buffer[i*WINDOW_W + j]))) {
			mVec3 Normal = Cobj.Normal(intersectPoint);
			Normal = (normalM * mVec4(Normal, 0.0)).tomVec3();
			framebuffer[i*WINDOW_W + j] = ADSshading_shadow(LightSourse, MyCamera->position, intersectPoint, Normal, (Cobj.Material()), uobjlist, objid, 64);
			depth_buffer[i*WINDOW_W + j] = distance;
			ObjsBuffer[i*WINDOW_W + j] = objid;
		}
	}
	else if (Cobjptype == XD_PLANE) {
		Plane Cobj = *((Plane*)Cobjp);
		Matrix4 M = *Cobj.pIModleMatrix;
		Matrix4 normalM = *Cobj.pNormalMatrix;
		Ray tRay = Ray((M*mVec4(cRay.ori, 1.0f)).tomVec3(), (M*mVec4(cRay.dir, 0.0f)).tomVec3());
		float time = intersection(Cobj, tRay);
		mVec3  intersectPoint = tRay.point(time);
		float distance = intersectPoint.z - MyCamera->position.z;
		distance = distance > 0 ? distance : -distance;
		if (time > 0 && (depth_buffer[i*WINDOW_W + j] == 0 || (depth_buffer[i*WINDOW_W + j] > 0 && distance <= depth_buffer[i*WINDOW_W + j]))) {
			mVec3 Normal = Cobj.n;
			Normal = (normalM * mVec4(Normal, 0.0)).tomVec3();
			framebuffer[i*WINDOW_W + j] = ADSshading_shadow(LightSourse, MyCamera->position, intersectPoint, Normal, (Cobj.Material()), uobjlist, objid, 64);
			depth_buffer[i*WINDOW_W + j] = distance;
			ObjsBuffer[i*WINDOW_W + j] = objid;
		}
	}
	else if (Cobjptype == XD_PARAMRID) {
		Triangular_pyramid Cobj = *((Triangular_pyramid*)Cobjp);
		Matrix4 M = *Cobj.pIModleMatrix;
		Matrix4 normalM = *Cobj.pNormalMatrix;
		Ray tRay = Ray((M*mVec4(cRay.ori, 1.0f)).tomVec3(), (M*mVec4(cRay.dir, 0.0f)).tomVec3());

		intersectInfo time_n = intersection(Cobj, tRay);
		float time = time_n.time;
		mVec3  intersectPoint = tRay.point(time);
		float distance = intersectPoint.z - MyCamera->position.z;
		distance = distance > 0 ? distance : -distance;
		if (time > 0 && (depth_buffer[i*WINDOW_W + j] == 0 || (depth_buffer[i*WINDOW_W + j] > 0 && distance <= depth_buffer[i*WINDOW_W + j]))) {
			mVec3 Normal = time_n.n;
			Normal = (normalM * mVec4(Normal, 0.0)).tomVec3();
			framebuffer[i*WINDOW_W + j] = ADSshading_shadow(LightSourse, MyCamera->position, intersectPoint, Normal, (Cobj.Material()), uobjlist, objid, 64);
			depth_buffer[i*WINDOW_W + j] = distance;
			ObjsBuffer[i*WINDOW_W + j] = objid;
		}
	}
	else if (Cobjptype == XD_AABB) {
		Cubid Cobj = *((Cubid*)Cobjp);
		Matrix4 M = *Cobj.pIModleMatrix;
		Matrix4 normalM = *Cobj.pNormalMatrix;
		Ray tRay = Ray((M*mVec4(cRay.ori, 1.0f)).tomVec3(), (M*mVec4(cRay.dir, 0.0f)).tomVec3());

		intersectInfo time_n = intersection(Cobj, tRay);
		float time = time_n.time;
		mVec3  intersectPoint = tRay.point(time);
		float distance = intersectPoint.z - MyCamera->position.z;
		distance = distance > 0 ? distance : -distance;
		if (time > 0 && (depth_buffer[i*WINDOW_W + j] == 0 || (depth_buffer[i*WINDOW_W + j] > 0 && distance <= depth_buffer[i*WINDOW_W + j]))) {
			mVec3 Normal = time_n.n;
			Normal = (normalM * mVec4(Normal, 0.0)).tomVec3();
			framebuffer[i*WINDOW_W + j] = ADSshading_shadow(LightSourse, MyCamera->position, intersectPoint, Normal, (Cobj.Material()), uobjlist, objid, 64);
			depth_buffer[i*WINDOW_W + j] = distance;
			ObjsBuffer[i*WINDOW_W + j] = objid;
		}
	}


}

//===============multi thread simple_ray_tracing===============

void RenderCols(int startcol, int endcol, ViewPlane vp, Object_List objlist, mVec3 LightSourse, mVec3 Eye, float* depth_buffer, FloatRGB* framebuffer, bool perspective)
{
	int objs_len = objlist.Len();
	for (int i = startcol; i < endcol; i++)
		for (int j = 0; j < WINDOW_H; j++) {
			mVec3  uv_pixel = vp.Cordinate_Translate_cr2uv(j, i);
			Ray  cRay;
			if (perspective == true) {
			//cRay = Ray(Eye, mVec3(uv_pixel, Eye));
			cRay = MyCamera->generatePerpectiveRay(j, i);
			}
			//for parallel projction
			else {
				/*mVec3 uv_ = mVec3(uv_pixel.x, uv_pixel.y, 0);
				cRay = Ray(Eye + uv_, mVec3(0, 0, 1));*/
				cRay = MyCamera->generateParallelRay(j, i);
			}
			for (int objid = 0; objid < objs_len; objid++) {
				parameter opt = { objlist, objid, j,i, cRay, LightSourse, Eye, depth_buffer, framebuffer };
				RenderObject(&opt);
			}//for obj
		}//for pixel

}

void setupSence() {
	
	mVec3  Up = mVec3(-1, 0, 0);
	mVec3  Down = mVec3(1, 0, 0);

	

	Sphere* spoRed =new Sphere( 0.4f, 4.0f, 1.0f, 1.0f);
	//Sphere* spoRed2 = new Sphere(0.4f, 4.0f, 1.0f, 1.0f);
	//Sphere* spoRed3 = new Sphere(0.4f, 4.0f, 1.0f, 1.0f);

	Sphere* spoGreen = new Sphere(mVec3(-4, 4, 4 ),0.3);//head

	Sphere* spoblue = new Sphere(mVec3 (-1, 2, 6),1);//budy
	Sphere* spoyellow = new Sphere(mVec3(-4, 2, 5), 0.3);

	Plane* Base_plane = new Plane(mVec3(1.4, 0, 0), Up);
	//Cone* coo = new Cone(mVec3(-4.6, 1.7, 6), Down, 6.0, 0.985);
	Cone* coo = new Cone(mVec3(-3.7, 4, 4), Down, 6.0, 0.985);
	Triangular_pyramid* pyramid = new Triangular_pyramid(mVec3(-1, -1 , 1.4 ), mVec3(1.4, -2 , 1 ), mVec3(1.4, -4.4 , 2 ), mVec3(1.4, -0.2 , 2));
	Cubid*pcubid= new Cubid(mVec3(1, 1, 2), 1.0f, 1.0f, 1.0f);

	pcubid->boudleMaterial(FloatRGB(0.3, 0.9, 1) * 0.01, FloatRGB(1, 0.9, 0.2) * 1, FloatRGB(1, 1, 1) *1.0);//6.3
	pyramid->boudleMaterial(FloatRGB(1, 0.9, 0) * 0.01, FloatRGB(1, 0.9, 0.2) * 1, FloatRGB(1, 1, 1) *1.0);//6.3
	spoRed->boudleMaterial(FloatRGB(1, 0, 0) * 0.01, FloatRGB(1, 0, 0)* 1.5, FloatRGB(1, 1, 1) *1.0);//6.3
	spoGreen->boudleMaterial(FloatRGB(0, 1, 0) * 0.01, FloatRGB(0, 1, 0)* 1.5, FloatRGB(1, 1, 1)*1.0);//8
	coo->boudleMaterial(FloatRGB(0, 0, 1) * 0.01, FloatRGB(0, 0, 1)*1.2, FloatRGB(1, 1, 1) * 1.0);//3.2
	Base_plane->boudleMaterial(FloatRGB(0.5, 0.5, 0.5) * 0.01, FloatRGB(0.5, 0.5, 0.5) * 1, FloatRGB(0.5, 0.5, 0.5) * 1);
	spoyellow->boudleMaterial(FloatRGB(1, 0.9, 0) * 0.01, FloatRGB(0.9, 0.9, 0.3) *1, FloatRGB(1, 1, 1) *1.0);
	spoblue->boudleMaterial(FloatRGB(0, 0.2, 1) * 0.01, FloatRGB(0, 0.2, 1)*1, FloatRGB(1, 1, 1) * 1.0);
	//spoRed2->boudleMaterial(FloatRGB(1, 0, 0) * 0.01, FloatRGB(1, 0, 0)* 1.5, FloatRGB(1, 1, 1) *1.0);//6.3
	
	objlist.AddObj(*spoGreen);
	objlist.AddObj(*spoblue);
	objlist.AddObj(*coo);
	objlist.AddObj(*pyramid);
	objlist.AddObj(*Base_plane);
	objlist.AddObj(*spoRed);
	objlist.AddObj(*spoyellow);
	objlist.AddObj(*pcubid);

	MyCamera = new Camera(eyepos, mVec3 (0,0,1), Down, 45.0f);
}

void multithread_simple_ray_tracing(FloatRGB* framebuffer, bool perspective, mVec3 Lightpos) {
	MyTimer timer;
	timer.begin();
	float* depth_buffer = new float[WINDOW_W*WINDOW_H]();
	//x,y,z -》y为列数 left and right（0-640），x为行数up and down（0-480）,x越大越向下； z越大越远
	//mVec  Up = mVec (-1, 0, 0)mVec  right = mVec(0, 1, 0);
	mVec3 Eye = eyepos;
	mVec3 LightSourse = Lightpos;
	ViewPlane vp = ViewPlane();
	//for selection
#ifdef CALSELECTION
	if (MousePosition_x != -1 && MousePosition_y != -1) {
		mVec3  uv_pixel = vp.Cordinate_Translate_cr2uv(MousePosition_x, MousePosition_y);
		Ray  cRay;
		if (perspective == true)
			cRay = Ray(Eye, mVec3(uv_pixel, Eye));
		//for parallel projction
		else {
			mVec3 uv_ = mVec3(uv_pixel.x, uv_pixel.y, 0);
			cRay = Ray(Eye + uv_, mVec3(0, 0, 1));
		}
			
		int object_count = objlist.Len();
		if (MouseSelectObjID != -1) {
			auto Cobjpre = objlist.GetPointer(MouseSelectObjID);
			Surface* Cobjprep = ((Surface*)Cobjpre);
			Cobjprep->ResetMaterial();
			MouseSelectObjID = -1;
		}
		float minZ = 10000000;
		
		for (int u = 0; u < object_count; u++)
		{
			float time;
			auto Cobjp =objlist.GetPointer(u);
			unsigned char Cobjptype =objlist.GetType(u);
			switch (Cobjptype) {
			case (XD_SPHERE):
			{Sphere Cobj = *((Sphere*)Cobjp);
			time = intersection(Cobj, cRay);
			break; }
			case (XD_PLANE):
			{Plane Cobj = *((Plane*)Cobjp);
			time = intersection(Cobj, cRay);
			break; }
			case (XD_CONE):
			{Cone Cobj = *((Cone*)Cobjp);
			time = intersection(Cobj, cRay);
			break; }
			case (XD_PARAMRID):
			{Triangular_pyramid Cobj = *((Triangular_pyramid*)Cobjp);
			intersectInfo time_n = intersection(Cobj, cRay);
			time = time_n.time;
			break; }
			}
			if (time > 0) {
				mVec3  intersectPoint = cRay.point(time);
				if (intersectPoint.z < minZ) {
					minZ = intersectPoint.z;
					MouseSelectObjID =u;
				}
			}
		}
		if (MouseSelectObjID != -1) {
			auto Cobjp = objlist.GetPointer(MouseSelectObjID);
			Surface* Cobj = ((Surface*)Cobjp);
			Cobj->updateMaterial(FloatRGB(1, 0.02, 1) * 0.01, FloatRGB(1, 0.02, 1)*1.2, FloatRGB(1, 1, 1) * 1.0);
	  }
		MousePosition_x = -1;
		MousePosition_y = -1;
	}
#endif
	if (MousePosition_x != -1 && MousePosition_y != -1) {
		if (MouseSelectObjID != -1) {
			auto Cobjpre = objlist.GetPointer(MouseSelectObjID);
			Surface* Cobjprep = ((Surface*)Cobjpre); 
			Cobjprep->ResetMaterial();
			MouseSelectObjID = -1;
		}
		MouseSelectObjID = ObjsBuffer[MousePosition_x*WINDOW_W+MousePosition_y];
		if (MouseSelectObjID != -1) {
			auto Cobjp = objlist.GetPointer(MouseSelectObjID);
			Surface* Cobj = ((Surface*)Cobjp);
			Cobj->updateMaterial(FloatRGB(1, 0.02, 1) * 0.01, FloatRGB(1, 0.02, 1)*1.2, FloatRGB(1, 1, 1) * 1.0);
		}
		MousePosition_x = -1;
		MousePosition_y = -1;
	}


	int threadcount = 8;
	int interval_count = vp.nx / threadcount;
	std::thread th0(RenderCols, 0, interval_count, vp, objlist, LightSourse, Eye, depth_buffer, framebuffer, perspective);
	std::thread th1(RenderCols, interval_count, 2 * interval_count, vp, objlist, LightSourse, Eye, depth_buffer, framebuffer, perspective);
	std::thread th2(RenderCols, 2 * interval_count, 3 * interval_count, vp, objlist, LightSourse, Eye, depth_buffer, framebuffer, perspective);
	std::thread th3(RenderCols, 3 * interval_count, 4 * interval_count, vp, objlist, LightSourse, Eye, depth_buffer, framebuffer, perspective);
	std::thread th4(RenderCols, 4 * interval_count, 5 * interval_count, vp, objlist, LightSourse, Eye, depth_buffer, framebuffer, perspective);
	std::thread th5(RenderCols, 5 * interval_count, 6 * interval_count, vp, objlist, LightSourse, Eye, depth_buffer, framebuffer, perspective);
	std::thread th6(RenderCols, 6 * interval_count, 7 * interval_count, vp, objlist, LightSourse, Eye, depth_buffer, framebuffer, perspective);
	std::thread th7(RenderCols, 7 * interval_count, 8 * interval_count, vp, objlist, LightSourse, Eye, depth_buffer, framebuffer, perspective);
	th0.join();
	th1.join();
	th2.join();
	th3.join();
	th4.join();
	th5.join();
	th6.join();
	th7.join();
	timer.end();
	cout << "Multi-thread Ray tracing on intel i5-6300 spends " << timer.interval << " s." << endl;
	delete depth_buffer;
}
//===============single thread simple_ray_tracing===============
void simple_ray_tracing(FloatRGB* framebuffer, bool perspective, mVec3 Lightpos) {
	MyTimer timer;
	timer.begin();
	float* depth_buffer = new float[WINDOW_W*WINDOW_H]();
	//x,y,z -》y为列数 left and right（0-640），x为行数up and down（0-480）,x越大越向下； z越大越远
	//mVec  Up = mVec (-1, 0, 0)mVec  right = mVec(0, 1, 0);
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
				parameter opt = { objlist, objid, i, j, cRay, LightSourse, Eye, depth_buffer, framebuffer };
				RenderObject(&opt);
			}//for obj
		}//for pixel
	timer.end();
	cout << "Ray tracing on intel i5-6300 spends " << timer.interval << " s." << endl;
	delete depth_buffer;
}






#ifdef DEBUG
void Test_intersection() {
	//测试intersection：	
	mVec3  Orig_P(0, 0, 0);
	mVec3  s(0, 9, 0);
	mVec U(s, Orig_P);
	auto ooo = Ray(Orig_P, U);
	//不相交  （0，0，0）+t（0,9,0）与 （1,0,0）（0,0，2）（9，0,4） 、
	mVec3  PA(1, 0, 0);
	mVec3  PB(0, 0, 2);
	mVec3  PC(9, 0, 4);
	Triangle  tri(PA, PB, PC);
	float kd = intersection(tri, ooo);
	//相交于边 （0，0，0）+t（0,9,0） 与 （1,2,0）（0,2，2）（0，2,-1） 、
	mVec3  PA1(1, 2, 0);
	mVec3  PB1(0, 2, 2);
	mVec3  PC1(0, 2, -1);
	Triangle  tri1(PA1, PB1, PC1);
	float kd1 = intersection(tri1, ooo);
	//相交于顶点（0，0，0）+t（0,9,0） 与 （1,2,0）（0,2，2）（0，4,0）、（1,2,0）（0,2，2）（0，2,0）
	mVec3  PA2(1, 2, 0);
	mVec3  PB2(0, 2, 2);
	mVec3  PC2(0, 4, 0);
	Triangle  tri2(PA2, PB2, PC2);
	float kd2 = intersection(tri2, ooo);
	//相交于面内（0，0，0）+t（0,9,0） 与 （1,0,0）（0,0，2）（-1，5,-2）
	mVec3  PA3(1, 0, 0);
	mVec3  PB3(0, 0, 2);
	mVec3  PC3(-1, 5, -2);
	Triangle  tri3(PA3, PB3, PC3);
	float kd3 = intersection(tri3, ooo);
	//光线从三角形出发向外射（0，0，0）+t（0,9,0）与 （1,0,0）（-1,0，-2）（-1，0,4）
	mVec3  PA4(1, 0, 0);
	mVec3  PB4(-1, 0, -2);
	mVec3  PC4(-1, 0, 4);
	Triangle  tri4(PA4, PB4, PC4);
	float kd4 = intersection(tri4, ooo);
}
void runTestFunction() {
	float a[3] = { 4,3,0 };
	mVec A(3, a);
	float b[3] = { 1,0,2 };
	mVec B(3, b);
	mVec  C = A + B;
	mVec  D = A - B;
	D -= A;
	D *= (int)2;
	D.normalize();
	double dotsum = A * B;
	auto u = A.getAngle(B);
	double x = sin(u / 180 * PI);
	auto R = A.cross_product(B);

	auto u2 = A.getAngle(A);
	auto u3 = B.getAngle(R);
	mVec3  PA(0, 2, 20);
	mVec3  PB(10, 1, -2);
	mVec3  PC(-10, 1, -2);
	mVec3  Orig_P(0, 0, 0);
	mVec3  s(0, 2, 0);
	mVec U(s, Orig_P);
	Triangle  tri(PA, PB, PC);
	/*tri.Square();*/
	auto ooo = Ray(Orig_P, U);
	float kd = intersection(tri, ooo);
	Test_intersection();
	mVec3  Eye(0, 0, -1);
	mVec3  Light(3, -3, 3);
	mVec3  intr(2, 2, 2);
	/*ShadingMaterial mty(0.3,0.3,0.9);
	ADSshading(Light,Eye, intr,U, mty);*/
	cout << "begin shading!" << endl;
	//simple_ray_tracing();
	cout << "finish shading!" << endl;
}
#endif // DEBUG
