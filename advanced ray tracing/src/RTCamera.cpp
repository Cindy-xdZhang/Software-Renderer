#include "RTCamera.h"

using namespace std;
//float rotateTime = 0.0f;
//eyepos should deprecate in future
//mVec3 eyepos = mVec3(0, 0, -1);
//int MousePosition_x=-1;
//int MousePosition_y=-1;
//mVec3 lightPos = mVec3(-12, 0, -12);
//Object_List objlist;
//int MouseSelectObjID = -1;
//int ObjsBuffer[WINDOW_H*WINDOW_W] = {-1};
//Camera* MyCamera;

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







