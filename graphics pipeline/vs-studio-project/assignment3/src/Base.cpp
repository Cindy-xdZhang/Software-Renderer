#include "Base.h"
#include "timer.hpp"
#include <thread>
using namespace std;


Camera::Camera(mVec3 _position, mVec3 _front, mVec3 _up)
{


	position = _position;

	//保存相机前方向
	front = normalize(_front);
	//	//计算相机上方向
	up = normalize(_up);
	//计算相机右方向

	right = front.cross_product(up);
	right.normalize();
	
	//nearPlane =  ( tan(fov * 0.5 *PI / 180)) ;
	fov = 90.0f;
	yaw = 90.0f;
	pitch = 0;
}

void Camera::UpdateFov(float dfov) {
	 this->fov = dfov+ fov;
	 if (this->fov > 160.0f) this->fov = 160.0f;
	 if (this->fov < 12.50f)	this->fov = 12.50f;


}
//fov only for vertical direction
void Camera::SetFrustm(float r, float l, float b, float t, float dfov, float farplane)
{   
	this->fardistance = farplane;
	this->r = (r);
	this->l = (l);
	this->b=(b);
	this->t=(t);
	this->fov = dfov;
	if (this->fov > 160.0f) this->fov = 160.0f;
	if (this->fov < 12.50f)	this->fov = 12.50f;
	
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
	direction.z = -sin(Radians(yaw))*cos(Radians(pitch));
	this->front = direction;
	//计算相机right方向
	right = front.cross_product(up);
	right.normalize();
	
}

//p的x y z分别带表 相机上下 左右 前后
void Camera::UpdatePos(mVec3 p) {
	position += (up*p.y+ right * p.x + front * p.z);
}


Matrix4 Camera::genViewMat() {
   return	ViewMatrix(this->position, front, up);
}

Matrix4 Camera::genPerspectiveMat() {
	float Flength = abs(b) / tan(Radians((this->fov / 2)));
	float n = ( Flength  * this->front.z);
	float f = (fardistance * this->front.z);
	return	PerspectiveMatrix(r,l,t,b,n,f);
}


mVec3 ArcBallControler::GetArcBallPositionVector(int x, int y) {
	float rx = (float(2 * x) / float(w)) - 1;
	float ry = (float(2 * y) / float(h)) - 1;
	float square = rx * rx + ry * ry;
	float z;
	if (square > 1.0f) {
		float dis = sqrtf(square);
		rx = rx / dis;
		rx = rx / dis;
		z = 0.0f;
	}
	else {
		z = sqrtf(1 - square);
	}

	return mVec3(rx, ry, z);
}

Matrix4 ArcBallControler::GetArcBallrotateMatrix(mVec3 a, mVec3 b) {
	float ElmentDotproduct = a * b;
	float aValue = a.x*a.x + a.y*a.y + a.z*a.z;
	float bValue = b.x*b.x + b.y*b.y + b.z*b.z;
	float cosTheta = ElmentDotproduct/(aValue*bValue);
	cosTheta = cosTheta > 1 ? 1: cosTheta;
	cosTheta = cosTheta < -1 ? -1 : cosTheta;
	float Theta = -acos(cosTheta);//in radians form
	//Theta = 180 * Theta / PI;
	mVec3 axis = a.cross_product(b);
	axis.normalize();
	Matrix4 tmp=rotateMatrix(axis, Theta);
	return tmp;

}

float Camera::getNearPlane() {
	float Flength = abs(b) / tan(Radians((this->fov / 2)));
	return (Flength  * this->front.z);
	
}
float Camera::getFarPlane() {
	return (fardistance * this->front.z);

}