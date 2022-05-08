#include "Camera.h"
#include <assert.h>
#include "timer.hpp"
#include <thread>
using namespace std;


Camera::Camera(mVec3f _position, mVec3f _front, mVec3f _up)
{


	position = _position;

	front = normalize(_front);
	up = normalize(_up);

	right = front.cross_product(up);
	right.normalize();
	
	fov = 60.0f;
	yaw = 90.0f;
	pitch = 0;
}

void Camera::UpdateFov(float dfov) {
	 this->fov = dfov+ fov;
	 if (this->fov > 160.0f) this->fov = 160.0f;
	 if (this->fov < 12.50f)	this->fov = 12.50f;
	 float FocalLength = abs(b) / tan(Radians((this->fov / 2)));
	 this->nearplane = FocalLength;

}
//fov only for vertical direction
void Camera::SetFrustm(float r, float l, float b, float t, float dfov, float farplane)
{   
	this->farPlane = farplane;
	this->r = (r);
	this->l = (l);
	this->b=(b);
	this->t=(t);
	this->fov = dfov;
	if (this->fov > 160.0f) this->fov = 160.0f;
	if (this->fov < 12.50f)	this->fov = 12.50f;
	float FocalLength = abs(b) / tan(Radians((this->fov / 2)));
	this->nearplane = FocalLength;

}

void Camera::UpdatePitchAngle(float dpitch) {
	pitch += dpitch;
	if (pitch > 89.0f)pitch = 89.0f;
	if (pitch < -89.0f)pitch = -89.0f;

	mVec3f direction(1.0f);
	direction.x = cos(Radians(yaw)) * cos(Radians(pitch));
	direction.y = sin(Radians(pitch));
	direction.z = -sin(Radians(yaw)) * cos(Radians(pitch));
	//front in word space
	this->front =normalize(direction) ;
	

	up = right.cross_product(front);
	up.normalize();
	//right = front.cross_product(up);
	//right.normalize();
}

void Camera::UpdateYawAngle( float dyaw) {
	yaw += dyaw;

	mVec3f direction(1.0f);
	direction.x = cos(Radians(yaw))*cos(Radians(pitch));
	direction.y = sin(Radians(pitch));
	direction.z = -sin(Radians(yaw))*cos(Radians(pitch));
	this->front = direction;
	right = front.cross_product(up);
	right.normalize();
	
}

void Camera::UpdatePos(mVec3f p) {
	position += (up*p.y+ right * p.x + front * p.z);

}


Matrix4 Camera::genViewMat() const {
   return	ViewMatrix(this->position, front, up);
}

Matrix4 Camera::genPerspectiveMat() const{
	/*float n = (nearplane * this->front.z);
	this is wrong as perseptiveMatrix is considered in eye space, so n,f always in -z direction*/
	float n = (nearplane * -1);
	float f = (farPlane * -1);

	return	PerspectiveMatrix(r,l,t,b,n,f);
}


mVec3f ArcBallControler::GetArcBallPositionVector(int x, int y) const{
	assert(w > 0 && h > 0);
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

	return mVec3f(rx, ry, z);
}

Matrix4 ArcBallControler::GetArcBallrotateMatrix(mVec3f a, mVec3f b) const {
	if (a==b)
	{
		return eye(4);
	}
	float ElmentDotproduct = a * b;
	float aValue = a.x*a.x + a.y*a.y + a.z*a.z;
	float bValue = b.x*b.x + b.y*b.y + b.z*b.z;
	float cosTheta = ElmentDotproduct/(aValue*bValue);
	cosTheta = cosTheta > 1 ? 1: cosTheta;
	cosTheta = cosTheta < -1 ? -1 : cosTheta;
	float Theta = -acos(cosTheta);//in radians form
	Theta *= sensitivity;
	//Theta = 180 * Theta / PI;
	mVec3f axis = a.cross_product(b);
	axis.normalize();
	Matrix4 tmp=rotateMatrix(axis, Theta);
	assert(!isnan(tmp.p[0][0]) && !isnan(tmp.p[1][1]) && !isnan(tmp.p[2][2]) && !isnan(tmp.p[3][3]));
	return tmp;

}

float Camera::getNearPlane() const {
	return (nearplane  * -1);
	
}
float Camera::getFarPlane() const {
	return (farPlane * -1);

}