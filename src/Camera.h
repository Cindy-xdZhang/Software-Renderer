#pragma   once  
#define WINDOW_H 480*1.5
#define WINDOW_W 640*1.5

#include<cmath>
#include "la.h"
#include "timer.hpp"
#include "matrix.h"

class Camera
{
private:

	float fardistance{};
	float yaw{};
	float pitch{};
	
public:
	float fov{};
	float r{};
	float l{};
	float b{};
	float t{};
	Camera()=default;
	Camera(mVec3f _position, mVec3f _front, mVec3f _up);
	float getNearPlane() const;
	float getFarPlane() const;
	mVec3f front{};
	mVec3f up{};
	mVec3f right{};
	mVec3f position{};


	//void Reset();
	void Camera::SetFrustm(float r, float l, float b, float t, float fov, float farplane);
	void Camera::UpdatePitchAngle(float dpitch);
	void Camera::UpdateYawAngle(float dyaw);
	void Camera::UpdatePos(mVec3f z);
	void Camera::UpdateFov(float fov);
	//void Camera::UpdateTarget(mVec3f _position); 
	Matrix4 Camera::genViewMat();
	Matrix4 Camera::genPerspectiveMat();
	//~Camera();
};
class ArcBallControler {
	int w,h;
public:
	ArcBallControler() = default;
	ArcBallControler(int w, int h): w(w), h(h){
	}

	mVec3f GetArcBallPositionVector(int x, int y) const;
	Matrix4 GetArcBallrotateMatrix(mVec3f a, mVec3f b) const;
};