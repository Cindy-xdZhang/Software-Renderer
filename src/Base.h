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

	float fardistance;
	float yaw;
	float pitch;
	
public:
	float fov;
	float r;
	float l;
	float b;
	float t;
	Camera()=default;
	Camera(mVec3 _position, mVec3 _front, mVec3 _up);
	//前方
	mVec3 front;
	//上方
	mVec3 up;
	//右方
	mVec3 right;
	//相机位置
	mVec3 position;


	//void Reset();
	void Camera::SetFrustm(float r, float l, float b, float t, float fov, float farplane);
	void Camera::UpdatePitchAngle(float dpitch);
	void Camera::UpdateYawAngle(float dyaw);
	void Camera::UpdatePos(mVec3 z);
	void Camera::UpdateFov(float fov);
	//void Camera::UpdateTarget(mVec3 _position); 
	Matrix4 Camera::genViewMat();
	Matrix4 Camera::genPerspectiveMat();
	float Camera::getNearPlane();
	float Camera::getFarPlane();
	//~Camera();
};
class ArcBallControler {
	int w,h;
public:
	ArcBallControler() = default;
	ArcBallControler(int w, int h): w(w), h(h){
	}
	mVec3 GetArcBallPositionVector(int x, int y);

	Matrix4 GetArcBallrotateMatrix(mVec3 x, mVec3 y);
};