#pragma   once  
#include"Parameters.h"
#include "la.h"
#include "timer.hpp"
#include "obj.h"
#include "matrix.h"
#include <random>
class Camera
{

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
	//Camera();
	Camera(mVec3 _position, mVec3 _front, mVec3 _up, float fov/*field of view 视场*/);
	//前方
	mVec3 front, Reset_front;
	//上方
	mVec3 up, Reset_up;
	//右方
	mVec3 right;
	//相机位置
	mVec3 position, Reset_Pos;
	float fov;
	float scale;
	float pitch;
	float yaw;
	//生成光线
	std::default_random_engine e; //引擎
	std::normal_distribution<double> n ; 
	Ray generatePerpectiveRay(int i, int j);
	Ray generateParallelRay(int j, int i);
	Ray Camera::generateDistributedPerpectiveRay(int j, int i, mVec3 pqn);


	void Reset();
	void InitViewPlane(int nx = WINDOW_W, int ny = WINDOW_H, float z = 0, float l = -2, float r = 2, float b = -1.5, float t = 1.5);
	void UpdateFov(float z);
	void Camera::UpdatePitchAngle(float dpitch);
	void Camera::UpdateYawAngle(float dyaw);
	void Camera::UpdatePos(mVec3 z);
	void Camera::UpdateTarget(mVec3 _position); 
	//~Camera();
};
void setupSence();
void simple_ray_tracing(mVec3* framebuffer, bool perspective,mVec3 Lightpos);
void multithread_simple_ray_tracing(mVec3* framebuffer);//no parallel view anymore