#pragma   once  
#include "la.h"
#include "matrix.h"
#include"Parameters.h"
#include "obj.h"
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
	Camera(mVec3 _position, mVec3 _front, mVec3 _up, float fov/*field of view �ӳ�*/);
	//ǰ��
	mVec3 front, Reset_front;
	//�Ϸ�
	mVec3 up, Reset_up;
	//�ҷ�
	mVec3 right;
	//���λ��
	mVec3 position, Reset_Pos;
	float fov;
	float scale;
	float pitch;
	float yaw;
	//���ɹ���
	Ray generatePerpectiveRay(int i, int j);
	Ray generateParallelRay(int j, int i);
	void Reset();
	void InitViewPlane(int nx = WINDOW_W, int ny = WINDOW_H, float z = 0, float l = -2, float r = 2, float b = -1.5, float t = 1.5);
	void UpdateFov(float z);
	void Camera::UpdatePitchAngle(float dpitch);
	void Camera::UpdateYawAngle(float dyaw);
	void Camera::UpdatePos(mVec3 z);
	void Camera::UpdateTarget(mVec3 _position); 
	//~Camera();
};
//void setupSence();
