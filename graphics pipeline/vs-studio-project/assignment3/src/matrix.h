#pragma once
#include <cmath>
#include <iostream>
#include<vector>
#include "la.h"

class Matrix4{
public:
	float p[4][4];
	//----------------------------
	//----------���캯��----------
	//----------------------------
	Matrix4();//Ĭ������Ϊ0
	Matrix4(float);//����ռ䣬Ĭ������Ϊ�������
	Matrix4(const Matrix4& r);
	Matrix4::Matrix4(std::vector<std::vector<float> > dvec);
	//-----------------------------
   //----------����������----------
   //------------------------------
	Matrix4& operator=(const Matrix4&);//����ĸ���
	Matrix4& operator+=(const Matrix4&);//�����+=����
	Matrix4& operator-=(const Matrix4&);//-=
	Matrix4& operator*=(const Matrix4&);//*=
	Matrix4 operator*(const Matrix4 & m)const;
	Matrix4 operator*(float scale)const;
	mVec4 Matrix4::operator*(mVec4 V)const;


	//-----------------------------
   //----------���ݻ�ȡ------------
   //------------------------------
	void Show() const;//������ʾ
	friend Matrix4 eye(int);//����һ����λ����
	//-----------------------------
   //----------�߼�����------------
   //------------------------------
	friend Matrix4 inv(Matrix4);//�����������
	float det();//����������ʽ
	Matrix4 Transpose();//����ת�õ�ʵ��,�Ҳ��ı����

};

//for transformation
Matrix4 scaleMatrix(float scalei);
Matrix4 scaleMatrix(float scalex, float scaley, float scalez);
Matrix4 translateMatrix(mVec3 translation);
Matrix4 rotateMatrix(int axis, float angle);
Matrix4 ViewPortMatrix(int Nx, int Ny);
Matrix4 OrthoMatrix(float r, float l, float t, float b, float n, float f);
Matrix4 PerspectiveMatrix(float r, float l, float t, float b, float n, float f);
Matrix4 ViewMatrix(mVec3 eyePos, mVec3 GazeDirection, mVec3 TopDirection);
Matrix4 rotateMatrix(mVec3 n, float theta);
Matrix4 squashMatrix(float n, float f);
