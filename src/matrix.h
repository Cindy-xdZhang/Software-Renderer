#pragma once
#include <cmath>
#include <iostream>
#include<vector>
#include "la.h"

class Matrix4{
public:
	float p[4][4];

	Matrix4();
	Matrix4(float);
	Matrix4(const Matrix4& r);
	Matrix4::Matrix4(std::vector<std::vector<float> > dvec);
	//-----------------------------
   //------------------------------
	Matrix4& operator=(const Matrix4&);
	Matrix4& operator+=(const Matrix4&);
	Matrix4& operator-=(const Matrix4&);//-=
	Matrix4& operator*=(const Matrix4&);//*=
	Matrix4 operator*(const Matrix4 & m)const;
	Matrix4 operator*(float scale)const;
	mVec4f operator*(mVec4f V)const;


	//-----------------------------
   //------------------------------
	void Show() const;
	friend Matrix4 eye(int);
	//-----------------------------
   //------------------------------
	friend Matrix4 inv(Matrix4);
	float det();
	Matrix4 Transpose();

};

//for transformation
Matrix4 scaleMatrix(float scalei);
Matrix4 scaleMatrix(float scalex, float scaley, float scalez);
Matrix4 translateMatrix(mVec3f translation);
Matrix4 rotateMatrix(int axis, float angle);
Matrix4 ViewPortMatrix(int Nx, int Ny);
Matrix4 OrthoMatrix(float r, float l, float t, float b, float n, float f);
Matrix4 PerspectiveMatrix(float r, float l, float t, float b, float n, float f);
Matrix4 ViewMatrix(mVec3f eyePos, mVec3f GazeDirection, mVec3f TopDirection);
Matrix4 rotateMatrix(mVec3f n, float theta);
Matrix4 squashMatrix(float n, float f);
