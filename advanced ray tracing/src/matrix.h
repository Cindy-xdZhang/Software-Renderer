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
	Matrix4(std::vector<std::vector<float> > dvec);

	Matrix4& operator=(const Matrix4&);//
	Matrix4& operator+=(const Matrix4&);//
	Matrix4& operator-=(const Matrix4&);//-=
	Matrix4& operator*=(const Matrix4&);//*=
	Matrix4 operator*(const Matrix4 & m)const;
	Matrix4 operator*(float scale)const;
	mVec4 operator*(mVec4 V)const;



	void Show() const;


	friend Matrix4 inv(Matrix4);
	friend Matrix4 inv2(const Matrix4 A);
	float det();
	Matrix4 Transpose();

};


Matrix4 eye(int n);

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
