#pragma once
#include <cmath>
#include <iostream>
#include<vector>
#include "la.h"

class Matrix4{
public:
	float p[4][4];
	//----------------------------
	//----------构造函数----------
	//----------------------------
	Matrix4();//默认数据为0
	Matrix4(float);//分配空间，默认数据为输入参数
	Matrix4(const Matrix4& r);
	Matrix4::Matrix4(std::vector<std::vector<float> > dvec);
	//-----------------------------
   //----------操作符重载----------
   //------------------------------
	Matrix4& operator=(const Matrix4&);//矩阵的复制
	Matrix4& operator+=(const Matrix4&);//矩阵的+=操作
	Matrix4& operator-=(const Matrix4&);//-=
	Matrix4& operator*=(const Matrix4&);//*=
	Matrix4 operator*(const Matrix4 & m)const;
	Matrix4 operator*(float scale)const;
	mVec4 Matrix4::operator*(mVec4 V)const;


	//-----------------------------
   //----------数据获取------------
   //------------------------------
	void Show() const;//矩阵显示
	friend Matrix4 eye(int);//制造一个单位矩阵
	//-----------------------------
   //----------高级运算------------
   //------------------------------
	friend Matrix4 inv(Matrix4);//求矩阵的逆矩阵
	float det();//求矩阵的行列式
	Matrix4 Transpose();//矩阵转置的实现,且不改变矩阵

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
