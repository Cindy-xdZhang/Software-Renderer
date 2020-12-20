#pragma once
#include <cmath>
#include <iostream>
#include<vector>
#include "la.h"
#ifdef VECTORBASE_Matrix
class Matrix {
private:
	int rows_num, cols_num;
	
	void allocMem(int, int);//初始化矩阵
	void swapRows(int, int);
public:
	std::vector< std::vector<float>>  p;
	//----------------------------
	//----------构造函数----------
	//----------------------------
	Matrix() = default;
	Matrix(int, int);//分配空间，默认数据为0
	Matrix(int, int, float);//分配空间，默认数据为输入参数
	Matrix(const Matrix& r);
	friend Matrix inv(Matrix);//求矩阵的逆矩阵
	//~Matrix();//析构函数应当是虚函数，除非此类不用做基类
	Matrix::Matrix(std::vector<std::vector<float> > dvec);
	//-----------------------------
   //----------操作符重载----------
   //------------------------------
	Matrix& operator=(const Matrix&);//矩阵的复制
	Matrix& operator+=(const Matrix&);//矩阵的+=操作
	Matrix& operator-=(const Matrix&);//-=
	Matrix& operator*=(const Matrix&);//*=
	Matrix operator*(const Matrix & m)const;
	Matrix operator*(float scale)const;
	mVec4 Matrix::operator*(mVec4 V)const;


	//-----------------------------
   //----------数据获取------------
   //------------------------------
	std::vector<std::vector<float> >&  getData();
	const std::vector<std::vector<float >>& getData() const;
	int row() const;
	int col() const;
	void Show() const;//矩阵显示
	float At(int i, int j) const;
	friend Matrix eye(int);//制造一个单位矩阵
	//-----------------------------
   //----------高级运算------------
   //------------------------------
	static Matrix Solve(const Matrix&, const Matrix&);//求解线性方程组Ax=b
	float det();//求矩阵的行列式
	 Matrix Transpose();//矩阵转置的实现,且不改变矩阵
	Matrix gaussianEliminate();//高斯消元法
	//for transformation
	friend Matrix scaleMatrix(float scalei);
	friend Matrix scaleMatrix(float scalex, float scaley, float scalez);
	//#TODO:
	friend   Matrix translateMatrix(mVec3 translation);
	friend Matrix rotateMatrix(int axis, float angle);
};
#endif

class Matrix4{
private:
	//void allocMem(int, int);//初始化矩阵
	//void swapRows(int, int);
public:
	float p[4][4];
	//----------------------------
	//----------构造函数----------
	//----------------------------
	Matrix4();//分配空间，默认数据为0
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
	//for transformation
	friend Matrix4 scaleMatrix(float scalei);
	friend Matrix4 scaleMatrix(float scalex, float scaley, float scalez);
	//#TODO:
	friend   Matrix4 translateMatrix(mVec3 translation);
	friend  Matrix4 rotateMatrix(int axis, float angle);
};

