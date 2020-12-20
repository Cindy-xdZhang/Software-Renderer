#pragma once
#include <cmath>
#include <iostream>
#include<vector>
#include "la.h"
#ifdef VECTORBASE_Matrix
class Matrix {
private:
	int rows_num, cols_num;
	
	void allocMem(int, int);//��ʼ������
	void swapRows(int, int);
public:
	std::vector< std::vector<float>>  p;
	//----------------------------
	//----------���캯��----------
	//----------------------------
	Matrix() = default;
	Matrix(int, int);//����ռ䣬Ĭ������Ϊ0
	Matrix(int, int, float);//����ռ䣬Ĭ������Ϊ�������
	Matrix(const Matrix& r);
	friend Matrix inv(Matrix);//�����������
	//~Matrix();//��������Ӧ�����麯�������Ǵ��಻��������
	Matrix::Matrix(std::vector<std::vector<float> > dvec);
	//-----------------------------
   //----------����������----------
   //------------------------------
	Matrix& operator=(const Matrix&);//����ĸ���
	Matrix& operator+=(const Matrix&);//�����+=����
	Matrix& operator-=(const Matrix&);//-=
	Matrix& operator*=(const Matrix&);//*=
	Matrix operator*(const Matrix & m)const;
	Matrix operator*(float scale)const;
	mVec4 Matrix::operator*(mVec4 V)const;


	//-----------------------------
   //----------���ݻ�ȡ------------
   //------------------------------
	std::vector<std::vector<float> >&  getData();
	const std::vector<std::vector<float >>& getData() const;
	int row() const;
	int col() const;
	void Show() const;//������ʾ
	float At(int i, int j) const;
	friend Matrix eye(int);//����һ����λ����
	//-----------------------------
   //----------�߼�����------------
   //------------------------------
	static Matrix Solve(const Matrix&, const Matrix&);//������Է�����Ax=b
	float det();//����������ʽ
	 Matrix Transpose();//����ת�õ�ʵ��,�Ҳ��ı����
	Matrix gaussianEliminate();//��˹��Ԫ��
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
	//void allocMem(int, int);//��ʼ������
	//void swapRows(int, int);
public:
	float p[4][4];
	//----------------------------
	//----------���캯��----------
	//----------------------------
	Matrix4();//����ռ䣬Ĭ������Ϊ0
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
	//for transformation
	friend Matrix4 scaleMatrix(float scalei);
	friend Matrix4 scaleMatrix(float scalex, float scaley, float scalez);
	//#TODO:
	friend   Matrix4 translateMatrix(mVec3 translation);
	friend  Matrix4 rotateMatrix(int axis, float angle);
};

