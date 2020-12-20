#pragma once
#define MAX_DIMENTION  10
#define MAIN_DATA_PRECISION float
#define PI 3.141592653589793238
#define FLOAT_ZERO 10e-15
#define MAX(a,b) ((a)>(b))? (a):(b)
#include<cmath>
#include <iostream>
#ifdef OldVersion
class Point3D {
public:
	MAIN_DATA_PRECISION x;
	MAIN_DATA_PRECISION y;
	MAIN_DATA_PRECISION z;
	Point3D() :x(0), y(0), z(0) {
	}
	Point3D(MAIN_DATA_PRECISION ix, MAIN_DATA_PRECISION iy, MAIN_DATA_PRECISION iz) :x(ix), y(iy), z(iz) {
	}
	Point3D operator *(float scale) const {
		return  Point3D(this->x*scale, this->y*scale, this->z*scale);
	}
	Point3D operator +(const Point3D scale) const {
		return  Point3D(this->x + scale.x, this->y + scale.y, this->z + scale.z);
	}
	Point3D operator -(const Point3D scale) const {
		return Point3D(this->x - scale.x, this->y - scale.y, this->z - scale.z);
	}
};
class mVec {
public:
	int di;
	MAIN_DATA_PRECISION data[MAX_DIMENTION];
	//datatype_eum datatype;
	MAIN_DATA_PRECISION getEuclideannNorms()const {
		MAIN_DATA_PRECISION sum = 0.0;
		for (int i = 0; i < di; i++)
		{
			sum += data[i] * data[i];
		}return sqrt(sum);

	}

	mVec(Point3D A, Point3D B) {
		this->di = 3;
		//this->data = new  MAIN_DATA_PRECISION[3];
		this->data[0] = A.x - B.x;
		this->data[1] = A.y - B.y;
		this->data[2] = A.z - B.z;
	}
	mVec(float a, float b, float c) {
		this->di = 3;
		this->data[0] = a;
		this->data[1] = b;
		this->data[2] = c;
	}
	mVec(float a, float b, float c, float cd) {
		this->di = 4;
		this->data[0] = a;
		this->data[1] = b;
		this->data[2] = c;
		this->data[3] = cd;
	}
     mVec(Point3D p,float w):mVec( p.x, p.y, p.z,w){
	}
	Point3D  toPoint3D() {
		return Point3D(this->data[0], this->data[1], this->data[2]);
	}
	mVec() {
		this->di = 0;
		//this->data = NULL;
	}
	mVec(float d) {
		this->di = d;
		this->data[0] = 0;
		this->data[1] = 0;
		this->data[2] = 0;
		this->data[3] = 0;
	}
	mVec(int d, const MAIN_DATA_PRECISION* array) {
		this->di = d;
		//this->data = new  float[d];
		memcpy_s(this->data, d * sizeof array, array, d * sizeof array);
	}
	~mVec() {
	}
	mVec operator +(const mVec& right) const {
		if (right.di != this->di)throw ("E: vetor + with different dimention.");
		mVec c = mVec(this->di, this->data);
		for (int i = 0; i < this->di; i++) {
			c.data[i] += right.data[i];
		}
		return c;
	}
	void operator +=(const mVec& right) {
		if (right.di != this->di)throw ("E: vetor += with different dimention.");
		for (int i = 0; i < this->di; i++) {
			this->data[i] += right.data[i];
		}
	}
	void operator -=(const mVec& right) {
		if (right.di != this->di)throw ("E: vetor -=  with different dimention.");
		for (int i = 0; i < this->di; i++) {
			this->data[i] -= right.data[i];
		}
	}
	mVec operator -(const mVec& right) const {
		if (right.di != this->di)throw ("E: vetor -=  with different dimention.");
		mVec c = mVec(this->di, this->data);
		for (int i = 0; i < this->di; i++) {
			c.data[i] -= right.data[i];
		}
		return c;
	}
	void normalize() {
		MAIN_DATA_PRECISION ma = (MAIN_DATA_PRECISION)this->getEuclideannNorms();
		(*this) *= (1 / ma);
	}
	template <typename T>
	void operator *=(T scale) {
		for (int i = 0; i < this->di; i++) {
			this->data[i] *= scale;
		}
	}
	template <typename T>
	mVec operator *(T scale) const {
		mVec c = mVec(this->di, this->data);
		for (int i = 0; i < this->di; i++) {
			c.data[i] *= scale;
		}return c;
	}
	MAIN_DATA_PRECISION operator *(const mVec& right) const { //向量相乘，结果为一个数，而不是向量
		if (right.di != this->di)throw ("E: vetor dot-product with different dimention.");
		MAIN_DATA_PRECISION ans = 0;
		for (int i = 0; i < this->di; i++) {
			ans += this->data[i] * right.data[i];
		}
		return ans;
	}

	mVec cross_product(const mVec & right) const {
		MAIN_DATA_PRECISION cro[3];
		cro[0] = (this->data[1] * right.data[2] - this->data[2] * right.data[1]);
		cro[1] = (this->data[2] * right.data[0] - this->data[0] * right.data[2]);
		cro[2] = (this->data[0] * right.data[1] - this->data[1] * right.data[0]);
		return mVec(3, cro);
	}

	//以角度制表示
	MAIN_DATA_PRECISION getAngle(const mVec& right)const {
		MAIN_DATA_PRECISION dot_product = (*this)*right;
		auto ma = this->getEuclideannNorms();
		auto mb = right.getEuclideannNorms();
		dot_product /= (ma*mb);
		return acos(dot_product) *180.0 / PI;
	}
};
#endif

class mVec3 {
public:
	MAIN_DATA_PRECISION x;
	MAIN_DATA_PRECISION y;
	MAIN_DATA_PRECISION z;
	//datatype_eum datatype;
	MAIN_DATA_PRECISION getEuclideannNorms()const;
	mVec3(mVec3 A, mVec3 B);
	mVec3(float a, float b, float c);

	mVec3();
	mVec3(float d);
	~mVec3();
	mVec3 operator +(const mVec3& right) const;
	void operator +=(const mVec3& right);
	void operator -=(const mVec3& right);
	mVec3 operator -(const mVec3& right) const;
	void normalize();
	void operator *=(MAIN_DATA_PRECISION scale);

	mVec3 operator *(MAIN_DATA_PRECISION scale) const;
	MAIN_DATA_PRECISION operator *(const mVec3& right) const;

	mVec3 cross_product(const mVec3 & right) const;
	
	//以角度制表示
	MAIN_DATA_PRECISION getAngle(const mVec3& right)const;
};
mVec3  normalize(mVec3 Iu);
class mVec4 {
public:
	MAIN_DATA_PRECISION x;
	MAIN_DATA_PRECISION y;
	MAIN_DATA_PRECISION z;
	MAIN_DATA_PRECISION w;
	//datatype_eum datatype;
	MAIN_DATA_PRECISION getEuclideannNorms()const;
	mVec3 tomVec3();
	mVec4(mVec4 A, mVec4 B);
	mVec4(float a, float b, float c, float cd);
	mVec4();
	mVec4(float d);
	mVec4(mVec3 r, float w);
	~mVec4();
	mVec4 operator +(const mVec4& right) const;
	void operator +=(const mVec4& right);
	void operator -=(const mVec4& right);
	mVec4 operator -(const mVec4& right) const;
	void normalize();

	void operator *=(MAIN_DATA_PRECISION scale);
	mVec4 operator *(MAIN_DATA_PRECISION scale) const;
	MAIN_DATA_PRECISION operator *(const mVec4& right) const;

	//mVec4 cross_product(const mVec4 & right) const;

	//以角度制表示
	MAIN_DATA_PRECISION getAngle(const mVec4& right)const;
};

