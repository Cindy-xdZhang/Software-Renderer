#pragma once
#define PI 3.141592653589793238
#define FLOAT_ZERO 10e-15
#define Radians(x) ((x)*PI/180.0)
#define MAX(a,b)    (((a)>(b)) ?  (a):(b))
#define MIN(a,b)    (((a)>(b)) ?  (b):(a))
#include<cmath>
#include <iostream>

template<typename  T>
struct mVec2{
	T x;
	T y;
	mVec2<T> operator +(mVec2<T>& right) const {
		return { right.x + x,right.y + y };
	}
	mVec2<T> operator -(mVec2<T>& right) const{
		return { x- right.x ,y-right.y };
	}
	mVec2<T> operator *(float right) const {
		return { x * right ,y * right};
	}
};


struct mVec3i {
	int x;
	int y;
	int z;

};
class mVec3 {
public:
	float x;
	float y;
	float z;
	//datatype_eum datatype;
	float getEuclideannNorms()const;
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
	void operator *=(float scale);

	mVec3 operator *(float scale) const;
	mVec3 operator /(float denominator) const;
	float operator *(const mVec3& right) const;

	mVec3 cross_product(const mVec3 & right) const;
	
	//以角度制表示
	float getAngle(const mVec3& right)const;
	void ColorClamp();
};
mVec3  normalize(mVec3 Iu);


class mVec4 {
public:
	float x;
	float y;
	float z;
	float w;
	//datatype_eum datatype;
	float getEuclideannNorms()const;
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

	void operator *=(float scale);
	mVec4 operator *(float scale) const;
	mVec4 operator /(float scale) const;
	float operator *(const mVec4& right) const;

	//mVec4 cross_product(const mVec4 & right) const;

	//以角度制表示
	float getAngle(const mVec4& right)const;
};

class Triangle {
public:
	mVec3 a, b, c;
	inline Triangle(mVec3 a, mVec3 b, mVec3 c) :a(a), b(b), c(c) {
	}
	Triangle() = default;
	 mVec2<float> Triangle::PointIsInTriangle(mVec3 p);
	 mVec3 Triangle::Normal();
	 mVec3 Triangle::Normal(mVec3 dir);
	 mVec3 Triangle::Center();

};

//
//class Plane 
//{
//public:
//	mVec3 n;		// normal
//	mVec3 px;	// solution to dot(n,p-p')=0  d
//	
//	Plane() = default;
//	Plane(mVec3 pm, mVec3 n); 
//};
//int clip_with_plane(Plane c_plane, std::vector<mVec3> vert_list, std::vector<mVec3>& in_list);