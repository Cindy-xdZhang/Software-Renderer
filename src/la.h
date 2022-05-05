#pragma once
#include"Macros.h"
#include<cmath>
#include <iostream>
#include <vector>

template<typename  T>
struct mVec2{
	T x;
	T y;
	mVec2 < T >() = default;
	inline mVec2 < T >(T x, T y) :x(x), y(y) {
	}

	mVec2 < T >(const mVec2 < T >&other) :x(other.x), y(other.y) {

	};
	mVec2 < T >& operator=(const mVec2 < T >&rhs) {
		this->x = rhs.x;
		this->y = rhs.y;
		return *this;
	}

	mVec2<T>(mVec2<T> && other) :x(std::move(other.x)), y(std::move(other.y))  {
	}

	mVec2<T>& operator=(mVec2<T> && rhs)  {
		this->x = rhs.x;
		this->y = rhs.y;
		return *this;
	}

	inline mVec2<T> operator +(mVec2<T>& right) const {
		return { right.x + x,right.y + y };
	}
	inline mVec2<T> operator -(mVec2<T>& right) const{
		return { x- right.x ,y-right.y };
	}
	inline mVec2<T> operator *(float right) const {
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
	MAIN_DATA_PRECISION x;
	MAIN_DATA_PRECISION y;
	MAIN_DATA_PRECISION z;
	//datatype_eum datatype;
	MAIN_DATA_PRECISION getEuclideannNorms()const;
	mVec3(mVec3 A, mVec3 B);
	mVec3(float a, float b, float c);
	inline mVec3(mVec3&& other):x(other.x), y(other.y),z(other.z) {
	}
	 mVec3(const mVec3& other) : x(other.x), y(other.y), z(other.z) {
	}
	 mVec3& operator =(const mVec3&) = default; 
	 mVec3& operator =( mVec3&&) = default;

	mVec3()=default;
	mVec3(float d);
	~mVec3();
	mVec3 operator +(const mVec3& right) const;
	void operator +=(const mVec3& right);
	void operator -=(const mVec3& right);
	mVec3 operator -(const mVec3& right) const;
	void normalize();
	void operator *=(MAIN_DATA_PRECISION scale);

	mVec3 operator *(MAIN_DATA_PRECISION scale) const;
	mVec3 operator /(MAIN_DATA_PRECISION denominator) const;
	MAIN_DATA_PRECISION operator *(const mVec3& right) const;

	mVec3 cross_product(const mVec3 & right) const;
	
	//以角度制表示
	MAIN_DATA_PRECISION getAngle(const mVec3& right)const;
	void ColorClamp();
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
	mVec4 operator /(MAIN_DATA_PRECISION scale) const;
	MAIN_DATA_PRECISION operator *(const mVec4& right) const;

	//mVec4 cross_product(const mVec4 & right) const;

	//以角度制表示
	MAIN_DATA_PRECISION getAngle(const mVec4& right)const;
};

class Triangle {
public:
	mVec3 a, b, c;
	__forceinline Triangle(mVec3 a, mVec3 b, mVec3 c) :a(a), b(b), c(c) {
	}
	Triangle() = default;
	__forceinline  mVec2<float> Triangle::PointIsInTriangle(const mVec3& p) const {
		//P = A + u * (C - A) + v * (B - A)       // Original equation
		//	(P - A) = u * (C - A) + v * (B - A)     // Subtract A from both sides
		//	v2 = u * v0 + v * v1                    // Substitute v0, v1, v2 for less writing
		//assert(p.z==0);
		mVec3 ta(a.x, a.y, 0.0f);
		mVec3 tb(b.x, b.y, 0.0f);
		mVec3 tc(c.x, c.y, 0.0f);
		mVec3 v0 = mVec3(tc, ta);
		mVec3 v1 = mVec3(tb, ta);
		mVec3 v2 = mVec3(p, ta);

		float deno = 1.0f / ((v0 * v0) * (v1 * v1) - (v0 * v1) * (v1 * v0));
		float u = (((v1 * v1) * (v2 * v0) - (v1 * v0) * (v2 * v1))) * deno;
		float v = (((v0 * v0) * (v2 * v1) - (v0 * v1) * (v2 * v0))) * deno;

		return { u,v };

	}
	 mVec3 Triangle::Normal();
	 mVec3 Triangle::Normal(mVec3 dir);
};

class Plane 
{
public:
	mVec3 n;		// normal
	mVec3 px;	// solution to dot(n,p-p')=0  d
	
	Plane() = default;
	Plane(mVec3 pm, mVec3 n); 
};
//int clip_with_plane(Plane c_plane, std::vector<mVec3> vert_list, std::vector<mVec3>& in_list);