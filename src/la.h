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

template<typename T = float>
class mVec3 {
public:
	T x;
	T y;
	T z;
	//datatype_eum datatype;
	T getEuclideannNorms()const{
		T sum = 0.0;
		sum += x * x;
		sum += y * y;
		sum += z * z;
		return sqrt(sum);

	}

	mVec3(mVec3 A, mVec3 B);
	mVec3(float a, float b, float c);
	inline mVec3(mVec3&& other):x(other.x), y(other.y),z(other.z) {
	}
	inline mVec3(const mVec3& other) : x(other.x), y(other.y), z(other.z) {
	}
	inline mVec3& operator =(const mVec3&) = default;
	inline mVec3& operator =( mVec3&&) = default;

	mVec3()=default;
	inline mVec3(float d);
	~mVec3() = default;
	inline mVec3 operator +(const mVec3& right) const;
	inline void operator +=(const mVec3& right);
	inline void operator -=(const mVec3& right);
	inline mVec3 operator -(const mVec3& right) const;

	inline void mVec3<T>::operator *=(T scale);
	inline mVec3<T> mVec3<T>::operator *(T scale) const;
	inline T operator *(const mVec3<T>& right) const;
	inline mVec3<T> operator/(T denominator) const;




	
	inline void normalize() {
		T ma = static_cast<T>(this->getEuclideannNorms());
		(*this) *= (1 / ma);
	}



	mVec3 cross_product(const mVec3 & right) const;
	
	T getAngle(const mVec3& right)const;
	void ColorClamp();
	
};




template<typename T = float>
inline  mVec3<T>  normalize(mVec3<T> Iu) {
	T ma = static_cast<T>(Iu.getEuclideannNorms());
	mVec3<T> tmp = Iu;
	tmp *= (1 / ma);
	return tmp;
}






template<typename T = float>
class mVec4 {
public:
	T x;
	T y;
	T z;
	T w;
	//datatype_eum datatype
	T getEuclideannNorms()const {
		T sum = 0.0;
		sum += x * x;
		sum += y * y;
		sum += z * z;
		sum += w * w;
		return sqrt(sum);

	}

	mVec3<T> tomVec3();
	mVec4(mVec4 A, mVec4 B);
	mVec4(float a, float b, float c, float cd);
	mVec4();
	mVec4(float d);
	mVec4(mVec3<T> r, float w);
	~mVec4();
	mVec4 operator +(const mVec4& right) const;
	void operator +=(const mVec4& right);
	void operator -=(const mVec4& right);
	mVec4 operator -(const mVec4& right) const;
	void normalize();

	void operator *=(T scale);
	mVec4 operator *(T scale) const;
	mVec4 operator /(T scale) const;
	T operator *(const mVec4& right) const;

	//mVec4 cross_product(const mVec4 & right) const;

	T getAngle(const mVec4& right)const{
		T dot_product = (*this) * right;
		auto ma = this->getEuclideannNorms();
		auto mb = right.getEuclideannNorms();
		dot_product /= (ma * mb);
		return acos(dot_product) * 180.0 / PI;
	}


};



using mVec3f = mVec3<float>;
using mVec4f = mVec4<float>;


class Triangle {
public:
	mVec3f a, b, c;
	__forceinline Triangle(mVec3f a, mVec3f b, mVec3f c) :a(a), b(b), c(c) {
	}
	Triangle() = default;
	__forceinline  mVec2<float> Triangle::PointIsInTriangle(const mVec3f & p) const {
		//P = A + u * (C - A) + v * (B - A)       // Original equation
		//	(P - A) = u * (C - A) + v * (B - A)     // Subtract A from both sides
		//	v2 = u * v0 + v * v1                    // Substitute v0, v1, v2 for less writing
		//assert(p.z==0);
		mVec3f ta(a.x, a.y, 0.0f);
		mVec3f tb(b.x, b.y, 0.0f);
		mVec3f tc(c.x, c.y, 0.0f);
		mVec3f v0 = mVec3f(tc, ta);
		mVec3f v1 = mVec3f(tb, ta);
		mVec3f v2 = mVec3f(p, ta);

		float deno = 1.0f / ((v0 * v0) * (v1 * v1) - (v0 * v1) * (v1 * v0));
		float u = (((v1 * v1) * (v2 * v0) - (v1 * v0) * (v2 * v1))) * deno;
		float v = (((v0 * v0) * (v2 * v1) - (v0 * v1) * (v2 * v0))) * deno;

		return { u,v };

	}
	mVec3f Triangle::Normal();
	mVec3f Triangle::Normal(mVec3f dir);
};

class Plane
{
public:
	mVec3f n;		// normal
	mVec3f px;	// solution to dot(n,p-p')=0  d

	Plane() = default;
	Plane(mVec3f pm, mVec3f n);
};
//int clip_with_plane(Plane c_plane, std::vector<mVec3> vert_list, std::vector<mVec3>& in_list);