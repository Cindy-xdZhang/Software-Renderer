#include "la.h"

template<typename T>
mVec3<T> mVec4<T>::tomVec3(){
	return mVec3<T>(this->x, this->y, this->z);
}

template<typename T>
mVec3<T>::mVec3<T>(mVec3<T> A, mVec3<T> B) :x(A.x - B.x), y(A.y - B.y), z(A.z - B.z) {
	//this->data = new  MAIN_DATA_PRECISION[3];
}

template<typename T>
mVec3<T>::mVec3<T>(float a, float b, float c):x(a),y(b),z(c) {
}

template<typename T>
 void mVec3<T>:: ColorClamp() {
	this->x = this->x > 1.0f ? 1.0f : this->x;
	this->y = this->y > 1.0f ? 1.0f : this->y;
	this->z= this->z > 1.0f ? 1.0f : this->z;
	this->x = this->x < 0.00f ? 0.00f : this->x;
	this->y = this->y < 0.00f ? 0.00f : this->y;
	this->z = this->z < 0.00f ? 0.00f : this->z;
}

 template<typename T>
 mVec3<T>::mVec3<T>(float d) : x(d), y(d), z(d) {
}

 template<typename T>
 mVec3<T> mVec3<T>::operator +(const mVec3<T>& right) const {
	return mVec3<T>(this->x + right.x, this->y + right.y, this->z + right.z);
}
template<typename T>
void mVec3<T>::operator +=(const mVec3<T>& right) {
	this->x += right.x;
	this->y += right.y;
	this->z += right.z;
}
template<typename T>
void mVec3<T>::operator -=(const mVec3<T>& right) {
	this->x -= right.x;
	this->y -= right.y;
	this->z -= right.z;
}
template<typename T>
mVec3<T> mVec3<T>::operator -(const mVec3<T>& right) const {
	return mVec3(this->x - right.x, this->y - right.y, this->z - right.z);
}

template<typename T>
mVec3<T> mVec3<T>::cross_product(const mVec3<T>& right) const {
	mVec3<T> cro(0);
	cro.x = (this->y * right.z - this->z * right.y);
	cro.y = (this->z * right.x - this->x * right.z);
	cro.z = (this->x * right.y - this->y * right.x);
	return cro;
}



template<typename T>
mVec4<T>::mVec4<T>(mVec4<T>A, mVec4<T> B) {
	//this->data = new  MAIN_DATA_PRECISION[3];
	this->x = A.x - B.x;
	this->y = A.y - B.y;
	this->z = A.z - B.z;
	this->z = A.w - B.w;
}

template<typename T>
mVec4<T>::mVec4<T>(float a, float b, float c, float cd) {
	this->x = a;
	this->y = b;
	this->z = c;
	this->w = cd;
}
template<typename T>
mVec4<T>::mVec4<T>() {
	this->x = 0;
	this->y = 0;
	this->z = 0;
	this->w = 0;
	//this->data = NULL;
}
template<typename T>
mVec4<T>::mVec4<T>(float d) {
	this->x = d;
	this->y = d;
	this->z = d;
	this->w = d;
}
template<typename T>
mVec4<T>::mVec4<T>(mVec3<T> r, float w) {
	this->x = r.x;
	this->y = r.y;
	this->z = r.z;
	this->w = w;
}
template<typename T>
mVec4<T>::~mVec4<T>() {
}
template<typename T>
mVec4<T> mVec4<T>:: operator +(const mVec4<T>& right) const {
	return mVec4<T>(this->x + right.x, this->y + right.y, this->z + right.z, this->w + right.w);
}
template<typename T>
mVec4<T> mVec4<T>::operator -(const mVec4<T>& right) const {
	return mVec4<T>(this->x - right.x, this->y - right.y, this->z - right.z, this->w - right.w);
}
template<typename T>
void mVec4<T>:: operator +=(const mVec4<T>& right) {
	this->x += right.x;
	this->y += right.y;
	this->z += right.z;
	this->w += right.w;
}
template<typename T>
void mVec4<T>:: operator -=(const mVec4<T>& right) {
	this->x -= right.x;
	this->y -= right.y;
	this->z -= right.z;
	this->w -= right.w;
}
template<typename T>
void mVec4<T>::normalize() {
	T ma = (T)this->getEuclideannNorms();
	(*this) *= (1 / ma);
}
template<typename T>
void mVec4<T>::operator *=(T scale) {
	this->x *= scale;
	this->y *= scale;
	this->z *= scale;
	this->w *= scale;
}
template<typename T>
mVec4<T> mVec4<T>::operator *(T scale) const {
	return mVec4<T>(this->x *scale, this->y *scale, this->z *scale,this->w*scale);
}

template<typename T>
T mVec4<T>::operator *(const mVec4<T>& right) const { //向量相乘，结果为一个数，而不是向量
	T ans = 0;
	ans += this->x* right.x;
	ans += this->y* right.y;
	ans += this->z* right.z;
	ans += this->w* right.w;
	return ans;
}




template<typename T>
mVec4<T> mVec4<T>::operator /(T scale) const {
	return mVec4<T>(this->x/ scale, this->y / scale, this->z / scale, this->w / scale);

}


mVec3f Triangle::Normal(mVec3f dir) {
	mVec3f L1(a, c);
	mVec3f L2(b, c);
	mVec3f out = L2.cross_product(L1);
	mVec3f out2 = L1.cross_product(L2);
	if (dir*out > 0)return out2;
	return  out;

}
mVec3f Triangle::Normal() {
	mVec3f L1(a, c);
	mVec3f L2(b, c);
	mVec3f out = L2.cross_product(L1);
	out.normalize();
	/*	mVec3f out2 = L1.cross_product(L2);
		if (dir*out > 0)return out2;*/
	return  out;

}

 

 Plane::Plane(mVec3f pm, mVec3f n) :px(pm), n(n) {
}


 template mVec3<float>;
 template mVec4<float>;