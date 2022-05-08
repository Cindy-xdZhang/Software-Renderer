#include "la.h"

template<typename T>
mVec3<T> mVec4<T>::tomVec3() const{
	return mVec3<T>(this->x, this->y, this->z);
}
template<typename T /*= float*/>
 mVec3<T> mVec4<T>::HomoCordinates2InHomoVec3()
{
	mVec4<T>tmp = (* this) / this->w;
	return {tmp.x, tmp.y, tmp.z};
}



template<typename T>
mVec3<T>::mVec3(mVec3<T> A, mVec3<T> B) :x(A.x - B.x), y(A.y - B.y), z(A.z - B.z) {
	//this->data = new  MAIN_DATA_PRECISION[3];
}

template<typename T>
mVec3<T>::mVec3(T a, T b, T c):x(a),y(b),z(c) {
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
 mVec3<T>::mVec3(T d) : x(d), y(d), z(d) {
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
void mVec3<T>::operator *=(T scale) {
	this->x *= scale;
	this->y *= scale;
	this->z *= scale;
}
template<typename T>
mVec3<T> mVec3<T>::operator *(T scale) const {
	return mVec3<T>(this->x * scale, this->y * scale, this->z * scale);
}

template<typename T>
T mVec3<T>::operator *(const mVec3<T>& right) const {
	T ans = 0;
	ans += this->x * right.x;
	ans += this->y * right.y;
	ans += this->z * right.z;
	return ans;
}


template<typename T /*= float*/>
mVec3<T> mVec3<T>::operator/(T denominator) const
{
	T scale = static_cast<T>(1.0f / denominator);
	return mVec3<T>(this->x * scale, this->y * scale, this->z * scale);
}










template<typename T>
mVec4<T>::mVec4(mVec4<T>A, mVec4<T> B): 
	x(A.x - B.x), y(A.y - B.y), z(A.z - B.z), w(A.w - B.w)
{
}

template<typename T>
mVec4<T>::mVec4(float a, float b, float c, float cd) {
	this->x = a;
	this->y = b;
	this->z = c;
	this->w = cd;
}
template<typename T>
mVec4<T>::mVec4():x(0),y(0),z(0),w(0) {
}
template<typename T>
mVec4<T>::mVec4(float d) {
	this->x = d;
	this->y = d;
	this->z = d;
	this->w = d;
}
template<typename T>
mVec4<T>::mVec4(mVec3<T> r, float w) {
	this->x = r.x;
	this->y = r.y;
	this->z = r.z;
	this->w = w;
}
template<typename T>
mVec4<T>::~mVec4() {
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


template <typename T/*=float*/>
T Plane<T>::cal_project_distance(mVec3<T> point) const
{
	return (point - this->ptx) * this->n;
}
template <typename T/*=float*/>
T Plane<T>::cal_intersectRatio(mVec3f ptx1, mVec3f ptx2) const{
	T d1 = cal_project_distance(ptx1);
	T d2 = cal_project_distance(ptx2);

	return d1 /(d1 -d2);
};
 
 template mVec3<int>;
 template mVec3<float>;
 template mVec4<float>;

 template Plane<float>;