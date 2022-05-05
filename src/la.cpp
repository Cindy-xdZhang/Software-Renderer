#include "la.h"

MAIN_DATA_PRECISION mVec3:: getEuclideannNorms()const {
	MAIN_DATA_PRECISION sum = 0.0;
	sum += x * x;
	sum += y * y;
	sum += z * z;
	return sqrt(sum);

}
mVec3 mVec4::tomVec3(){
	return mVec3(this->x, this->y, this->z);
}
mVec3::mVec3(mVec3 A, mVec3 B) :x(A.x - B.x), y(A.y - B.y), z(A.z - B.z) {
	//this->data = new  MAIN_DATA_PRECISION[3];
}
mVec3::mVec3(float a, float b, float c):x(a),y(b),z(c) {
}
 void mVec3:: ColorClamp() {
	this->x = this->x > 1.0f ? 1.0f : this->x;
	this->y = this->y > 1.0f ? 1.0f : this->y;
	this->z= this->z > 1.0f ? 1.0f : this->z;
	this->x = this->x < 0.00f ? 0.00f : this->x;
	this->y = this->y < 0.00f ? 0.00f : this->y;
	this->z = this->z < 0.00f ? 0.00f : this->z;
}

mVec3::mVec3(float d) : x(d), y(d), z(d) {
}
mVec3::~mVec3() {
}
mVec3 mVec3::operator +(const mVec3& right) const {
	return mVec3(this->x + right.x, this->y + right.y, this->z + right.z);
}
void mVec3::operator +=(const mVec3& right) {
	this->x += right.x;
	this->y += right.y;
	this->z += right.z;
}
void mVec3::operator -=(const mVec3& right) {
	this->x -= right.x;
	this->y -= right.y;
	this->z -= right.z;
}
mVec3 mVec3::operator -(const mVec3& right) const {
	return mVec3(this->x - right.x, this->y - right.y, this->z - right.z);
}
void mVec3::normalize() {
	MAIN_DATA_PRECISION ma = (MAIN_DATA_PRECISION)this->getEuclideannNorms();
	(*this) *= (1 / ma);
}
 mVec3  normalize(mVec3 I) {
	MAIN_DATA_PRECISION ma = (MAIN_DATA_PRECISION)I.getEuclideannNorms();
	return I * (1 / ma);
}
void mVec3::operator *=(MAIN_DATA_PRECISION scale) {
	this->x *= scale;
	this->y *= scale;
	this->z *= scale;
}

mVec3 mVec3::operator *(MAIN_DATA_PRECISION scale) const {
	return mVec3(this->x *scale, this->y *scale, this->z *scale);
}
MAIN_DATA_PRECISION mVec3::operator *(const mVec3& right) const { //向量相乘，结果为一个数，而不是向量
	MAIN_DATA_PRECISION ans = 0;
	ans += this->x* right.x;
	ans += this->y* right.y;
	ans += this->z* right.z;
	return ans;
}

mVec3 mVec3::cross_product(const mVec3 & right) const {
	mVec3 cro(0);
	cro.x = (this->y * right.z - this->z * right.y);
	cro.y = (this->z * right.x - this->x * right.z);
	cro.z = (this->x * right.y - this->y * right.x);
	return cro;
}
mVec3 mVec3::operator /(MAIN_DATA_PRECISION denominator) const {
	float scale = 1.0f / denominator;
	return mVec3(this->x *scale, this->y *scale, this->z *scale);

}
//以角度制表示
MAIN_DATA_PRECISION mVec3::getAngle(const mVec3& right)const {
	MAIN_DATA_PRECISION dot_product = (*this)*right;
	auto ma = this->getEuclideannNorms();
	auto mb = right.getEuclideannNorms();
	dot_product /= (ma*mb);
	return acos(dot_product) *180.0 / PI;
}
//datatype_eum datatype;
MAIN_DATA_PRECISION mVec4::getEuclideannNorms()const {
	MAIN_DATA_PRECISION sum = 0.0;
	sum += x * x;
	sum += y * y;
	sum += z * z;
	sum += w * w;
	return sqrt(sum);

}




mVec4::mVec4(mVec4 A, mVec4 B) {
	//this->data = new  MAIN_DATA_PRECISION[3];
	this->x = A.x - B.x;
	this->y = A.y - B.y;
	this->z = A.z - B.z;
	this->z = A.w - B.w;
}
mVec4::mVec4(float a, float b, float c, float cd) {
	this->x = a;
	this->y = b;
	this->z = c;
	this->w = cd;
}
mVec4::mVec4() {
	this->x = 0;
	this->y = 0;
	this->z = 0;
	this->w = 0;
	//this->data = NULL;
}
mVec4::mVec4(float d) {
	this->x = d;
	this->y = d;
	this->z = d;
	this->w = d;
}
mVec4::mVec4(mVec3 r, float w) {
	this->x = r.x;
	this->y = r.y;
	this->z = r.z;
	this->w = w;
}
mVec4::~mVec4() {
}
mVec4 mVec4:: operator +(const mVec4& right) const {
	return mVec4(this->x + right.x, this->y + right.y, this->z + right.z, this->w + right.w);
}
mVec4 mVec4::operator -(const mVec4& right) const {
	return mVec4(this->x - right.x, this->y - right.y, this->z - right.z, this->w - right.w);
}
void mVec4:: operator +=(const mVec4& right) {
	this->x += right.x;
	this->y += right.y;
	this->z += right.z;
	this->w += right.w;
}
void mVec4:: operator -=(const mVec4& right) {
	this->x -= right.x;
	this->y -= right.y;
	this->z -= right.z;
	this->w -= right.w;
}

void mVec4::normalize() {
	MAIN_DATA_PRECISION ma = (MAIN_DATA_PRECISION)this->getEuclideannNorms();
	(*this) *= (1 / ma);
}
void mVec4::operator *=(MAIN_DATA_PRECISION scale) {
	this->x *= scale;
	this->y *= scale;
	this->z *= scale;
	this->w *= scale;
}

mVec4 mVec4::operator *(MAIN_DATA_PRECISION scale) const {
	return mVec4(this->x *scale, this->y *scale, this->z *scale,this->w*scale);
}
MAIN_DATA_PRECISION mVec4::operator *(const mVec4& right) const { //向量相乘，结果为一个数，而不是向量
	MAIN_DATA_PRECISION ans = 0;
	ans += this->x* right.x;
	ans += this->y* right.y;
	ans += this->z* right.z;
	ans += this->w* right.w;
	return ans;
}



//以角度制表示
MAIN_DATA_PRECISION mVec4::getAngle(const mVec4& right)const {
	MAIN_DATA_PRECISION dot_product = (*this)*right;
	auto ma = this->getEuclideannNorms();
	auto mb = right.getEuclideannNorms();
	dot_product /= (ma*mb);
	return acos(dot_product) *180.0 / PI;
}
mVec4 mVec4::operator /(MAIN_DATA_PRECISION scale) const {
	return mVec4(this->x/ scale, this->y / scale, this->z / scale, this->w / scale);

}
mVec3 Triangle::Normal(mVec3 dir) {
	mVec3 L1(a, c);
	mVec3 L2(b, c);
	mVec3 out = L2.cross_product(L1);
	mVec3 out2 = L1.cross_product(L2);
	if (dir*out > 0)return out2;
	return  out;

}
mVec3 Triangle::Normal() {
	mVec3 L1(a, c);
	mVec3 L2(b, c);
	mVec3 out = L2.cross_product(L1);
	out.normalize();
	/*	mVec3 out2 = L1.cross_product(L2);
		if (dir*out > 0)return out2;*/
	return  out;

}

 

 Plane::Plane(mVec3 pm, mVec3 n) :px(pm), n(n) {
}

