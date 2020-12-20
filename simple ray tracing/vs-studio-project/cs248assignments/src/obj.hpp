#pragma   once  
#ifndef OBJ
#define OBJ
#include "la.h"
#include "matrix.h"
#define MAX(a,b) (a)>(b)?(a):(b)
#define MIN(a,b) (a)>(b)?(b):(a)
#define MAX_OBJECTS 10
enum datatype_eum {
	XD_INT = 1,
	XD_FLOAT = 2,
	XD_DOUBLE = 3,
	XD_SPHERE = 4,
	XD_CONE = 5,
	XD_PLANE = 6,
	XD_PARAMRID = 7,
	XD_AABB = 8
};
struct intersectInfo {
	float time;
	mVec3 n;
};
class Ray {
public:
	mVec3 ori;
	mVec3 dir;//ori+t*dir;
	inline Ray(mVec3 ori, mVec3 Veco);
	inline mVec3 point(float time) const;
	inline Ray();
	inline ~Ray();
};
inline Ray::Ray(mVec3 ori, mVec3 Veco) {
	this->ori = ori;
	this->dir = Veco;
	this->dir.normalize();
}
inline mVec3 Ray::point(float time) const {
	mVec3 bt = mVec3(ori.x + time * this->dir.x, \
		ori.y + time * this->dir.y, ori.z + time * this->dir.z);
	return bt;
}
inline Ray::~Ray() {
}
inline Ray::Ray() {
}
typedef struct FloatRGB {
	float r;
	float g;
	float b;
	inline FloatRGB(float r = 0, float g = 0, float b = 0) :r(r), g(g), b(b) {

	}
	template < typename T>
	inline FloatRGB operator *(T scale) const {
		FloatRGB c = FloatRGB();
		c.r = this->r*scale;
		c.g = this->g*scale;
		c.b = this->b*scale;
		return c;
	}
	inline FloatRGB operator +(const FloatRGB scale) const {
		FloatRGB c = FloatRGB(0, 0, 0);
		c.r = this->r + scale.r;
		c.g = this->g + scale.g;
		c.b = this->b + scale.b;
		return c;
	}
	inline bool operator >(float scale) const {
		return (this->r > scale || this->b > scale || this->g > scale);
	}
	inline bool operator >=(float scale) const {
		return (this->r >= scale || this->b >= scale || this->g >= scale);
	}
	inline void Clamp() {
		this->r = this->r > 1.000000f ? 1.00f : this->r;
		this->g = this->g > 1.000000f ? 1.00f : this->g;
		this->b = this->b > 1.000000f ? 1.00f : this->b;
		this->r = this->r < 0.00f ? 0.00f : this->r;
		this->g = this->g < 0.00f ? 0.00f : this->g;
		this->b = this->b < 0.00f ? 0.00f : this->b;
	}
} FloatRGB;
class ShadingMaterial {
public:
	FloatRGB Ka;
	FloatRGB Kd;
	FloatRGB Ks;
	inline ShadingMaterial() {
	} 
	inline ShadingMaterial(FloatRGB a, FloatRGB d, FloatRGB s) {
		this->Ka = a;
		this->Kd = d;
		this->Ks = s;
	}
	inline ~ShadingMaterial() {
	}
};

//卡了一个bug 如果在surface 析构函数里调用delete pModelMatrix
//会卡进这个问题里：https://community.intel.com/t5/Intel-C-Compiler/Window-x64-17-0-1-compiler-Visual-Studio-2015-and-std-list/td-p/1111088?profile.language=zh-CN
class Surface {
protected:
	ShadingMaterial shma;	// material
	ShadingMaterial reset_shma;	// material
	void updateIModelMatrix() {
		*this->pIModleMatrix = InverseModelMatrix();
	}
	void updateNormalMatrix() {
		*this->pNormalMatrix = this->pIModleMatrix->Transpose();
	}
	Matrix4 InverseModelMatrix() {
		Matrix4 InvScale = eye(4);
		InvScale.p[0][0] = 1.0f / (pScaleMatrix->p[0][0]);
		InvScale.p[1][1] = 1.0f / (pScaleMatrix->p[1][1]);
		InvScale.p[2][2] = 1.0f / (pScaleMatrix->p[2][2]);
		Matrix4 InvR = (*pRotateMatrix).Transpose();

		Matrix4 InvT = eye(4);
		InvT.p[0][3] = -pTransMatrix->p[0][3];
		InvT.p[1][3] = -pTransMatrix->p[1][3];
		InvT.p[2][3] = -pTransMatrix->p[2][3];
		return InvScale * InvR * InvT;
	}
public:
	Matrix4* pTransMatrix=NULL;
	Matrix4* pScaleMatrix = NULL;
	Matrix4* pRotateMatrix = NULL;
	Matrix4* pIModleMatrix = NULL;
	Matrix4* pNormalMatrix = NULL;
	inline void updateMaterial(FloatRGB a, FloatRGB d, FloatRGB s) {
		this->shma = ShadingMaterial(a, d, s);
	}
	inline void boudleMaterial(FloatRGB a, FloatRGB d, FloatRGB s) {
		this->reset_shma = ShadingMaterial(a, d, s);
		this->shma = ShadingMaterial(a, d, s);
	}
	inline void ResetMaterial() {
		this->shma = ShadingMaterial(reset_shma.Ka, reset_shma.Kd, reset_shma.Ks);
	}
	inline ShadingMaterial Material() {
		return this->shma;
	}
	void initModelMatrix() {
		pTransMatrix = new Matrix4(eye(4));
		pScaleMatrix = new Matrix4(eye(4));
		pRotateMatrix = new Matrix4(eye(4));
		pIModleMatrix = new Matrix4(eye(4));
		pNormalMatrix = new Matrix4(eye(4));
		//pNormalMatrix = new Matrix4(eye(4));
	}
	void ResetModelMatrix() {
		*pTransMatrix =  Matrix4(eye(4));
		*pScaleMatrix =  Matrix4(eye(4));
		*pRotateMatrix =  Matrix4(eye(4));
		*pIModleMatrix =  Matrix4(eye(4));
		*pNormalMatrix =  Matrix4(eye(4));
	}
	void updateTranslateMatrix(Matrix4& right) {
		*pTransMatrix = right  * (*pTransMatrix);
		updateIModelMatrix();
		updateNormalMatrix();
	}
	void updateRotateMatrix(Matrix4& right) {
		*pRotateMatrix = right * (*pRotateMatrix);
		updateIModelMatrix();
		updateNormalMatrix();
	}
	void updateScaleMatrix(Matrix4& right) {
		*pScaleMatrix = right * (*pScaleMatrix);
		updateIModelMatrix();
		updateNormalMatrix();
	}
	Matrix4 IverseModelMatrix() { return *(this->pIModleMatrix); }
	void deleteModelMatrix() {
		delete pTransMatrix;
		delete pScaleMatrix;
		delete pRotateMatrix;
	}
	Surface() = default;
};
class Plane :public Surface
{
public:
	mVec3 pm;	// solution to dot(n,p-p')=0  d
	mVec3 n;		// normal
	Plane() = default;
	inline Plane(mVec3 pm, mVec3 n) :pm(pm), n(n) {
		initModelMatrix();
	}
};
////aabb:axis aligned box
class Cubid :public Surface
{
public:
	mVec3 Center;
	Plane Faces[6];
	float h, w, l;
	Cubid() = default;
	inline Cubid(mVec3 pm, float h = 1.0f, float w = 1.0f, float l = 1.0f) :Center(pm), h(h), w(w), l(l) {
		Faces[0] = Plane(Center + (mVec3(0, 0, -1))*0.5*l, mVec3(0, 0, -1));
		Faces[1] = Plane(Center + (mVec3(0, 0, 1))*0.5*l, mVec3(0, 0, 1));
		Faces[2] = Plane(Center + (mVec3(0, 1, 0))*0.5*w, mVec3(0, 1, 0));
		Faces[3] = Plane(Center + (mVec3(0, -1, 0))*0.5*w, mVec3(0, -1, 0));
		Faces[4] = Plane(Center + (mVec3(1, 0, 0))*0.5*h, mVec3(1, 0, 0));
		Faces[5] = Plane(Center + (mVec3(-1, 0, 0))*0.5*h, mVec3(-1, 0, 0));
		initModelMatrix();
	}
	inline Cubid(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax) {
		Center = mVec3((xmax + xmin) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2);
		float h = xmax - xmin;
		float w = ymax - ymin;
		float l = zmax - zmin;
		Faces[0] = Plane(Center + (mVec3(0, 0, -1))*0.5*l, mVec3(0, 0, -1));
		Faces[1] = Plane(Center + (mVec3(0, 0, 1))*0.5*l, mVec3(0, 0, 1));
		Faces[2] = Plane(Center + (mVec3(0, 1, 0))*0.5*w, mVec3(0, 1, 0));
		Faces[3] = Plane(Center + (mVec3(0, -1, 0))*0.5*w, mVec3(0, -1, 0));
		Faces[4] = Plane(Center + (mVec3(1, 0, 0))*0.5*h, mVec3(1, 0, 0));
		Faces[5] = Plane(Center + (mVec3(-1, 0, 0))*0.5*h, mVec3(-1, 0, 0));
		initModelMatrix();
	}
};

class Triangle {
public:
	mVec3 Points[3];
	inline Triangle(mVec3 a, mVec3 b, mVec3 c) {
		Points[0] = a;
		Points[1] = b;
		Points[2] = c;
	}
	inline Triangle() {
	}
	inline float Square() {
		mVec3 L1(Points[0], Points[1]);
		mVec3 L2(Points[0], Points[2]);
		float S1 = 0.5 *L1.getEuclideannNorms()*L2.getEuclideannNorms()*sin(L1.getAngle(L2) / 180 * PI);
		return S1;
	}

	inline mVec3 Normal(mVec3 dir) {
		mVec3 L1(Points[0], Points[2]);
		mVec3 L2(Points[1], Points[2]);
		mVec3 out = L2.cross_product(L1);
		mVec3 out2 = L1.cross_product(L2);
		if (dir*out > 0)return out2;
		return  out;

	}


};
class Sphere :public Surface {
public:
	mVec3 center;
	float R;
	float scalex, scaley, scalez;//a2 b2  c2
	//Cubid AABB;
	inline Sphere( float R = 1, float scalex = 1.0f, float scaley = 1.0f, float scalez = 1.0f) : scalex(scalex), scaley(scaley), scalez(scalez) {
		this->center = mVec3(0,0,0);
		this->R = R;
		initModelMatrix();

		//for aabb
	/*	float xmin = center.x-scalex* R;
		float xmax = center.x + scalex * R;
		float ymin = center.y - scaley * R;
		float ymax = center.y +scaley * R;
		float zmin = center.z - scalez * R;
		float zmax = center.z +scalez * R;
		AABB = Cubid(xmin, xmax, ymin, ymax, zmin, zmax);*/
	}

	inline Sphere(mVec3 ori, float R = 1, float scalex = 1.0f, float scaley = 1.0f, float scalez = 1.0f) :scalex(scalex), scaley(scaley), scalez(scalez) {
		this->center = ori;
		this->R = R;
		initModelMatrix();
	}
	inline mVec3 Normal(mVec3 ptx) {
		mVec3 bt = mVec3(ptx, this->center);
		bt.x /= scalex;
		bt.y /= scaley;
		bt.z /= scalez;
		return bt;
	}
};
class Cone :public Surface {
public:
	float cosa;	// half cone angle
	float h;	// height
	mVec3  c;		// tip position
	mVec3   v;		// axis
	inline Cone() {
		c= mVec3(0, 0, 0);
	    v = mVec3(1, 0, 0);
		h = 1;
		cosa = 0.95;
		initModelMatrix();
	}
	//deprecate in the future
	inline Cone(mVec3  c, mVec3   iv, float h = 2, float cosa = 0.95) : c(c), cosa(cosa), h(h) {
		iv.normalize();
		v = iv;
		initModelMatrix();
	}
	inline mVec3  Normal(mVec3  ptx) {
		mVec3  cp = mVec3(ptx, this->c);
		float h = (cp*this->v);
		//if (h < 0. || h > s.h) return -100;
		mVec3  n = cp * ((this->v* cp) / (cp*cp)) - this->v;
		n.normalize();
		return n;
	}
};


class Triangular_pyramid :public Surface {
public:
	mVec3 vertex[4];
	Triangle tris[4];
	//Cubid AABB;
	inline Triangular_pyramid(mVec3 A, mVec3 B, mVec3 C, mVec3 D){
		vertex[0] = A;
		vertex[1] = B;
		vertex[2] = C;
		vertex[3] = D;
		//float xmin = MIN(MIN(C.x,MIN(A.x, B.x)), D.x);
		//float xmax = MAX(MAX(C.x, MAX(A.x, B.x)), D.x);
		//float ymin = MIN(MIN(C.y, MIN(A.y, B.y)), D.y);
		//float ymax = MAX(MAX(C.y, MAX(A.y, B.y)), D.y);
		//float zmin = MIN(MIN(C.z, MIN(A.z, B.z)), D.z);
		//float zmax = MAX(MAX(C.z, MAX(A.z, B.z)), D.z);
		//AABB = Cubid(xmin, xmax, ymin, ymax, zmin, zmax);

		tris[0] = Triangle(A, B, C);
		tris[1] = Triangle(A, B, D);
		tris[2] = Triangle(A, C, D);
		tris[3] = Triangle(B, C, D);

		initModelMatrix();

	}
};



inline float intersection(Plane& p, Ray& Ry)
{
	float dotn_d = p.n*Ry.dir;
	if (dotn_d == 0.)return  -100;
	//if (dotn_d > 0.) return -100;//!!!!理解为啥是大于0就不要了，大于零代表光从面下方穿过故不要了，小于0则留下

	float t = mVec3(p.pm, Ry.ori)* p.n / dotn_d;
	/*if (t > 0.)	return t;
	return  -100;*/
	return t;//t<0 外部函数自动判定不相交
}
//aabb和光线相交
inline intersectInfo intersection(Cubid AABB, Ray& cRay) {
	float time[6] = { 0 };
	time[0] = intersection(AABB.Faces[0], cRay);
	time[1] = intersection(AABB.Faces[1], cRay);
	time[2] = intersection(AABB.Faces[2], cRay);
	time[3] = intersection(AABB.Faces[3], cRay);
	time[4] = intersection(AABB.Faces[4], cRay);
	time[5] = intersection(AABB.Faces[5], cRay);


	float tmin1 = MIN(time[0], time[1]);
	float tmax1 = MAX(time[0], time[1]);
	float tmin2 = MIN(time[2], time[3]);
	float tmax2 = MAX(time[2], time[3]);
	float tmin3 = MIN(time[4], time[5]);
	float tmax3 = MAX(time[4], time[5]);
	float tmin = MAX(tmin3, MAX(tmin1, tmin2));
	float tmax = MIN(tmax1, MIN(tmax2, tmax3));

	if (tmin < tmax) {
		for (int i = 0; i < 6; i++)
			if (tmin == time[i])
				return { tmin, AABB.Faces[i].n };
	}
	else return{ -100,mVec3(0) };
	//float tmin1 = MIN(time[0], time[1]); 
	//float tmin2 = MIN(time[2], time[3]);
	//float tmin3 = MIN(time[4], time[5]);
	//float tmin = -100;
	//int faceId;
	//if (tmin < tmin1) {
	//	tmin = tmin1;
	//}
	//if (tmin < tmin1) {
	//	tmin = tmin1;
	//}
	//if (tmin < tmin1) {
	//	tmin = tmin1;
	//}
	//
}

inline float intersection(Triangle& tri, Ray& Ray) {
	//相交则返回t大于零，否则返回-100
	//如果平行则省略，视为不相交
	//相交于边 （0，0，0）+t（0,9,0） 与 （1,2,0）（0,2，2）（0，2,-1） ,视为相交
	//光线从三角形出发向外射（0，0，0）+t（0,9,0）与 （1,0,0）（-1,0，-2）（-1，0,4）,或相较于顶点视为不相交
	mVec3 L1(tri.Points[0], tri.Points[1]);
	mVec3 L2(tri.Points[0], tri.Points[2]);
	mVec3 Normal = L1.cross_product(L2);
	if (Normal*Ray.dir < FLOAT_ZERO && Normal*Ray.dir*-1 < FLOAT_ZERO)return -100;
	//x+t * xd=xa+B(xb-xa)+G(xc-xa)
	//y+t * yd=ya+B(yb-ya)+G(yc-ya)
	//z+t * zd=za+B(zb-za)+G(zc-za)
	//=>D:                              Dt:                  Dg:                 Db:
	// B(xb-xa)+G(xc-xa)-t * xd =x-xa  (xb-xa) (xc-xa) x-xa  (xb-xa) (x-xa) -xd  (x-xa) (xc-xa) -xd 
	//B(yb-ya)+G(yc-ya)-t * yd=y-ya    (yb-ya) (yc-ya) y-ya  (yb-ya) (y-ya) -yd  (y-ya) (yc-ya) -yd
	//B(zb-za)+G(zc-za)-t * zd =z-za   (zb-za) (zc-za) z-za  (zb-za) (z-za) -zd  (z-za) (zc-za) -zd    
	float xa = (float)tri.Points[0].x;
	float xb = (float)tri.Points[1].x;
	float xc = (float)tri.Points[2].x;
	float ya = (float)tri.Points[0].y;
	float yb = (float)tri.Points[1].y;
	float yc = (float)tri.Points[2].y;
	float za = (float)tri.Points[0].z;
	float zb = (float)tri.Points[1].z;
	float zc = (float)tri.Points[2].z;
	float xd = (float)Ray.dir.x;
	float yd = (float)Ray.dir.y;
	float zd = (float)Ray.dir.z;
	float x = (float)Ray.ori.x;
	float y = (float)Ray.ori.y;
	float z = (float)Ray.ori.z;
	//compute D(系数矩阵行列式)=	
   //(xb - xa)[(yc-ya) * (-zd) +yd*(zc-za)]-(xc - xa)[ (yb-ya) * -zd +yd*(zb-za)]- xd*[ (yb-ya) * (zc-za) -(zb-za)*(yc-ya)]
	float D = (xb - xa)*((yc - ya) * (-zd) + yd * (zc - za)) - \
		(xc - xa)*((yb - ya) * (-zd) + yd * (zb - za)) - xd * ((yb - ya) * (zc - za) - (zb - za)*(yc - ya));
	//compute t=Dt/D
	float Dt = (xb - xa)*((yc - ya) * (z - za) - (y - ya) * (zc - za)) - \
		(xc - xa)*((yb - ya) * (z - za) - (y - ya) * (zb - za)) + (x - xa) * ((yb - ya) * (zc - za) - (zb - za)*(yc - ya));
	float t = Dt / D;
	if (t < 0.)return -100;
	//compute gama=Dg/D
	float Dg = (xb - xa)*((y - ya) * (-zd) + yd * (z - za)) - \
		(x - xa)*((yb - ya) * (-zd) + yd * (zb - za)) - xd * ((yb - ya) * (z - za) - (zb - za)*(y - ya));
	float gama = Dg / D;
	if (gama < 0. || gama >1)return -100;
	//compute beta=Db/D
	float Db = (x - xa)*((yc - ya) * (-zd) + yd * (zc - za)) - \
		(xc - xa)*((y - ya) * (-zd) + yd * (z - za)) - xd * ((y - ya) * (zc - za) - (z - za)*(yc - ya));
	float beta = Db / D;
	if (beta < 0. || beta >1)return -100;
	if (beta + gama > 1)return -100;
	return t;


}

inline intersectInfo intersection(Triangular_pyramid& tri, Ray& Ray) {

	//intersectInfo AABBtime= intersection(tri.AABB, Ray);
	//if (AABBtime.time < 0) return{ -100.0f ,mVec3(0)};

	float time = (float)INT_MAX;
	int id = -1;
	for (int i = 0; i < 4; i++) {
		float t = intersection(tri.tris[i], Ray); 
		if (t > 0 && time > t) {
			time = t;
			id = i;
		}
	}
	return { time!= (float)INT_MAX? time:-100.0f,tri.tris[id].Normal(Ray.dir) };

}


inline float intersection(Cone& s, Ray& Ray) {


		mVec3 co = mVec3(Ray.ori, s.c);
		float a = (Ray.dir* s.v)*(Ray.dir* s.v) - s.cosa*s.cosa;

		float b = 2. *((Ray.dir* s.v)*(co*s.v) - (Ray.dir* co)*  s.cosa*s.cosa);
		float c = (co*s.v)*(co* s.v) - (co* co)*s.cosa*s.cosa;

		float det = b * b - 4.*a*c;
		if (det < 0.) return -100;

		det = sqrt(det);
		float t1 = (-b - det) / (2. * a);
		float t2 = (-b + det) / (2. * a);
		
		float t= t1 > t2 ? t2 : t1;
		mVec3 cp = mVec3(Ray.point(t), s.c);
		float h = (cp*s.v);
		if (h < 0. || h > s.h) return -100;

		//float m1 = t1 > t2 ? t2 : t1;
		//float m2 = t1 > t2 ? t1 : t2;//大根
		//float t;
		//if (m2 <= 0)return -100;
		//else if (m1 < 0) t = m2;//光源在球内部
		//else t = m1;//光在求外部，优先返回小根

		
		return t;

	
}

inline float intersection(Sphere& sphere, Ray& Ray) {

	//intersectInfo AABBtime = intersection(sphere.AABB, Ray);
	//if (AABBtime.time < 0) return -100.0f ;
	Ray.dir.x /= sphere.scalex;
	Ray.dir.y /= sphere.scaley;
	Ray.dir.z /= sphere.scalez;
	Ray.ori.x /= sphere.scalex;
	Ray.ori.y /= sphere.scaley;
	Ray.ori.z /= sphere.scalez;
	//只返回近的交点，优先返回小根,相切则省略
	float a = Ray.dir* Ray.dir;
	mVec3 E_C(Ray.ori, sphere.center);
	float b = Ray.dir*E_C;//*2.0f
	float c = E_C * E_C - sphere.R*sphere.R;

	float delta = b * b - a * c;
	if (delta < 0.)return -100;
	float t1 = (sqrt(delta) - b) / a;
	float t2 = (- b - sqrt(delta)) / a;

	return t1 > t2 ? t2 : t1;
	//float m1 = t1 > t2 ? t2 : t1;
	//float m2 = t1 > t2 ? t1 : t2;//大根

	//if (m2 <=0)return -100;
	//else if (m1 < 0) return m2;//光源在球内部
	//else return m1;//光在求外部，优先返回小根


}


class Object_List {
	int object_count;
	unsigned char* pointers[MAX_OBJECTS];
	unsigned char type[MAX_OBJECTS];
public:
	Object_List() {
		object_count = 0;
	}
	int Len() {
		return object_count;
	}
	unsigned char* GetPointer(int id) {
		return pointers[id];
	}
	unsigned char  GetType(int id) {
		return type[id];
	}
	void AddObj(Sphere& spo) {
		if (object_count == MAX_OBJECTS)throw("AddObj Fail for out of momeory.");
		else {
			pointers[object_count] = (unsigned char*)&spo;
			type[object_count] = XD_SPHERE;
			object_count++;
		}
	}
	void AddObj(Triangular_pyramid& spo) {
		if (object_count == MAX_OBJECTS)throw("AddObj Fail for out of momeory.");
		else {
			pointers[object_count] = (unsigned char*)&spo;
			type[object_count] = XD_PARAMRID;
			object_count++;
		}
	}
	void AddObj(Plane& spo) {
		if (object_count == MAX_OBJECTS)throw("AddObj Fail for out of momeory.");
		else {
			pointers[object_count] = (unsigned char*)&spo;
			type[object_count] = XD_PLANE;
			object_count++;
		}
	}
	void AddObj(Cone& spo) {
		if (object_count == MAX_OBJECTS)throw("AddObj Fail for out of momeory.");
		else {
			pointers[object_count] = (unsigned char*)&spo;
			type[object_count] = XD_CONE;
			object_count++;
		}
	}

	void AddObj(Cubid& spo) {
		if (object_count == MAX_OBJECTS)throw("AddObj Fail for out of momeory.");
		else {
			pointers[object_count] = (unsigned char*)&spo;
			type[object_count] = XD_AABB;
			object_count++;
		}
	}
	void DelObj(int id) {
		if (id >= object_count - 1)throw("DelObj Fail.");
		else {
			delete pointers[object_count];
			for (int p = id; p < object_count - 1; p++) {
				pointers[p] = pointers[p + 1];
				type[p] = type[p + 1];
			}
			object_count--;
		}
	}
	bool OccludeByOther(int skip, Ray  cRay) {
		for (int u = 0; u < object_count; u++)
		{
			if (u == skip)continue;
			float time = -100;
			auto Cobjp = pointers[u];
			switch (type[u]) {
			case (XD_SPHERE):
			{Sphere Cobj = *((Sphere*)Cobjp);
			Matrix4 M = *Cobj.pIModleMatrix;
			Ray tRay = Ray((M*mVec4(cRay.ori, 1.0f)).tomVec3(), (M*mVec4(cRay.dir, 0.0f)).tomVec3());
			time = intersection(Cobj, tRay);
			break; }
			case (XD_PLANE):
			{Plane Cobj = *((Plane*)Cobjp);
			Matrix4 M = *Cobj.pIModleMatrix;
			Ray tRay = Ray((M*mVec4(cRay.ori, 1.0f)).tomVec3(), (M*mVec4(cRay.dir, 0.0f)).tomVec3());
			time = intersection(Cobj, tRay);
			break; }
			case (XD_CONE):
			{Cone Cobj = *((Cone*)Cobjp);
			Matrix4 M = *Cobj.pIModleMatrix;
			Ray tRay = Ray((M*mVec4(cRay.ori, 1.0f)).tomVec3(), (M*mVec4(cRay.dir, 0.0f)).tomVec3());
			time = intersection(Cobj, tRay);
			break; }
			case (XD_PARAMRID):
			{Triangular_pyramid Cobj = *((Triangular_pyramid*)Cobjp);
			Matrix4 M = *Cobj.pIModleMatrix;
			Ray tRay = Ray((M*mVec4(cRay.ori, 1.0f)).tomVec3(), (M*mVec4(cRay.dir, 0.0f)).tomVec3());
			intersectInfo time_n = intersection(Cobj, tRay);
			time = time_n.time;
			break; }
			}
			if (time > 0)return true;
		}
		return false;
	}

};













#endif