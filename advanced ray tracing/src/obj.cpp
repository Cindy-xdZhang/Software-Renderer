#include "obj.h"


bool TureisUsePredefinedNormalFalseIsUseAverageNormal = false;

 Ray::Ray(mVec3 ori, mVec3 Veco) {
	this->ori = ori;
	this->dir = Veco;
	//this->dir.normalize();
}

 mVec3 Ray::point(float time) const {
	mVec3 bt = mVec3(ori.x + time * this->dir.x, \
		ori.y + time * this->dir.y, ori.z + time * this->dir.z);
	return bt;
}
 Ray::~Ray() {
}
 Ray::Ray() {
}

////this is not aabb
//class Cubid :public Surface
//{
//public:
//	mVec3 Center;
//	Plane Faces[6];
//	float h, w, l;
//	Cubid() = default;
//	inline Cubid(mVec3 pm, float h = 1.0f, float w = 1.0f, float l = 1.0f) :Center(pm), h(h), w(w), l(l) {
//		Faces[0] = Plane(Center + (mVec3(0, 0, -1))*0.5*l, mVec3(0, 0, -1));
//		Faces[1] = Plane(Center + (mVec3(0, 0, 1))*0.5*l, mVec3(0, 0, 1));
//		Faces[2] = Plane(Center + (mVec3(0, 1, 0))*0.5*w, mVec3(0, 1, 0));
//		Faces[3] = Plane(Center + (mVec3(0, -1, 0))*0.5*w, mVec3(0, -1, 0));
//		Faces[4] = Plane(Center + (mVec3(1, 0, 0))*0.5*h, mVec3(1, 0, 0));
//		Faces[5] = Plane(Center + (mVec3(-1, 0, 0))*0.5*h, mVec3(-1, 0, 0));
//		initModelMatrix();
//	}
//	inline Cubid(float xmin, float xmax, float ymin, float ymax, float zmin, float zmax) {
//		Center = mVec3((xmax + xmin) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2);
//		float h = xmax - xmin;
//		float w = ymax - ymin;
//		float l = zmax - zmin;
//		Faces[0] = Plane(Center + (mVec3(0, 0, -1))*0.5*l, mVec3(0, 0, -1));
//		Faces[1] = Plane(Center + (mVec3(0, 0, 1))*0.5*l, mVec3(0, 0, 1));
//		Faces[2] = Plane(Center + (mVec3(0, 1, 0))*0.5*w, mVec3(0, 1, 0));
//		Faces[3] = Plane(Center + (mVec3(0, -1, 0))*0.5*w, mVec3(0, -1, 0));
//		Faces[4] = Plane(Center + (mVec3(1, 0, 0))*0.5*h, mVec3(1, 0, 0));
//		Faces[5] = Plane(Center + (mVec3(-1, 0, 0))*0.5*h, mVec3(-1, 0, 0));
//		initModelMatrix();
//	}
//};

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
	}
	inline Sphere(mVec3 ori, float R = 1, float scalex = 1.0f, float scaley = 1.0f, float scalez = 1.0f) :scalex(scalex), scaley(scaley), scalez(scalez) {
		this->center = ori;
		this->R = R;
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

Surface::Surface() {
	initModelMatrix();

}

//class Triangular_pyramid :public Surface {
//public:
//	mVec3 vertex[4];
//	Triangle tris[4];
//	//Cubid AABB;
//	inline Triangular_pyramid(mVec3 A, mVec3 B, mVec3 C, mVec3 D){
//		vertex[0] = A;
//		vertex[1] = B;
//		vertex[2] = C;
//		vertex[3] = D;
//		//float xmin = MIN(MIN(C.x,MIN(A.x, B.x)), D.x);
//		//float xmax = MAX(MAX(C.x, MAX(A.x, B.x)), D.x);
//		//float ymin = MIN(MIN(C.y, MIN(A.y, B.y)), D.y);
//		//float ymax = MAX(MAX(C.y, MAX(A.y, B.y)), D.y);
//		//float zmin = MIN(MIN(C.z, MIN(A.z, B.z)), D.z);
//		//float zmax = MAX(MAX(C.z, MAX(A.z, B.z)), D.z);
//		//AABB = Cubid(xmin, xmax, ymin, ymax, zmin, zmax);
//
//		tris[0] = Triangle(A, B, C);
//		tris[1] = Triangle(A, B, D);
//		tris[2] = Triangle(A, C, D);
//		tris[3] = Triangle(B, C, D);
//
//		initModelMatrix();
//
//	}
//};
 Plane::Plane(mVec3 pm, mVec3 n) :pm(pm), n(n) {

}
 void Surface::updateMaterial(mVec3 a, mVec3 d, mVec3 s, float f ) {
	this->shma = { a, d, s,f };
}
 void  Surface::boudleMaterial(mVec3 a, mVec3 d, mVec3 s, float  f) {
	this->reset_shma = { a, d, s,f };
	this->shma = { a, d, s,f };
}
 void  Surface::ResetMaterial() {
	this->shma = this->reset_shma;
}
 ShadingMaterial  Surface::Material() {
	return this->shma;
}
void  Surface::initModelMatrix() {
	pModleMatrix = Matrix4(eye(4));
	pIModleMatrix = Matrix4(eye(4));
	pNormalMatrix = Matrix4(eye(4));
}
void  Surface::updateModelMatrix(Matrix4& right) {

	pModleMatrix = right * (pModleMatrix);
	pIModleMatrix = inv(pModleMatrix);
	pNormalMatrix = pIModleMatrix.Transpose();


}
















 BVHnode* CreateBVHtree(std::vector<TriangleInMesh>triangleList, const ShareVertexMesh*  const pReferMesh) {
	int TriangleNumber = triangleList.size();
	//init a node
	BVHnode* NewNode = (BVHnode*)malloc(sizeof(BVHnode));
	NewNode->left = NULL;
	NewNode->right = NULL;
	NewNode->LocalTriangleIds = NULL;
	NewNode->LocalTriangleNumber = -1;
	NewNode->LocalMeshIds = -999;
	auto CombineBoundingBoxes = [](BoundingBox& Box1, BoundingBox& Box2)-> BoundingBox {
		float Xmax = MAX(Box1.Xmax, Box2.Xmax);
		float Ymax = MAX(Box1.Ymax, Box2.Ymax);
		float Zmax = MAX(Box1.Zmax, Box2.Zmax);
		float Xmin = MIN(Box1.Xmin, Box2.Xmin);
		float Ymin = MIN(Box1.Ymin, Box2.Ymin);
		float Zmin = MIN(Box1.Zmin, Box2.Zmin);
		return{ Xmax,Xmin, Ymax,Ymin, Zmax ,Zmin };
	};
	auto FindArrayBoxes = [&](std::vector<TriangleInMesh>&localtriangleList)-> BoundingBox {
		float xmin = (float)INT_MAX, xmax = (float)-INT_MAX;
		float ymin = (float)INT_MAX, ymax = (float)-INT_MAX;
		float zmin = (float)INT_MAX, zmax = (float)-INT_MAX;
		for (auto vid : localtriangleList) {
			Triangle me = { pReferMesh->Vertexes[vid.v1],pReferMesh->Vertexes[vid.v2],pReferMesh->Vertexes[vid.v3] };
			mVec3 v = me.a;
			if (v.x < xmin) {
				xmin = v.x;
			}
			if (v.y < ymin) {
				ymin = v.y;
			}
			if (v.z < zmin) {
				zmin = v.z;
			}

			if (v.x > xmax) {
				xmax = v.x;
			}
			if (v.y > ymax) {
				ymax = v.y;
			}
			if (v.z > zmax) {
				zmax = v.z;
			}
			v = me.b;
			if (v.x < xmin) {
				xmin = v.x;
			}
			if (v.y < ymin) {
				ymin = v.y;
			}
			if (v.z < zmin) {
				zmin = v.z;
			}

			if (v.x > xmax) {
				xmax = v.x;
			}
			if (v.y > ymax) {
				ymax = v.y;
			}
			if (v.z > zmax) {
				zmax = v.z;
			}
			v = me.c;
			if (v.x < xmin) {
				xmin = v.x;
			}
			if (v.y < ymin) {
				ymin = v.y;
			}
			if (v.z < zmin) {
				zmin = v.z;
			}

			if (v.x > xmax) {
				xmax = v.x;
			}
			if (v.y > ymax) {
				ymax = v.y;
			}
			if (v.z > zmax) {
				zmax = v.z;
			}

		}

		return { xmax,xmin, ymax,ymin, zmax,zmin };
	};


	if (TriangleNumber <= 20) {
		NewNode->left = NULL;
		NewNode->right = NULL;
		NewNode->LocalTriangleNumber = TriangleNumber;
		NewNode->LocalTriangleIds = (TriangleInMesh*)malloc(sizeof(TriangleInMesh) * TriangleNumber);
		memcpy(NewNode->LocalTriangleIds, &triangleList[0], TriangleNumber * sizeof(TriangleInMesh));
		NewNode->LocalBoundingBox = FindArrayBoxes(triangleList);
		return NewNode;
	}
	else {
		auto b = FindArrayBoxes(triangleList);
		float xlen = b.Xmax - b.Xmin;
		float ylen = b.Ymax - b.Ymin;
		float zlen = b.Zmax - b.Zmin;

		std::vector<TriangleInMesh>LeftArray; LeftArray.reserve(TriangleNumber);
		std::vector<TriangleInMesh>RightArray; RightArray.reserve(TriangleNumber);

		if (xlen >= ylen && xlen >= zlen) {//x
			float divideX = (b.Xmax + b.Xmin)*0.5;
			for (auto triangle : triangleList) {
				Triangle RealTriangle = Triangle(pReferMesh->Vertexes[triangle.v1], pReferMesh->Vertexes[triangle.v2]\
					, pReferMesh->Vertexes[triangle.v3]);
				if (RealTriangle.Center().x < divideX)
					LeftArray.emplace_back(triangle);
				else RightArray.emplace_back(triangle);
			}

		}
		else if (ylen >= xlen && ylen >= zlen) {//y
			float divideY = (b.Ymax + b.Ymin)*0.5;
			for (auto triangle : triangleList) {
				Triangle RealTriangle = Triangle(pReferMesh->Vertexes[triangle.v1], pReferMesh->Vertexes[triangle.v2]\
					, pReferMesh->Vertexes[triangle.v3]);
				if (RealTriangle.Center().y < divideY)
					LeftArray.emplace_back(triangle);
				else RightArray.emplace_back(triangle);
			}
		}
		else {//z
			float divideZ = (b.Zmax + b.Zmin)*0.5;
			for (auto triangle : triangleList) {
				Triangle RealTriangle = Triangle(pReferMesh->Vertexes[triangle.v1], pReferMesh->Vertexes[triangle.v2]\
					, pReferMesh->Vertexes[triangle.v3]);
				if (RealTriangle.Center().z < divideZ)
					LeftArray.emplace_back(triangle);
				else RightArray.emplace_back(triangle);
			}
		}
		if (LeftArray.size() == triangleList.size()) {
			auto last = LeftArray.back();
			RightArray.push_back(last);
			LeftArray.pop_back();
		}
		if (RightArray.size() == triangleList.size()) {
			auto last = RightArray.back();
			LeftArray.push_back(last);
			RightArray.pop_back();
		}

		NewNode->left = CreateBVHtree(LeftArray, pReferMesh);
		NewNode->right = CreateBVHtree(RightArray, pReferMesh);
		BoundingBox Lb = NewNode->left->LocalBoundingBox;
		BoundingBox Rb = NewNode->right->LocalBoundingBox;
		NewNode->LocalBoundingBox = CombineBoundingBoxes(Lb, Rb);
		return NewNode;
	}
}



MeshObjects::MeshObjects(ShareVertexMesh p)  {
		this->VertexPerFaceNormal = p.VertexPerFaceNormal;
		this->Vertexes = p.Vertexes;
		this->TriangleIds = p.TriangleIds;
		this->VertexPredefinedNormal = p.VertexPredefinedNormal;
		initModelMatrix();
		//initBV();
		pBVHmesh=CreateBVHtree(p.TriangleIds,&p);
	}













//inline float intersection(const Sphere& sphere,  Ray& Ray) {
//
//	//intersectInfo AABBtime = intersection(sphere.AABB, Ray);
//	//if (AABBtime.time < 0) return -100.0f ;
//	Ray.dir.x /= sphere.scalex;
//	Ray.dir.y /= sphere.scaley;
//	Ray.dir.z /= sphere.scalez;
//	Ray.ori.x /= sphere.scalex;
//	Ray.ori.y /= sphere.scaley;
//	Ray.ori.z /= sphere.scalez;
//	//只返回近的交点，优先返回小根,相切则省略
//	float a = Ray.dir* Ray.dir;
//	mVec3 E_C(Ray.ori, sphere.center);
//	float b = Ray.dir*E_C;//*2.0f
//	float c = E_C * E_C - sphere.R*sphere.R;
//
//	float delta = b * b - a * c;
//	if (delta < 0.)return -100;
//	float t1 = (sqrt(delta) - b) / a;
//	float t2 = (- b - sqrt(delta)) / a;
//
//	return t1 > t2 ? t2 : t1;
//	//float m1 = t1 > t2 ? t2 : t1;
//	//float m2 = t1 > t2 ? t1 : t2;//大根
//
//	//if (m2 <=0)return -100;
//	//else if (m1 < 0) return m2;//光源在球内部
//	//else return m1;//光在求外部，优先返回小根
//
//
//}




//inline intersectInfo intersectionwithShareMeshNOBVH(const MeshObjects& Mesh, const Ray& Ray) {
//	/*if(Mesh.rootBV.intersection(Ray)==-1)
//		return { -100.0f, };*/
//	float time = (float)INT_MAX;
//	int id = -1;
//	int intersecTriangleId = -1;
//	intersectInfo rInfo;
//	for (auto triangle : Mesh.TriangleIds) {
//		id += 1;
//		Triangle RealTriangle = Triangle(Mesh.Vertexes[triangle.v1], Mesh.Vertexes[triangle.v2]\
//			,Mesh.Vertexes[triangle.v3]);
//		//TODO:三角形相交还要做重心插值
//		intersectInfo Info = intersectionTri(RealTriangle, Ray);
//		float t = Info.time;
//		if (t > 0 && time > t) {
//			time = t;
//			intersecTriangleId = id;
//			rInfo = Info;
//		
//		}
//	}
//	if (intersecTriangleId == -1)
//		return { -100.0f, };
//	else {
//		TriangleInMesh Ids = Mesh.TriangleIds[intersecTriangleId];
//		mVec3 interpolatedNormal=Mesh.VertexPerFaceNormal[Ids.v1] * (1 - rInfo.n.x - rInfo.n.y) +
//		Mesh.VertexPerFaceNormal[Ids.v2] * (rInfo.n.x ) +
//		Mesh.VertexPerFaceNormal[Ids.v3] * ( rInfo.n.y);
//		interpolatedNormal.normalize();
//		return { time,interpolatedNormal };
//	}
//}



unsigned char* Object_List::GetPointer(int id) const {
	if (id >= this->Len() ) {
		return  pointers[this->Len()-1];
	}
	if (id <0) {
		return  pointers[0];
	}
	return pointers[id];
}
unsigned char  Object_List::GetType(int id) const {
	if (id >= this->Len()) {
		return  type[this->Len() - 1];
	}
	if (id < 0) {
		return  type[0];
	}
	return type[id];
}
void Object_List::AddObj(Plane& spo) {
	if (object_count == MAX_OBJECTS)throw("AddObj Fail for out of momeory.");
	else {
		pointers[object_count] = (unsigned char*)&spo;
		type[object_count] = XD_PLANE;
		object_count++;
	}
}
void Object_List::AddObj(BvhTopLevelStructure& spo) {
	if (object_count == MAX_OBJECTS)throw("AddObj Fail for out of momeory.");
	else {
		pointers[object_count] = (unsigned char*)&spo;
		type[object_count] = XD_TOPLEVELMESH;
		object_count++;
	}
}
void Object_List::AddObj(MeshObjects& spo) {
	if (object_count == MAX_OBJECTS)throw("AddObj Fail for out of momeory.");
	else {
		pointers[object_count] = (unsigned char*)&spo;
		type[object_count] = XD_MESH;
		object_count++;
	}
}
bool Object_List::OccludeByOther ( int skipid,Ray  cRay) const {
	cRay = Ray(cRay.ori + cRay.dir*0.0001, cRay.dir);
	for (int u = 0; u < object_count; u++)
	{
		
		float time = -100;
		auto Cobjp = pointers[u];
		switch (type[u]) {
		case (XD_PLANE):
		{Plane Cobj = *((Plane*)Cobjp);
		Matrix4 M = Cobj.pIModleMatrix;
		Ray tRay = Ray((M*mVec4(cRay.ori, 1.0f)).tomVec3(), (M*mVec4(cRay.dir, 0.0f)).tomVec3());
		time = intersection(Cobj, tRay);
		break; }
		case (XD_MESH):
		{MeshObjects* Cobj = ((MeshObjects*)Cobjp);
		Matrix4 M = Cobj->pIModleMatrix;
		Ray tRay = Ray((M*mVec4(cRay.ori, 1.0f)).tomVec3(), (M*mVec4(cRay.dir, 0.0f)).tomVec3());
		auto In= Cobj->intersection(Cobj->pBVHmesh,tRay);
		time = In.time;
		break; }
		case (XD_TOPLEVELMESH):
		{BvhTopLevelStructure* Cobj = ((BvhTopLevelStructure*)Cobjp);
		auto In = Cobj->intersection(Cobj->pBVHmesh, cRay);
		time = In.time;
		break; }
	}

	if (time > 0)return true;
	}
	return false;
}













