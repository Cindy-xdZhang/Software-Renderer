#pragma   once  
#include"Parameters.h"
#include "la.h"
#include "matrix.h"
#include "Mesh.h"
#include <stack>
#include <iostream>
#include <fstream>
#include <iomanip>
extern bool TureisUsePredefinedNormalFalseIsUseAverageNormal;

enum datatype_eum {
	XD_INT = 1,
	XD_FLOAT = 2,
	XD_DOUBLE = 3,
	XD_SPHERE = 4,
	XD_TOPLEVELMESH = 4,
	XD_PLANE = 6,
	XD_PARAMRID = 7,
	XD_AABB = 8,
	XD_MESH = 9,
};

struct ShadingMaterial {
	mVec3 Ka;
	mVec3 Kd;
	mVec3 Ks;
	float f;
};

struct intersectInfo {
	float time;
	mVec3 n;
	ShadingMaterial shm;
};

class Ray {
public:
	mVec3 ori;
	mVec3 dir;//ori+t*dir;
	 Ray(mVec3 ori, mVec3 Veco);
	 mVec3 point(float time) const;
	 Ray();
	 ~Ray();
};

class Surface {
protected:
	ShadingMaterial shma;	// material
	ShadingMaterial reset_shma;	// material
public:
	Matrix4 pModleMatrix;
	Matrix4 pIModleMatrix;
	Matrix4 pNormalMatrix;
	 void updateMaterial(mVec3 a, mVec3 d, mVec3 s, float f = 16);
	 void boudleMaterial(mVec3 a, mVec3 d, mVec3 s, float  f = 16);
	 void ResetMaterial();
	 ShadingMaterial Material();
	void initModelMatrix();
	void updateModelMatrix(Matrix4& right);
	Surface();
};

class Plane :public Surface
{
public:
	mVec3 pm;	// solution to dot(n,p-p')=0  d
	mVec3 n;		// normal
	Plane() = default;
	 Plane(mVec3 pm, mVec3 n);
};


////aabb:axis aligned box
struct BoundingBox {
	float Xmax, Xmin;
	float Ymax, Ymin;
	float Zmax, Zmin;
	inline mVec3 Center() {
		return { (Xmax + Xmin)*0.5f ,(Ymax + Ymin)*0.5f ,(Zmax + Zmin)*0.5f };
	}
	inline bool intersection(const Ray& Ray) const {
		float t_xmin, t_xmax, t_ymin, t_zmin, t_ymax, t_zmax;
		float ax = 1.0f / Ray.dir.x;
		if (Ray.dir.x >= 0) {
			t_xmin = (this->Xmin - Ray.ori.x) * ax;
			t_xmax = (this->Xmax - Ray.ori.x) * ax;
		}
		else {
			t_xmax = (this->Xmin - Ray.ori.x) * ax;
			t_xmin = (this->Xmax - Ray.ori.x) * ax;
		}
		float ay = 1.0f / Ray.dir.y;
		if (Ray.dir.y >= 0) {
			t_ymin = (this->Ymin - Ray.ori.y) * ay;
			t_ymax = (this->Ymax - Ray.ori.y) * ay;
		}
		else {
			t_ymax = (this->Ymin - Ray.ori.y) * ay;
			t_ymin = (this->Ymax - Ray.ori.y) * ay;
		}
		float az = 1.0f / Ray.dir.z;
		if (Ray.dir.z >= 0) {
			t_zmin = (this->Zmin - Ray.ori.z) * az;
			t_zmax = (this->Zmax - Ray.ori.z) *  az;
		}
		else {
			t_zmax = (this->Zmin - Ray.ori.z) *  az;
			t_zmin = (this->Zmax - Ray.ori.z) *  az;
		}
		float t_enter = MAX(t_zmin, MAX(t_xmin, t_ymin));
		float t_exit = MIN(t_zmax, MIN(t_xmax, t_ymax));
		if (t_exit > 0 && t_exit > t_enter)return true;
		else return false;
	}

	inline BoundingBox operator*(const Matrix4 ModelMatix) const {
		float xmin = (float)INT_MAX, xmax = (float)-INT_MAX;
		float ymin = (float)INT_MAX, ymax = (float)-INT_MAX;
		float zmin = (float)INT_MAX, zmax = (float)-INT_MAX;
		std::vector<mVec3>P8; P8.reserve(8);
		P8.emplace_back((ModelMatix * mVec4(Xmin, Ymin, Zmin, 1.0f)).tomVec3());
		P8.emplace_back((ModelMatix * mVec4(Xmin, Ymin, Zmax, 1.0f)).tomVec3());
		P8.emplace_back((ModelMatix * mVec4(Xmin, Ymax, Zmin, 1.0f)).tomVec3());
		P8.emplace_back((ModelMatix * mVec4(Xmin, Ymax, Zmax, 1.0f)).tomVec3());
		P8.emplace_back((ModelMatix * mVec4(Xmax, Ymin, Zmin, 1.0f)).tomVec3());
		P8.emplace_back((ModelMatix * mVec4(Xmax, Ymin, Zmax, 1.0f)).tomVec3());
		P8.emplace_back((ModelMatix * mVec4(Xmax, Ymax, Zmin, 1.0f)).tomVec3());
		P8.emplace_back((ModelMatix * mVec4(Xmax, Ymax, Zmax, 1.0f)).tomVec3());
		for (auto v : P8) {
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
		
		return { xmax,xmin,ymax,ymin,zmax,zmin };
	}
};

struct BVHnode {

	//for internal node
	BVHnode* left;
	BVHnode* right;
	BoundingBox LocalBoundingBox;

	//for leaf node
	TriangleInMesh* LocalTriangleIds;
	int LocalMeshIds=-1;//-1 toplevel leaf;  mesh tree internol ;mesh tree leaf 
	int LocalTriangleNumber;

};

class MeshObjects :public Surface, public ShareVertexMesh {

public:
	BVHnode* pBVHmesh;
	MeshObjects() = default;
	MeshObjects(ShareVertexMesh p);
	intersectInfo intersection(const BVHnode * const ROOT, const Ray& Ray);
};

class BvhTopLevelStructure {
private:
	BVHnode* BuildTree(std::vector<int>&tmpMeshids) {
		auto CombineBoundingBoxes = [](BoundingBox& Box1, BoundingBox& Box2)-> BoundingBox {
			float Xmax = MAX(Box1.Xmax, Box2.Xmax);
			float Ymax = MAX(Box1.Ymax, Box2.Ymax);
			float Zmax = MAX(Box1.Zmax, Box2.Zmax);
			float Xmin = MIN(Box1.Xmin, Box2.Xmin);
			float Ymin = MIN(Box1.Ymin, Box2.Ymin);
			float Zmin = MIN(Box1.Zmin, Box2.Zmin);
			return{ Xmax,Xmin, Ymax,Ymin, Zmax ,Zmin };
		};
		auto FindArrayBoxes = [&](std::vector<int>&localMeshidsList)-> BoundingBox {
			float xmin = (float)INT_MAX, xmax = (float)-INT_MAX;
			float ymin = (float)INT_MAX, ymax = (float)-INT_MAX;
			float zmin = (float)INT_MAX, zmax = (float)-INT_MAX;
			BoundingBox restul = { xmax,xmin, ymax,ymin, zmax,zmin };
			for (auto vid : localMeshidsList) {
				BoundingBox temp = Meshes[vid].pBVHmesh->LocalBoundingBox * Meshes[vid].pModleMatrix;
				restul = CombineBoundingBoxes(restul, temp);
			}
			return restul;
		};
		int MeshesNumber = tmpMeshids.size();
		//init a node
		BVHnode* NewNode = (BVHnode*)malloc(sizeof(BVHnode));
		NewNode->left = NULL;
		NewNode->right = NULL;
		NewNode->LocalTriangleIds = NULL;
		NewNode->LocalTriangleNumber = -1;
		NewNode->LocalMeshIds = NULL;
	 if (MeshesNumber== 1) {
			NewNode->left = NULL;
			NewNode->right = NULL;
			NewNode->LocalTriangleNumber = MeshesNumber;
			NewNode->LocalMeshIds = tmpMeshids[0];
			NewNode->LocalBoundingBox = FindArrayBoxes(tmpMeshids);
			return NewNode;
		}
		else {
			auto b = FindArrayBoxes(tmpMeshids);
			float xlen = b.Xmax - b.Xmin;
			float ylen = b.Ymax - b.Ymin;
			float zlen = b.Zmax - b.Zmin;

			std::vector<int>LeftArray; LeftArray.reserve(MeshesNumber);
			std::vector<int>RightArray; RightArray.reserve(MeshesNumber);

			if (xlen >= ylen && xlen >= zlen) {
				float divideX = (b.Xmax + b.Xmin)*0.5;
				for (auto triangle : tmpMeshids) {
					BoundingBox Bbbo = Meshes[triangle].pBVHmesh->LocalBoundingBox;
					if (Bbbo.Center().x < divideX)
						LeftArray.emplace_back(triangle);
					else RightArray.emplace_back(triangle);
				}

			}
			else if (ylen >= xlen && ylen >= zlen) {//y
				float divideY = (b.Ymax + b.Ymin)*0.5;
				for (auto triangle : tmpMeshids) {
					BoundingBox Bbbo = Meshes[triangle].pBVHmesh->LocalBoundingBox;
					if (Bbbo.Center().y < divideY)
						LeftArray.emplace_back(triangle);
					else RightArray.emplace_back(triangle);
				}
			}
			else {//z
				float divideZ = (b.Zmax + b.Zmin)*0.5;
				for (auto triangle : tmpMeshids) {
					BoundingBox Bbbo = Meshes[triangle].pBVHmesh->LocalBoundingBox;
					if (Bbbo.Center().z < divideZ)
						LeftArray.emplace_back(triangle);
					else RightArray.emplace_back(triangle);
				}
			}
			if (LeftArray.size() == tmpMeshids.size()) {
				auto last = LeftArray.back();
				RightArray.push_back(last);
				LeftArray.pop_back();
			}
			if (RightArray.size() == tmpMeshids.size()) {
				auto last = RightArray.back();
				LeftArray.push_back(last);
				RightArray.pop_back();
			}

			NewNode->left = BuildTree(LeftArray);
			NewNode->right = BuildTree(RightArray);
			BoundingBox Lb = NewNode->left->LocalBoundingBox;
			BoundingBox Rb = NewNode->right->LocalBoundingBox;
			NewNode->LocalBoundingBox = CombineBoundingBoxes(Lb, Rb);
			return NewNode;
		}
	}

	void Traverse(BVHnode* root) {
		if (root != NULL) {
			if (root->left = NULL)//is leaf for top level tree?
			{
				return;
			}
			else {
				nodes.push(root);
				Traverse(root->left);
				Traverse(root->right);
			}
		}
		
	}
public:
	std::vector<MeshObjects>Meshes;
	BVHnode* pBVHmesh;
	std::stack<BVHnode*> nodes;
	std::vector<int> MeshIDs;
	void deleteTree(){
		
		Traverse(pBVHmesh);
		while (nodes.size() > 0) {
			free(nodes.top());
			nodes.pop();
		}
		
		


	}

	inline void UpdateBoundingBOx( ) {
		deleteTree();
		pBVHmesh = BuildTree(MeshIDs);
	}
	inline BvhTopLevelStructure(const std::vector<MeshObjects>&Meshes){
		this->Meshes = Meshes;
		for (int i = 0; i < Meshes.size(); i++)
			MeshIDs.emplace_back(i);
		pBVHmesh=BuildTree(MeshIDs);

		


	}
	inline intersectInfo intersection(BVHnode*ROOT,const Ray& cRay) {
		if (ROOT->LocalBoundingBox.intersection(cRay) == false) return { -100.0f, };

		if (ROOT->left == NULL) {// top level mesh tree LEAF->go into local mesh tree root
			MeshObjects* pPossiblehitMesh = &Meshes[ROOT->LocalMeshIds];
			Matrix4 M = (pPossiblehitMesh->pIModleMatrix);
			Matrix4 normalM = (pPossiblehitMesh->pNormalMatrix);
			Ray tRay = Ray((M*mVec4(cRay.ori, 1.0f)).tomVec3(), (M*mVec4(cRay.dir, 0.0f)).tomVec3());
			intersectInfo rInfo = pPossiblehitMesh->intersection(pPossiblehitMesh->pBVHmesh, tRay);
			rInfo.n = (normalM * mVec4(rInfo.n, 0.0)).tomVec3();
			rInfo.shm = pPossiblehitMesh->Material();
			return rInfo;
		}
		intersectInfo h1 = intersection(ROOT->left, cRay);
		intersectInfo h2 = intersection(ROOT->right, cRay);
		if (h1.time > 0 && h2.time > 0)
			return h1.time > h2.time ? h2 : h1;
		else if (h1.time > 0)return h1;
		else return h2;
	}
	void storeModelMatrix() {
		std::ofstream outFile;
		outFile.open("test.txt", std::ios::out);
		for (auto& mesh : Meshes) {
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 4; j++) {
					outFile << std::setw(7) << std::setfill('0') << std::setiosflags(std::ios::fixed) << std::setprecision(6) << mesh.pModleMatrix.p[i][j] << " ";
				}
			outFile << std::endl;
		}
		outFile.close();
	}
	void LoadModelMatrix() {
		std::fstream inFile;
		inFile.open("test.txt", std::ios::in);
		if (!inFile.is_open())
		{
			std::cout << "Error opening file"; exit(1);
		}
		for (auto& mesh : Meshes) {
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 4; j++) {
					inFile >> std::setw(7) >> std::setfill('0') >> std::setiosflags(std::ios::fixed) >> std::setprecision(6) >> mesh.pModleMatrix.p[i][j];
				}
			mesh.pIModleMatrix = inv(mesh.pModleMatrix);
			mesh.pNormalMatrix = mesh.pIModleMatrix.Transpose();
		}
		inFile.close();
		UpdateBoundingBOx();
	}
	
};


class Object_List {
	int object_count;
	unsigned char* pointers[MAX_OBJECTS];
	unsigned char type[MAX_OBJECTS];
public:
	Object_List() {
		object_count = 0;
	}
	int Len() const {
		return object_count;
	}
	unsigned char* GetPointer(int id) const;
	unsigned char  GetType(int id) const;
	void AddObj(Plane& spo);
	void AddObj(MeshObjects& spo);
	void AddObj(BvhTopLevelStructure& spo);
	bool OccludeByOther(int skipid, Ray  cRay)const;
};

struct VisibleLightSource {

	
};


inline intersectInfo intersectionTri(const Triangle& tri, const Ray& Ray) {
	//相交则返回t大于零，否则返回-100
	//如果平行则省略，视为不相交
	//相交于边 （0，0，0）+t（0,9,0） 与 （1,2,0）（0,2，2）（0，2,-1） ,视为相交
	//光线从三角形出发向外射（0，0，0）+t（0,9,0）与 （1,0,0）（-1,0，-2）（-1，0,4）,或相较于顶点视为不相交
	//mVec3 L1(tri.a, tri.b);
	//mVec3 L2(tri.a, tri.c);
	//mVec3 Normal = L1.cross_product(L2);
	//if (Normal*Ray.dir < FLOAT_ZERO && Normal*Ray.dir*-1 < FLOAT_ZERO)return -100;

	//x+t * xd=xa+B(xb-xa)+G(xc-xa)
	//y+t * yd=ya+B(yb-ya)+G(yc-ya)
	//z+t * zd=za+B(zb-za)+G(zc-za)
	//=>D:                              Dt:                  Dg:                 Db:
	// B(xb-xa)+G(xc-xa)-t * xd =x-xa  (xb-xa) (xc-xa) x-xa  (xb-xa) (x-xa) -xd  (x-xa) (xc-xa) -xd 
	//B(yb-ya)+G(yc-ya)-t * yd=y-ya    (yb-ya) (yc-ya) y-ya  (yb-ya) (y-ya) -yd  (y-ya) (yc-ya) -yd
	//B(zb-za)+G(zc-za)-t * zd =z-za   (zb-za) (zc-za) z-za  (zb-za) (z-za) -zd  (z-za) (zc-za) -zd    
	float xa = (float)tri.a.x;
	float xb = (float)tri.b.x;
	float xc = (float)tri.c.x;
	float ya = (float)tri.a.y;
	float yb = (float)tri.b.y;
	float yc = (float)tri.c.y;
	float za = (float)tri.a.z;
	float zb = (float)tri.b.z;
	float zc = (float)tri.c.z;
	float xd = (float)Ray.dir.x;
	float yd = (float)Ray.dir.y;
	float zd = (float)Ray.dir.z;
	float x = (float)Ray.ori.x;
	float y = (float)Ray.ori.y;
	float z = (float)Ray.ori.z;
	//compute D(系数矩阵行列式)=(xb - xa)[(yc-ya) * (-zd) +yd*(zc-za)]-(xc - xa)[ (yb-ya) * -zd +yd*(zb-za)]- xd*[ (yb-ya) * (zc-za) -(zb-za)*(yc-ya)]
	float D = (xb - xa)*((yc - ya) * (-zd) + yd * (zc - za)) - \
		(xc - xa)*((yb - ya) * (-zd) + yd * (zb - za)) - xd * ((yb - ya) * (zc - za) - (zb - za)*(yc - ya));
	float iD = 1.0 / D;
	//compute t=Dt/D
	float Dt = (xb - xa)*((yc - ya) * (z - za) - (y - ya) * (zc - za)) - \
		(xc - xa)*((yb - ya) * (z - za) - (y - ya) * (zb - za)) + (x - xa) * ((yb - ya) * (zc - za) - (zb - za)*(yc - ya));
	float t = Dt * iD;
	if (t < 0.)return { -100 };
	//compute gama=Dg/D
	float Dg = (xb - xa)*((y - ya) * (-zd) + yd * (z - za)) - \
		(x - xa)*((yb - ya) * (-zd) + yd * (zb - za)) - xd * ((yb - ya) * (z - za) - (zb - za)*(y - ya));
	float gama = Dg * iD;
	if (gama < 0. || gama >1)return { -100 };
	//compute beta=Db/D
	float Db = (x - xa)*((yc - ya) * (-zd) + yd * (zc - za)) - \
		(xc - xa)*((y - ya) * (-zd) + yd * (z - za)) - xd * ((y - ya) * (zc - za) - (z - za)*(yc - ya));
	float beta = Db * iD;
	if (beta < 0. || beta >1 || beta + gama > 1)return { -100 };

	return { t, mVec3(beta,gama,0) };


}

inline float intersection(const Plane& p, const Ray& Ry)
{

	float dotn_d = p.n*Ry.dir;
	if (dotn_d == 0.)return  -100;
	//if (dotn_d > 0.) return -100;//!!!!理解为啥是大于0就不要了，大于零代表光从面下方穿过故不要了，小于0则留下

	float t = mVec3(p.pm, Ry.ori)* p.n / dotn_d;
	/*if (t > 0.)	return t;
	return  -100;*/
	return t;//t<0 外部函数自动判定不相交
}

inline intersectInfo MeshObjects::intersection(const BVHnode * const ROOT, const Ray& Ray) {
	if (ROOT->LocalBoundingBox.intersection(Ray) == false) return { -100.0f, };

	if (ROOT->left == NULL) {//LEAF
		float time = (float)INT_MAX;
		TriangleInMesh Ids = { -1, };
		intersectInfo rInfo;
		for (int localID = 0; localID < ROOT->LocalTriangleNumber; localID++) {
			auto triangle = ROOT->LocalTriangleIds[localID];
			Triangle RealTriangle = Triangle(this->Vertexes[triangle.v1], this->Vertexes[triangle.v2]\
				, this->Vertexes[triangle.v3]);

			intersectInfo Info = intersectionTri(RealTriangle, Ray);
			float t = Info.time;
			if (t > 0 && time > t) {
				time = t;
				Ids = triangle;
				rInfo = Info;
			}
		}
		if (Ids.v1 == -1)
			return { -100.0f, };
		else {
			mVec3 interpolatedNormal;
			if (TureisUsePredefinedNormalFalseIsUseAverageNormal)
				interpolatedNormal = this->VertexPredefinedNormal[Ids.vn1] * (1 - rInfo.n.x - rInfo.n.y) +
				this->VertexPredefinedNormal[Ids.vn2] * (rInfo.n.x) +
				this->VertexPredefinedNormal[Ids.vn3] * (rInfo.n.y);
			else
				interpolatedNormal = this->VertexPerFaceNormal[Ids.v1] * (1 - rInfo.n.x - rInfo.n.y) +
				this->VertexPerFaceNormal[Ids.v2] * (rInfo.n.x) +
				this->VertexPerFaceNormal[Ids.v3] * (rInfo.n.y);
			interpolatedNormal.normalize();
			return { time,interpolatedNormal };
		}
	}

	intersectInfo h1 = intersection(ROOT->left, Ray);
	intersectInfo h2 = intersection(ROOT->right, Ray);
	if (h1.time > 0 && h2.time > 0)
		return h1.time > h2.time ? h2 : h1;
	else if (h1.time > 0)return h1;
	else return h2;
}

inline BVHnode* CreateBVHtree(std::vector<TriangleInMesh>triangleList, \
	const ShareVertexMesh*  const pReferMesh);


