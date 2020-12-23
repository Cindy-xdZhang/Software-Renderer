#pragma once
#include"la.h"
#include<vector>
typedef mVec3  VertexInMesh ;

struct TriangleInMesh {
	int v1,v2,v3;
	int vt1,vt2,vt3;
	int vn1,vn2,vn3;
};


class ShareVertexMesh {
public:
	
	std::vector<TriangleInMesh>TriangleIds;
	std::vector<VertexInMesh>Vertexes;
	std::vector<mVec3>VertexPredefinedNormal;
	std::vector<mVec3>VertexPerFaceNormal;
	ShareVertexMesh() = default;
	ShareVertexMesh(std::vector<TriangleInMesh>AllTriangles, std::vector<VertexInMesh> AllVertexes, std::vector<mVec3> VertexPredefinedNormal,
		std::vector<mVec3> VertexPerFaceNormal);

};



//#TODO
class SimpleMesh {
public:
	std::vector<Triangle>AllTriangles;
	std::vector<VertexInMesh> AllVertexes;
	SimpleMesh() = default;
	SimpleMesh(std::vector<Triangle>AllTriangles, std::vector<VertexInMesh> AllVertexes);

};