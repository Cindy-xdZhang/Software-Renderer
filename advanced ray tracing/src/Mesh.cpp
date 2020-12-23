#include"Mesh.h"
ShareVertexMesh::ShareVertexMesh(std::vector<TriangleInMesh>AllTriangles, 
	std::vector<VertexInMesh> AllVertexes, std::vector<mVec3> VertexPredefinedNormal,
	std::vector<mVec3> VertexPerFaceNormal):
	TriangleIds(AllTriangles), Vertexes(AllVertexes), VertexPredefinedNormal(VertexPredefinedNormal), VertexPerFaceNormal(VertexPerFaceNormal)
{

}









SimpleMesh::SimpleMesh(std::vector<Triangle>AllTriangles, std::vector<VertexInMesh> AllVertexes):
	AllTriangles(AllTriangles), AllVertexes(AllVertexes) {

}
