#pragma once
#include <vector>
#include <string>
#include <iostream>
#include<fstream>
#include"la.h"
#include"matrix.h"
#include "Mesh.h"



class OBjReader {
private:
	void readVertexs(std::string Filepath);
	void CalculatePerfaceAverageNormal();
public:
	//for object order rendering
	std::vector<mVec3> Vertexes;
	std::vector<mVec2<float>> VertexesTexture;
	std::vector<mVec3> VertexesNormal;
	std::vector<mVec3i> TrianglesIdx;


	//for ray-tracing
	std::vector<mVec3> VertexesPerfaceAverageNormal;
	std::vector<VertexInMesh> VertexInMeshs;
	std::vector<TriangleInMesh> TriangleswithAttibIdx;

	void readTriangleObjFile(std::string Filepath);
	//f a b c d可能是四个点也可能是三个点
	ShareVertexMesh readObj2ShareMesh(std::string Filepath);

	void Load_CalculateNormal_Store(std::string inputFilepath, std::string outputFilepath);
	void Load_CalculateTextureCoordinate_Store(std::string inputFilepath, std::string outputFilepath);
};