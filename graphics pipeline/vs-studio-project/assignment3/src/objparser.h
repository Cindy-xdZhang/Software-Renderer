#pragma once
#include <vector>
#include <string>
#include <iostream>
#include<fstream>
#include"la.h"
#include"matrix.h"
class OBjReader {
private:
	void readVertexs(std::string Filepath);
public:
	std::vector<mVec3> Vertexes;
	std::vector<mVec2<float>> VertexesTexture;
	std::vector<mVec3> VertexesNormal;
	std::vector<mVec3i> TrianglesIdx;
	std::vector<mVec3i>& read(std::string Filepath);
	void Load_CalculateNormal_Store(std::string inputFilepath, std::string outputFilepath);
	void Load_CalculateTextureCoordinate_Store(std::string inputFilepath, std::string outputFilepath);
};