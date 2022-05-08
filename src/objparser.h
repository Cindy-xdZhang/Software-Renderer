#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <filesystem>
#include<fstream>
#include"la.h"
#include"matrix.h"
class OBjReader {

public:
	std::vector<mVec3f> Vertexes;
	std::vector<mVec2<float>> VertexesTexture;
	std::vector<mVec3f> VertexesNormal;
	std::vector<mVec3i> TrianglesIdx;

	static std::vector<mVec2<float>> AssignTextureCoordinates(const std::vector<mVec3f>& vtxs);

	void readObjFile(const std::string& Filepath);
	void OBjReader::Load_CalculateNormal_Store(std::string inputFilepath, std::string outputFilepath);
	void OBjReader::Load_CalculateTextureCoordinate_Store(std::string inputFilepath, std::string outputFilepath);
	inline void clear() {
		Vertexes.clear();
		VertexesNormal.clear();
		VertexesTexture.clear();
		TrianglesIdx.clear();
	}
};