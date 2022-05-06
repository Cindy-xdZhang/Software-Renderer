#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <filesystem>
#include<fstream>
#include"la.h"
#include"matrix.h"
class OBjReader {
private:
	template<typename FilePathType = std::filesystem::path>
	void readVertexs(const FilePathType& Filepath);
public:
	std::vector<mVec3f> Vertexes;
	std::vector<mVec2<float>> VertexesTexture;
	std::vector<mVec3f> VertexesNormal;
	std::vector<mVec3i> TrianglesIdx;
	std::vector<mVec3i>& read(std::filesystem::path Filepath, int numberOfTriangles = -1);
	std::vector<mVec3i>& OBjReader::read(std::string Filepath, int numberOfTriangles = -1);
	void OBjReader::Load_CalculateNormal_Store(std::string inputFilepath, std::string outputFilepath);
	void OBjReader::Load_CalculateTextureCoordinate_Store(std::string inputFilepath, std::string outputFilepath);
};