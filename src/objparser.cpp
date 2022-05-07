#include "objparser.h"
#include <xutility>



void OBjReader::readObjFile(const std::string& Filepath) {
	std::fstream f(Filepath);
	std::string  line; 
	Vertexes.reserve(32332);
	TrianglesIdx.reserve(67580);
	VertexesNormal.reserve(32332);
	while (getline(f, line))
	{
		if (line.size() > 0) {
			char beging = line.at(0);
			if (beging == 'v') {
				char beging2 = line.at(1);
				if (beging2 == 'n') {
					const char* Cstr = line.c_str();
					float Vx = 0, Vy = 0, Vz = 0;
					std::ignore=sscanf(Cstr, "vn %f %f %f", &Vx, &Vy, &Vz);
					VertexesNormal.emplace_back(mVec3f(Vx, Vy, Vz));
				}
				else if(beging2 == 't'){
					const char* Cstr = line.c_str();
					float Vx = 0, Vy = 0;
					std::ignore=sscanf(Cstr, "vt %f %f", &Vx, &Vy);
					VertexesTexture.emplace_back(mVec2<float>{Vx, Vy});
				}
				else {
					const char* Cstr = line.c_str();
					float Vx = 0, Vy = 0, Vz = 0;
					std::ignore = sscanf(Cstr, "v %f %f %f", &Vx, &Vy, &Vz);
					Vertexes.emplace_back(mVec3f(Vx, Vy, Vz));

				}
			}
			if (beging == 'f') {
				int count_slash = std::count(line.begin(), line.end(), '/');;
				const char* Cstr = line.c_str();
				int Vx = 0, Vy = 0, Vz = 0;
				if (count_slash==6)
				{
					int not_used=0;
					std::ignore = sscanf(Cstr, "f %d/%d/%d %d/%d/%d %d/%d/%d", &Vx, &not_used, &not_used, &Vy, &not_used, &not_used, &Vz, &not_used, &not_used);
				}
				else {
					std::ignore = sscanf(Cstr, "f %d %d %d", &Vx, &Vy, &Vz);
				}
				mVec3i tmp = { Vx - 1, Vy - 1, Vz - 1 };
				TrianglesIdx.emplace_back(tmp);
			}
	
		}
	
	}
	f.close();
	return ;
}


struct vertex2triangleMap {
	int vertexid=-1;
	std::vector<int>belongTriangles;
};

int CopyFileFunc(std:: string path1, std::string path2) {
	std::ifstream ifs(path1, std::ios::binary);
	if (!ifs.is_open())
	{
		std::ignore = printf("open %s failed\n", path1.c_str());
		return(2);
	}

	std::ofstream ofs(path2, std::ios::binary | std::ios::trunc);
	if (!ofs.is_open())
	{
		std::ignore = printf("create %s failed\n", path2.c_str());
		return(3);
	}

	ofs << ifs.rdbuf();

	ofs.close();
	ifs.close();
	return EXIT_SUCCESS;
}
//LoadObjFiles_withoutNormal_storeObjFileswith_normal
void OBjReader::Load_CalculateNormal_Store(std::string inputFilepath, std::string outputFilepath){
	//load in
	readObjFile(inputFilepath);
	//build key-value map:   key: vertexId, value: trianglesIds(1,2,3,4,..etc)
	std::vector<vertex2triangleMap>Kmaps;
	Kmaps.resize(Vertexes.size());
	for (int tid = 0; tid < TrianglesIdx.size();tid++) {
		auto idxs = TrianglesIdx[tid];
		Kmaps[idxs.x].belongTriangles.push_back(tid);
		Kmaps[idxs.y].belongTriangles.push_back(tid);
		Kmaps[idxs.z].belongTriangles.push_back(tid);
	}
	//calculate all triangle normal
	std::vector<mVec3f> TargetRenderTrianglesNormal;
	TargetRenderTrianglesNormal.reserve(TrianglesIdx.size());
	for (auto idxs : TrianglesIdx) {
		Triangle  tri(Vertexes[idxs.x], Vertexes[idxs.y], Vertexes[idxs.z]);
		mVec3f normal = tri.Normal();
		TargetRenderTrianglesNormal.emplace_back(normal);
	}
	//for vertex :normal= average(triangles normal)
	std::vector<mVec3f> VertexNormal;
	VertexNormal.reserve(Vertexes.size());
	for (auto VtMAP : Kmaps) {
		int id = VtMAP.vertexid;
		mVec3f avgnormal = { 0,0,0 };
		for (auto tri : VtMAP.belongTriangles) {
			avgnormal += TargetRenderTrianglesNormal[tri];
		}
		avgnormal = avgnormal / (float)(VtMAP.belongTriangles.size());
		VertexNormal.emplace_back(avgnormal);
	}
	//store back
	CopyFileFunc(inputFilepath, outputFilepath);
	std::fstream f(outputFilepath, std::ios::app);//创建一个fstream文件流对象
	if (f.fail()) {
		std::cout << "Fail to open file." << std::endl;
	}
	for (auto nor : VertexNormal) {
		char buff[50];
		sprintf_s(buff, "vn %.6f  %.6f  %.6f\n", nor.x, nor.y, nor.z);
		std::string  line = buff;
		f <<line;
	}
	f.close();
	

}


void OBjReader::Load_CalculateTextureCoordinate_Store(std::string inputFilepath, std::string outputFilepath) {
	//load in
	readObjFile(inputFilepath);
	CopyFileFunc(inputFilepath, outputFilepath);
	auto STs=AssignTextureCoordinates(Vertexes);

	std::fstream f(outputFilepath, std::ios::app);
	if (f.fail()) {
		std::cout << "Fail to open file." << std::endl;
	}
	for (auto vt : STs) {
		char buff[50];
		sprintf_s(buff, "vt %.6f  %.6f\n", vt.x, vt.y);
		std::string  line = buff;
		f << line;
	}
	f.close();
}


 std::vector<mVec2<float>> OBjReader::AssignTextureCoordinates(const std::vector<mVec3f>& vtxs) {
	std::vector<mVec2<float>>STs;
	STs.reserve(vtxs.size());

	float co = float(1.0f / 125.0f);
	Matrix4 sm = scaleMatrix(co);

	//r2=x2+y2+zw
	//thta=arctan(  squrt(x2+y2)/z)
	//fi = arctan(y /x)
	for (auto vertex : vtxs) {
		auto HomoLocal =  mVec4f(vertex, 1);
		HomoLocal = sm * HomoLocal;

		float  fi = atan(HomoLocal.y / HomoLocal.x);

		float  thta = atan(sqrt(HomoLocal.x * HomoLocal.x + HomoLocal.y * HomoLocal.y) / HomoLocal.z);


		if (fi < 0)fi += 2 * PI;
		if (thta < 0)thta += 2 * PI;
		fi /= 2 * PI;
		thta /= 2 * PI;
		STs.emplace_back(mVec2<float>{fi, thta });

	}

	return STs;
}