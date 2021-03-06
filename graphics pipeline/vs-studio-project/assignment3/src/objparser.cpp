#include "objparser.h"

std::vector<mVec3i>& OBjReader::read(std::string Filepath,int numberOfTriangles) {
	readVertexs(Filepath);
	if (numberOfTriangles == -1) {
		return TrianglesIdx;
	}
	else {
		TrianglesIdx.resize(numberOfTriangles);
		return TrianglesIdx;
	}
}

void OBjReader::readVertexs(std::string Filepath) {
	std::fstream f(Filepath);//创建一个fstream文件流对象
	std::string  line; //保存读入的每一行
	Vertexes.reserve(32332);
	VertexesNormal.reserve(32332);
	TrianglesIdx.reserve(67580);
	while (getline(f, line))//会自动把\n换行符去掉 
	{
		if (line.size() > 0) {
			char beging = line.at(0);
			if (beging == 'v') {
				char beging2 = line.at(1);
				if (beging2 == 'n') {
					const char* Cstr = line.c_str();
					float Vx = 0, Vy = 0, Vz = 0;
					sscanf(Cstr, "vn %f %f %f", &Vx, &Vy, &Vz);
					VertexesNormal.emplace_back(mVec3(Vx, Vy, Vz));
				}
				else if(beging2 == 't'){
					const char* Cstr = line.c_str();
					float Vx = 0, Vy = 0;
					sscanf(Cstr, "vt %f %f ", &Vx, &Vy);
					VertexesTexture.emplace_back(mVec2<float>{Vx, Vy});
				}
				else {
					const char* Cstr = line.c_str();
					float Vx = 0, Vy = 0, Vz = 0;
					sscanf(Cstr, "v %f %f %f", &Vx, &Vy, &Vz);
					Vertexes.emplace_back(mVec3(Vx, Vy, Vz));

				}
			}
			if (beging == 'f') {
				const char* Cstr = line.c_str();
				int Vx = 0, Vy = 0, Vz = 0;
				sscanf(Cstr, "f %d %d %d", &Vx, &Vy, &Vz);
				mVec3i tmp = { Vx - 1, Vy - 1, Vz - 1 };
				TrianglesIdx.emplace_back(tmp);
			}
	
		}
	
	}
	f.close();
	return ;
}
struct vertex2triangleMap {
	int vertexid;
	std::vector<int>belongTriangles;
};
int CopyFileFunc(std:: string path1, std::string path2) {
	std::ifstream ifs(path1, std::ios::binary);
	if (!ifs.is_open())
	{
		printf("open %s failed\n", path1);
		return(2);
	}

	std::ofstream ofs(path2, std::ios::binary | std::ios::trunc);
	if (!ofs.is_open())
	{
		printf("create %s failed\n", path2);
		return(3);
	}

	ofs << ifs.rdbuf();

	ofs.close();
	ifs.close();
}
//LoadObjFiles_withoutNormal_storeObjFileswith_normal
void OBjReader::Load_CalculateNormal_Store(std::string inputFilepath, std::string outputFilepath){
	//load in
	read(inputFilepath);
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
	std::vector<mVec3> TargetRenderTrianglesNormal;
	TargetRenderTrianglesNormal.reserve(TrianglesIdx.size());
	for (auto idxs : TrianglesIdx) {
		Triangle  tri(Vertexes[idxs.x], Vertexes[idxs.y], Vertexes[idxs.z]);
		mVec3 normal = tri.Normal();
		TargetRenderTrianglesNormal.emplace_back(normal);
	}
	//for vertex :normal= average(triangles normal)
	std::vector<mVec3> VertexNormal;
	VertexNormal.reserve(Vertexes.size());
	for (auto VtMAP : Kmaps) {
		int id = VtMAP.vertexid;
		mVec3 avgnormal = { 0,0,0 };
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
	read(inputFilepath);
	Matrix4  ModleMatrix = translateMatrix(mVec3(-125, -125, -125));
	float co = float(1.0f / 125.0f);
	Matrix4 sm = scaleMatrix(co);
	std::vector<mVec2<float>>STs;
	STs.reserve(Vertexes.size());
	//r2=x2+y2+zw
	//thta=arctan(  squrt(x2+y2)/z)
	//fi = arctan(y /x)
	for (auto vertex : Vertexes) {
		auto HomoLocal=   ModleMatrix*mVec4(vertex, 1);
		HomoLocal = sm * HomoLocal;

		float  fi= atan(HomoLocal.y/ HomoLocal.x); 

		float  thta = atan(sqrt(HomoLocal.x*HomoLocal.x+ HomoLocal.y* HomoLocal.y )/ HomoLocal.z  );
		
	
		if (fi<0)fi+= 2 * PI;
		if (thta < 0)thta += 2 * PI;
		fi /= 2 * PI;
		thta /= 2 * PI;
		STs.emplace_back(mVec2<float>{fi ,thta });

	}
	CopyFileFunc(inputFilepath, outputFilepath);

	std::fstream f(outputFilepath, std::ios::app);//创建一个fstream文件流对象
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