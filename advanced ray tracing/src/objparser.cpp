#include "objparser.h"

int countSpaceInCstr(const char* Cstr) {
	const char* ch = Cstr;
	int spaces = 0;
	while (*ch != '\0')
	{
		if (*ch == ' ')
			++spaces;
		ch++;
	}
	return spaces;
}

struct vertex2triangleMap {
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




void OBjReader::Load_CalculateTextureCoordinate_Store(std::string inputFilepath, std::string outputFilepath) {
	//load in
	readTriangleObjFile(inputFilepath);
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


ShareVertexMesh OBjReader::readObj2ShareMesh(std::string Filepath) {
	std::fstream f(Filepath);//创建一个fstream文件流对象
	std::string  line; //保存读入的每一行
	VertexInMeshs.reserve(32332);
	VertexesNormal.reserve(32332);
	VertexesTexture.reserve(32332);
	//TrianglesIdx.reserve(67580);
	TriangleswithAttibIdx.reserve(67580);
	while (getline(f, line))//会自动把\n换行符去掉 
	{
		if (line.size() > 0) {
			char beging = line.at(0);
			if (beging == 'v') {
				char beging2 = line.at(1);
				if (beging2 == 'n') {//vn
					const char* Cstr = line.c_str();
					float Vx = 0, Vy = 0, Vz = 0;
					sscanf(Cstr, "vn %f %f %f", &Vx, &Vy, &Vz);
					VertexesNormal.emplace_back(mVec3(Vx, Vy, Vz));
				}
				else if (beging2 == 't') {//vt
					const char* Cstr = line.c_str();
					float Vx = 0, Vy = 0;
					sscanf(Cstr, "vt %f %f ", &Vx, &Vy);
					VertexesTexture.emplace_back(mVec2<float>{Vx, Vy});
				}
				else {//v
					const char* Cstr = line.c_str();
					float Vx = 0, Vy = 0, Vz = 0;
					sscanf(Cstr, "v %f %f %f", &Vx, &Vy, &Vz);
					VertexInMeshs.emplace_back(mVec3(Vx, Vy, Vz)*2);

				}
			}
			if (beging == 'f') {
				const char* Cstr = line.c_str();
				//for 3 vertex face
			   /*int Vx = 0, Vy = 0, Vz = 0;
				sscanf(Cstr, "f %d %d %d", &Vx, &Vy, &Vz);*/
				//for 4 vertex face
				mVec3i a, b, c, d;
				TriangleInMesh tmp1;
				TriangleInMesh tmp2;
				d.z = -1;
				int spaces = countSpaceInCstr(Cstr);
				if (spaces == 5) {
					sscanf(Cstr, "f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d", &a.x, &a.y, &a.z, \
						&b.x, &b.y, &b.z, &c.x, &c.y, &c.z, \
						&d.x, &d.y, &d.z);
					tmp1.v1 = a.x - 1; tmp1.vt1 = a.y - 1; tmp1.vn1 = a.z - 1;
					tmp1.v2 = b.x - 1; tmp1.vt2 = b.y - 1; tmp1.vn2 = b.z - 1;
					tmp1.v3 = d.x - 1; tmp1.vt3 = d.y - 1; tmp1.vn3 = d.z - 1;

					tmp2.v1 = c.x - 1; tmp2.vt1 = c.y - 1; tmp2.vn1 = c.z - 1;
					tmp2.v2 = d.x - 1; tmp2.vt2 = d.y - 1; tmp2.vn2 = d.z - 1;
					tmp2.v3 = b.x - 1; tmp2.vt3 = b.y - 1; tmp2.vn3 = b.z - 1;
					TriangleswithAttibIdx.emplace_back(tmp1);
					TriangleswithAttibIdx.emplace_back(tmp2);
				}
				else {
					sscanf(Cstr, "f %d/%d/%d %d/%d/%d %d/%d/%d", &a.x, &a.y, &a.z, \
						&b.x, &b.y, &b.z, &c.x, &c.y, &c.z);
					tmp1.v1 = a.x - 1; tmp1.vt1 = a.y - 1; tmp1.vn1 = a.z - 1;
					tmp1.v2 = b.x - 1; tmp1.vt2 = b.y - 1; tmp1.vn2 = b.z - 1;
					tmp1.v3 = c.x - 1; tmp1.vt3 = c.y - 1; tmp1.vn3 = c.z - 1;
					TriangleswithAttibIdx.emplace_back(tmp1);
		
				}

		

			}

		}

	}
	f.close();
	CalculatePerfaceAverageNormal();
	return 	ShareVertexMesh(TriangleswithAttibIdx, VertexInMeshs, VertexesNormal,VertexesPerfaceAverageNormal);
}


void OBjReader::CalculatePerfaceAverageNormal() {
	
	//build key-value map:   key: vertexId, value: trianglesIds(1,2,3,4,..etc)
	std::vector<vertex2triangleMap>Kmaps;
	Kmaps.resize(VertexInMeshs.size());
	for (int tid = 0; tid < TriangleswithAttibIdx.size(); tid++) {
		auto idxs = TriangleswithAttibIdx[tid];
		Kmaps[idxs.v1].belongTriangles.push_back(tid);
		Kmaps[idxs.v2].belongTriangles.push_back(tid);
		Kmaps[idxs.v3].belongTriangles.push_back(tid);
	}
	//calculate all triangle normal
	std::vector<mVec3> TargetRenderTrianglesNormal;
	TargetRenderTrianglesNormal.reserve(TrianglesIdx.size());
	for (auto idxs : TriangleswithAttibIdx) {
		Triangle  tri(VertexInMeshs[idxs.v1], VertexInMeshs[idxs.v2], VertexInMeshs[idxs.v3]);
		mVec3 normal = tri.Normal();
	
		if (normal* tri.a < 0|| normal * tri.b < 0|| normal * tri.c < 0)normal *= -1;
		normal.normalize();
		
		TargetRenderTrianglesNormal.emplace_back(normal);
	}
	//for vertex :normal= average(triangles normal)
	VertexesPerfaceAverageNormal.reserve(VertexInMeshs.size());
	for (auto VtMAP : Kmaps) {
		mVec3 avgnormal = { 0,0,0 };
		for (auto tri : VtMAP.belongTriangles) {
			avgnormal += TargetRenderTrianglesNormal[tri];
		}
		avgnormal = avgnormal / (float)(VtMAP.belongTriangles.size());
		avgnormal.normalize();
		VertexesPerfaceAverageNormal.emplace_back(avgnormal);
	}
}


//deprecate:
void OBjReader::Load_CalculateNormal_Store(std::string inputFilepath, std::string outputFilepath) {
	//load in
	readTriangleObjFile(inputFilepath);
	//build key-value map:   key: vertexId, value: trianglesIds(1,2,3,4,..etc)
	std::vector<vertex2triangleMap>Kmaps;
	Kmaps.resize(Vertexes.size());
	for (int tid = 0; tid < TrianglesIdx.size(); tid++) {
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
		f << line;
	}
	f.close();


}
//deprecate:
void OBjReader::readTriangleObjFile(std::string Filepath) {
	readVertexs(Filepath);
}

//deprecate:
void OBjReader::readVertexs(std::string Filepath) {
	std::fstream f(Filepath);//创建一个fstream文件流对象
	std::string  line; //保存读入的每一行
	Vertexes.reserve(32332);
	VertexesNormal.reserve(32332);
	VertexesTexture.reserve(32332);
	//TrianglesIdx.reserve(67580);
	TriangleswithAttibIdx.reserve(67580);
	while (getline(f, line))//会自动把\n换行符去掉 
	{
		if (line.size() > 0) {
			char beging = line.at(0);
			if (beging == 'v') {
				char beging2 = line.at(1);
				if (beging2 == 'n') {//vn
					const char* Cstr = line.c_str();
					float Vx = 0, Vy = 0, Vz = 0;
					sscanf(Cstr, "vn %f %f %f", &Vx, &Vy, &Vz);
					VertexesNormal.emplace_back(mVec3(Vx, Vy, Vz));
				}
				else if (beging2 == 't') {//vt
					const char* Cstr = line.c_str();
					float Vx = 0, Vy = 0;
					sscanf(Cstr, "vt %f %f ", &Vx, &Vy);
					VertexesTexture.emplace_back(mVec2<float>{Vx, Vy});
				}
				else {//v
					const char* Cstr = line.c_str();
					float Vx = 0, Vy = 0, Vz = 0;
					sscanf(Cstr, "v %f %f %f", &Vx, &Vy, &Vz);
					Vertexes.emplace_back(mVec3(Vx, Vy, Vz));

				}
			}
			if (beging == 'f') {
				const char* Cstr = line.c_str();
				//for 3 vertex face
			   /*int Vx = 0, Vy = 0, Vz = 0;
				sscanf(Cstr, "f %d %d %d", &Vx, &Vy, &Vz);*/
				//for 4 vertex face
				mVec3i a, b, c, d;
			/*	FaceTriangle tmp1;
				FaceTriangle tmp2;
				sscanf(Cstr, "f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d", &a.x, &a.y, &a.z, \
					&b.x, &b.y, &b.z, &b.x, &c.y, &c.z, \
					&d.x, &d.y, &d.z);
				tmp1.p1.x = a.x - 1; tmp1.p1.y = a.y - 1; tmp1.p1.z = a.z - 1;
				tmp1.p2.x = b.x - 1; tmp1.p2.y = b.y - 1; tmp1.p2.z = b.z - 1;
				tmp1.p3.x = d.x - 1; tmp1.p3.y = d.y - 1; tmp1.p3.z = d.z - 1;

				tmp2.p1.x = c.x - 1; tmp2.p1.y = c.y - 1; tmp2.p1.z = c.z - 1;
				tmp2.p2.x = d.x - 1; tmp2.p2.y = d.y - 1; tmp2.p2.z = d.z - 1;
				tmp2.p3.x = b.x - 1; tmp2.p3.y = b.y - 1; tmp2.p3.z = b.z - 1;
				TriangleswithAttibIdx.emplace_back(tmp1);
				TriangleswithAttibIdx.emplace_back(tmp2);*/
			
			}

		}

	}
	f.close();
	return;
}
