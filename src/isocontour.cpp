#include "isocontour.h"
#include<assert.h>
#include <random>

double Epsilon = 0.001;


MarchingSquaresDrawer::MarchingSquaresDrawer() {
	this->hasdata = false;
}

bool MarchingSquaresDrawer::SetScalarFieldData(int dimx, int dimy, int dimz, unsigned short * dataptr) {
	
	if (dataptr == NULL) {
		std::cout << "No Scalar data." << std::endl;
		return false;
	}
	else {
		//this->data = dataptr;
		this->dimx = dimx;
		this->dimy = dimy;
		this->dimz = dimz;
		if (hasdata) {
			delete this->data;		
		}

		this->data = new double[dimx * dimy * dimz];
		double max = -INFINITY;
		double min = INFINITY;
		for (int z = 0; z < dimz; z++) {
			for (int y = 0; y < dimy; y++) {
				for (int x = 0; x < dimx; x++) {
					
					int flatindex = x + (y * dimx) + (z * dimy * dimx);
					this->data[flatindex] =(double) (dataptr[flatindex]) ;

					if (this->data[flatindex] > max) {
						max = this->data[flatindex];
					}
					if (this->data[flatindex] < min) {
						min = this->data[flatindex];
					}

				}
			}
		}
		for (int z = 0; z < dimz; z++) {
			for (int y = 0; y < dimy; y++) {
				for (int x = 0; x < dimx; x++) {
					int flatindex = x + (y * dimx) + (z * dimy * dimx);
					this->data[flatindex] = (this->data[flatindex]-min)/ (max-min);
				}
			}
		}

		return true;

	}
}

MarchingSquaresDrawer::~MarchingSquaresDrawer() {
	if (this->hasdata)delete data;
}






int  CodingCell(const double a, const double b, const double c, const double d, const double iso_value) {
	/*unsigned char n = 0;
	if (dll > iso_value) n += 1;
	if (dul > iso_value) n += 4;
	if (dur > iso_value) n += 8;

	if (dlr > iso_value) n += 2;
	return n;*/
	double epsilo = Epsilon;
	int n = 0;
	if (a > iso_value- epsilo) n += 1;
	if (b > iso_value - epsilo) n += 8;
	if (c > iso_value - epsilo) n += 4;
	if (d > iso_value - epsilo) n += 2;
	return n;

}



#if defined(OPENGL)

void MarchingSquaresDrawer::MydrawIsoCurve(int CurrentSlice, int CurrentAxis, double X_MIN,  double Y_MIN, double X_MAX, double Y_MAX,  double iso_value) {

	int N_X, N_Y;
	if (CurrentAxis == 0) {
		N_X = dimy;
		N_Y = dimz;
	}
	if (CurrentAxis == 1) {
		N_X = dimx;
		N_Y = dimz;
	}
	if (CurrentAxis == 2) {
		N_X = dimx;
		N_Y = dimy;
	}

	auto draw_one = [iso_value, X_MAX, X_MIN, Y_MAX, Y_MIN, N_X, N_Y](int num, int i, int j, double a, double b, double c, double d)
	{
		double x1, y1, x2, y2;
		double ox, oy;
		double dx, dy;
		dx = (X_MAX - X_MIN) / (N_X - 1.0);
		dy = (Y_MAX - Y_MIN) / (N_Y - 1.0);
		ox = X_MIN + i * dx;
		oy = Y_MIN + j * dy;
		switch (num)
		{
		case 1: case 14:
			x1 = ox;
			y1 = oy + dy * (iso_value - a) / (d - a);
			x2 = ox + dx * (iso_value - a) / (b - a);
			y2 = oy;
			break;
		case 2: case 13:
			x1 = ox;
			y1 = oy + dy * (iso_value - a) / (d - a);
			x2 = ox + dx * (iso_value - d) / (c - d);
			y2 = oy + dy;
			break;
		case 4: case 11:
			x1 = ox + dx * (iso_value - d) / (c - d);
			y1 = oy + dy;
			x2 = ox + dx;
			y2 = oy + dy * (iso_value - b) / (c - b);
			break;
		case 7: case 8:
			x1 = ox + dx * (iso_value - a) / (b - a);
			y1 = oy;
			x2 = ox + dx;
			y2 = oy + dy * (iso_value - b) / (c - b);
			break;
		}
	
		glBegin(GL_LINES);
		
		glVertex2d(x1, y1);
		glVertex2d(x2, y2);
		glEnd();
	};
	auto draw_adjacent = [iso_value, X_MAX, X_MIN, Y_MAX, Y_MIN, N_X, N_Y](int num, int i, int j, double a, double b, double c
		, double d)
	{
		double x1, y1, x2, y2;
		double ox, oy;
		double dx, dy;
		dx = (X_MAX - X_MIN) / (N_X - 1.0);
		dy = (Y_MAX - Y_MIN) / (N_Y - 1.0);
		ox = X_MIN + i * (X_MAX - X_MIN) / (N_X - 1.0);
		oy = Y_MIN + j * (Y_MAX - Y_MIN) / (N_Y - 1.0);
		switch (num)
		{
		case 3: case 12:
			x1 = ox + dx * (iso_value - a) / (b - a);
			y1 = oy;
			x2 = ox + dx * (iso_value - d) / (c - d);
			y2 = oy + dy;
			break;
		case 6: case 9:
			x1 = ox;
			y1 = oy + dy * (iso_value - a) / (d - a);
			x2 = ox + dx;
			y2 = oy + dy * (iso_value - b) / (c - b);
			break;
		}

		glBegin(GL_LINES);
		
		glVertex2d(x1, y1);
		glVertex2d(x2, y2);
		glEnd();
	};
	auto draw_opposite = [iso_value, X_MAX, X_MIN, Y_MAX, Y_MIN, N_X, N_Y](int num, int i, int j, double a, double b, double c
		, double d)
	{
		double x1, y1, x2, y2, x3, y3, x4, y4;
		double ox, oy;
		double dx, dy;
		dx = (X_MAX - X_MIN) / (N_X - 1.0);
		dy = (Y_MAX - Y_MIN) / (N_Y - 1.0);
		ox = X_MIN + i * (X_MAX - X_MIN) / (N_X - 1.0);
		oy = Y_MIN + j * (Y_MAX - Y_MIN) / (N_Y - 1.0);
		switch (num)
		{
		case 5:
			x1 = ox;
			y1 = oy + dy * (iso_value - a) / (d - a);
			x2 = ox + dx * (iso_value - a) / (b - a);
			y2 = oy;
			x3 = ox + dx * (iso_value - d) / (c - d);
			y3 = oy + dy;
			x4 = ox + dx;
			y4 = oy + dy * (iso_value - b) / (c - b);
			break;

		case 10:
			x1 = ox;
			y1 = oy + dy * (iso_value - a) / (d - a);
			x2 = ox + dx * (iso_value - d) / (c - d);
			y2 = oy + dy;
			x3 = ox + dx * (iso_value - d) / (c - d);
			y3 = oy;
			x4 = ox + dx;
			y4 = oy + dy * (iso_value - b) / (c - b);
			break;
		}
	
		glBegin(GL_LINES);
		glVertex2d(x1, y1);
		glVertex2d(x2, y2);
		glVertex2d(x3, y3);
		glVertex2d(x4, y4);
		glEnd();
	};
	///* draw line segments for each case */
	auto drawlines = [&draw_opposite, &draw_one, draw_adjacent](int num, int i, int j, double a, double b, double c, double d)
	{

		switch (num)
		{
		case 1: case 2: case 4: case 7: case 8: case 11: case 13: case
		14:
			draw_one(num, i, j, a, b, c, d);
			break;
		case 3:  case 6:  case 9:  case 12:
			draw_adjacent(num, i, j, a, b, c, d);
			break;
		case 5:  case 10:
			draw_opposite(num, i, j, a, b, c, d);
			break;
		case 0:  case 15:
			break;
		}
	};

	auto SliceGRIDData = [this, CurrentSlice, CurrentAxis]( int coordinates1, int coordinates2) {

		if (CurrentAxis == 0) {
			return this->data[CurrentSlice + (coordinates1 * dimx) + (coordinates2 * dimy * dimx)];
		}
		if (CurrentAxis ==1) {
			return this->data[ coordinates1 + (CurrentSlice * dimx) + (coordinates2 * dimy * dimx)];
		}
		if (CurrentAxis == 2) {
			return this->data[coordinates1 + (coordinates2 * dimx) + (CurrentSlice * dimy * dimx)];
		}
	};
	
	glDisable(GL_TEXTURE_3D);
	glLineWidth(3.0);
	glColor3f(0.0,1.0, 0.1);



	for (int i = 0; i < N_X-1; i++)
		for (int j = 0; j < N_Y-1; j++)
		{

			double a = SliceGRIDData(i, j);
			double b= SliceGRIDData(i+1, j); //go up a row
			double c = SliceGRIDData(i+1, j+1);
			double d = SliceGRIDData(i, j+1);
			int code = CodingCell( a,b,c,d, iso_value);
			
			drawlines(code, i, j, a, b, c, d);
		}
	glEnable(GL_TEXTURE_3D);
}
#endif




//marching cubes



MarchingCubesDrawer::MarchingCubesDrawer() {
	this->hasdata=false; 
	this->normals = std::make_unique<std::vector <Vertex>>();
	this->vertices = std::make_unique<std::vector <Vertex>>();
	this->normals->reserve(10000);
	this->vertices->reserve(10000);
}

void MarchingCubesDrawer::randomGeneratVolumeData(double Noise, mVec3i GridXYZ) {

	//get mean value by random 
	std::uniform_real_distribution<double>r(-1, 1);
	std::default_random_engine e(time(NULL));
	//get normal distribution using random mean
	std::seed_seq seed2{ r(e), r(e), r(e), r(e), r(e), r(e), r(e), r(e) };
	std::mt19937 e2(seed2);
	std::normal_distribution<> normal_dist(0, Noise*Noise);

	dimx = GridXYZ.x;
	dimy = GridXYZ.y;
	dimz = GridXYZ.z;
	assert(dimz>0&&dimx>0&&dimy>0);
	this->volumeData = new double[dimx * dimy * dimz];
	assert(volumeData);
	for (int z = 0; z < dimz; z++) {
		for (int y = 0; y < dimy; y++) {
			for (int x = 0; x < dimx; x++) {

				int flatindex = x + (y * dimx) + (z * dimy * dimx);
				this->volumeData[flatindex] = normal_dist(e2);
			}
		}
	}

}

//procedual landmass generation
void MarchingCubesDrawer::randomGenrateLandmass(mVec3i GridXYZ, double Noise) {

	constexpr double iso_value = 20;
	//randomGeneratVolumeData( Noise, GridXYZ);
	dimx = GridXYZ.x;
	dimy = GridXYZ.y;
	dimz = GridXYZ.z;
	assert(dimz > 0 && dimx > 0 && dimy > 0);
	this->volumeData = new double[dimx * dimy * dimz];
	assert(volumeData);
	memset(volumeData,0, dimx * dimy * dimz*sizeof(float));

	std::uniform_real_distribution<double>r(-1, 1);
	std::default_random_engine e(time(NULL));
	std::seed_seq seed2{ r(e), r(e), r(e), r(e), r(e), r(e), r(e), r(e) };
	std::mt19937 e2(seed2);
	std::normal_distribution<> normal_dist(0, Noise * Noise);

	constexpr int MainLandZ0 = 1;
	constexpr int MainLandZ1 = 2;
	constexpr int MainLandZ2 = 3;

	assert(MainLandZ2 <GridXYZ.z-1);
	for (int y = 0; y < dimy; y++) {
		for (int x = 0; x < dimx; x++) {
			double randomV = normal_dist(e2);
			int thisz = MainLandZ1;
			if (randomV > 0.2) thisz= MainLandZ2;
			if (randomV <- 0.2) thisz = MainLandZ1;
			
			int flatindex = x + (y * dimx) + (thisz* dimy * dimx);
			this->volumeData[flatindex] += iso_value;
		}
	}


	MeshBuiding_MarchingCubes(iso_value);

}


bool MarchingCubesDrawer::SetScalarFieldData(int dimx, int dimy, int dimz, unsigned short* dataptr) {

	if (dataptr == NULL) {
		std::cout << "No Scalar data." << std::endl;
		return false;
	}
	else {
		//this->data = dataptr;
		this->dimx = dimx;
		this->dimy = dimy;
		this->dimz = dimz;
		if (hasdata) {
			delete this->volumeData;
		}

		this->volumeData = new double[dimx * dimy * dimz];
		double max = -INFINITY;
		double min = INFINITY;
		for (int z = 0; z < dimz; z++) {
			for (int y = 0; y < dimy; y++) {
				for (int x = 0; x < dimx; x++) {

					int flatindex = x + (y * dimx) + (z * dimy * dimx);
					this->volumeData[flatindex] = (double)(dataptr[flatindex]);

					if (this->volumeData[flatindex] > max) {
						max = this->volumeData[flatindex];
					}
					if (this->volumeData[flatindex] < min) {
						min = this->volumeData[flatindex];
					}

				}
			}
		}
		for (int z = 0; z < dimz; z++) {
			for (int y = 0; y < dimy; y++) {
				for (int x = 0; x < dimx; x++) {
					int flatindex = x + (y * dimx) + (z * dimy * dimx);
					this->volumeData[flatindex] = (this->volumeData[flatindex] - min) / (max - min);
				}
			}
		}

		return true;

	}
}

double MarchingCubesDrawer::Griddata(int x, int y, int z) {
	x = x <= -1 ? 0 : x >= this->dimx ? x - 1 : x;
	y = y <= -1 ? 0 : y >= this->dimy ? y - 1 : y;
	z = z <= -1 ? 0 : z >= this->dimz ? z - 1 : z;
	int flatindex = x + (y * dimx) + (z * dimy * dimx);
	return	this->volumeData[flatindex];
}

void MarchingCubesDrawer::MeshBuiding_MarchingCubes(const double iso_value) {
	auto CalculateDifference = [this](int x, int y, int z) {
		Vertex n;
		n.x = 0.5 * (Griddata(x + 1, y, z) - Griddata(x - 1, y, z));
		n.y = 0.5 * (Griddata(x, y + 1, z) - Griddata(x, y - 1, z));
		n.z = 0.5 * (Griddata(x, y, z + 1) - Griddata(x, y, z - 1));
		return n;
	};
	const int  EdgesPointIdx[12][2] = {
			{0, 1},
			{1, 2},
			{2, 3},
			{3, 0},
			{4, 5},
			{5, 6},
			{6, 7},
			{7, 4},
			{0, 4},
			{1, 5},
			{2, 6},
			{3, 7}
	};
	const	Vertex VertexInCubeCoordinatesOffset[8] = {
			{0, 0, 0},
			{0, 1, 0},
			{0, 1, 1},
			{0, 0, 1},
			{1, 0, 0},
			{1, 1, 0},
			{1, 1, 1},
			{1, 0, 1},
	};

	vertices->clear();
	normals->clear();

	for (int x = 0; x < dimx - 1; x++) {
		for (int y = 0; y < dimy - 1; y++) {
			for (int z = 0; z < dimz - 1; z++) {
				Vertex p[8];
				double Scalarval[8];
				uint8_t cubeIndex = 0;
				//get cube data
				for (uint8_t i = 0; i < 8; i++) {
					p[i].x = x + VertexInCubeCoordinatesOffset[i].x;
					p[i].y = y + VertexInCubeCoordinatesOffset[i].y;
					p[i].z = z + VertexInCubeCoordinatesOffset[i].z;
					Scalarval[i] = Griddata(p[i].x, p[i].y, p[i].z);

					//coding this cube
					if (Scalarval[i] < iso_value - Epsilon) {
						cubeIndex |= 1u << i;
					}
				}

				//build triangle
				Vertex vertlist[12];
				Vertex normaList[12];
				for (int i = 0; i < 12; i++) {

					uint16_t bit = 1u << i;

					if (edgeTable3D[cubeIndex] & bit) {
						// if this edge has contour on it: edgeTable3D is a a list, using cubeIndex,
						//e.g. if an element of edgeTable3D  is like 0000 0000 1000 0001:means the 0 edge, 7-th edge has contour
						int a = EdgesPointIdx[i][0];
						int b = EdgesPointIdx[i][1];

						Vertex& p1 = p[a];
						Vertex& p2 = p[b];
						double v1 = Scalarval[a];
						double v2 = Scalarval[b];
						Vertex n1 = CalculateDifference(p1.x, p1.y, p1.z);
						Vertex n2 = CalculateDifference(p2.x, p2.y, p2.z);
						double mu = (iso_value - v1) / (v2 - v1);

						Vertex& vertex = vertlist[i];
						vertex.x = p1.x + mu * (p2.x - p1.x);
						vertex.y = p1.y + mu * (p2.y - p1.y);
						vertex.z = p1.z + mu * (p2.z - p1.z);

						Vertex& normal = normaList[i];
						normal.x = n1.x + mu * (n2.x - n1.x);
						normal.y = n1.y + mu * (n2.y - n1.y);
						normal.z = n1.z + mu * (n2.z - n1.z);
					}
					else {
						vertlist[i] = {};
					}
				}
				for (int i = 0; triangleTable[cubeIndex][i] != -1; i++) {

					//e.g. if an element of triangleTable  is like {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
					//means first vertex is point on edge 0, then the second vertex is on 8-th edge...(the order how triangles are connected together)
					Vertex v1 = vertlist[triangleTable[cubeIndex][i]];
					v1.x = v1.x * 2 / dimx - 1;
					v1.y = v1.y * 2 / dimy - 1;
					v1.z = v1.z * 2 / dimz - 1;
					vertices->emplace_back(v1);
					normals->emplace_back(normaList[triangleTable[cubeIndex][i]]);
				}
			}
		}
	}

#ifdef OPENGL
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
#endif // OPENGL

}




void MarchingCubesDrawer::DrawSurfaceMesh() {
#ifdef OPENGL
	glBegin(GL_TRIANGLES);
	for (int i = 0; i < vertices->size(); i++) {
		Vertex& normal = normals->at(i);
		Vertex& vertex = vertices->at(i);
		glNormal3f(normal.x, normal.y, normal.z);
		glVertex3f(vertex.x, vertex.y, vertex.z);
	}
	glEnd();
#endif // OPENGL


}