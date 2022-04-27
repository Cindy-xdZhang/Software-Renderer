#include "matrix.h"

using std::endl;
using std::cout;
using std::istream;
const double EPS = 1e-20;



//声明一个全0矩阵
Matrix4::Matrix4()
{

	p[0][0] = 0; p[0][1] = 0; p[0][2] = 0; p[0][3] = 0;
	p[1][0] = 0; p[1][1] = 0; p[1][2] = 0; p[1][3] = 0;
	p[2][0] = 0; p[2][1] = 0; p[2][2] = 0; p[2][3] = 0;
	p[3][0] = 0; p[3][1] = 0; p[3][2] = 0; p[3][3] = 0;
}
//声明一个值全部为value的矩阵
Matrix4::Matrix4( float value)
{
	p[0][0] = value; p[0][1] = value; p[0][2] = value; p[0][3] = value;
	p[1][0] = value; p[1][1] = value; p[1][2] = value; p[1][3] = value;
	p[2][0] = value; p[2][1] = value; p[2][2] = value; p[2][3] = value;
	p[3][0] = value; p[3][1] = value; p[3][2] = value; p[3][3] = value;
}
Matrix4::Matrix4(const Matrix4& r)   // 拷贝构造函数
{
	p[0][0] = r.p[0][0]; p[0][1] = r.p[0][1]; p[0][2] = r.p[0][2]; p[0][3] = r.p[0][3];
	p[1][0] = r.p[1][0]; p[1][1] = r.p[1][1]; p[1][2] = r.p[1][2]; p[1][3] = r.p[1][3];
	p[2][0] = r.p[2][0]; p[2][1] = r.p[2][1]; p[2][2] = r.p[2][2]; p[2][3] = r.p[2][3];
	p[3][0] = r.p[3][0]; p[3][1] = r.p[3][1]; p[3][2] = r.p[3][2]; p[3][3] = r.p[3][3];

}

Matrix4::Matrix4(std::vector<std::vector<float> > dvec) {
	p[0][0] = dvec[0][0]; p[0][1] = dvec[0][1]; p[0][2] = dvec[0][2]; p[0][3] = dvec[0][3];
	p[1][0] = dvec[1][0]; p[1][1] = dvec[1][1]; p[1][2] = dvec[1][2]; p[1][3] = dvec[1][3];
	p[2][0] = dvec[2][0]; p[2][1] = dvec[2][1]; p[2][2] = dvec[2][2]; p[2][3] = dvec[2][3];
	p[3][0] = dvec[3][0]; p[3][1] = dvec[3][1]; p[3][2] = dvec[3][2]; p[3][3] = dvec[3][3];
}

//实现矩阵的复制
Matrix4& Matrix4::operator=(const Matrix4& m)
{
	if (this == &m) {
		return *this;
	}
	p[0][0] = m.p[0][0]; p[0][1] = m.p[0][1]; p[0][2] = m.p[0][2]; p[0][3] = m.p[0][3];
	p[1][0] = m.p[1][0]; p[1][1] = m.p[1][1]; p[1][2] = m.p[1][2]; p[1][3] = m.p[1][3];
	p[2][0] = m.p[2][0]; p[2][1] = m.p[2][1]; p[2][2] = m.p[2][2]; p[2][3] = m.p[2][3];
	p[3][0] = m.p[3][0]; p[3][1] = m.p[3][1]; p[3][2] = m.p[3][2]; p[3][3] = m.p[3][3];

	return *this;
}

//+=操作
Matrix4& Matrix4::operator+=(const Matrix4& m)
{
	p[0][0] += m.p[0][0]; p[0][1] += m.p[0][1]; p[0][2] += m.p[0][2]; p[0][3] += m.p[0][3];
	p[1][0] += m.p[1][0]; p[1][1] += m.p[1][1]; p[1][2] += m.p[1][2]; p[1][3] += m.p[1][3];
	p[2][0] += m.p[2][0]; p[2][1] += m.p[2][1]; p[2][2] += m.p[2][2]; p[2][3] += m.p[2][3];
	p[3][0] += m.p[3][0]; p[3][1] += m.p[3][1]; p[3][2] += m.p[3][2]; p[3][3] += m.p[3][3];
	return *this;
}
//实现-=
Matrix4& Matrix4::operator-=(const Matrix4& m)
{
	p[0][0] -= m.p[0][0]; p[0][1] -= m.p[0][1]; p[0][2] -= m.p[0][2]; p[0][3] -= m.p[0][3];
	p[1][0] -= m.p[1][0]; p[1][1] -= m.p[1][1]; p[1][2] -= m.p[1][2]; p[1][3] -= m.p[1][3];
	p[2][0] -= m.p[2][0]; p[2][1] -= m.p[2][1]; p[2][2] -= m.p[2][2]; p[2][3] -= m.p[2][3];
	p[3][0] -= m.p[3][0]; p[3][1] -= m.p[3][1]; p[3][2] -= m.p[3][2]; p[3][3] -= m.p[3][3];
	return *this;
}




//矩阵显示
void Matrix4::Show() const {
	//cout << rows_num <<" "<<cols_num<< endl;//显示矩阵的行数和列数
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << p[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

float Matrix4::det() {
	//为计算行列式做一个备份
	float ** back_up;
	back_up = new float *[4];
	for (int i = 0; i < 4; i++) {
		back_up[i] = new float[4];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			back_up[i][j] = p[i][j];
		}
	}
	float ans = 1;
	for (int i = 0; i < 4; i++) {
		//通过行变化的形式，使得矩阵对角线上的主元素不为0
		if (abs(p[i][i]) <= EPS) {
			bool flag = false;
			for (int j = 0; (j < 4) && (!flag); j++) {
				//若矩阵的一个对角线上的元素接近于0且能够通过行变换使得矩阵对角线上的元素不为0
				if ((abs(p[i][j]) > EPS) && (abs(p[j][i]) > EPS)) {
					flag = true;
					//注：进行互换后,p[i][j]变为p[j][j]，p[j][i]变为p[i][i]
					//对矩阵进行行变换
					float temp;
					for (int k = 0; k < 4; k++) {
						temp = p[i][k];
						p[i][k] = p[j][k];
						p[j][k] = temp;
					}
				}
			}
			if (flag)
				return 0;
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = i + 1; j < 4; j++) {
			for (int k = i + 1; k < 4; k++) {
				p[j][k] -= p[i][k] * (p[j][i] * p[i][i]);
			}
		}
	}
	for (int i = 0; i < 4; i++) {
		ans *= p[i][i];
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			p[i][j] = back_up[i][j];
		}
	}
	return ans;
}

//求矩阵的逆矩阵-> 这种方法看了好几个博客，甚至有教学官方使用，但实际上在一些情况下会有问题，导致物体突然畸变（用M-1乘以光线）
//Matrix4 inv(Matrix4 A) {
//	double temp;
//	Matrix4 A_B = Matrix4();
//	A_B = A;//为矩阵A做一个备份
//	Matrix4 B = eye(4);
//	//将小于EPS的数全部置0
//	for (int i = 0; i < 4; i++) {
//		for (int j = 0; j < 4; j++) {
//			if (abs(A.p[i][j]) <= EPS) {
//				A.p[i][j] = 0;
//			}
//		}
//	}
//	//选择需要互换的两行选主元
//	for (int i = 0; i < 4; i++) {
//		if (abs(A.p[i][i]) <= EPS) {
//			bool flag = false;
//			for (int j = 0; (j < 4) && (!flag); j++) {
//				if ((abs(A.p[i][j]) > EPS) && (abs(A.p[j][i]) > EPS)) {
//					flag = true;
//					for (int k = 0; k < 4; k++) {
//						temp = A.p[i][k];
//						A.p[i][k] = A.p[j][k];
//						A.p[j][k] = temp;
//						temp = B.p[i][k];
//						B.p[i][k] = B.p[j][k];
//						B.p[j][k] = temp;
//					}
//				}
//			}
//			if (!flag) {
//				std::cout << "逆矩阵不存在\n";
//				std::abort();
//			}
//		}
//	}
//	//通过初等行变换将A变为上三角矩阵
//	float temp_rate;
//	for (int i = 0; i < 4; i++) {
//		for (int j = i + 1; j < 4; j++) {
//			temp_rate = A.p[j][i] / A.p[i][i];
//			for (int k = 0; k < 4; k++) {
//				A.p[j][k] -= A.p[i][k] * temp_rate;
//				B.p[j][k] -= B.p[i][k] * temp_rate;
//			}
//			A.p[j][i] = 0;
//		}
//	}
//	//使对角元素均为1
//	for (int i = 0; i < 4; i++) {
//		temp = A.p[i][i];
//		for (int j = 0; j < 4; j++) {
//			A.p[i][j] /= temp;
//			B.p[i][j] /= temp;
//		}
//	}
//	//std::cout<<"算法可靠性检测，若可靠，输出上三角矩阵"<<std::endl;
//	//将已经变为上三角矩阵的A，变为单位矩阵
//	for (int i = 4 - 1; i >= 1; i--) {
//		for (int j = i - 1; j >= 0; j--) {
//			temp = A.p[j][i];
//			for (int k = 0; k <4; k++) {
//				A.p[j][k] -= A.p[i][k] * temp;
//				B.p[j][k] -= B.p[i][k] * temp;
//			}
//		}
//	}
//	A = A_B;//还原A
//	return B;//返回该矩阵的逆矩阵
//}

Matrix4 inv(const Matrix4 AB)
{
	float m[16];
	float Output[16];
	memcpy(m, &AB.p[0], 16*sizeof (float));

	double inv[16], det;
	int i;
		inv[0] = m[5] * m[10] * m[15] -
			m[5] * m[11] * m[14] -
			m[9] * m[6] * m[15] +
			m[9] * m[7] * m[14] +
			m[13] * m[6] * m[11] -
			m[13] * m[7] * m[10];

		inv[4] = -m[4] * m[10] * m[15] +
			m[4] * m[11] * m[14] +
			m[8] * m[6] * m[15] -
			m[8] * m[7] * m[14] -
			m[12] * m[6] * m[11] +
			m[12] * m[7] * m[10];

		inv[8] = m[4] * m[9] * m[15] -
			m[4] * m[11] * m[13] -
			m[8] * m[5] * m[15] +
			m[8] * m[7] * m[13] +
			m[12] * m[5] * m[11] -
			m[12] * m[7] * m[9];

		inv[12] = -m[4] * m[9] * m[14] +
			m[4] * m[10] * m[13] +
			m[8] * m[5] * m[14] -
			m[8] * m[6] * m[13] -
			m[12] * m[5] * m[10] +
			m[12] * m[6] * m[9];

		inv[1] = -m[1] * m[10] * m[15] +
			m[1] * m[11] * m[14] +
			m[9] * m[2] * m[15] -
			m[9] * m[3] * m[14] -
			m[13] * m[2] * m[11] +
			m[13] * m[3] * m[10];

		inv[5] = m[0] * m[10] * m[15] -
			m[0] * m[11] * m[14] -
			m[8] * m[2] * m[15] +
			m[8] * m[3] * m[14] +
			m[12] * m[2] * m[11] -
			m[12] * m[3] * m[10];

		inv[9] = -m[0] * m[9] * m[15] +
			m[0] * m[11] * m[13] +
			m[8] * m[1] * m[15] -
			m[8] * m[3] * m[13] -
			m[12] * m[1] * m[11] +
			m[12] * m[3] * m[9];

		inv[13] = m[0] * m[9] * m[14] -
			m[0] * m[10] * m[13] -
			m[8] * m[1] * m[14] +
			m[8] * m[2] * m[13] +
			m[12] * m[1] * m[10] -
			m[12] * m[2] * m[9];

		inv[2] = m[1] * m[6] * m[15] -
			m[1] * m[7] * m[14] -
			m[5] * m[2] * m[15] +
			m[5] * m[3] * m[14] +
			m[13] * m[2] * m[7] -
			m[13] * m[3] * m[6];

		inv[6] = -m[0] * m[6] * m[15] +
			m[0] * m[7] * m[14] +
			m[4] * m[2] * m[15] -
			m[4] * m[3] * m[14] -
			m[12] * m[2] * m[7] +
			m[12] * m[3] * m[6];

		inv[10] = m[0] * m[5] * m[15] -
			m[0] * m[7] * m[13] -
			m[4] * m[1] * m[15] +
			m[4] * m[3] * m[13] +
			m[12] * m[1] * m[7] -
			m[12] * m[3] * m[5];

		inv[14] = -m[0] * m[5] * m[14] +
			m[0] * m[6] * m[13] +
			m[4] * m[1] * m[14] -
			m[4] * m[2] * m[13] -
			m[12] * m[1] * m[6] +
			m[12] * m[2] * m[5];

		inv[3] = -m[1] * m[6] * m[11] +
			m[1] * m[7] * m[10] +
			m[5] * m[2] * m[11] -
			m[5] * m[3] * m[10] -
			m[9] * m[2] * m[7] +
			m[9] * m[3] * m[6];

		inv[7] = m[0] * m[6] * m[11] -
			m[0] * m[7] * m[10] -
			m[4] * m[2] * m[11] +
			m[4] * m[3] * m[10] +
			m[8] * m[2] * m[7] -
			m[8] * m[3] * m[6];

		inv[11] = -m[0] * m[5] * m[11] +
			m[0] * m[7] * m[9] +
			m[4] * m[1] * m[11] -
			m[4] * m[3] * m[9] -
			m[8] * m[1] * m[7] +
			m[8] * m[3] * m[5];
		inv[15] = m[0] * m[5] * m[10] -
			m[0] * m[6] * m[9] -
			m[4] * m[1] * m[10] +
			m[4] * m[2] * m[9] +
			m[8] * m[1] * m[6] -
			m[8] * m[2] * m[5];

		det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

		if (det == 0)
			return false;

		det = 1.0 / det;

		for (i = 0; i < 16; i++)
			Output[i] = inv[i] * det;

		Matrix4 B=eye(0);
		memcpy(&B.p[0], Output, 16 * sizeof(float));
		return B;
}






//实现矩阵的转置
Matrix4 Matrix4::Transpose()
{
	Matrix4 mt;
	mt.p[0][0] = p[0][0]; mt.p[0][1] = p[1][0]; mt.p[0][2] = p[2][0]; mt.p[0][3] = p[3][0];
	mt.p[1][0] = p[0][1]; mt.p[1][1] = p[1][1]; mt.p[1][2] = p[2][1]; mt.p[1][3] = p[3][1];
	mt.p[2][0] = p[0][2]; mt.p[2][1] = p[1][2]; mt.p[2][2] = p[2][2]; mt.p[2][3] = p[3][2];
	mt.p[3][0] = p[0][0]; mt.p[3][1] = p[1][3]; mt.p[3][2] = p[2][3]; mt.p[3][3] = p[3][3];
	return mt;
}



//实现*=
Matrix4& Matrix4::operator*=(const Matrix4& m)
{
	Matrix4 temp;//若C=AB,则矩阵C的行数等于矩阵A的行数，C的列数等于B的列数。
	for (int i = 0; i <4; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				temp.p[i][j] += (p[i][k] * m.p[k][j]);
			}
		}
	}
	*this = temp;
	return *this;
}
//实现矩阵的乘法
Matrix4 Matrix4::operator*(const Matrix4 & m)const {
	Matrix4 ba_M(0.0f);
	//optimization from :https://zhuanlan.zhihu.com/p/146250334?from_voters_page=true
	for (int i = 0; i < 4; i++) {
		for (int k = 0; k < 4; k++) {
			float s = p[i][k];
			for (int j = 0; j < 4; j++) {
				ba_M.p[i][j] += ( s* m.p[k][j]);
			}
		}
	}
	return ba_M;
}

Matrix4 Matrix4::operator*(float value)const {
	Matrix4 ba_M(*this);
	ba_M.p[0][0] *= value; ba_M.p[0][1] *= value; ba_M.p[0][2] *= value; ba_M.p[0][3] *= value;
	ba_M.p[1][0] *= value; ba_M.p[1][1] *= value; ba_M.p[1][2] *= value; ba_M.p[1][3] *= value;
	ba_M.p[2][0] *= value; ba_M.p[2][1] *= value; ba_M.p[2][2] *= value; ba_M.p[2][3] *= value;
	ba_M.p[3][0] *= value; ba_M.p[3][1] *= value; ba_M.p[3][2] *= value; ba_M.p[3][3] *= value;
	return ba_M;

}
mVec4 Matrix4::operator*(mVec4 V)const {
	float a=0.0f, b = 0.0f, c = 0.0f, d = 0.0f;

	a += this->p[0][0] * V.x;
	a += this->p[0][1] * V.y;
	a += this->p[0][2] * V.z;
	a += this->p[0][3] * V.w;


	b += this->p[1][0] * V.x;
	b += this->p[1][1] * V.y;
	b += this->p[1][2] * V.z;  
	b += this->p[1][3] * V.w;


	c += this->p[2][0] * V.x;
	c += this->p[2][1] * V.y;
	c += this->p[2][2] * V.z;
	c += this->p[2][3] * V.w;


	d += this->p[3][0] * V.x;
	d += this->p[3][1] * V.y;
	d += this->p[3][2] * V.z;
	d += this->p[3][3] * V.w;
	return mVec4(a,b,c,d);
}



Matrix4  eye(int n) {
	Matrix4 A(0.0f);
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j) {
				A.p[i][j] = 1;
			}
			else {
				A.p[i][j] = 0;
			}
		}
	}
	return A;
}


Matrix4 scaleMatrix(float scalei) {
	return scaleMatrix(scalei, scalei, scalei);
}
Matrix4 scaleMatrix(float scalex, float scaley, float scalez) {
	Matrix4 tmp = eye(4);
	tmp.p[0][0] *= scalex;
	tmp.p[1][1] *= scaley;
	tmp.p[2][2] *= scalez;
	return tmp;
}
//@translation: mVec: xtrans, ytrans,ztrans
Matrix4 translateMatrix(mVec3 translation) {
	Matrix4 tmp = eye(4);
	tmp.p[0][3] = translation.x;
	tmp.p[1][3] = translation.y;
	tmp.p[2][3] = translation.z;
	return tmp;
}
//axis=0,1,2->x,y,z angle in degree,but cos sin function works for radiant 
Matrix4 rotateMatrix(int axis, float angle) {
	Matrix4 tmp;
	float RotateAngle = angle * PI / 180.0f;
	if (axis == 0) {
		auto datap = std::vector< std::vector<float>>({ {1., 0.,0.,0.} ,{0.0f,cos(RotateAngle),-sin(RotateAngle),0.} ,{ 0.,sin(RotateAngle),cos(RotateAngle),0. } ,{0., 0.,0.,1.} });
		tmp = Matrix4(datap);

	}
	else if (axis == 1) {
		auto datap = std::vector< std::vector<float>>({ {cos(RotateAngle),0.0f,sin(RotateAngle),0.}  ,{0., 1.,0.,0.} ,{ -sin(RotateAngle),0.,cos(RotateAngle),0. } ,{0., 0.,0.,1.} });
		tmp = Matrix4(datap);

	}
	else if (axis == 2) {
		auto datap = std::vector< std::vector<float>>({ {cos(RotateAngle),-sin(RotateAngle),0.,0.} ,{ sin(RotateAngle),cos(RotateAngle),0.,0} ,{0., 0.,1.,0.} ,{0., 0.,0.,1.} });
		tmp = Matrix4(datap);

	}
	return tmp;
}


Matrix4 ViewPortMatrix(int Nx, int Ny) {

	auto datap = std::vector< std::vector<float>>({ \
	{ float(Nx) / 2, 0.,0.,float(Nx - 1) / 2} ,\
	{0.0f,float(Ny) / 2,0.0f, float(Ny - 1) / 2} ,\
	{ 0,0,1,0 } ,\
	{0, 0,0,1} });
	return Matrix4(datap);
}
Matrix4 OrthoMatrix(float r, float l, float t, float b, float n, float f) {
	auto datap = std::vector< std::vector<float>>({ \
	{ float(2 / (r - l)) , 0, 0, float(-(r + l) / (r - l))} ,\
	{ 0,float(2 / (t - b)), 0, float(-(t + b) / (t - b)) },\
	{ 0,0,float(2 / (n - f)), float(-(n + f) / (n - f)) } ,\
	{0, 0,0,1} });

	return Matrix4(datap);

}

Matrix4 PerspectiveMatrix(float r, float l, float t, float b, float n, float f) {

	auto datap = std::vector< std::vector<float>>({ \
	{ float((2 * n) / (r - l)), 0,  float((r + l) / (l - r)), 0} ,\
	{ 0,float((2 * n) / (t - b)),  float((t + b) / (b - t)),0 },\
	{ 0,0, float((n + f) / (n - f)),float((2 * f*n) / (f - n)) } ,\
	{0, 0,1,0} });
	return Matrix4(datap);



}

Matrix4 ViewMatrix(mVec3 eyePos, mVec3 GazeDirection, mVec3 TopDirection) {
	mVec3 GxT = GazeDirection.cross_product(TopDirection);
	mVec3 minus_G = GazeDirection * -1.0f;
	mVec3 minus_eye = eyePos * -1.0f;
	auto dataRotate_View = std::vector< std::vector<float>>({ \
	{ GxT.x, GxT.y, GxT.z,0} ,\
	{TopDirection.x, TopDirection.y, TopDirection.z,0} ,\
	{ minus_G.x, minus_G.y, minus_G.z,0 } ,\
	{0, 0,0,1} });
	Matrix4 R_view(dataRotate_View);
	auto dataRotate_T = std::vector< std::vector<float>>({ \
	{1, 0, 0, minus_eye.x}, \
	{0, 1, 0, minus_eye.y}, \
	{0, 0, 1, minus_eye.z}, \
	{0, 0, 0, 1} });
	Matrix4 T_view(dataRotate_T);

	return R_view * T_view;
}


Matrix4 rotateMatrix(mVec3 n, float theta) {

	float SinTheta = sin(theta);
	float CosTheta = cos(theta);
	auto datap = std::vector< std::vector<float>>({ \
{ n.x*n.x*(1 - CosTheta) + CosTheta,\
  n.x*n.y*(1 - CosTheta) + n.z*SinTheta,\
  n.x*n.z*(1 - CosTheta) - n.y*SinTheta, \
0 } ,\
{ n.x*n.y*(1 - CosTheta) - n.z*SinTheta, \
  n.y*n.y*(1 - CosTheta) + CosTheta, \
  n.y*n.z*(1 - CosTheta) + n.x*SinTheta, \
0 }, \
{ n.z*n.x*(1 - CosTheta) + n.y*SinTheta, \
  n.z*n.y*(1 - CosTheta) - n.x*SinTheta, \
  n.z*n.z*(1 - CosTheta) + CosTheta, \
		0 }, \
	{0, 0, 0, 1	}
});

	return Matrix4(datap);


}

Matrix4 squashMatrix(float n, float f) {
	auto datap = std::vector< std::vector<float>>({ \
	{ n, 0,  0, 0} ,\
	{ 0,n,  0,0 },\
	{ 0,0, n + f, -f*n } ,\
	{0, 0,1,0} });
	return Matrix4(datap);

}