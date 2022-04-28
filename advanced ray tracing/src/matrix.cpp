#include "matrix.h"

using std::endl;
using std::cout;
using std::istream;
const double EPS = 1e-20;



//����һ��ȫ0����
Matrix4::Matrix4()
{

	p[0][0] = 0; p[0][1] = 0; p[0][2] = 0; p[0][3] = 0;
	p[1][0] = 0; p[1][1] = 0; p[1][2] = 0; p[1][3] = 0;
	p[2][0] = 0; p[2][1] = 0; p[2][2] = 0; p[2][3] = 0;
	p[3][0] = 0; p[3][1] = 0; p[3][2] = 0; p[3][3] = 0;
}
//����һ��ֵȫ��Ϊvalue�ľ���
Matrix4::Matrix4( float value)
{
	p[0][0] = value; p[0][1] = value; p[0][2] = value; p[0][3] = value;
	p[1][0] = value; p[1][1] = value; p[1][2] = value; p[1][3] = value;
	p[2][0] = value; p[2][1] = value; p[2][2] = value; p[2][3] = value;
	p[3][0] = value; p[3][1] = value; p[3][2] = value; p[3][3] = value;
}
Matrix4::Matrix4(const Matrix4& r)   // �������캯��
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

//ʵ�־���ĸ���
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

//+=����
Matrix4& Matrix4::operator+=(const Matrix4& m)
{
	p[0][0] += m.p[0][0]; p[0][1] += m.p[0][1]; p[0][2] += m.p[0][2]; p[0][3] += m.p[0][3];
	p[1][0] += m.p[1][0]; p[1][1] += m.p[1][1]; p[1][2] += m.p[1][2]; p[1][3] += m.p[1][3];
	p[2][0] += m.p[2][0]; p[2][1] += m.p[2][1]; p[2][2] += m.p[2][2]; p[2][3] += m.p[2][3];
	p[3][0] += m.p[3][0]; p[3][1] += m.p[3][1]; p[3][2] += m.p[3][2]; p[3][3] += m.p[3][3];
	return *this;
}
//ʵ��-=
Matrix4& Matrix4::operator-=(const Matrix4& m)
{
	p[0][0] -= m.p[0][0]; p[0][1] -= m.p[0][1]; p[0][2] -= m.p[0][2]; p[0][3] -= m.p[0][3];
	p[1][0] -= m.p[1][0]; p[1][1] -= m.p[1][1]; p[1][2] -= m.p[1][2]; p[1][3] -= m.p[1][3];
	p[2][0] -= m.p[2][0]; p[2][1] -= m.p[2][1]; p[2][2] -= m.p[2][2]; p[2][3] -= m.p[2][3];
	p[3][0] -= m.p[3][0]; p[3][1] -= m.p[3][1]; p[3][2] -= m.p[3][2]; p[3][3] -= m.p[3][3];
	return *this;
}




//������ʾ
void Matrix4::Show() const {
	//cout << rows_num <<" "<<cols_num<< endl;//��ʾ���������������
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << p[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

float Matrix4::det() {
	//Ϊ��������ʽ��һ������
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
		//ͨ���б仯����ʽ��ʹ�þ���Խ����ϵ���Ԫ�ز�Ϊ0
		if (abs(p[i][i]) <= EPS) {
			bool flag = false;
			for (int j = 0; (j < 4) && (!flag); j++) {
				//�������һ���Խ����ϵ�Ԫ�ؽӽ���0���ܹ�ͨ���б任ʹ�þ���Խ����ϵ�Ԫ�ز�Ϊ0
				if ((abs(p[i][j]) > EPS) && (abs(p[j][i]) > EPS)) {
					flag = true;
					//ע�����л�����,p[i][j]��Ϊp[j][j]��p[j][i]��Ϊp[i][i]
					//�Ծ�������б任
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

//�����������-> ���ַ������˺ü������ͣ������н�ѧ�ٷ�ʹ�ã���ʵ������һЩ����»������⣬��������ͻȻ���䣨��M-1���Թ��ߣ�
//Matrix4 inv(Matrix4 A) {
//	double temp;
//	Matrix4 A_B = Matrix4();
//	A_B = A;//Ϊ����A��һ������
//	Matrix4 B = eye(4);
//	//��С��EPS����ȫ����0
//	for (int i = 0; i < 4; i++) {
//		for (int j = 0; j < 4; j++) {
//			if (abs(A.p[i][j]) <= EPS) {
//				A.p[i][j] = 0;
//			}
//		}
//	}
//	//ѡ����Ҫ����������ѡ��Ԫ
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
//				std::cout << "����󲻴���\n";
//				std::abort();
//			}
//		}
//	}
//	//ͨ�������б任��A��Ϊ�����Ǿ���
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
//	//ʹ�Խ�Ԫ�ؾ�Ϊ1
//	for (int i = 0; i < 4; i++) {
//		temp = A.p[i][i];
//		for (int j = 0; j < 4; j++) {
//			A.p[i][j] /= temp;
//			B.p[i][j] /= temp;
//		}
//	}
//	//std::cout<<"�㷨�ɿ��Լ�⣬���ɿ�����������Ǿ���"<<std::endl;
//	//���Ѿ���Ϊ�����Ǿ����A����Ϊ��λ����
//	for (int i = 4 - 1; i >= 1; i--) {
//		for (int j = i - 1; j >= 0; j--) {
//			temp = A.p[j][i];
//			for (int k = 0; k <4; k++) {
//				A.p[j][k] -= A.p[i][k] * temp;
//				B.p[j][k] -= B.p[i][k] * temp;
//			}
//		}
//	}
//	A = A_B;//��ԭA
//	return B;//���ظþ���������
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






//ʵ�־����ת��
Matrix4 Matrix4::Transpose()
{
	Matrix4 mt;
	mt.p[0][0] = p[0][0]; mt.p[0][1] = p[1][0]; mt.p[0][2] = p[2][0]; mt.p[0][3] = p[3][0];
	mt.p[1][0] = p[0][1]; mt.p[1][1] = p[1][1]; mt.p[1][2] = p[2][1]; mt.p[1][3] = p[3][1];
	mt.p[2][0] = p[0][2]; mt.p[2][1] = p[1][2]; mt.p[2][2] = p[2][2]; mt.p[2][3] = p[3][2];
	mt.p[3][0] = p[0][0]; mt.p[3][1] = p[1][3]; mt.p[3][2] = p[2][3]; mt.p[3][3] = p[3][3];
	return mt;
}



//ʵ��*=
Matrix4& Matrix4::operator*=(const Matrix4& m)
{
	Matrix4 temp;//��C=AB,�����C���������ھ���A��������C����������B��������
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
//ʵ�־���ĳ˷�
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