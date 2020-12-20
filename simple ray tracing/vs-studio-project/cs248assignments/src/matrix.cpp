#include "matrix.h"

using std::endl;
using std::cout;
using std::istream;
const float EPS = 1e-10;
#ifdef VECTORBASE_Matrix
void Matrix::allocMem(int _x, int _y) {//初始化矩阵大小
	std::vector< std::vector<float> > temp(_x, std::vector<float>(_y));
	this->p = temp;
}
//声明一个全0矩阵
Matrix::Matrix(int rows, int cols)
{
	rows_num = rows;
	cols_num = cols;
	allocMem(rows, cols);
	for (int i = 0; i < rows_num; i++) {
		for (int j = 0; j < cols_num; j++) {
			p[i][j] = 0;
		}
	}
}
//声明一个值全部为value的矩阵
Matrix::Matrix(int rows, int cols, float value)
{
	rows_num = rows;
	cols_num = cols;
	allocMem(rows_num, cols_num);
	for (int i = 0; i < rows_num; i++) {
		for (int j = 0; j < cols_num; j++) {
			p[i][j] = value;
		}
	}
}
Matrix::Matrix(const Matrix& r)   // 拷贝构造函数
{
	rows_num = r.rows_num;
	cols_num = r.cols_num;
	this->allocMem(rows_num, cols_num);
	this->p = r.getData();

}

Matrix::Matrix(std::vector<std::vector<float> > dvec) {
	rows_num = dvec.size();
	cols_num = dvec[0].size();
	p = dvec;
}

std::vector<std::vector<float> >& Matrix::getData() {
	return p;
}

/**
 * 以 const 方式获取成员，可用于安全读
 */
const std::vector<std::vector<float> >& Matrix::getData() const {
	return p;
}
//实现矩阵的复制
Matrix& Matrix::operator=(const Matrix& m)
{
	if (this == &m) {
		return *this;
	}

	if (rows_num != m.rows_num || cols_num != m.cols_num) {
		rows_num = m.rows_num;
		cols_num = m.cols_num;
		allocMem(rows_num, cols_num);
	}
	p = m.getData();
	return *this;
}

//+=操作
Matrix& Matrix::operator+=(const Matrix& m)
{
	for (int i = 0; i < rows_num; i++) {
		for (int j = 0; j < cols_num; j++) {
			p[i][j] += m.p[i][j];
		}
	}
	return *this;
}
//实现-=
Matrix& Matrix::operator-=(const Matrix& m)
{
	for (int i = 0; i < rows_num; i++) {
		for (int j = 0; j < cols_num; j++) {
			p[i][j] -= m.p[i][j];
		}
	}
	return *this;
}

//解方程Ax=b
Matrix Matrix::Solve(const Matrix &Ain, const Matrix &bin)
{
	Matrix A(Ain);
	Matrix b(bin);
	//高斯消去法实现Ax=b的方程求解
	for (int i = 0; i < A.rows_num; i++) {
		if (A.p[i][i] == 0) {

			cout << "请重新输入" << endl;
		}
		for (int j = i + 1; j < A.rows_num; j++) {
			for (int k = i + 1; k < A.cols_num; k++) {
				A.p[j][k] -= A.p[i][k] * (A.p[j][i] / A.p[i][i]);
				if (abs(A.p[j][k]) < EPS)
					A.p[j][k] = 0;
			}
			b.p[j][0] -= b.p[i][0] * (A.p[j][i] / A.p[i][i]);
			if (abs(A.p[j][0]) < EPS)
				A.p[j][0] = 0;
			A.p[j][i] = 0;
		}
	}

	// 反向代换
	Matrix x(b.rows_num, 1);
	x.p[x.rows_num - 1][0] = b.p[x.rows_num - 1][0] / A.p[x.rows_num - 1][x.rows_num - 1];
	if (abs(x.p[x.rows_num - 1][0]) < EPS)
		x.p[x.rows_num - 1][0] = 0;
	for (int i = x.rows_num - 2; i >= 0; i--) {
		float sum = 0;
		for (int j = i + 1; j < x.rows_num; j++) {
			sum += A.p[i][j] * x.p[j][0];
		}
		x.p[i][0] = (b.p[i][0] - sum) / A.p[i][i];
		if (abs(x.p[i][0]) < EPS)
			x.p[i][0] = 0;
	}

	return x;
}

//矩阵显示
void Matrix::Show() const {
	//cout << rows_num <<" "<<cols_num<< endl;//显示矩阵的行数和列数
	for (int i = 0; i < rows_num; i++) {
		for (int j = 0; j < cols_num; j++) {
			cout << p[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}
//实现行变换
void Matrix::swapRows(int a, int b)
{
	a--;
	b--;
	auto temp = p[a];
	p[a] = p[b];
	p[b] = temp;
}
//计算矩阵行列式的值
float Matrix::det() {
	//为计算行列式做一个备份
	float ** back_up;
	back_up = new float *[rows_num];
	for (int i = 0; i < rows_num; i++) {
		back_up[i] = new float[cols_num];
	}
	for (int i = 0; i < rows_num; i++) {
		for (int j = 0; j < cols_num; j++) {
			back_up[i][j] = p[i][j];
		}
	}
	if (rows_num != cols_num) {
		std::abort();//只有方阵才能计算行列式，否则调用中断强制停止程序
	}
	float ans = 1;
	for (int i = 0; i < rows_num; i++) {
		//通过行变化的形式，使得矩阵对角线上的主元素不为0
		if (abs(p[i][i]) <= EPS) {
			bool flag = false;
			for (int j = 0; (j < cols_num) && (!flag); j++) {
				//若矩阵的一个对角线上的元素接近于0且能够通过行变换使得矩阵对角线上的元素不为0
				if ((abs(p[i][j]) > EPS) && (abs(p[j][i]) > EPS)) {
					flag = true;
					//注：进行互换后,p[i][j]变为p[j][j]，p[j][i]变为p[i][i]
					//对矩阵进行行变换
					float temp;
					for (int k = 0; k < cols_num; k++) {
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
	for (int i = 0; i < rows_num; i++) {
		for (int j = i + 1; j < rows_num; j++) {
			for (int k = i + 1; k < cols_num; k++) {
				p[j][k] -= p[i][k] * (p[j][i] * p[i][i]);
			}
		}
	}
	for (int i = 0; i < rows_num; i++) {
		ans *= p[i][i];
	}
	for (int i = 0; i < rows_num; i++) {
		for (int j = 0; j < cols_num; j++) {
			p[i][j] = back_up[i][j];
		}
	}
	return ans;
}
//返回矩阵第i行第j列的数
float Matrix::At(int i, int j) const {
	return this->p[i][j];
}
//求矩阵的逆矩阵
Matrix inv(Matrix A) {
	if (A.rows_num != A.cols_num) {
		std::cout << "只有方阵能求逆矩阵" << std::endl;
		std::abort();//只有方阵能求逆矩阵
	}
	float temp;
	Matrix A_B = Matrix(A.rows_num, A.cols_num);
	A_B = A;//为矩阵A做一个备份
	Matrix B = eye(A.rows_num);
	//将小于EPS的数全部置0
	for (int i = 0; i < A.rows_num; i++) {
		for (int j = 0; j < A.cols_num; j++) {
			if (abs(A.p[i][j]) <= EPS) {
				A.p[i][j] = 0;
			}
		}
	}
	//选择需要互换的两行选主元
	for (int i = 0; i < A.rows_num; i++) {
		if (abs(A.p[i][i]) <= EPS) {
			bool flag = false;
			for (int j = 0; (j < A.rows_num) && (!flag); j++) {
				if ((abs(A.p[i][j]) > EPS) && (abs(A.p[j][i]) > EPS)) {
					flag = true;
					for (int k = 0; k < A.cols_num; k++) {
						temp = A.p[i][k];
						A.p[i][k] = A.p[j][k];
						A.p[j][k] = temp;
						temp = B.p[i][k];
						B.p[i][k] = B.p[j][k];
						B.p[j][k] = temp;
					}
				}
			}
			if (!flag) {
				std::cout << "逆矩阵不存在\n";
				std::abort();
			}
		}
	}
	//通过初等行变换将A变为上三角矩阵
	float temp_rate;
	for (int i = 0; i < A.rows_num; i++) {
		for (int j = i + 1; j < A.rows_num; j++) {
			temp_rate = A.p[j][i] / A.p[i][i];
			for (int k = 0; k < A.cols_num; k++) {
				A.p[j][k] -= A.p[i][k] * temp_rate;
				B.p[j][k] -= B.p[i][k] * temp_rate;
			}
			A.p[j][i] = 0;
		}
	}
	//使对角元素均为1
	for (int i = 0; i < A.rows_num; i++) {
		temp = A.p[i][i];
		for (int j = 0; j < A.cols_num; j++) {
			A.p[i][j] /= temp;
			B.p[i][j] /= temp;
		}
	}
	//std::cout<<"算法可靠性检测，若可靠，输出上三角矩阵"<<std::endl;
	//将已经变为上三角矩阵的A，变为单位矩阵
	for (int i = A.rows_num - 1; i >= 1; i--) {
		for (int j = i - 1; j >= 0; j--) {
			temp = A.p[j][i];
			for (int k = 0; k < A.cols_num; k++) {
				A.p[j][k] -= A.p[i][k] * temp;
				B.p[j][k] -= B.p[i][k] * temp;
			}
		}
	}
	//std::cout << "算法可靠性检测，若可靠，输出单位矩阵" << std::endl;
	//for (int i = 0; i < A.rows_num; i++) {
	//	for (int j = 0; j < A.cols_num; j++) {
	//		printf("%7.4lf\t\t", A.p[i][j]);
	//	}
	//	cout << endl;
	//}
	A = A_B;//还原A
	return B;//返回该矩阵的逆矩阵
}

//读取矩阵行列数
int Matrix::row() const {
	return rows_num;
}
int Matrix::col() const {
	return cols_num;
}
//实现矩阵的转置
Matrix Matrix::Transpose()
{
	int col_size = this->cols_num;
	int row_size = this->rows_num;
	Matrix mt(col_size, row_size);
	for (int i = 0; i < row_size; i++) {
		for (int j = 0; j < col_size; j++) {
			mt.p[j][i] = this->p[i][j];
		}
	}
	return mt;
}
//高斯消元法
Matrix Matrix::gaussianEliminate()
{
	Matrix Ab(*this);
	int rows = Ab.rows_num;
	int cols = Ab.cols_num;
	int Acols = cols - 1;

	int i = 0; //跟踪行
	int j = 0; //跟踪列
	while (i < rows)
	{
		bool flag = false;
		while (j < Acols && !flag)
		{
			if (Ab.p[i][j] != 0) {
				flag = true;
			}
			else {
				int max_row = i;
				float max_val = 0;
				for (int k = i + 1; k < rows; ++k)
				{
					float cur_abs = Ab.p[k][j] >= 0 ? Ab.p[k][j] : -1 * Ab.p[k][j];
					if (cur_abs > max_val)
					{
						max_row = k;
						max_val = cur_abs;
					}
				}
				if (max_row != i) {
					Ab.swapRows(max_row, i);
					flag = true;
				}
				else {
					j++;
				}
			}
		}
		if (flag)
		{
			for (int t = i + 1; t < rows; t++) {
				for (int s = j + 1; s < cols; s++) {
					Ab.p[t][s] = Ab.p[t][s] - Ab.p[i][s] * (Ab.p[t][j] / Ab.p[i][j]);
					if (abs(Ab.p[t][s]) < EPS)
						Ab.p[t][s] = 0;
				}
				Ab.p[t][j] = 0;
			}
		}
		i++;
		j++;
	}
	return Ab;
}

//实现*=
Matrix& Matrix::operator*=(const Matrix& m)
{
	Matrix temp(rows_num, m.cols_num);//若C=AB,则矩阵C的行数等于矩阵A的行数，C的列数等于B的列数。
	for (int i = 0; i < temp.rows_num; i++) {
		for (int j = 0; j < temp.cols_num; j++) {
			for (int k = 0; k < cols_num; k++) {
				temp.p[i][j] += (p[i][k] * m.p[k][j]);
			}
		}
	}
	*this = temp;
	return *this;
}
//实现矩阵的乘法
Matrix Matrix::operator*(const Matrix & m)const {
	Matrix ba_M(rows_num, m.cols_num, 0.0);
	for (int i = 0; i < rows_num; i++) {
		for (int j = 0; j < m.cols_num; j++) {
			for (int k = 0; k < cols_num; k++) {
				ba_M.p[i][j] += (p[i][k] * m.p[k][j]);
			}
		}
	}
	return ba_M;
}

Matrix Matrix::operator*(float scale)const {
	Matrix ba_M(*this);
	for (int i = 0; i < rows_num; i++) {
		for (int j = 0; j < cols_num; j++) {
			ba_M.p[i][j] *= scale;
		}
	}
	return ba_M;

}
mVec4 Matrix::operator*(mVec4 V)const {
	mVec4 tmp(0.0f);

	tmp.x += this->p[0][0] * V.x;
	tmp.x += this->p[0][1] * V.y;
	tmp.x += this->p[0][2] * V.z;
	tmp.x += this->p[0][3] * V.w;


	tmp.y += this->p[1][0] * V.x;
	tmp.y += this->p[1][1] * V.y;
	tmp.y += this->p[1][2] * V.z;
	tmp.y += this->p[1][3] * V.w;


	tmp.z += this->p[2][0] * V.x;
	tmp.z += this->p[2][1] * V.y;
	tmp.z += this->p[2][2] * V.z;
	tmp.z += this->p[2][3] * V.w;


	tmp.w += this->p[3][0] * V.x;
	tmp.w += this->p[3][1] * V.y;
	tmp.w += this->p[3][2] * V.z;
	tmp.w += this->p[3][3] * V.w;
	return tmp;
}


//制造一个单位矩阵
Matrix  eye(int n) {
	Matrix A(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
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


Matrix scaleMatrix(float scalei) {
	return scaleMatrix(scalei, scalei, scalei);
}
Matrix scaleMatrix(float scalex, float scaley, float scalez) {
	Matrix tmp = eye(4);
	tmp.p[0][0] *= scalex;
	tmp.p[1][1] *= scaley;
	tmp.p[2][2] *= scalez;
	return tmp;
}
//@translation: mVec: xtrans, ytrans,ztrans
Matrix translateMatrix(mVec3 translation) {
	Matrix tmp = eye(4);
	tmp.p[0][3] = translation.x;
	tmp.p[1][3] = translation.y;
	tmp.p[2][3] = translation.z;
	return tmp;
}
//axis=0,1,2->x,y,z angle in degree,but cos sin function works for radiant 
Matrix rotateMatrix(int axis, float angle) {
	Matrix tmp;
	float RotateAngle = angle * PI / 180.0f;
	if (axis == 0) {
		auto datap = std::vector< std::vector<float>>({ {1., 0.,0.,0.} ,{0.0f,cos(RotateAngle),-sin(RotateAngle),0.} ,{ 0.,sin(RotateAngle),cos(RotateAngle),0. } ,{0., 0.,0.,1.} });
		tmp = Matrix(datap);

	}
	else if (axis == 1) {
		auto datap = std::vector< std::vector<float>>({ {cos(RotateAngle),0.0f,sin(RotateAngle),0.}  ,{0., 1.,0.,0.} ,{ -sin(RotateAngle),0.,cos(RotateAngle),0. } ,{0., 0.,0.,1.} });
		tmp = Matrix(datap);

	}
	else if (axis == 2) {
		auto datap = std::vector< std::vector<float>>({ {cos(RotateAngle),-sin(RotateAngle),0.,0.} ,{ sin(RotateAngle),cos(RotateAngle),0.,0} ,{0., 0.,1.,0.} ,{0., 0.,0.,1.} });
		tmp = Matrix(datap);

	}
	return tmp;
}

#endif








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
//实现行变换
//void Matrix4::swapRows(int a, int b)
//{
//	a--;
//	b--;
//	auto temp = p[a];
//	p[a] = p[b];
//	p[b] = temp;
//}
//计算矩阵行列式的值
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

//求矩阵的逆矩阵
Matrix4 inv(Matrix4 A) {
	float temp;
	Matrix4 A_B = Matrix4();
	A_B = A;//为矩阵A做一个备份
	Matrix4 B = eye(4);
	//将小于EPS的数全部置0
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (abs(A.p[i][j]) <= EPS) {
				A.p[i][j] = 0;
			}
		}
	}
	//选择需要互换的两行选主元
	for (int i = 0; i < 4; i++) {
		if (abs(A.p[i][i]) <= EPS) {
			bool flag = false;
			for (int j = 0; (j < 4) && (!flag); j++) {
				if ((abs(A.p[i][j]) > EPS) && (abs(A.p[j][i]) > EPS)) {
					flag = true;
					for (int k = 0; k < 4; k++) {
						temp = A.p[i][k];
						A.p[i][k] = A.p[j][k];
						A.p[j][k] = temp;
						temp = B.p[i][k];
						B.p[i][k] = B.p[j][k];
						B.p[j][k] = temp;
					}
				}
			}
			if (!flag) {
				std::cout << "逆矩阵不存在\n";
				std::abort();
			}
		}
	}
	//通过初等行变换将A变为上三角矩阵
	float temp_rate;
	for (int i = 0; i < 4; i++) {
		for (int j = i + 1; j < 4; j++) {
			temp_rate = A.p[j][i] / A.p[i][i];
			for (int k = 0; k < 4; k++) {
				A.p[j][k] -= A.p[i][k] * temp_rate;
				B.p[j][k] -= B.p[i][k] * temp_rate;
			}
			A.p[j][i] = 0;
		}
	}
	//使对角元素均为1
	for (int i = 0; i < 4; i++) {
		temp = A.p[i][i];
		for (int j = 0; j < 4; j++) {
			A.p[i][j] /= temp;
			B.p[i][j] /= temp;
		}
	}
	//std::cout<<"算法可靠性检测，若可靠，输出上三角矩阵"<<std::endl;
	//将已经变为上三角矩阵的A，变为单位矩阵
	for (int i = 4 - 1; i >= 1; i--) {
		for (int j = i - 1; j >= 0; j--) {
			temp = A.p[j][i];
			for (int k = 0; k <4; k++) {
				A.p[j][k] -= A.p[i][k] * temp;
				B.p[j][k] -= B.p[i][k] * temp;
			}
		}
	}
	//std::cout << "算法可靠性检测，若可靠，输出单位矩阵" << std::endl;
	//for (int i = 0; i < A.rows_num; i++) {
	//	for (int j = 0; j < A.cols_num; j++) {
	//		printf("%7.4lf\t\t", A.p[i][j]);
	//	}
	//	cout << endl;
	//}
	A = A_B;//还原A
	return B;//返回该矩阵的逆矩阵
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
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				ba_M.p[i][j] += (p[i][k] * m.p[k][j]);
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


//制造一个单位矩阵
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