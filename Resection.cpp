#include "Resection.h"

Resection::Resection(double* data,int row,int col, double mX0, double mY0, double mF, double mM)
{
	DataMatrix=Mat(row, col, CV_64FC1, data);
	X0 = mX0; Y0 = mY0; F = mF; m = mM;
	Xs = 0.0, Ys = 0.0, Zs = 0.0, Q = 0.0, W = 0.0, K = 0.0;
	R = Mat(3, 3, CV_64FC1);
	L = Mat(8, 1, CV_64FC1);
	A = Mat(8, 6, CV_64FC1);
	AXG = Mat(8, 1, CV_64FC1);
	V = Mat(8, 1, CV_64FC1);
	VT = Mat(1, 8, CV_64FC1);
	VTV = Mat(1, 1, CV_64FC1);
	D = Mat(6, 6, CV_64FC1);
	XYZtemp = Mat(3, 1, CV_64FC1);
}

Resection::Resection(double** data, int row, int col, double mX0, double mY0, double mF,double mM)
{
	DataMatrix=Mat(row, col, CV_64FC1, data);
	X0 = mX0; Y0 = mY0; F = mF; m = mM;
	Xs = 0.0, Ys = 0.0, Zs = 0.0, Q = 0.0, W = 0.0, K = 0.0;
	R=Mat(3, 3, CV_64FC1);
	L=Mat(8, 1, CV_64FC1);
	A=Mat(8, 6, CV_64FC1);
	AXG = Mat(8, 1, CV_64FC1);
	V = Mat(8, 1, CV_64FC1);
	VT = Mat(1, 8, CV_64FC1);
	VTV = Mat(1, 1, CV_64FC1);
	D = Mat(6, 6, CV_64FC1);
	XYZtemp=Mat(3, 1, CV_64FC1);
}

Mat Resection::�Y��(double deltaDis, double deltaAng)
{
	//δ֪���γ�����Q�����ޤ���
	for (int i = 0; i < 4; i++) {
		Xs += DataMatrix.at<double>(i, 2);
		Ys += DataMatrix.at<double>(i, 3);
		Zs += DataMatrix.at<double>(i, 4);
	}
	Xs /= 4; Ys /= 4; Zs /= 4 + m * F;
	int count = 0;					//���ͥ�����
	do									//Ӌ����R�귵��
	{
		������ܞ����(R, Q, W, K);
		//�S�����Фȶ���헤�Ӌ��

		for (int i = 0, j = 0, k = 0; i < 4; i++, k++, j++) {
			//�¤���XYZ�΂���ȡ�ä��ޤ���
			�������ˤ�ä�(XYZtemp, DataMatrix, i, Xs, Ys, Zs);
			Mat temp = �خ������Ё\��(R, XYZtemp);
			X = temp.at<double>(0, 0); Y = temp.at<double>(1, 0); Z = temp.at<double>(2, 0);
			//����ѥ��`������L
			L.at<double>(j, 0) = DataMatrix.at<double>(i, 0) - (X0 - F * X / Z);
			L.at<double>(j + 1, 0) = DataMatrix.at<double>(i, 1) - (Y0 - F * Y / Z);
			j++;
			�ѥ��`���ޥȥ�å�����ȡ��(A, R, DataMatrix, k, i, F, K, W, Z, X0, Y0);
			k++;
		}
		transpose(A, AT);
		
		//�g��ƽ����Ф���δ֪����⤤�Ƥ���������
		ATA = AT * A;
		ATA = ATA.inv();
		ATL = AT * L;
		XG = ATA * ATL;
		Xs = Xs + XG.at<double>(0, 0);
		Ys = Ys + XG.at<double>(1, 0);
		Zs = Zs + XG.at<double>(2, 0);
		Q = Q + XG.at<double>(3, 0);
		W = W + XG.at<double>(4, 0);
		K = K + XG.at<double>(5, 0);
		count++;
	} while (XG.at<double>(3, 0) >= deltaAng || XG.at<double>(4, 0) >= deltaAng || XG.at<double>(5, 0) >= deltaAng||
		XG.at<double>(0, 0) >= deltaDis|| XG.at<double>(1,0) >= deltaDis|| XG.at<double>(2,0) >= deltaDis);//���ȗʳ����Ф�
	double* myResult = new double[6];
	double mytempResult[6] = {Xs,Ys,Zs,Q,W,K};
	for (int i = 0; i < 6; i++) {
		myResult[i] = mytempResult[i];
	}
	���^�Y���ξ���();
	return Mat(6,1,CV_64FC1,myResult);
}

Mat Resection::��ܞ�ޥȥ�å������ä��ޤ�()
{
	return this->R;
}

Resection::~Resection()
{
	DataMatrix.release();
	L.release(); A.release(); XYZtemp.release();
	AT.release(); ATA.release(); ATL.release(); XG.release();
	R.release();
	AXG.release();
	V.release();
	VT.release();
	VTV.release();
	D.release();
}

void Resection::������ܞ����(Mat& R, double Q, double W, double K)
{
	R.at<double>(0, 0) = cos(Q) * cos(K) - sin(Q) * sin(W) * sin(K);
	R.at<double>(0, 1) = -cos(Q) * sin(K) - sin(Q) * sin(W) * cos(K);
	R.at<double>(0, 2) = -sin(Q) * cos(W);
	R.at<double>(1, 0) = cos(W) * sin(K);
	R.at<double>(1, 1) = cos(W) * cos(K);
	R.at<double>(1, 2) = -sin(W);
	R.at<double>(2, 0) = sin(Q) * cos(K) + cos(Q) * sin(W) * sin(K);
	R.at<double>(2, 1) = -sin(Q) * sin(K) + cos(Q) * sin(W) * cos(K);
	R.at<double>(2, 2) = cos(Q) * cos(W);
}

void Resection::�������ˤ�ä�(Mat& XYZtemp, Mat& DataMatrix, int p, double Xs, double Ys, double Zs)
{
	XYZtemp.at<double>(0, 0) = DataMatrix.at<double>(p, 2) - Xs;
	XYZtemp.at<double>(1, 0) = DataMatrix.at<double>(p, 3) - Ys;
	XYZtemp.at<double>(2, 0) = DataMatrix.at<double>(p, 4) - Zs;
}

void Resection::�ѥ��`���ޥȥ�å�����ȡ��(Mat& A, Mat& R, Mat& DataMatrix, int k, int i, double F, double K, double W, double Z, double X0, double Y0)
{
	A.at<double>(k, 0) = (R.at<double>(0, 0) * F + R.at<double>(0, 2) * (DataMatrix.at<double>(i, 0) - X0)) / Z;
	A.at<double>(k, 1) = (R.at<double>(1, 0) * F + R.at<double>(1, 2) * (DataMatrix.at<double>(i, 0) - X0)) / Z;
	A.at<double>(k, 2) = (R.at<double>(2, 0) * F + R.at<double>(2, 2) * (DataMatrix.at<double>(i, 0) - X0)) / Z;
	A.at<double>(k, 3) = (DataMatrix.at<double>(i, 1) - Y0) * sin(W) - ((DataMatrix.at<double>(i, 0) - X0) * ((DataMatrix.at<double>(i, 0) - X0) * cos(K) - (DataMatrix.at<double>(i, 1) - Y0) * sin(K)) / F + F * cos(K)) * cos(W);
	A.at<double>(k, 4) = -F * sin(K) - (DataMatrix.at<double>(i, 0) - X0) * ((DataMatrix.at<double>(i, 0) - X0) * sin(K) + (DataMatrix.at<double>(i, 1) - Y0) * cos(K)) / F;
	A.at<double>(k, 5) = DataMatrix.at<double>(i, 1) - Y0;
	A.at<double>(k + 1, 0) = (R.at<double>(0, 1) * F + R.at<double>(0, 2) * (DataMatrix.at<double>(i, 1) - Y0)) / Z;
	A.at<double>(k + 1, 1) = (R.at<double>(1, 1) * F + R.at<double>(1, 2) * (DataMatrix.at<double>(i, 1) - Y0)) / Z;
	A.at<double>(k + 1, 2) = (R.at<double>(2, 1) * F + R.at<double>(2, 2) * (DataMatrix.at<double>(i, 1) - Y0)) / Z;
	A.at<double>(k + 1, 3) = -(DataMatrix.at<double>(i, 0) - X0) * sin(W) - ((DataMatrix.at<double>(i, 1) - Y0) * ((DataMatrix.at<double>(i, 0) - X0) * cos(K) - (DataMatrix.at<double>(i, 1) - Y0) * sin(K)) / F - F * sin(K)) * cos(W);
	A.at<double>(k + 1, 4) = -F * cos(K) - (DataMatrix.at<double>(i, 1) - Y0) * ((DataMatrix.at<double>(i, 0) - X0) * sin(K) + (DataMatrix.at<double>(i, 1) - Y0) * cos(K)) / F;
	A.at<double>(k + 1, 5) = -DataMatrix.at<double>(i, 0) - X0;
}

Mat Resection::�خ������Ё\��(Mat& R, Mat& DataMatrix)
{
	Mat temp(3, 1, CV_64FC1);
	temp.at<double>(0, 0) = R.at<double>(0, 0) * XYZtemp.at<double>(0, 0) + R.at<double>(1, 0) * XYZtemp.at<double>(1, 0) + R.at<double>(2, 0) * XYZtemp.at<double>(2, 0);
	temp.at<double>(1, 0) = R.at<double>(0, 1) * XYZtemp.at<double>(0, 0) + R.at<double>(1, 1) * XYZtemp.at<double>(1, 0) + R.at<double>(2, 1) * XYZtemp.at<double>(2, 0);
	temp.at<double>(2, 0) = R.at<double>(0, 2) * XYZtemp.at<double>(0, 0) + R.at<double>(1, 2) * XYZtemp.at<double>(1, 0) + R.at<double>(2, 2) * XYZtemp.at<double>(2, 0);
	return temp;
}

void Resection::���^�Y���ξ���()
{
	double m0;
	AXG = A * XG;
	for (int i = 0; i < 8; i++) {
		V.at<double>(i, 0) = AXG.at<double>(i, 0) - L.at<double>(i, 0);
	}
	transpose(V, VT);
	VTV = VT * V;
	m0 =sqrt( VTV.at<double>(0, 0) / 2);
	
	cout << endl << endl << ATA << endl << endl;
	for (int i = 0; i < 6; i++) {
		D.at<double>(i, i) = m0 * sqrt(ATA.at<double>(i, i));
	}
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			if (i != j)D.at<double>(i, j) = 0.0;
		}
	}
	cout<<"����" << m0<<endl;
	cout << "����X��������" <<endl<<D << endl;

}
