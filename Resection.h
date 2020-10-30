#pragma once
#include<iostream>
#include<opencv2/opencv.hpp>
#include<opencv2/core.hpp>
using namespace cv;
using namespace std;
class Resection
{
public:
	
	
private:
	double Xs, Ys, Zs, Q, W, K;//���Ԫ��
	double X0, Y0, F;//�ڷ�λԪ��
	double m;//д��s��
	Mat R;//��ܞ����
	double X, Y, Z;//�������
	Mat DataMatrix;//���ڥǩ`������
	Mat L;//�g��ƽ��ζ����
	Mat A;//�g��ƽ��΂S���
	Mat XYZtemp;//X��Y��Z��һ�r�A����
	Mat AT, ATA, ATL, XG;//�g��ƽ������g����
	Mat AXG, V, VT, VTV, D;
	
public://���饹�Ό��⥤�󥿩`�ե��`��
	Resection(double *data, int row, int col, double mX0, double mY0, double mF, double mM);
	Resection(double **data, int row, int col, double mX0, double mY0, double mF, double mM);
	Mat �Y��(double deltaDis=1e-6,double deltaAng=1e-3);
	Mat ��ܞ�ޥȥ�å������ä��ޤ�();
	~Resection();
private://���饹�Υ��󥿩`�ե��`��
	//�v�B����Ӌ���v��
	void ������ܞ����(Mat& R, double Q, double W, double K);
	void �������ˤ�ä�(Mat& XYZtemp, Mat& DataMatrix, int p, double Xs, double Ys, double Zs);
	void �ѥ��`���ޥȥ�å�����ȡ��(Mat& A, Mat& R, Mat& DataMatrix, int k, int i, double F, double K, double W, double Z, double X0, double Y0);
	Mat �خ������Ё\��(Mat& R, Mat& DataMatrix);
	void ���^�Y���ξ���();
};

