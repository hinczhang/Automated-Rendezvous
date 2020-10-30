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
	double Xs, Ys, Zs, Q, W, K;//外角元素
	double X0, Y0, F;//内方位元素
	double m;//写真縮尺
	Mat R;//回転行列
	double X, Y, Z;//物点座標
	Mat DataMatrix;//初期データ行列
	Mat L;//間接平差の定数項
	Mat A;//間接平差の係数陣
	Mat XYZtemp;//X，Y，Zの一時預かり
	Mat AT, ATA, ATL, XG;//間接平差の中間行列
	Mat AXG, V, VT, VTV, D;
	
public://クラスの対外インターフェース
	Resection(double *data, int row, int col, double mX0, double mY0, double mF, double mM);
	Resection(double **data, int row, int col, double mX0, double mY0, double mF, double mM);
	Mat 結果(double deltaDis=1e-6,double deltaAng=1e-3);
	Mat 回転マトリックスが得られます();
	~Resection();
private://クラスのインターフェース
	//関連する計算関数
	void 良い回転行列(Mat& R, double Q, double W, double K);
	void 対物座標を得る(Mat& XYZtemp, Mat& DataMatrix, int p, double Xs, double Ys, double Zs);
	void パラメータマトリックスを取得(Mat& A, Mat& R, Mat& DataMatrix, int k, int i, double F, double K, double W, double Z, double X0, double Y0);
	Mat 特異な行列乗算(Mat& R, Mat& DataMatrix);
	void 検証結果の精度();
};

