#pragma once

#include "VxApi.h"
#include "Eigen/Core"
#include "Eigen/Dense"

int FittingInit(int format, void* data);
int FittingExit(void);
int Fitting(int idData, void* echo);
int FittingEO(int idData, void* echo);
int FittingTerm(int idData);
//绘制多段线echo函数
int DrawPolyLineEO(int idData);
//绘制多段线
int DrawPolyLine(svxPoint* ptrPntList, int nPntCount);
//插值型拟合方法
//多项式拟合
void PloynomialFitting(svxPoint* ptrPntList, int nPntCount);
//高斯拟合
void GaussFitting(svxPoint* ptrPntList, int nPntCount);

//Lagrange插值多项式

//Newton插值多项式
// 
//计算多项式拟合表达式结果
void CalculatePolynominalExpress(double dStart, double dEnd);

//计算Gauss拟合曲线上采样点的值
void CalculateGaussExpress(double dStart, double dEnd, svxPoint* ptrPntList);

//逼近型拟合方法

//成员变量
Eigen::VectorXd m_vecdCoe;
const int m_iSampleCount = 100;
svxPoint* m_ptrFittingPntList = nullptr;