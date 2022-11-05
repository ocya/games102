#pragma once

#include "VxApi.h"
#include "Eigen/Core"
#include "Eigen/Dense"

int FittingInit(int format, void* data);
int FittingExit(void);
int Fitting(int idData, void* echo);
int FittingEO(int idData, void* echo);
int FittingTerm(int idData);
//���ƶ����echo����
int DrawPolyLineEO(int idData);
//���ƶ����
int DrawPolyLine(svxPoint* ptrPntList, int nPntCount);
//��ֵ����Ϸ���
//����ʽ���
void PloynomialFitting(svxPoint* ptrPntList, int nPntCount);
//��˹���
void GaussFitting(svxPoint* ptrPntList, int nPntCount);

//Lagrange��ֵ����ʽ

//Newton��ֵ����ʽ
// 
//�������ʽ��ϱ��ʽ���
void CalculatePolynominalExpress(double dStart, double dEnd);

//����Gauss��������ϲ������ֵ
void CalculateGaussExpress(double dStart, double dEnd, svxPoint* ptrPntList);

//�ƽ�����Ϸ���

//��Ա����
Eigen::VectorXd m_vecdCoe;
const int m_iSampleCount = 100;
svxPoint* m_ptrFittingPntList = nullptr;