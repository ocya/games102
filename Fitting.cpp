#include "Fitting.h"
#include <cmath>
#define	FITTING_UI "Fitting"

int FittingInit(int format, void* data)
{
	cvxCmdFunc("Fitting", (void*)Fitting, VX_CODE_GENERAL);

	cvxCmdFunc("FittingEO", (void*)FittingEO, VX_CODE_GENERAL);
	cvxCmdFunc("FittingTerm", (void*)FittingTerm, VX_CODE_GENERAL);
	return 0;
}


int FittingExit(void)
{
	cvxCmdFuncUnload("Fitting");
	return 0;
}

//int cmp(const void* pntFirst, const void* pntSecond) 
//{
//	return (*((svxPoint*)pntSecond)).x - (*((svxPoint*)pntFirst)).x;
//}

int Fitting(int idData, void *echo)
{
	//ʹģ�Ϳռ���ӽǱ�Ϊxyƽ��
	cvxCmdBuffer("^viewTop", 0);
	//acquire user input points list
	int nRet = 0;
	int nPntCount = 0;
	svxPoint* ptrPntList = nullptr;
	nRet = cvxDataGetPnts(idData, 2, &nPntCount, &ptrPntList);
	if (1 == nRet)
	{
		cvxMemFree((void**)&ptrPntList);
		return nRet;
	}

	//qsort(ptrPntList, nPntCount, sizeof(svxPoint), cmp);
	//���û�����ĵ㼯�Ͻ�������	
	//int nListLength = sizeof(ptrPntList) / sizeof(svxPoint);
	std::sort(ptrPntList, ptrPntList + nPntCount, [](svxPoint pntFirst, svxPoint pntSecond) {
		return pntFirst.x < pntSecond.x;
		});

	//�����û�ѡ���
	nRet = cvxPartPnts(nPntCount, ptrPntList);
	if (1 == nRet)
	{
		cvxMemFree((void**)&ptrPntList);
		return nRet;
	}
	//���ƶ���ʽ�������
	//nRet = DrawPolyLine(ptrPntList, nPntCount);
	PloynomialFitting(ptrPntList, nPntCount);
	if (nullptr != m_ptrFittingPntList)
	{
		nRet = DrawPolyLine(m_ptrFittingPntList, m_iSampleCount);
	}

	if (1 == nRet)
	{
		cvxMemFree((void**)&ptrPntList); //�����û��������ڴ�
		return nRet;
	}

	//����Gauss�������
	GaussFitting(ptrPntList, nPntCount);
	if (nullptr != m_ptrFittingPntList)
	{
		nRet = DrawPolyLine(m_ptrFittingPntList, m_iSampleCount);
	}

	//�����ڴ�
	cvxMsgDisp("Fitting successed");
	cvxMemFree((void**)&ptrPntList); //�����û��������ڴ�
	//cvxMemFree((void**)&m_ptrFittingPntList); //������Ϸ����ϲ�������ڴ�


	return nRet;
}


int FittingEO(int idData, void* echo)
{
	cvxEchoStart();
	DrawPolyLineEO(idData);
	cvxEchoEnd();
	return 0;
}


int FittingTerm(int idData)
{
	return 0;
}


//���ƶ����echo����
int DrawPolyLineEO(int idData)
{
	
	int nPolylineIdData = 0;	
	svxData tempData;
	cvxDataInit("CdLnPoly", &nPolylineIdData);
	cvxDataZero(&tempData);
	tempData.isPoint = 1;
	tempData.PntType = VX_PNT3_ABS;
	
	int nRet = 0;
	int nPntCount = 0;
	svxPoint* ptrPntList = nullptr;

	//update acquired points list
	cvxDataGetPnts(idData, 2, &nPntCount, &ptrPntList);

	for (int i = 0; i < nPntCount; i++)
	{
		tempData.Pnt = ptrPntList[i];
		cvxDataSet(nPolylineIdData, 1, &tempData);
	}
	cvxCmdExec(nPolylineIdData);


	return 0;
}

//���ƶ����
int DrawPolyLine(svxPoint* ptrPntList, int nPntCount)
{
	//��ʼ������߻��Ƶ�ui
	int nRet = 0;
	int nPolylineIdData = 0;
	svxData tempData;
	nRet = cvxDataInit("CdLnPoly", &nPolylineIdData);
	if (1 == nRet)
	{
		return nRet;
	}

	//���ݸ�ֵ
	for (int i = 0; i < nPntCount; i++)
	{
		tempData.Pnt = ptrPntList[i];
		//����PointList����������cvxDataSet����Ϊ��List�����Point
		nRet = cvxDataSet(nPolylineIdData, 1, &tempData);
		if (1 == nRet)
		{
			return nRet;
		}
	}

	//ִ��ָ��
	nRet = cvxCmdExec(nPolylineIdData);


	return nRet;
}

//����ʽ���
void PloynomialFitting(svxPoint* ptrPntList, int nPntCount)
{
	//���ݵ㹹��x�����y����
	Eigen::MatrixXd matdX(nPntCount, nPntCount);
	Eigen::VectorXd vecdY(nPntCount);
	for (int i = 0; i < nPntCount; i++)
	{
		vecdY(i) = ptrPntList[i].y;
		matdX(i, 0) = 1;
		for (int j = 1; j < nPntCount; j++)
		{
			matdX(i, j) = matdX(i, j - 1) * ptrPntList[i].x;
		}
	}
	//Eigen::VectorXd vecdCoe = matdX.inverse() * vecdY;
	m_vecdCoe.resize(nPntCount);// n�� 1��
	m_vecdCoe = matdX.inverse() * vecdY; //inverse()ת�þ���
	
	
	//����Ϸ��̽��в�������û�ͼʱ����Ҫ�����ݵ�
	double dStart = ptrPntList[0].x;
	double dEnd = ptrPntList[nPntCount - 1].x;
	CalculatePolynominalExpress(dStart, dEnd);
}


//��˹���
void GaussFitting(svxPoint* ptrPntList, int nPntCount)
{
	Eigen::MatrixXd matdX(nPntCount+1, nPntCount + 1);// n+1�� n+1��
	Eigen::VectorXd vecdY(nPntCount+1);// 1�� n+1��
	//���ݵ㹹��x�����y����
	//����δ֪������Ϊn+1�����̸���Ϊn������һ��Լ���㣬�����������ͨ��X��ֵ�㣬Y��ֵ�㡣
	double dConstraintX = 0;
	double dConstraintY = 0;
	for (int i = 0; i < nPntCount; i++)
	{
		vecdY(i) = ptrPntList[i].y;
		dConstraintY += ptrPntList[i].y / nPntCount;//��ֹ�������������
		dConstraintX += ptrPntList[i].x / nPntCount;
		matdX(i, 0) = 1;
		double dX = ptrPntList[i].x;
		double dXi = 0;
		for (int j = 1; j < nPntCount + 1; j++)
		{
			dXi = ptrPntList[j - 1].x;
			matdX(i, j) = exp((dX - dXi) * (dX - dXi) * (-0.5));
		}
	}
	vecdY(nPntCount) = dConstraintY;
	matdX(nPntCount, 0) = 1;
	for (int i = 1; i < nPntCount + 1; i++)
	{
		matdX(nPntCount, i) = exp((dConstraintX - ptrPntList[i - 1].x) * (dConstraintX - ptrPntList[i - 1].x) * (-0.5));
	}

	m_vecdCoe.resize(nPntCount+1);//ϵ������1�� n+1��
	m_vecdCoe = matdX.inverse() * vecdY;

	//����Ϸ��̽��в�������û�ͼʱ����Ҫ�����ݵ�
	double dStart = ptrPntList[0].x;
	double dEnd = ptrPntList[nPntCount - 1].x;
	CalculateGaussExpress(dStart, dEnd, ptrPntList);

}

//�������ʽ��������ϲ������ֵ
void CalculatePolynominalExpress(double dStart, double dEnd)
{
	
	double x = 0;
	double dStepLength = (dEnd - dStart) / (m_iSampleCount - 1);
	int n = m_vecdCoe.size();//����ʽ����ָ��
	m_ptrFittingPntList = new svxPoint[m_iSampleCount];
	
	for (int i = 0; i < m_iSampleCount; i++)
	{
		x = dStart + i * dStepLength;
		m_ptrFittingPntList[i].x = x;
		double dXExpTemp = 1;
		double y = 0;
		for (int j = 0; j < n; j++)
		{
			if(0 != j)
			{
				dXExpTemp *= x;//��x����ָ��
			}		
			y += m_vecdCoe[j] * dXExpTemp;
		}
		m_ptrFittingPntList[i].y = y;
		m_ptrFittingPntList[i].z = 0;
	}
}

//����Gauss��������ϲ������ֵ
void CalculateGaussExpress(double dStart, double dEnd, svxPoint* ptrPntList)
{
	double dX = 0;
	double dStepLength = (dEnd - dStart) / (m_iSampleCount - 1);
	int n = m_vecdCoe.size();//����ʽ����ָ��
	m_ptrFittingPntList = new svxPoint[m_iSampleCount];

	for (int i = 0; i < m_iSampleCount; i++)
	{
		dX = dStart + i * dStepLength;
		m_ptrFittingPntList[i].x = dX;

		double dGXi = 1;
		double y = 0;
		double dXi = 0;
		for (int j = 0; j < n; j++)
		{
			if (0 != j)
			{
				dXi = ptrPntList[j - 1].x;
				dGXi = exp((dX - dXi) * (dX - dXi) * (-0.5));
			}
			y += m_vecdCoe[j] * dGXi;
		}
		m_ptrFittingPntList[i].y = y;
		m_ptrFittingPntList[i].z = 0;
	}
}

//void memoryClear()
//{
//	cvxMemFree((void**)&m_ptrFittingPntList);
//}
