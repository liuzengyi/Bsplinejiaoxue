// BsplineCurve.h: interface for the BsplineCurve class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BSPLINECURVE_H__6070F1A8_773E_4016_AC04_CD095008D5B5__INCLUDED_)
#define AFX_BSPLINECURVE_H__6070F1A8_773E_4016_AC04_CD095008D5B5__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include "Vec.h"
#include <vector>
using namespace std;

class BsplineCurve  
{
public:
	BsplineCurve();
	virtual ~BsplineCurve();
	BsplineCurve(vector<vec> ctlps);

	vector<vec2>            pts;
	void                    GetPts();

	void		            SetCP(vector<vec>  p);
	vector<vec>	            GetCP();
	vector<double>			GetU();
	void					SetU(vector<double> knots);
	int						Getk();
	void					Setk(int d);
	//
	void					chordLengthParam(vector<vec> pt, vector<double> &ub);			//�����ҳ�������
	void					Hartley_Judd_Param(vector<vec> pt);											//������-�ֵ²�����
	int						findSpan(double u) ;
	vec					    Dbordv( double u);	
	vector<vec>             globalInterp(vector<vec> datapoints, int d);
protected:
  vector<double>            U;						//�ڵ�ʸ��
  vector<vec>	            CP;						//���ƶ���
  int			            k;						//���ߴ���
};

#endif // !defined(AFX_BSPLINECURVE_H__6070F1A8_773E_4016_AC04_CD095008D5B5__INCLUDED_)
