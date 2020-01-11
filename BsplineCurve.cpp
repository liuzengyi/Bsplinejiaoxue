// BsplineCurve.cpp: implementation of the BsplineCurve class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Bsplinejiaoxue.h"
#include "BsplineCurve.h"
#include "FYMatrixMN.h"
typedef FYNAMESPACE::FYMatrixMN<double> FyMatrix;
#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BsplineCurve::BsplineCurve()
{
	U.clear();
	CP.clear();
	pts.clear();
}

BsplineCurve::~BsplineCurve()
{
	U.clear();
	CP.clear();
	pts.clear();
}
BsplineCurve::BsplineCurve(vector<vec> ctlps)
{
	pts.clear();
	CP=ctlps;
	int n=CP.size();
	if (n<2)
	{
		return;
	}
	else if (n<4)
	{
		k=n-1;
	}
	else
	{
		k=3;
	}

	Hartley_Judd_Param(ctlps);
}
void BsplineCurve::Hartley_Judd_Param(vector<vec> pt)
{
	U.clear();
	int i;
	int n=pt.size();
	if (n<2)
	{
		return;
	}
	U.resize(n+k+1,0);

	vector<double> ll(n,0);
	vector<double> delta(n+1,0);
	for (int i=1; i<n;i++)
	{
		vec li=pt[i]-pt[i-1];
		ll[i]=len(li);
	}
	double fenmu=0.0;
	for ( i=k+1; i<n+1; i++)
	{
		for (int j=i-k; j<i; j++)
		{
			fenmu+=ll[j];
		}
	}
	for ( i=k+1; i<n+1; i++)
	{
		double fenzi=0.0;
		for (int j=i-k; j<i; j++)
		{
			fenzi+=ll[j];
		}
		delta[i]=fenzi/fenmu;
	}

	for ( i=k+1; i<n; i++)
	{
		for (int j=k+1; j<=i; j++)
		{
			U[i]+=delta[j];
		}
	}
	for ( i=n; i<n+k+1; i++)
	{
		U[i]=1.0;
	}
}

int BsplineCurve::findSpan(double u)
{
	int ctlpsize=CP.size();
	if(u>=U[ctlpsize]) 
		return ctlpsize-1 ;
	if(u<=U[k])
		return k ;

	int low  = 0 ;
	int high = ctlpsize+1 ; 
	int mid = (low+high)/2 ;

	while(u<U[mid] || u>= U[mid+1])
	{
		if(u<U[mid])
			high = mid ;
		else
			low = mid ;
		mid = (low+high)/2 ;
	}
	return mid ;
}

vec BsplineCurve::Dbordv( double u)
{
	int span=findSpan(u);
	vector<point> dt;
	dt.resize(CP.size());
	for (int i=0; i<CP.size(); i++)
	{
		dt[i]=CP[i];
	}
	for (int ll=1; ll<k+1; ll++)
	{
		for (int j=span-k; j<span-ll+1; j++)
		{
			//generate_Geodesic(dt[j],dt[j+1]);
			float alpha,denom;
			denom=U[j+k+1]-U[j+ll];
			if (fabs(denom)<1e-5)
			{
				alpha=0.0;
			}
			else
			{
				alpha=(u-U[j+ll])/denom;
			}
			dt[j]=(1-alpha)*dt[j]+alpha*dt[j+1];
		}
	}
	return dt[span-k];
}

void BsplineCurve::chordLengthParam(vector<vec> pt, vector<double> &ub)
{
	//应首先定义ub的大小，否则出错
	int n=pt.size();
	double s=0;
	ub[0]=0;
	for (int i=0; i<n-1;i++)
	{
		vec dl=pt[i+1]-pt[i];
		s=s+len(dl);
	}
	if (s>0)
	{
		for (int i=1; i<(ub.size())-1; i++)
		{
			ub[i] = ub[i-1] + len(pt[i]-pt[i-1])/s ;
		}
		ub[ub.size()-1]=1.0;
	}
	else
	{
		for (int i=1; i<(ub.size()-1); i++)
		{
			ub[i] = (double)i/(double)(ub.size()-1);
		}
		ub[ub.size()-1]=1.0;
	}
}

void BsplineCurve::GetPts()
{
	pts.clear();
	double t=0;
	int i,n;
	n=101;
	for (i=0;i<n;i++)
	{
		t=(double)i/(double)(n-1);
		vec p=Dbordv(t);
		vec2 pt(p[0],p[1]);
		pts.push_back(pt);
	}
}

void BsplineCurve::SetCP(vector<vec> p)
{
	CP=p;
}
//
vector<vec> BsplineCurve::globalInterp(vector<vec> datapoints, int d)
{
	int qsize=datapoints.size();
	vector<double> ub(qsize);
	int n=qsize+d-1;

	k = d ;
	U.resize(n+k+1);

	chordLengthParam(datapoints,ub);

	vector<vec> P(n);
	int bc=1;
	int ctlpsize=n;
	int kotsize=n+k+1;
	int i;
	for ( i=0;i<k;i++)
	{
		U[i]=0;
	}
	for ( i=kotsize-d;i<kotsize;i++)
	{
		U[i]=1.0;
	}
	for ( i=d;i<kotsize-d;i++)
	{
		U[i]=ub[i-d];
	}
	vector<double> dlt(kotsize);
	for ( i=0;i<kotsize-1;i++)
	{
		dlt[i]=U[i+1]-U[i];
	}
	vector<double> a(ctlpsize);
	vector<double> b(ctlpsize);
	vector<double> c(ctlpsize);
	vector<vec> e(ctlpsize-2);
	for ( i=2;i<ctlpsize-2;i++)
	{
		a[i]=(dlt[i+2]*dlt[i+2])/(dlt[i]+dlt[i+1]+dlt[i+2]);
		b[i]=(dlt[i+2]*(dlt[i]+dlt[i+1]))/(dlt[i]+dlt[i+1]+dlt[i+2])+(dlt[i+1]*(dlt[i+2]+dlt[i+3]))/(dlt[i+1]+dlt[i+2]+dlt[i+3]);
		c[i]=(dlt[i+1]*dlt[i+1])/(dlt[i+1]+dlt[i+2]+dlt[i+3]);
		e[i-1]=datapoints[i-1]*(float)(dlt[i+1]+dlt[i+2]);//e[]的下标往前移动了一位，便于数组操作
	}
	//
	FyMatrix A(ctlpsize-2,ctlpsize-2);
	A.AllZero();

	if (bc)//抛物线边界条件
	{
		A(0,0)=1.0-dlt[3]*dlt[4]/((dlt[3]+dlt[4])*(dlt[3]+dlt[4]));
		e[0]=(datapoints[0]+datapoints[1]*2.0f)/3.0f;
		A(ctlpsize-3,ctlpsize-3)=dlt[ctlpsize-2]*dlt[ctlpsize-1]/((dlt[ctlpsize-2]+dlt[ctlpsize-1])*(dlt[ctlpsize-2]+dlt[ctlpsize-1]))-1.0;
		e[ctlpsize-3]=(datapoints[ctlpsize-d]+datapoints[ctlpsize-d-1]*2.0f)*(-1.0f/3.0f);
	}
	else//自由端点条件
	{
		A(0,0)=2.0-dlt[3]*dlt[4]/((dlt[3]+dlt[4])*(dlt[3]+dlt[4]));
		e[0]=datapoints[0]+datapoints[1];
		A(ctlpsize-3,ctlpsize-3)=dlt[ctlpsize-2]*dlt[ctlpsize-1]/((dlt[ctlpsize-2]+dlt[ctlpsize-1])*(dlt[ctlpsize-2]+dlt[ctlpsize-1]))-2.0;
		e[ctlpsize-3]=datapoints[ctlpsize-d]*(-1.0f)-datapoints[ctlpsize-d-1];
	}
	A(0,1)=(dlt[3]/(dlt[3]+dlt[4]))*(dlt[4]/(dlt[3]+dlt[4])-dlt[3]/(dlt[3]+dlt[4]+dlt[5]));
	A(0,2)=dlt[3]*dlt[3]/((dlt[3]+dlt[4])*(dlt[3]+dlt[4]+dlt[5]));

	A(ctlpsize-3,ctlpsize-5)=-dlt[ctlpsize-1]*dlt[ctlpsize-1]/((dlt[ctlpsize-2]+dlt[ctlpsize-1])*(dlt[ctlpsize-2]+dlt[ctlpsize-1]+dlt[ctlpsize-3]));
	A(ctlpsize-3,ctlpsize-4)=dlt[ctlpsize-1]/(dlt[ctlpsize-2]+dlt[ctlpsize-1])*(dlt[ctlpsize-1]/(dlt[ctlpsize-2]+dlt[ctlpsize-1]+dlt[ctlpsize-3])-dlt[ctlpsize-2]/(dlt[ctlpsize-2]+dlt[ctlpsize-1]));

	FyMatrix xx(ctlpsize-2,3);
	FyMatrix ee(ctlpsize-2,3);
	for ( i=0; i<ctlpsize-2; i++)
	{
		for (int j=0; j<3; j++)
		{
			ee(i,j)=e[i][j];
		}
	}

	for( i=1;i<ctlpsize-3;i++)
		for (int j=0;j<ctlpsize-2;j++)
		{
			if(i==j)
				A(i,j)=b[i+1];
			else if(i==j+1)
				A(i,j)=a[i+1];
			else if(i==j-1)
				A(i,j)=c[i+1];
		}

		A=A.Inverse();
		xx=A*ee;

		P[0]=datapoints[0];
		P[ctlpsize-1]=datapoints[ctlpsize-k];

		for( i=1;i<ctlpsize-1;i++)
		{
			P[i][0] = xx[i-1][0];
			P[i][1] = xx[i-1][1];
			P[i][2] = xx[i-1][2];
		}

	return P;
}
