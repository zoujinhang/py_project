#include<stdio.h>
#include<math.h>
/*
这里需要用c99进行编译。因为有自定义长度的数组
编译命令：

gcc -std=c99 WT_weight.c -shared -o WT_weight.so.6
gcc -std=c99 WT_weight.c -shared -o WT_weight.so
*/


void work_WT(double wt[], double sigma,double dt[],double t[],int len_t,int len_dt);
//float sum(float x[],int len_x);
//void weight(float x[],float rr[],int len_x, float x0, float sigma);

void work_WT(double wt[], double sigma,double dt[],double t[],int len_t,int len_dt)
{
	int i,j;
	double w_sum,wj,w_dt,exp_in;

	for(i=0;i<len_t;i++)
	{
		w_sum = 0;
		w_dt = 0;
		for(j=0;j<len_dt;j++)
		{
		    if(j<i){
		        exp_in = -1*(t[j]-t[i])*(t[j]-t[i])*sigma;
		    }
		    else{
		        exp_in = -1*(t[j+1]-t[i])*(t[j+1]-t[i])*sigma;
		    }
		    wj = exp(exp_in);
		    w_sum += wj;
		    w_dt += dt[j]*wj;
		}
		wt[i] = w_dt/w_sum;
	}

}
/*
float sum(float x[],int len_x)
{
	float sum_=0.0;
	int i;
	for(i=0;i<len_x;i++)
	{
		sum_ += x[i];		
	}
	return sum_;
}


void weight(float x[],float rr[],int len_x, float x0, float sigma)
{
	int i;
	//float rr[len_x];
	printf("e = %f\n",exp(1));
	for(i=0;i<len_x;i++)
	{
		rr[i] = exp(-(x[i]-x0)*(x[i]-x0)*sigma);
		printf("rr[i] = %f\n",rr[i]);
	}
	//return rr;	

}
*/

