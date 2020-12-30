#include<stdio.h>
/*
spec[]是输入的光谱，e_lo[]能量小，e_hi[]能量大，
retur[]返回值，lspec 光谱的长度，
nn 合并数，e_n retur的长度
gcc -std=c99 c_baseline.c -shared -o c_baseline.so.6

*/

void c_baseline_kernel(double spec[] ,int len_spec,double xx[] ,int lefts[], int rights[] ,int len_xx, int w[], int it);

int min_int3(int i,int j,int k);
double min_int2(double i,double j);

void c_baseline_kernel(double spec[] ,int len_spec,double xx[] ,int lefts[], int rights[] ,int len_xx, int w[], int it)
{
	int i,j,k,h,v;
	double a;
	
	
	for(i=0;i<len_xx;i++)
	{
		a = 0.;
		h = 0;
		for(j=lefts[i];j<=rights[i];j++)
		{
			//printf("spec %f \n",spec[j]);
			a = a + spec[j];
			h++;
		};
		xx[i] = a/(h*1.0);
		//printf("xx[i] %f\n h %d \n a %f\n",xx[i],h,a);
	};

	i = 0;
	while (i < it)
	{
		for(j=1;j<len_xx-1;j++)
		{
			v = min_int3(j,w[i],len_xx-j-1);
			a = 0.;
			for(k = j-v;k<=j+v;k++)
			{
				a = a + xx[k];
			}
			a = a/(2.*v+1.);
			xx[j] = min_int2(xx[j],a);			
		};
		for(j=1;j<len_xx-1;j++)
		{
			h = len_xx-1-j;
			v = min_int3(j,w[i],len_xx-j-1);
			a = 0.;
			for(k = h-v;k<=h+v;k++)
			{
				a = a + xx[k];
			}
			a = a/(2.*v+1.);
			xx[h] = min_int2(xx[h],a);
		};
		i++;
	};
		
}

int min_int3(int i,int j,int k)
{
	int  r;
	r = i;
	if(j<r){r = j;};
	if(k<r){r = k;};
	
	return r;
}


double min_int2(double i,double j)
{
	double  r;
	r = i;
	if(j<r){r = j;};	
	return r;
}




