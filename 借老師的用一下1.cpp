#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
int Group_2(double *x, double *y, int N);
int Group_3(double *x, double *y, int N);
int Group_p(double *x, double *y, int N,int P);

int main()
{
	int i, N, m, n, p;
	double *x, *y, *t;
	
	printf("Input N:");
	scanf("%d", &N);
	x = (double *) malloc(N*sizeof(double));
	y = (double *) malloc(N*sizeof(double));
	
	
	for(i=0;i<N;++i)
	{
		x[i] = i;
	}
	for(i=0;i<N;++i)
		printf("%f\n", x[i]);
		
		printf("\n");
		
	n = N;
	while(n>1)
	{
		m = 0;
		while(m < N)
		{
			if(n%2 == 0)   
			{
				p = 2;
			}
			else if (n%3 == 0)
			{
				p = 3;
			}
			else if (n%5==0)
			{
                 p=5;
            }
            
            
			Group_p(x+m,y+m,n,p);
			m = m + n;
		}
		t = x; x = y; y = t;
		n = n / p;
	}
	for(i=0;i<N;++i)
		printf("%f\n", y[i]);	
		
		system("pause");
}


int Group_p(double *x, double *y, int N, int p)
{
	int i;
	for(i=0;i<N;++i)
	{
		y[(i/p) + (N/p)*(i%p)] = x[i];
	}	
	return 0;
}
