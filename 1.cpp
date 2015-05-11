#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
int bit_reverse(int N,int B);
int fft_for_2r(double *x_r, double *x_i, double *y_r, double *y_i, int N);
int fft_for_3r(double *x_r, double *x_i, double *y_r, double *y_i, int N);
int fft_for_5r(double *x_r, double *x_i, double *y_r, double *y_i, int N);


int main()
{   double  *x_r,*x_i,*y_r,*y_i; 
	int i,N;
	N=25;
	x_r = (double *) malloc(N*sizeof(double));
	x_i = (double *) malloc(N*sizeof(double));
	y_r = (double *) malloc(N*sizeof(double));
	y_i = (double *) malloc(N*sizeof(double));

	
	
	for(i=0;i<N;i++)
	{
	x_r[i]=i;
	x_i[i]=0;
	
	}

		



	
	
	fft_for_5r(x_r,x_i,y_r,y_i,N);

	
	for(i=0;i<N;i++)
		printf("%f+%fi\n",y_r[i],y_i[i]); 	
		
	return 0;
}





int bit_reverse(int N,int B)//n個數用B進位法的BIT_REVERSE 
{
	int i,j,m;
	i=j=0;
	m=N/B;
	
	while(i<N)
	{
		printf("%d=>%d\n",i,j);
	
		m=N/B;
		//進位問題
		//如果該BIT再加上一位就要進位的話要扣掉改成加下一位
		//  
		while(j>=(B-1)*m & m>0)
		{
			
			j=j-(B-1)*m; 
			m=m/B;      
			
			
			
			
		 } 
		
		i=i+1;	//i多1,j需要在最高為數多1=>+N/B 
		j=j+m;
	}
	
	
	
	return 0;
}

//2^r的FFT   
int fft_for_2r(double *x_r, double *x_i, double *y_r, double *y_i, int N)
{	
	int i;
	
	for(i=0;i<N;i++)
	{
		y_r[i]=x_r[i];
		y_i[i]=x_i[i];
	}
	
	

	
	int j,m;
	double tmp_r,tmp_i;
	
	
	i=j=0;
		while(i<N)
	{
		//swap y[i]<=>y[j]
		printf("%d <-> %d\n",i,j);
		if(i<j)
		{
						
		tmp_r=y_r[i];
		tmp_i=y_i[i];
		y_r[i]=y_r[j];
		y_i[i]=y_i[j];
		y_r[j]=tmp_r;
		y_i[j]=tmp_i;
		
		}
	
		m=N/2;
		  
		while(j>=m & m>0)
		{		
			j=j-m; 
			m=m/2;      	
				
		 } 
		
		i=i+1;	 
		j=j+m;
	}
	

	int n=1;
	double w_r,w_i;
	double wn_r,wn_i;
	w_r=1.0,w_i=0.0;
	
	
	
	while(n<=N/2)
	{   
		
		w_r=1.0,w_i=0.0; 
	    wn_r=cos(-M_PI/n);
	    wn_i=sin(-M_PI/n);
		for(i=0;i<n;i++)   // big group
		{   
			for(j=i;j<N;j=j+2*n)     //small group
			{  
				//yj=yj+w*y(j+n)
				//y(j+n)=yj-w*y(j+n)
			tmp_r=y_r[j+n]*w_r-y_i[j+n]*w_i;
			tmp_i=y_r[j+n]*w_i+y_i[j+n]*w_r;
			y_r[j+n]=y_r[j]-tmp_r;
			y_i[j+n]=y_i[j]-tmp_i;
			y_r[j]=y_r[j]+tmp_r;
			y_i[j]=y_i[j]+tmp_i;
				
				
			}
		tmp_r=w_r,tmp_i=w_i;
		w_r=tmp_r*wn_r-tmp_i*wn_i;
		w_i=tmp_r*wn_i+tmp_i*wn_r;  
		
		
		}
		
		n=2*n;
		
	}
			
		
	
return 0;	
}






int fft_for_3r(double *x_r, double *x_i, double *y_r, double *y_i, int N)
{	
	int i;
	
	for(i=0;i<N;i++)
	{
		y_r[i]=x_r[i];
		y_i[i]=x_i[i];
		
		
		
		
	}
	
	

	
	int j,m;
	double tmp_r,tmp_i;
	
	
	i=j=0;
		while(i<N)
	{
		//swap y[i]<=>y[j]
		
		if(i<j)
		{
						
		tmp_r=y_r[i];
		tmp_i=y_i[i];
		y_r[i]=y_r[j];
		y_i[i]=y_i[j];
		y_r[j]=tmp_r;
		y_i[j]=tmp_i;
		
		
		}
		m=N/3;
		  
		while(j>=2*m & m>0)
		{		
			j=j-2*m; 
			m=m/3;      	
				
		 } 
		
		i=i+1;	 
		j=j+m;
	}


	int n=1;
	double w_r,w_i,wsquare_r,wsquare_i;
	double wn_r,wn_i;
	double w3_r,w3_i,w3square_r,w3square_i;
	double tmp1_r,tmp1_i,tmp2_r,tmp2_i;
	w3_r=cos(-2.0*M_PI/3);
	w3_i=sin(-2.0*M_PI/3);
	w3square_r=w3_r*w3_r-w3_i*w3_i;
	w3square_i=2.0*w3_r*w3_i;
	
	while(n<N)      
	{   w_r=1.0,w_i=0.0;             //W從頭開始 
		wsquare_r=1.0;
		wsquare_i=0.0;	  
	    wn_r=cos(-2.0*M_PI/3/n);       //轉的角度跟N有關 
	    wn_i=sin(-2.0*M_PI/3/n);
	    
	   
	    
		for(i=0;i<n;i++)   // big group
		{   
		
			for(j=i;j<N;j=j+3*n)     //small group
			{   
				//yj=yj+w*y(j+n)+w^2*y(j+2n)
				//y(j+n)=yj+w3*w*y(j+n)+w3^2*w^2*y(j+2n)
				//y(j+2n)=yj+w3^2*w*y(j+n)+w3*w^2*y(j+2n)
				tmp1_r=w_r*y_r[j+n]-w_i*y_i[j+n];    //tmp1=w*y(j+n)
				tmp1_i=w_r*y_i[j+n]+w_i*y_r[j+n];     
				tmp2_r=wsquare_r*y_r[j+2*n]-wsquare_i*y_i[j+2*n];//tmp2=w^2*y(j+2n)
				tmp2_i=wsquare_r*y_i[j+2*n]+wsquare_i*y_r[j+2*n]; 
				
				y_r[j+2*n]=y_r[j]+(w3square_r*tmp1_r-w3square_i*tmp1_i)+(w3_r*tmp2_r-w3_i*tmp2_i);
				y_i[j+2*n]=y_i[j]+(w3square_r*tmp1_i+w3square_i*tmp1_r)+(w3_r*tmp2_i+w3_i*tmp2_r);
				y_r[j+n]=y_r[j]+(w3_r*tmp1_r-w3_i*tmp1_i)+(w3square_r*tmp2_r-w3square_i*tmp2_i);
				y_i[j+n]=y_i[j]+(w3_r*tmp1_i+w3_i*tmp1_r)+(w3square_r*tmp2_i+w3square_i*tmp2_r);
				y_r[j]=y_r[j]+tmp1_r+tmp2_r;
				y_i[j]=y_i[j]+tmp1_i+tmp2_i;
				
				
				
				
				
				
			}
				tmp_r=w_r,tmp_i=w_i;     //w轉 -2.0PI/3*n度 
				w_r=tmp_r*wn_r-tmp_i*wn_i;
				w_i=tmp_r*wn_i+tmp_i*wn_r;
				wsquare_r=w_r*w_r-w_i*w_i;
				wsquare_i=2.0*w_r*w_i;	  
		
		
		}
		
		n=3*n;
		
	}
			

	
return 0;	
}


int fft_for_5r(double *x_r, double *x_i, double *y_r, double *y_i, int N)
{	
	int i;
	
	for(i=0;i<N;i++)
	{
		y_r[i]=x_r[i];
		y_i[i]=x_i[i];
		
		
		
		
	}
	
	

	
	int j,m;
	double tmp_r,tmp_i;
	
	
	i=j=0;
		while(i<N)
	{
		//swap y[i]<=>y[j]
		
		if(i<j)
		{
						
		tmp_r=y_r[i];
		tmp_i=y_i[i];
		y_r[i]=y_r[j];
		y_i[i]=y_i[j];
		y_r[j]=tmp_r;
		y_i[j]=tmp_i;
		
		
		}
		m=N/5;
		  
		while(j>=5*m & m>0)
		{		
			j=j-5*m; 
			m=m/5;      	
				
		 } 
		
		i=i+1;	 
		j=j+m;
	}


	int n=1;
	double w_r[5],w_r[5];
	double wn_r,wn_i;
	double w5_r,w5_i[5];
	double tmp1_r,tmp1_i,tmp2_r,tmp2_i;
	w3_r=cos(-2.0*M_PI/3);
	w3_i=sin(-2.0*M_PI/3);
	w3square_r=w3_r*w3_r-w3_i*w3_i;
	w3square_i=2.0*w3_r*w3_i;
	
	while(n<N)      
	{   w_r=1.0,w_i=0.0;             //W從頭開始 
		wsquare_r=1.0;
		wsquare_i=0.0;	  
	    wn_r=cos(-2.0*M_PI/3/n);       //轉的角度跟N有關 
	    wn_i=sin(-2.0*M_PI/3/n);
	    
	   
	    
		for(i=0;i<n;i++)   // big group
		{   
		
			for(j=i;j<N;j=j+3*n)     //small group
			{   
				//yj=yj+w*y(j+n)+w^2*y(j+2n)
				//y(j+n)=yj+w3*w*y(j+n)+w3^2*w^2*y(j+2n)
				//y(j+2n)=yj+w3^2*w*y(j+n)+w3*w^2*y(j+2n)
				tmp1_r=w_r*y_r[j+n]-w_i*y_i[j+n];    //tmp1=w*y(j+n)
				tmp1_i=w_r*y_i[j+n]+w_i*y_r[j+n];     
				tmp2_r=wsquare_r*y_r[j+2*n]-wsquare_i*y_i[j+2*n];//tmp2=w^2*y(j+2n)
				tmp2_i=wsquare_r*y_i[j+2*n]+wsquare_i*y_r[j+2*n]; 
				
				y_r[j+2*n]=y_r[j]+(w3square_r*tmp1_r-w3square_i*tmp1_i)+(w3_r*tmp2_r-w3_i*tmp2_i);
				y_i[j+2*n]=y_i[j]+(w3square_r*tmp1_i+w3square_i*tmp1_r)+(w3_r*tmp2_i+w3_i*tmp2_r);
				y_r[j+n]=y_r[j]+(w3_r*tmp1_r-w3_i*tmp1_i)+(w3square_r*tmp2_r-w3square_i*tmp2_i);
				y_i[j+n]=y_i[j]+(w3_r*tmp1_i+w3_i*tmp1_r)+(w3square_r*tmp2_i+w3square_i*tmp2_r);
				y_r[j]=y_r[j]+tmp1_r+tmp2_r;
				y_i[j]=y_i[j]+tmp1_i+tmp2_i;
				
				
				
				
				
				
			}
				tmp_r=w_r,tmp_i=w_i;     //w轉 -2.0PI/3*n度 
				w_r=tmp_r*wn_r-tmp_i*wn_i;
				w_i=tmp_r*wn_i+tmp_i*wn_r;
				wsquare_r=w_r*w_r-w_i*w_i;
				wsquare_i=2.0*w_r*w_i;	  
		
		
		}
		
		n=3*n;
		
	}
			

	
return 0;	
}



