#include "nrutil.h"
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>

void printToFile(int p, int q, int r, double ***M);
void copy(int p, int q, int r, double ***source, double ***dest);

int main()
{
    clock_t bgn, end;
    bgn = clock();
    
    double C = 10;
    int n = 80;
    double *xx = dvector(1,n);
    for (int i = 1; i < n; i++) {
    	xx[i] = -2 + 4*(double)i/n;
    }
    int MAX_ITER = 1000;
    int nsteps = 100;
    double ***xold = d3tensor(1,n,1,n,1,n);
    for (int i = 1; i < n; i++) {
    	for (int j = 1; j < n; j++) {
	   for (int k = 1; k < n; k++) {
	   	xold[i][j][k]=exp(-2*pow(xx[i],2))*exp(-2*pow(xx[j],2))*exp(-2*pow(xx[k],2));
	   }
	}
    }
    double ***xnew = d3tensor(1,n,1,n,1,n);
    double ***x = d3tensor(1,n,1,n,1,n);
    
    for (int i=1; i<=n; i++) {   //boundary conditions are six planes
      for (int j=1; j<=n; j++) {
      	 for (int k=1; k<=n; k++) {
   		xnew[i][j][1]=0; xnew[i][j][n]=0;
   		xnew[i][1][k]=0; xnew[i][n][k]=0;
   		xnew[1][j][k]=0; xnew[n][j][k]=0;
      	 }
      }
    }
    
    for (int step=1; step<nsteps; step++) {

    	copy(n,n,n,xold,x);

	for (int m=1; m<MAX_ITER; m++) {
    	   
	   for (int i = 2; i < n; i++) {
    	      for (int j = 2; j < n; j++) {
	         for (int k = 2; k < n; k++) {
	   		xnew[i][j][k]=C/(6*C+1)*(x[i-1][j][k]+x[i+1][j][k]+x[i][j-1][k]
	   				+x[i][j+1][k]+x[i][j][k-1]+x[i][j][k+1])+
					(double)1/(6*C+1)*xold[i][j][k];
	         }
	      }
	   }
 
	   //ascertain if error is less than threshold in order to break
	   double sum=0;
	   for (int i=1; i<=n; i++) { 
	      for (int j=1; j<=n; j++) {
	      	for (int k=1; k<=n; k++) {
	  		sum+=fabs(x[i][j][k]-xnew[i][j][k]);
		}
	      }
	   }
	   double mean = sum/(double)(n*n*n);
	   //printf("Step: %d, Iter: %d, Sum: %lf, Mean: %lf\n",step,m,sum,mean);
	   if (mean < 0.000000001) break;

	   copy(n,n,n,xnew,x);
	}
	
	/*if (step%20==0) {
		printToFile(n,n,n,xnew);
	}*/

	copy(n,n,n,x,xold);
   }
   end = clock();
   printf("%d\n",end);
   double diff = (end - bgn) * 1000 / CLOCKS_PER_SEC;
   printf("Time elapsed for %d cubic size is %lf ms.\n", n, diff);
}

void printToFile(int p, int q, int r, double ***M) {
	FILE *fp;
	if ((fp = fopen("output", "a")) == NULL) {
		printf("Can't open file.\n");
		exit(1);
	}

	double * dataptr = &M[1][1][1];
	fwrite(dataptr, sizeof(double), p*q*r, fp);
	
	fclose(fp);
}

void copy(int p, int q, int r, double ***source, double ***dest) {
	for (int i=1;i<=p;i++) {
	   for (int j=1;j<=q;j++) {
	      for (int k=1;k<=r;k++) {
	      	  dest[i][j][k] = source[i][j][k];
	      }
	   }
	}
}
