#include<math.h>
#include<stdio.h>
#include "nrutil.h"
#include<stdlib.h>
#include<time.h>
void printToFile(int p, int q, int r, float ***M);
void upper_triangulate(float **A, float *b, int m);
void back_sub(float **A, float *x, float *b, int m);
void move_pivot(float **A, float * b, int j, int m);

int main(void) {
	int nx = 20;
	int ny = 20;
	int nz = 20;
	double L = 1;
	int nsteps = 10;
	double alpha = .001;

	double dx = (double)L/nx;
	double dy = (double)L/ny;
	double dz = (double)L/nz;
	double dt = .0005;

	double C = alpha*dt/(dx*dx);

	double *xx = dvector(1,nx);
	double *yy = dvector(1,ny);
	double *zz = dvector(1,nz);
	
	//empty the output file
	FILE * fp = fopen("3Doutput","w");
	fclose(fp);

	for (int i=1; i<=nx; i++) {
		xx[i] = (double)i/nx;
	}
	for (int i=1; i<=ny; i++) {
		yy[i] = (double)i/ny;
	}
	for (int i=1; i<=nz; i++) {
		zz[i] = (double)i/nz;
	}
	
	float ***T = f3tensor(1,nx,1,ny,1,nz);
	for (int i=1; i<=nx; i++) {
	   for (int j=1; j<=ny; j++) {
	      for (int k=1; k<=nz; k++) {
		double exponentX = -1*pow((5*xx[i]-2.5),2);
		double exponentY = -1*pow((5*yy[j]-2.5),2);
		double exponentZ = -1*pow((5*zz[k]-2.5),2);
	   	T[i][j][k] = exp(exponentX)*exp(exponentY)*exp(exponentZ);
	      }
	   }
	}
	float ***Tnew = f3tensor(1,nx,1,ny,1,nz);

	for (int i=1; i<=nx; i++) {   //boundary conditions are six planes
	   for (int j=1; j<=ny; j++) {
	   	T[i][j][1]=0; T[i][j][nz]=0;
	   }
	}
	for (int i=1; i<=nx; i++) {  
	   for (int k=1; k<=nz; k++) {
	   	T[i][1][k]=0; T[i][ny][k]=0;
	   }
	}
	for (int j=1; j<=ny; j++) { 
	   for (int k=1; k<=nz; k++) {
	   	T[1][j][k]=0; T[nx][j][k]=0;
	   }
	}
	
	//Crank-Nicolson method setup
   	int dims = nx*ny*nz;			//dimensions of T as vector
   	float *tnew = Tnew[1][1];
   
   	float **A = matrix(1,dims,1,dims);
   	for (int i=1;i<=dims;i++) {
   		for (int j=1;j<=dims;j++) {
	  		if (i==j)
	   			A[i][j] = 1+3*C;
	   		else if ((i-j)==1 || (j-i)==1 || (i-j)==nx || (j-i)==nx
						|| (i-j)==nx*ny || (j-i)==nx*ny)
	   			A[i][j] = -1*C/2;
	   		else 
	   			A[i][j] = 0;
		}
   	}

   	//direct solve Ax=b
   	float ***B = f3tensor(1,nx,1,ny,1,nz);
   	float *b;
	
	time_t bgn, end;	//timers for part 5
	bgn = clock();

	//compute Tnew with FTCS or Crank-Nicolson
	for (int n=1; n<=nsteps; n++) {
		/*
		//FTCS method
		for (int i=2; i<nx; i++) {
		   for (int j=2; j<ny; j++) {
		      for (int k=2; k<nz; k++) {
			   Tnew[i][j][k] = T[i][j][k] 
			   			+ C*(T[i][j-1][k]+T[i][j+1][k]
						+ T[i-1][j][k]+T[i+1][j][k]
						+T[i][j][k-1]+T[i][j][k+1]-6*T[i][j][k]);
		      }
		   }
		}
		*/
		//Crank-Nicolson method
		for (int i=2;i<nx;i++) {	//construct b for Ax=b
		   for (int j=2;j<ny;j++) {
	   		for (int k=2;k<nz;k++) {
	   			B[i][j][k] = (1-3*C)*T[i][j][k]+C/2*(T[i-1][j][k]+T[i+1][j][k]
	    					+T[i][j-1][k]+T[i][j+1][k]
						+T[i][j][k-1]+T[i][j][k+1]);
	   		}
	   	   }
		}
	   	
		for (int i=1; i<=nx; i++) {   //boundary conditions are six planes
   		   for (int j=1; j<=ny; j++) {
   	   		B[i][j][1]=0; B[i][j][nz]=0;
	   	   }
   		}
   		for (int i=1; i<=nx; i++) { 
	   	   for (int k=1; k<=nz; k++) {
	   		B[i][1][k]=0; B[i][ny][k]=0;
	   	   }
   		}
   		for (int j=1; j<=ny; j++) { 
	   	   for (int k=1; k<=nz; k++) {
	   		B[1][j][k]=0; B[nx][j][k]=0;
	   	   }
   		}
	
		b = B[1][1];
		upper_triangulate(A,b,dims);	//Tnew = A\b
		back_sub(A,tnew,b,dims);
		
		//arbitrary time-independent "source term"
		float ***S = f3tensor(1,nx,1,ny,1,nz);
		for (int i=2; i<nx; i++) {
		   for (int j=2; j<ny; j++) {
	     		for (int k=2; k<nz; k++) {
	      	            S[i][j][k]=0;		//can set to 0 for no source term
			    Tnew[i][j][k]+=S[i][j][k]*dt;
     	      		}
	   	   }
		}
		for (int i=1; i<=nx; i++) {   //boundary conditions are six planes
	   	    for (int j=1; j<=ny; j++) {
	   	         Tnew[i][j][1]=T[i][j][nz]; Tnew[i][j][nz]=T[i][j][2];	//can set to 0 for constant
		 								//at present, demonstrating
										//periodic boundary conditions
	   	    } 
	      	}
		for (int i=1; i<=nx; i++) { 
	   	    for (int k=1; k<=nz; k++) {
	   		Tnew[i][1][k]=T[i][ny][k]; Tnew[i][ny][k]=T[i][2][k];
	            }
		}
		for (int j=1; j<=ny; j++) { 
	      	    for (int k=1; k<=nz; k++) {
	   		Tnew[1][j][k]=T[nx][j][k]; Tnew[nx][j][k]=T[2][j][k];
	   	    }
		}
		
		
		if (n%10==0) {
			printToFile(nx,ny,nz,Tnew);
		}
		
		//reset T to Tnew
		for (int i=1; i<=nx; i++) {
		   for (int j=1;j<=ny;j++) {
		      for (int k=1;k<=nz;k++) {
			    T[i][j][k] = Tnew[i][j][k];
		      }
		   }
		}
		if (n==2) {
			end = clock();
			double diff = (end-bgn)*1000/CLOCKS_PER_SEC;
			printf("%lf ms",diff);
		}
	}
}

void print3d(int p, int q, int r, double M[p][q][r]) {
	for (int i=1;i<p;i++) {
	   for (int j=1;j<q;j++) {
	      	printf("%d ", M[i][j][5]);
	   }
	   printf("\n");
	}
	printf("----------------------------------------\n");
}

void printToFile(int p, int q, int r, float ***M) {
	FILE *fp;
	if ((fp = fopen("3Doutput", "a")) == NULL) {
		printf("Can't open file.\n");
		exit(1);
	}

	float * dataptr = &M[1][1][1];
	fwrite(dataptr, sizeof(float), p*q*r, fp);
	
	fclose(fp);
}
