#include "nrutil.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void upper_triangulate(float **A, float *b, int m);
void back_sub(float **A, float *x, float *b, int m);
void move_pivot(float **A, float * b, int j, int m);

/*
int main(int argc, char **argv){
  float **A, *b, *x;
  int n;
  if (argv[1][1]=='d')   //matrix dimensions
     n = 2;
  else n = get_int();
  A = dmatrix(1,n,1,n);
  b = dvector(1,n);
  x = dvector(1,n);
  if (argv[1][1]=='d') {
     A[1][1] = 0.00300; A[1][2]= 59.14;
     A[2][1] = 5.291; A[2][2]= -6.130;
  
     b[1] = 59.17; b[2] = 46.78;
  }  
  else {
     mfill(A,n);
     vfill(b,n);
  }

  mprint(A,n,"A original");
  vprint(b,n,"b original");

  upper_triangulate(A,b,n);
  back_sub(A,x,b,n);
  mprint(A,n,"A");
  vprint(b,n,"b");
  vprint(x,n,"x");
}

void mprint(float **matrix, int m, char *label){
  int i,j;
  printf("%s:\n",label);

  for (i = 1; i <= m; ++i){
    for (j = 1; j <= m; ++j){
      printf("%10.5f ", matrix[i][j]);
    }
    printf("\n");
  }
  printf("\n------------------------\n");
}

void vprint(float *vector, int m, char *label){
  int i;
  printf("%s:\n",label);
  
  for (i = 1; i <= m; ++i){
    printf("%10.2f ", vector[i]);
  }
  
  printf("\n------------------------\n");
}
*/

void upper_triangulate(float **A, float *b, int m){
  int i,j,k;
  
  float scale;
  for (j=1;j<m;++j){           /* loop over columns */
    move_pivot(A, b, j, m);    /* exchange row with largest pivot row */
    for (i=j+1;i<=m;++i){      /* loop over rows beneath pivot */
      if (A[i][j] != 0){       /* if entry not zero already */
	scale = A[i][j]/A[j][j];  /* zero out based on pivot */
	for (k=1;k<=m;++k){
	  A[i][k] = A[i][k] - A[j][k]*scale;
	}
	b[i] = b[i] - b[j]*scale; /* same for b */
      }
    }
  }
}

void back_sub(float **A, float *x, float *b, int m){
  int i,j;
  x[m] = b[m]/A[m][m];

  for (i=m-1;i>=1;--i){
    x[i] = b[i];
    for (j=i+1;j<=m;++j){
      x[i] -= A[i][j]*x[j];
    }
    x[i]/=A[i][i];
  }
}

void move_pivot(float **A, float * b, int j, int m) {
   int maxrow=j,h,k;
   for (h=j+1;h<=m;h++)      /* locate row with the biggest pivot */	
      if (abs(A[h][j]) > abs(A[maxrow][j]))
	maxrow = h;
   if (maxrow!=j) {          /* exchange that row with jth row */
        for (k=1;k<=m;k++) {
	   float tmp = A[j][k];
	   A[j][k]=A[maxrow][k];
	   A[maxrow][k]=tmp;
	}
	float tmp = b[j];   /* exchange values in b too */
	b[j] = b[maxrow];
	b[maxrow] = tmp;
   }
}

/*
int get_int() {
  int input;
  char ch;

  printf("For your n x n matrix, what is n? ");
  while (scanf("%d", &input) != 1 || input <= 0)
  {
     while ((ch=getchar()) != '\n')
        putchar(ch);  //dispose of bad input
     printf("Invalid n value. Please enter a positive integer value. ");
  }
  return input;
} 

void mfill(float **A, int m) {
  printf("Input matrix A (n by n floats): ");
  int i,j;
  for (i=1;i<=m;i++) {
     for (j=1;j<=m;j++)
        scanf("%lf", &A[i][j]);
  }
}

void vfill(float *b, int m) {
  printf("Input vector b (n floats): ");
  int i;
  for (i=1;i<=m;i++) {
    scanf("%lf", &b[i]);
  }
}
*/


