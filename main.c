#include "mg.h"
#include "nrutil.h"
#include <stdio.h>
#include<time.h>
int main(int argc, char **argv){
  int i,j,k;
  FILE *outfile;
  double ***f;
  int n=33;
  int ncycle=2;
  f = d3tensor(1,n,1,n,1,n);
  f[16][16][16]=1.0;
  //  for (i=2;i<n;++i)
  //  for (j=2;j<n;++j)
  //    f[i][j] = 2.0;
  time_t bgn, end;
  bgn = clock();
  mglin(f,n,ncycle);
  end = clock();
  int diff = (end - bgn)*1000/CLOCKS_PER_SEC;
  printf("Time in ms for n=%d is %d.\n", n, diff);
  outfile = fopen("soln.dat", "w");
  fwrite(&f[1][1][1],sizeof(double),n*n*n,outfile);
  fclose(outfile);
}
