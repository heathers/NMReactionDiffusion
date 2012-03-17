void resid(double ***res, double ***u, double ***rhs, int n)
/*
Returns minus the residual for the model problem. Input quantities are u[1..n][1..n] and
rhs[1..n][1..n], while res[1..n][1..n] is returned.
*/
{
  int i,j,k;
  //double h,h2i;
  //h=1.0/(n-1);
  //h2i=1.0/(h*h);
  double C = .5;
  /* Interior points.*/
  for (k=2;k<n;k++)
    for (j=2;j<n;j++) 
      for (i=2;i<n;i++)
         res[i][j][k] = -C*(u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+u[i][j-1][k]+
	 		u[i][j][k+1]+u[i][j][k-1]-6.0*u[i][j][k])+rhs[i][j][k];
  /* Boundary planes.*/
  for (i=1;i<=n;i++)
    for (j=1;j<=n;j++)
        res[i][j][1]=res[i][j][n]=res[1][i][j]=res[n][i][j]=res[i][1][j]=res[i][n][j]=0.0;
}
