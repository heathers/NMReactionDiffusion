void relax(double ***u, double ***rhs, int n)
/*
  Red-black Gauss-Seidel relaxation for model problem. Updates the current value of the solution
  u[1..n][1..n], using the right-hand side function rhs[1..n][1..n].
*/
{
  int i,j,k;
  //double h,h2;
  //h=1.0/(n-1);
  //h2=h*h;
  double C = 10;
  /* Red and black sweeps.*/
  /* jsw and isw toggle between 1 and 2 and
     determine starting row in each column
     for given pass 
  */
  for (i=2;i<n;i++)
    for (j=2;j<n;j++)
       for (k=2;k<n;k++)	/*Gauss-Seidel formula.*/ 
	    u[i][j][k]=C/(6*C+1)*(u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]
		        +u[i][j-1][k]+u[i][j][k+1]+u[i][j][k-1])
		        -(double)1/(6*C+1)*rhs[i][j][k];
}
