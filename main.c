
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include <lapacke.h>
#include "arnoldi.h"

int main(int argc, char* argv[])
{

  double** A;
  double** h;
  double** Q;

  double* q0;
  double* b;

  int n = 4, m=3;
  double tol =1e-7;

  // allocation

/*
  vector_alloc(q0, n);
  vector_alloc(b, n);


  matrix_alloc(A, n, n);
  matrix_alloc(h, m+1, m);
  matrix_alloc(Q, n, m+1); */


  q0 =  malloc(sizeof(double) * n);
  b  =  malloc(sizeof(double) * n);

  A  =  malloc(sizeof(double) * n);
  for (int i=0; i< n; i++)
  {
    A[i]=malloc(n*sizeof(double));
  }

  h  =  malloc((m+1)*sizeof(double));

  for (int i=0; i< m+1; i++)
  {
    h[i]=malloc(m*sizeof(double));
  }
  Q  = malloc(n*sizeof(double));
  for (int i=0; i< n; i++)
  {
    Q[i]=malloc((m+1)*sizeof(double));
  }

  // call the function
  init(A, n, n);
  initial_vector_basis(q0, n);

  arnoldi_projec(A, q0, h, Q,  n, m);

  printf("\nMatrix H (%d x %d) is:\n", m+1, m);
  print_matrix(h, m+1, m);
  printf("\nMatrix Q (%d x %d) is:\n", n, m+1);
  print_matrix(Q, n, m+1);

 /*  for(int i = 0; i < m+1; i++)
  {
     for(int j=0; j< m; j++)
     {
      printf("H[%d][%d]=%g\n",i,j, h[i][j]);
     }
  }
  for(int i = 0; i < n; i++)
  {
     for(int j=0; j< m+1; j++)
     {
      printf("Q[%d][%d]=%g\n",i,j, Q[i][j]);
     }
  } */

/*for(int i=0; i<n; i++)
{
  free(A[i]);
}
free_matrix(A);
/*
for(int i=0; i< m+1; i++)
{
  free(h[i]);
}
free(h);
for(int i=0; i<n; i++)
{
  free(Q[i]);
}*/
free_matrix(A,n);
free_matrix(h,n);
free_matrix(Q,m+1);
free_vector(q0);
free_vector(b);

return 0;

}
