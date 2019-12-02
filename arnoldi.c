#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "arnoldi.h"



void matrix_alloc(double** M, int size1, int size2)
{
  M  =  malloc(sizeof(double) * size1);
  for (int i=0; i< size1; i++)
  {
    M[i]=malloc(size2*sizeof(double));
    if( M[i] == NULL )
     {
          fprintf(stderr,"Allocation impossible");
          exit(EXIT_FAILURE);
     }

  }
}

void vector_alloc(double* v, int size1)
{
  v = malloc(size1*sizeof(double));
  if( v == NULL )
   {
        fprintf(stderr,"Allocation impossible");
        exit(EXIT_FAILURE);
   }

}

void get_matrix_columns(double** M, double* c, int j, int row)
{
  for(int i = 0; i< row; i++)
  {
    c[i] = M[i][j];
  }
}
double produit_scalaire(double* u, double* v, int n)
{
  double sum=0.;
  for (int i=0; i<n; i++){
      sum+=u[i]*v[i];
  }
  return sum;

}
double norm2(double* u, int n)
{
  double norm_2=0.;
  double norm=0.;

  for (int i=0; i<n; i++){

      norm_2 +=u[i]*u[i];
  }
  norm = sqrt(norm_2);
  return norm;

}
void mat_vec_mul(double** M, double* x, double* b, int n)
{
  for(int i=0; i<n; i++)
  {
    b[i]=0;
    for(int j=0; j<n; j++)
    {
      b[i] = b[i] + M[i][j]*x[j];
    }
  }
}

void init(double** M, int size1, int size2)
{
 for(int i = 0; i < size1; i++)
 {
    for(int j=0; j < size2; j++)
    {
     M[i][j] = 1.0 / (i + j + 1);
    }
 }
}
void initial_vector_basis(double* q0, int n)
{
  q0[0]=1;
  for(int i=1; i<n; i++)
  {
    q0[i] = 0;
  }
}



void arnoldi_projec(double** M, double* q0, double** h, double** Q,  int n, int m)
{
/* """Computes a basis of the (n + 1)-Krylov subspace of A: the space
    spanned by {b, Ab, ..., A^n b}. */
/* input */

/* A matrix n * n
q0 : initial vector (length m)
m: dimension of Krylov subspace, must be >=  1*/
/*Returns Q, h
      Q: n x (m + 1) array, the columns are an orthonormal basis of the
        Krylov subspace.
      h: (m + 1) x m array, A on basis Q. It is upper Hessenberg.*/

double* v;
double* q;
double* qj;
double* qk;




// Memory allocation

v  =  malloc(sizeof(double) * n);
q  =  malloc(sizeof(double) * n);
qj  =  malloc(sizeof(double) * n);
qk  =  malloc(sizeof(double) * n);

/*
vector_alloc(v, n);
vector_alloc(q, n);
vector_alloc(qj, n);
vector_alloc(qk, n);*/

// Initialization

  for (int i=0; i< n; i++)
  {
      q[i]= q0[i]/norm2(q0, n); //Normalize the input vector
      Q[i][0] = q[i];  // Use it as the first Krylov vector
  }
  for(int k=0; k<m; k++)
  {
     get_matrix_columns(Q, qk, k, n);
     mat_vec_mul(M, qk, v, n); // Generate a new candidate vector v
     for(int j=0; j<k+1; j++)
     {
       get_matrix_columns(Q, qj, j, n);
       h[j][k]= produit_scalaire(v, qj, n);
       for(int i=0; i<n; i++)
       {
         v[i] = v[i] - h[j][k]*qj[i];
       }// end i
     } // end j
     h[k+1][k]= norm2(v,n);


     for(int i=0; i<n; i++)
     {
       q[i] = v[i]/h[k+1][k];
       Q[i][k+1] = q[i];
      } // end i
  } // end k


free_vector(q);
free_vector(qj);
free_vector(qk);
free_vector(v);

}
void print_matrix(double** M, int size1, int size2)
{
  if(M != NULL)
  {
    for(int i = 0; i < size1; i++)
    {
      for(int j=0; j< size2; j++)
      {
        printf("[%8.4f]", M[i][j]);
      }
      printf ("\n" );
    }
  }
}

void free_matrix(double** M, int size1)
{
  for(int i=0; i<size1; i++)
  {
    free(M[i]);
  }
  free(M);
}

void free_vector(double* v)
{
  free(v);
}
