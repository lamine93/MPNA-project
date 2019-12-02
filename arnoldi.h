#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef arnoldi__H
#define arnoldi__H


void get_matrix_columns(double** M, double* c, int j, int row);

double produit_scalaire(double* u, double* v, int n);

double norm2(double* u, int n);

void mat_vec_mul(double** M, double* x, double* b, int n);

void init(double** M, int size1, int size2);

void initial_vector_basis(double* q0, int n);

void arnoldi_projec(double** A, double* q0, double** h, double** Q,  int n, int m);

void matrix_alloc(double** M, int size1, int size2);

void vector_alloc(double* v, int size1);

void free_matrix(double** M, int size1);

void free_vector(double* v);

void print_matrix(double** M, int size1, int size2);


#endif //arnoldi__H
