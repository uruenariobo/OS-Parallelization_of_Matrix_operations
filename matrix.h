#ifndef MATRIX_H
#define MATRIX_H

#include <stdint.h>

typedef struct Vector Vector;
struct Vector
{
    int size;
    double *elements;
};

typedef struct Matrix Matrix;
struct Matrix
{
    int rows;
    int cols;
    double **elements;
};

Vector *create_vector(int size);
Vector *create_vector_from_file(const char *file_path, int size);
Matrix *create_matrix(int rows, int cols);
Matrix *create_matrix_from_file(const char *file_path, int rows, int cols);

void init_matrix_rand(Matrix *M);
void init_vector_rand(Vector *V);
void copy_matrix(Matrix *dst, const Matrix *src);
void copy_vector(Vector *dst, const Vector *src);

void free_vector(Vector *v);
void free_matrix(Matrix *M);

Vector *add_vector(const Vector *a, const Vector *b);
Vector *dot_vector_matrix(const Vector *v, const Matrix *M);
Matrix *dot_matrix(const Matrix *M, const Matrix *N);


Vector* matrix_col_sum(const Matrix* M);

Vector *matrix_col_sum(const Matrix *M);


void scalar_vector(Vector *V, double k);

void print_vector(const Vector *v);
void print_matrix(const Matrix *M);

//1. Calcular la media de cada columna de una matriz
Vector *matrix_col_mean(const Matrix *M);
Vector *matrix_col_mean_parallel(const Matrix *M);

//2. Calcular la varianza de cada columna de una matriz
Vector *matrix_col_vrz(const Matrix *M);
Vector *matrix_col_vrz_parallel(const Matrix *M);

//3. Calcular la desviacion estandar de cada columna de una matriz
Vector *matrix_col_std(const Matrix *M);
Vector *matrix_col_std_parallel(const Matrix *M);

//5. Calcular la suma de dos Matrices
Matrix *add_matrix(const Matrix *M, const Matrix *N);
void *add_rows_thread(void *arg);
Matrix *add_matrix_parallel(const Matrix *M, const Matrix *N, int n_threads);

//7. Calcular la multiplicación de una matriz por un escalar
void scalar_matrix(Matrix *M, double k);
void *scalar_matrix_thread(void *arg);
void scalar_matrix_parallel(Matrix *M, double k);

// 4 Calcular el valor mínimo y el valor máximo de cada columna de una matriz
Vector* matrix_col_max(const Matrix* M);
void* max_cols_thread(void* arg);

Vector* matrix_col_min(const Matrix* M);
void* min_cols_thread(void* arg);

void min_max(Matrix* matrix);
void min_max_parallel(Matrix* matrix, int num_threads);
int min_max_by_columns(int rows, int cols, int num_threads);

//8. Normalizar una matriz columna por columna con el valor minimo y maximo
void normalize_matrix(Matrix* M, Vector* min, Vector* max);
void normalize_matrix_parallel(Matrix* matrix, Vector* max, Vector* min, int n);

//9. Normalizar una matriz columna por columna de acuerdo con la siguiente formula: x'=(x-u)/r, donde x’ es el nuevo valor que tomara cada elemento de la matriz, u es la media de cada columna y r es la desviacion estandar de cada columna.

Matrix* normalize_matrix_2(const Matrix* M);
Matrix* normalize_matrix_parallel_2(Matrix* matrix);

#endif
