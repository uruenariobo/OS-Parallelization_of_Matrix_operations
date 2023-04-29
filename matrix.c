#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <pthread.h>

Vector *create_vector(int size)
{
    Vector *v = calloc(1, sizeof(Vector));
    v->size = size;
    v->elements = calloc(size, sizeof(double));
    return v;
}

Vector *create_vector_from_file(const char *file_path, int size)
{
    Vector *v = create_vector(size);
    FILE *fp = fopen(file_path, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Failed to open \"%s\". Not possible to create vector from file.\n", file_path);
        return NULL;
    }

    double d;
    for (int i = 0; i < v->size; ++i)
    {
        const int r = fscanf(fp, "%lf", &d);
        if (r != 1)
        {
            fprintf(stderr, "fscanf failed.\n");
            fclose(fp);
            return NULL;
        }
        v->elements[i] = d;
    }

    fclose(fp);
    return v;
}

Matrix *create_matrix(int rows, int cols)
{
    Matrix *M = malloc(sizeof(Matrix));
    M->rows = rows;
    M->cols = cols;
    M->elements = calloc(rows, sizeof(double *));
    for (int i = 0; i < rows; ++i)
    {
        M->elements[i] = calloc(cols, sizeof(double));
    }

    return M;
}

Matrix *create_matrix_from_file(const char *file_path, int rows, int cols)
{
    Matrix *M = create_matrix(rows, cols);
    FILE *fp = fopen(file_path, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Failed to open \"%s\". Not possible to create matrix from file.\n", file_path);
        return NULL;
    }

    double d;
    for (int i = 0; i < M->rows; ++i)
    {
        for (int j = 0; j < M->cols; ++j)
        {
            const int r = fscanf(fp, "%lf", &d);
            if (r != 1)
            {
                fprintf(stderr, "fscanf failed.\n");
                fclose(fp);
                return NULL;
            }
            M->elements[i][j] = d;
        }
    }

    fclose(fp);
    return M;
}

void init_matrix_rand(Matrix *M)
{
    for (int i = 0; i < M->rows; ++i)
    {
        for (int j = 0; j < M->cols; ++j)
        {
            M->elements[i][j] = (double)rand() / (double)RAND_MAX;
        }
    }
}

void init_vector_rand(Vector *V)
{
    for (int i = 0; i < V->size; ++i)
    {
        V->elements[i] = (double)rand() / (double)RAND_MAX;
    }
}

void copy_matrix(Matrix *dst, const Matrix *src)
{
    for (int i = 0; i < src->rows; ++i)
    {
        for (int j = 0; j < src->cols; ++j)
        {
            dst->elements[i][j] = src->elements[i][j];
        }
    }
}

void copy_vector(Vector *dst, const Vector *src)
{
    for (int i = 0; i < src->size; ++i)
    {
        dst->elements[i] = src->elements[i];
    }
}

void free_vector(Vector *v)
{
    free(v->elements);
    free(v);
}

void free_matrix(Matrix *M)
{
    for (int i = 0; i < M->rows; ++i)
    {
        free(M->elements[i]);
    }
    free(M->elements);
    free(M);
}

Vector *add_vector(const Vector *a, const Vector *b)
{
    if (a->size != b->size)
    {
        fprintf(stderr, "Invalid size. %d and %d\n", a->size, b->size);
        return NULL;
    }

    Vector *r = create_vector(a->size);
    for (int i = 0; i < r->size; ++i)
    {
        r->elements[i] = a->elements[i] + b->elements[i];
    }

    return r;
}

Vector *dot_vector_matrix(const Vector *v, const Matrix *M)
{
    if (v->size != M->rows)
    {
        fprintf(stderr, "Invalid size. %d and (%d, %d)\n", v->size, M->rows, M->cols);
        return NULL;
    }

    Vector *r = create_vector(M->cols);
    for (int i = 0; i < r->size; ++i)
    {
        double d = 0.0;
        for (int j = 0; j < M->rows; ++j)
        {
            d += (v->elements[j] * M->elements[j][i]);
        }
        r->elements[i] = d;
    }

    return r;
}

Matrix *dot_matrix(const Matrix *M, const Matrix *N)
{
    if (M->cols != N->rows)
    {
        fprintf(stderr, "Invalid size. (%d, %d) and (%d, %d)\n", M->rows, M->cols, N->rows, N->cols);
        return NULL;
    }

    Matrix *A = create_matrix(M->rows, N->cols);
    for (int i = 0; i < M->rows; ++i)
    {
        for (int j = 0; j < N->cols; ++j)
        {
            double d = 0.0;
            for (int k = 0; k < M->cols; ++k)
            {
                d += M->elements[i][k] * N->elements[k][j];
            }
            A->elements[i][j] = d;
        }
    }

    return A;
}

Vector *matrix_col_mean(const Matrix *M)
{
    Vector *V = create_vector(M->cols);

    for (int i = 0; i < M->cols; ++i)
    {
        double sum = 0.0;
        for (int j = 0; j < M->rows; ++j)
        {
            sum += M->elements[j][i];
        }
        V->elements[i] = sum / M->rows;
    }

    return V;
}

Vector *matrix_col_sum(const Matrix *M)
{
    Vector *v = create_vector(M->cols);

    for (int i = 0; i < M->cols; ++i)
    {
        double sum = 0.0;
        for (int j = 0; j < M->rows; ++j)
        {
            sum += M->elements[j][i];
        }
        v->elements[i] = sum;
    }

    return v;
}

void scalar_matrix(Matrix *M, double k)
{
    for (int i = 0; i < M->rows; ++i)
    {
        for (int j = 0; j < M->cols; ++j)
        {
            M->elements[i][j] *= k;
        }
    }
}

void print_vector(const Vector *v)
{
    printf("size=%d, [", v->size);
    for (int i = 0; i < v->size; ++i)
    {
        printf("%lf ", v->elements[i]);
    }
    printf("]\n");
}

void print_matrix(const Matrix *M)
{
    printf("rows=%d,cols=%d,[", M->rows, M->cols);
    for (int i = 0; i < M->rows; ++i)
    {
        printf("[");
        for (int j = 0; j < M->cols; ++j)
        {
            printf("%lf ", M->elements[i][j]);
        }
        printf("]\n");
    }
    printf("]\n");
}

//1. Calcular la media de cada columna de una matriz
typedef struct ThreadData ThreadData;

struct ThreadData
{
    const Matrix *M;
    Vector *col_means;
    int col_index;
};

void *compute_col_mean(void *arg)
{
    ThreadData *data = (ThreadData *)arg;
    const Matrix *M = data->M;
    Vector *col_means = data->col_means;
    int col_index = data->col_index;

    double sum = 0.0;
    for (int i = 0; i < M->rows; i++)
    {
        sum += M->elements[i][col_index];
    }
    col_means->elements[col_index] = sum / M->rows;

    pthread_exit(NULL);
}

Vector *matrix_col_mean_parallel(const Matrix *M)
{
    Vector *col_means = create_vector(M->cols);

    pthread_t threads[M->cols];
    ThreadData thread_data[M->cols];

    for (int i = 0; i < M->cols; i++)
    {
        thread_data[i].M = M;
        thread_data[i].col_means = col_means;
        thread_data[i].col_index = i;
        pthread_create(&threads[i], NULL, compute_col_mean, (void *)&thread_data[i]);
    }

    for (int i = 0; i < M->cols; i++)
    {
        pthread_join(threads[i], NULL);
    }

    return col_means;
}

//2. Calcular la varianza de cada columna de una matriz
typedef struct ThreadData1 ThreadData1;
struct ThreadData1
{
    Matrix *M;
    int col;
    double result;
};

Vector *matrix_col_vrz(const Matrix *M)
{
    Vector *variance = create_vector(M->cols);

    for (int j = 0; j < M->cols; ++j)
    {
        double sum = 0.0;
        double sum_sq = 0.0;

        for (int i = 0; i < M->rows; ++i)
        {
            double x = M->elements[i][j];
            sum += x;
            sum_sq += x * x;
        }

        double mean = sum / M->rows;
        variance->elements[j] = sum_sq / M->rows - mean * mean;
    }
    return variance;
}

struct ThreadArgs
{
    Matrix *M;
    Vector *variances;
    int col;
};

void *compute_column_variance(void *arg)
{
    struct ThreadArgs *targs = (struct ThreadArgs *)arg;
    double sum = 0.0, sq_sum = 0.0;
    for (int i = 0; i < targs->M->rows; ++i)
    {
        double x = targs->M->elements[i][targs->col];
        sum += x;
        sq_sum += x * x;
    }
    double mean = sum / targs->M->rows;
    double variance = sq_sum / targs->M->rows - mean * mean;
    targs->variances->elements[targs->col] = variance;
    return NULL;
}

Vector *matrix_col_vrz_parallel(const Matrix *M)
{
    Vector *variances = create_vector(M->cols);

    pthread_t threads[M->cols];
    struct ThreadArgs targs[M->cols];
    for (int i = 0; i < M->cols; ++i)
    {
        targs[i].M = M;
        targs[i].variances = variances;
        targs[i].col = i;
        pthread_create(&threads[i], NULL, compute_column_variance, &targs[i]);
    }

    for (int i = 0; i < M->cols; ++i)
    {
        pthread_join(threads[i], NULL);
    }

    return variances;
}

//3. Calcular la desviacion estandar de cada columna de una matriz
Vector *matrix_col_std(const Matrix *matrix)
{
    Vector *variances = matrix_col_vrz(matrix);
    Vector *standard_deviations = create_vector(variances->size);
    for (int j = 0; j < matrix->cols; j++)
    {
        standard_deviations->elements[j] = sqrt(variances->elements[j]);
    }
    return standard_deviations;
}

typedef struct
{
    Matrix *matrix;
    int column;
    double *result;
} ColumnData;

void *calculate_standard_deviation(void *data)
{
    Vector *v = (Vector *)data;
    for (int i = 0; i < v->size; ++i)
    {
        v->elements[i] = sqrt(v->elements[i]);
    }
    pthread_exit(NULL);
}

Vector *matrix_col_std_parallel(const Matrix *matrix)
{
    Vector *variances = matrix_col_vrz_parallel(matrix);
    Vector *standard_deviations = create_vector(variances->size);
    copy_vector(standard_deviations, variances);
    pthread_t thread;
    pthread_create(&thread, NULL, calculate_standard_deviation, standard_deviations);
    pthread_join(thread, NULL);
    free_vector(variances);
    return standard_deviations;
}

// 4.Calcular el valor mínimo y el valor máximo de cada columna de una matriz 


// Minimos

Vector* matrix_col_min(const Matrix* M){
    Vector* v = create_vector(M->cols);
    for (int i = 0; i < M->cols; ++i) {
        double min = M->elements[0][i];
        for (int j = 0; j < M->rows; ++j) {
            if (min > M->elements[j][i]) {
                min = M->elements[j][i];
            }
        }
        v->elements[i] = min;
    }
    return v; 
}

typedef struct {
    struct Matrix* matrix;
    pthread_mutex_t* mutex;
    Vector* min_values;
    int start_col;
    int end_col;
} MinArgs;

typedef struct {
    pthread_mutex_t mutex;
    double min_value;
} Min;

void* min_cols_thread(void* arg){
    MinArgs* args = (MinArgs*) arg;

    for (int i = 0; i < args->matrix->cols; ++i) {
        double min = args->matrix->elements[0][i];
        for (int j = 0; j < args->matrix->rows; ++j) {
            if (min > args->matrix->elements[j][i]) {
                pthread_mutex_lock(args->mutex);
                min = args->matrix->elements[j][i];
                pthread_mutex_unlock(args->mutex);
            }
        }
        pthread_mutex_lock(args->mutex);
        args->min_values->elements[i] = min;
        pthread_mutex_unlock(args->mutex);
    }
    pthread_exit(NULL);
    return NULL;
}

// Maximos 

Vector* matrix_col_max(const Matrix* M){
    Vector* v = create_vector(M->cols);
    for (int i = 0; i < M->cols; ++i) {
        double max = M->elements[0][i];
        for (int j = 0; j < M->rows; ++j) {
            if (max < M->elements[j][i]) {
                max = M->elements[j][i];
            }
        }
        v->elements[i] = max;
    }
    return v; 
}

typedef struct {
    struct Matrix* matrix;
    pthread_mutex_t* mutex;
    Vector* max_values;
    int start_col;
    int end_col;
} MaxArgs;

typedef struct {
    pthread_mutex_t mutex;
    double max_value;
} Max;

void* max_cols_thread(void* arg){
    MaxArgs* args = (MaxArgs*) arg;

    for (int i = 0; i < args->matrix->cols; ++i) {
        double max = args->matrix->elements[0][i];
        for (int j = 0; j < args->matrix->rows; ++j) {
            if (max < args->matrix->elements[j][i]) {
                pthread_mutex_lock(args->mutex);
                max = args->matrix->elements[j][i];
                pthread_mutex_unlock(args->mutex);
            }
        }
        pthread_mutex_lock(args->mutex);
        args->max_values->elements[i] = max;
        pthread_mutex_unlock(args->mutex);
    }
    pthread_exit(NULL);
    return NULL;
}

void min_max(Matrix* matrix){
    Vector *max_col = matrix_col_max(matrix);
    Vector *min_col = matrix_col_min(matrix);
    printf("Los valores mínimos de la matriz: \n");
    print_vector(min_col);
    printf("Los valores maximos de la matriz: \n");
    print_vector(max_col);
}

// Calculo de valores minimos y valores maximos con paralelismo

void min_max_parallel(Matrix* matrix, int num_threads){

    pthread_t threads[num_threads];

    MinArgs min_args[num_threads];
    MaxArgs max_args[num_threads];

    Vector* min_values = create_vector(matrix->cols);
    Vector* max_values = create_vector(matrix->cols);

    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

    const int pthread_size = matrix->cols / num_threads;

    int start_col = 0;
    for (int i = 0; i < num_threads; ++i){
        int end_col = start_col + pthread_size;
        if (i == num_threads - 1){
            end_col = matrix->cols;
        }

        min_args[i] = (MinArgs){.start_col = start_col, .end_col = end_col, .matrix = matrix, .min_values = min_values, .mutex = &mutex};
        max_args[i] = (MaxArgs){.start_col = start_col, .end_col = end_col, .matrix = matrix, .max_values = max_values, .mutex = &mutex};

        pthread_create(&threads[i], NULL, min_cols_thread, &min_args[i]);
        start_col = end_col;
        pthread_create(&threads[i], NULL, max_cols_thread, &max_args[i]);
        start_col = end_col;
    }

    for (int j = 0; j < num_threads; j++) {
        pthread_join(threads[j], NULL);
    }

    pthread_mutex_destroy(&mutex);
    printf("Los valores mínimos de la matriz: \n");
    print_vector(min_values);
    printf("Los valores maximos de la matriz: \n");
    print_vector(max_values);
}

int min_max_by_columns(int rows, int cols, int num_threads) {
    Matrix* matrix = create_matrix(rows, cols);
    init_matrix_rand(matrix);
    print_matrix(matrix);
    min_max_parallel(matrix, num_threads);
    return 0;
}


// 8. Normalizar una matriz columna por columna x´ = (x - x(min) / (x(max) - x(min)))

typedef struct {
    struct Matrix* matrix;
    Vector* min;
    Vector* max;
    int start_row;
    int end_row;
    pthread_mutex_t* lock;
} ColumnNormalizeData;

void normalize_matrix(Matrix* M, Vector* min, Vector* max){
    for (int i = 0; i < M->rows; ++i) {
        for (int j = 0; j < M->cols; ++j) {
            M->elements[i][j] = (M->elements[i][j] - min->elements[j]) / (max->elements[j] - min->elements[j]);
        }
    }

    print_matrix(M);
}  

void* normalize(void* args) {
    ColumnNormalizeData* args_cnd = (ColumnNormalizeData*) args;
    for (int i = args_cnd->start_row; i < args_cnd->end_row; ++i) {
        for (int j = 0; j < args_cnd->matrix->cols; ++j) {
            pthread_mutex_lock(args_cnd->lock);
            args_cnd->matrix->elements[i][j] = (args_cnd->matrix->elements[i][j] - args_cnd->min->elements[j]) / (args_cnd->max->elements[j] - args_cnd->min->elements[j]);
            pthread_mutex_unlock(args_cnd->lock);
        }
    }

    return NULL;
}

void normalize_matrix_parallel(Matrix* matrix, Vector* min, Vector* max, int n) {

    pthread_mutex_t lock;
    pthread_mutex_init(&lock, NULL);
    pthread_t threads[n];

    ColumnNormalizeData args[n];

    int pthread_size = matrix->rows / n;
    int extra_rows = matrix->rows % n;
    
    for (int i = 0; i < n; ++i) {
        args[i].matrix = matrix;
        args[i].min = min;
        args[i].max = max;
        args[i].lock = &lock;

        int start_row = i * pthread_size;
        int end_row = (i + 1) * pthread_size;

        if (i == n - 1) {
            end_row += extra_rows;
        }
        args[i].start_row = start_row;
        args[i].end_row = end_row;

        pthread_create(&threads[i], NULL, normalize, &args[i]);
    }

    for (int i = 0; i < n; ++i) {
        pthread_join(threads[i], NULL);
    }

    pthread_mutex_destroy(&lock);
    print_matrix(matrix);
    free_matrix(matrix);
}

  
//5. Calcular la suma de dos matrices
typedef struct ThreadArgsMatrixSum ThreadArgsMatrixSum;
struct ThreadArgsMatrixSum
{
    const Matrix *M;
    const Matrix *N;
    Matrix *R;
    int *row_counter;
};

Matrix *add_matrix(const Matrix *M, const Matrix *N)
{
    if (M->rows != N->rows || M->cols != N->cols)
    {
        fprintf(stderr, "Invalid size. (%d, %d) and (%d, %d)\n", M->rows, M->cols, N->rows, N->cols);
        return NULL;
    }

    Matrix *R = create_matrix(M->rows, M->cols);
    for (int i = 0; i < M->rows; ++i)
    {
        for (int j = 0; j < M->cols; ++j)
        {
            R->elements[i][j] = M->elements[i][j] + N->elements[i][j];
        }
    }

    return R;
}

void *add_rows_thread(void *arg)
{
    ThreadArgsMatrixSum *args = (ThreadArgsMatrixSum *)arg;
    const Matrix *M = args->M;
    const Matrix *N = args->N;
    Matrix *R = args->R;
    int *row_counter = args->row_counter;

    int row;
    while (1)
    {
        row = (*row_counter)++;
        if (row >= M->rows)
        {
            break;
        }

        for (int j = 0; j < M->cols; ++j)
        {
            R->elements[row][j] = M->elements[row][j] + N->elements[row][j];
        }
    }

    pthread_exit(NULL);
}

Matrix *add_matrix_parallel(const Matrix *M, const Matrix *N, int n_threads)
{
    if (M->rows != N->rows || M->cols != N->cols)
    {
        fprintf(stderr, "Invalid size. (%d, %d) and (%d, %d)\n", M->rows, M->cols, N->rows, N->cols);
        return NULL;
    }

    Matrix *R = create_matrix(M->rows, M->cols);

    int row_counter = 0;
    ThreadArgsMatrixSum args[n_threads];

    pthread_t threads[n_threads];
    for (int i = 0; i < n_threads; ++i)
    {
        args[i].M = M;
        args[i].N = N;
        args[i].R = R;
        args[i].row_counter = &row_counter;

        pthread_create(&threads[i], NULL, add_rows_thread, (void *)&args[i]);
    }

    for (int i = 0; i < n_threads; ++i)
    {
        pthread_join(threads[i], NULL);
    }

    return R;
}

//7. Calcular la multiplicación de una matriz por un escalar
typedef struct ThreadDataScalarMatrix ThreadDataScalarMatrix;
struct ThreadDataScalarMatrix
{
    Matrix *M;
    double k;
    int row_start;
    int row_end;
};

void scalar_vector(Vector *V, double k)
{
    for (int i = 0; i < V->size; ++i)
    {
        V->elements[i] *= k;
    }
}

void *scalar_matrix_thread(void *arg)
{
    ThreadDataScalarMatrix *data = (ThreadDataScalarMatrix *)arg;
    Matrix *M = data->M;
    double k = data->k;
    int row_start = data->row_start;
    int row_end = data->row_end;

    for (int i = row_start; i < row_end; ++i)
    {
        for (int j = 0; j < M->cols; ++j)
        {
            M->elements[i][j] *= k;
        }
    }

    pthread_exit(NULL);
}

void scalar_matrix_parallel(Matrix *M, double k)
{
    const int num_threads = M->rows;
    pthread_t threads[num_threads];
    ThreadDataScalarMatrix thread_data[num_threads];

    for (int i = 0; i < num_threads; ++i)
    {
        thread_data[i].M = M;
        thread_data[i].k = k;
        thread_data[i].row_start = i;
        thread_data[i].row_end = i + 1;

        pthread_create(&threads[i], NULL, scalar_matrix_thread, &thread_data[i]);
    }

    for (int i = 0; i < num_threads; ++i)
    {
        pthread_join(threads[i], NULL);
    }
}

//9. Normalizar una matriz columna por columna de acuerdo con la siguiente formula: x'=(x-u)/r, donde x’ es el nuevo valor que tomara cada elemento de la matriz, u es la media de cada columna y r es la desviacion estandar de cada columna.
Matrix *normalize_matrix2(const Matrix *M)
{
    Matrix *result = create_matrix(M->rows, M->cols);
    Vector *means = matrix_col_mean(M);
    Vector *stds = matrix_col_std(M);

    for (int j = 0; j < M->cols; ++j)
    {
        double mean = means->elements[j];
        double std = stds->elements[j];

        for (int i = 0; i < M->rows; ++i)
        {
            double x = M->elements[i][j];
            result->elements[i][j] = (x - mean) / std;
        }
    }

    free_vector(means);
    free_vector(stds);

    return result;
}

  
typedef struct
{
    Matrix *M;
    Vector *means;
    Vector *stds;
    int col_idx;
} ColumnNormalizeData2;

  
void *normalize_column(void *data)
{
    ColumnNormalizeData2 *d = (ColumnNormalizeData2 *)data;
    Matrix *M = d->M;
    Vector *means = d->means;
    Vector *stds = d->stds;
  
    int col_idx = d->col_idx;

    for (int i = 0; i < M->rows; ++i)
    {
        double x = M->elements[i][col_idx];
        double u = means->elements[col_idx];
        double r = stds->elements[col_idx];
        M->elements[i][col_idx] = (x - u) / r;
    }

    pthread_exit(NULL);
}


Matrix *normalize_matrix_parallel2(Matrix *M)
{
    Vector *means = matrix_col_mean_parallel(M);
    Vector *stds = matrix_col_std_parallel(M);
    int num_cols = M->cols;
    pthread_t threads[num_cols];
    ColumnNormalizeData2 thread_data[num_cols];

    for (int i = 0; i < num_cols; ++i)
    {
        thread_data[i].M = M;
        thread_data[i].means = means;
        thread_data[i].stds = stds;
        thread_data[i].col_idx = i;
        int result = pthread_create(&threads[i], NULL, normalize_column, (void *)&thread_data[i]);
        if (result != 0)
        {
            fprintf(stderr, "Failed to create thread for column %d. Error code: %d\n", i, result);
            return NULL;
        }
    }

    for (int i = 0; i < num_cols; ++i)
    {
        int result = pthread_join(threads[i], NULL);
        if (result != 0)
        {
            fprintf(stderr, "Failed to join thread for column %d. Error code: %d\n", i, result);
            return NULL;
        }
    }

    return M;
}