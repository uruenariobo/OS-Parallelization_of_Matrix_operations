#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <pthread.h>

Vector* create_vector(int size) {
    Vector* v = calloc(1, sizeof(Vector));
    v->size = size;
    v->elements = calloc(size, sizeof(double));
    return v;
}

Vector* create_vector_from_file(const char* file_path, int size) {
    Vector* v = create_vector(size);
    FILE* fp = fopen(file_path, "r");
    if (fp == NULL) {
        fprintf(stderr, "Failed to open \"%s\". Not possible to create vector from file.\n", file_path);
        return NULL;
    }

    double d;
    for (int i = 0; i < v->size; ++i) {
        const int r = fscanf(fp, "%lf", &d);
        if (r != 1) {
            fprintf(stderr, "fscanf failed.\n");
            fclose(fp);
            return NULL;
        }
        v->elements[i] = d;
    }

    fclose(fp);
    return v;
}

Matrix* create_matrix(int rows, int cols) {
    Matrix* M = malloc(sizeof(Matrix));
    M->rows = rows;
    M->cols = cols;
    M->elements = calloc(rows, sizeof(double*));
    for (int i = 0; i < rows; ++i) {
        M->elements[i] = calloc(cols, sizeof(double));
    }

    return M;
}

Matrix* create_matrix_from_file(const char* file_path, int rows, int cols) {
    Matrix* M = create_matrix(rows, cols);
    FILE* fp = fopen(file_path, "r");
    if (fp == NULL) {
        fprintf(stderr, "Failed to open \"%s\". Not possible to create matrix from file.\n", file_path);
        return NULL;
    }

    double d;
    for (int i = 0; i < M->rows; ++i) {
        for (int j = 0; j < M->cols; ++j) {
            const int r = fscanf(fp, "%lf", &d);
            if (r != 1) {
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

void init_matrix_rand(Matrix* M) {
    for (int i = 0; i < M->rows; ++i) {
        for (int j = 0; j < M->cols; ++j) {
            M->elements[i][j] = (double)rand() / (double)RAND_MAX;
        }
    }
}

void init_vector_rand(Vector* V){
    for(int i = 0; i < V->size; ++i){
	V->elements[i] = (double)rand() / (double)RAND_MAX;
    }
}

void copy_matrix(Matrix* dst, const Matrix* src) {
    for (int i = 0; i < src->rows; ++i) {
        for (int j = 0; j < src->cols; ++j) {
            dst->elements[i][j] = src->elements[i][j];
        }
    }
}

void copy_vector(Vector* dst, const Vector* src) {
    for (int i = 0; i < src->size; ++i) {
        dst->elements[i] = src->elements[i];
    }
}

void free_vector(Vector* v) {
    free(v->elements);
    free(v);
}

void free_matrix(Matrix* M) {
    for (int i = 0; i < M->rows; ++i) {
        free(M->elements[i]);
    }
    free(M->elements);
    free(M);
}

Vector* add_vector(const Vector* a, const Vector* b) {
    if (a->size != b->size) {
        fprintf(stderr, "Invalid size. %d and %d\n", a->size, b->size);
        return NULL;
    }

    Vector* r = create_vector(a->size);
    for (int i = 0; i < r->size; ++i) {
        r->elements[i] = a->elements[i] + b->elements[i];
    }

    return r;
}

Vector* dot_vector_matrix(const Vector* v, const Matrix* M) {
    if (v->size != M->rows) {
        fprintf(stderr, "Invalid size. %d and (%d, %d)\n", v->size, M->rows, M->cols);
        return NULL;
    }

    Vector* r = create_vector(M->cols);
    for (int i = 0; i < r->size; ++i) {
        double d = 0.0;
        for (int j = 0; j < M->rows; ++j) {
            d += (v->elements[j] * M->elements[j][i]);
        }
        r->elements[i] = d;
    }

    return r;
}

Matrix* add_matrix(const Matrix* M, const Matrix* N) {
    if (M->rows != N->rows || M->cols != N->cols) {
        fprintf(stderr, "Invalid size. (%d, %d) and (%d, %d)\n", M->rows, M->cols, N->rows, N->cols);
        return NULL;
    }

    Matrix* R = create_matrix(M->rows, M->cols);
    for (int i = 0; i < M->rows; ++i) {
        for (int j = 0; j < M->cols; ++j) {
            R->elements[i][j] = M->elements[i][j] + N->elements[i][j];
        }
    }

    return R;
}

Matrix* dot_matrix(const Matrix* M, const Matrix* N) {
    if (M->cols != N->rows) {
        fprintf(stderr, "Invalid size. (%d, %d) and (%d, %d)\n", M->rows, M->cols, N->rows, N->cols);
        return NULL;
    }

    Matrix* A = create_matrix(M->rows, N->cols);
    for (int i = 0; i < M->rows; ++i) {
        for (int j = 0; j < N->cols; ++j) {
            double d = 0.0;
            for (int k = 0; k < M->cols; ++k) {
                d += M->elements[i][k] * N->elements[k][j];
            }
            A->elements[i][j] = d;
        }
    }

    return A;
}

Vector* matrix_col_mean(const Matrix* M) {
    Vector* V = create_vector(M->cols);

    for (int i = 0; i < M->cols; ++i) {
        double sum = 0.0;
        for (int j = 0; j < M->rows; ++j) {
            sum += M->elements[j][i];
        }
        V->elements[i] = sum / M->rows;
    }

    return V;
}

Vector* matrix_col_sum(const Matrix* M) {
    Vector* v = create_vector(M->cols);

    for (int i = 0; i < M->cols; ++i) {
        double sum = 0.0;
        for (int j = 0; j < M->rows; ++j) {
            sum += M->elements[j][i];
        }
        v->elements[i] = sum;
    }

    return v;
}

void scalar_matrix(Matrix* M, double k) {
    for (int i = 0; i < M->rows; ++i) {
       for (int j = 0; j < M->cols; ++j) {
            M->elements[i][j] *= k;
       }
    }
}

void scalar_vector(Vector* V, double k) {
    for (int i = 0; i < V->size; ++i) {
        V->elements[i] *= k;
    }
}

void print_vector(const Vector* v) {
    printf("size=%d, [", v->size);
    for (int i = 0; i < v->size; ++i) {
        printf("%lf ", v->elements[i]);
    }
    printf("]\n");
}

void print_matrix(const Matrix* M) {
    printf("rows=%d,cols=%d,[", M->rows, M->cols);
    for (int i = 0; i < M->rows; ++i) {
        printf("[");
        for (int j = 0; j < M->cols; ++j) {
            printf("%lf ", M->elements[i][j]);
        }
        printf("]\n");
    }
    printf("]\n");
}

typedef struct ThreadData ThreadData;

struct ThreadData {
    const Matrix* M;
    Vector* col_means;
    int col_index;
};

void* compute_col_mean(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    const Matrix* M = data->M;
    Vector* col_means = data->col_means;
    int col_index = data->col_index;

    double sum = 0.0;
    for (int i = 0; i < M->rows; i++) {
        sum += M->elements[i][col_index];
    }
    col_means->elements[col_index] = sum / M->rows;

    pthread_exit(NULL);
}

Vector* matrix_col_mean_parallel(const Matrix* M) {
    Vector* col_means = create_vector(M->cols);

    pthread_t threads[M->cols];
    ThreadData thread_data[M->cols];

    for (int i = 0; i < M->cols; i++) {
        thread_data[i].M = M;
        thread_data[i].col_means = col_means;
        thread_data[i].col_index = i;
        pthread_create(&threads[i], NULL, compute_col_mean, (void*)&thread_data[i]);
    }

    for (int i = 0; i < M->cols; i++) {
        pthread_join(threads[i], NULL);
    }

    return col_means;
}

typedef struct ThreadData1 ThreadData1;
struct ThreadData1 {
    Matrix* M;
    int col;
    double result;
};

Vector* matrix_col_vrz(const Matrix* M) {
    Vector* variance = create_vector(M->cols);

    for (int j = 0; j < M->cols; ++j) {
        double sum = 0.0;
        double sum_sq = 0.0;

        for (int i = 0; i < M->rows; ++i) {
            double x = M->elements[i][j];
            sum += x;
            sum_sq += x * x;
        }

        double mean = sum / M->rows;
        variance->elements[j] = sum_sq / M->rows - mean * mean;
    }
    return variance;
}

struct ThreadArgs {
    Matrix* M;
    Vector* variances;
    int col;
};

void* compute_column_variance(void* arg) {
    struct ThreadArgs* targs = (struct ThreadArgs*)arg;
    double sum = 0.0, sq_sum = 0.0;
    for (int i = 0; i < targs->M->rows; ++i) {
        double x = targs->M->elements[i][targs->col];
        sum += x;
        sq_sum += x * x;
    }
    double mean = sum / targs->M->rows;
    double variance = sq_sum / targs->M->rows - mean * mean;
    targs->variances->elements[targs->col] = variance;
    return NULL;
}

Vector* matrix_col_vrz_parallel(const Matrix* M) {
    Vector* variances = create_vector(M->cols);

    pthread_t threads[M->cols];
    struct ThreadArgs targs[M->cols];
    for (int i = 0; i < M->cols; ++i) {
        targs[i].M = M;
        targs[i].variances = variances;
        targs[i].col = i;
        pthread_create(&threads[i], NULL, compute_column_variance, &targs[i]);
    }

    for (int i = 0; i < M->cols; ++i) {
        pthread_join(threads[i], NULL);
    }

    return variances;
}