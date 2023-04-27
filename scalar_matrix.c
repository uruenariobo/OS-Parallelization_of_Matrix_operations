#include <pthread.h>

typedef struct ThreadData ThreadData;
struct ThreadData {
    Matrix* M;
    double k;
    int row_start;
    int row_end;
};

void* scalar_matrix_thread(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    Matrix* M = data->M;
    double k = data->k;
    int row_start = data->row_start;
    int row_end = data->row_end;
    
    for (int i = row_start; i < row_end; ++i) {
        for (int j = 0; j < M->cols; ++j) {
            M->elements[i][j] *= k;
        }
    }
    
    pthread_exit(NULL);
}

void scalar_matrix_concurrent(Matrix* M, double k) {
    const int num_threads = M->rows;
    pthread_t threads[num_threads];
    ThreadData thread_data[num_threads];
    
    for (int i = 0; i < num_threads; ++i) {
        thread_data[i].M = M;
        thread_data[i].k = k;
        thread_data[i].row_start = i;
        thread_data[i].row_end = i + 1;
        
        pthread_create(&threads[i], NULL, scalar_matrix_thread, &thread_data[i]);
    }
    
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(threads[i], NULL);
    }
}