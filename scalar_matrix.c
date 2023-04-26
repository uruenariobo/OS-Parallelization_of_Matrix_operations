#include <pthread.h>
#include "matrix.h"
#include "matrix.c"

typedef struct
{
  const Matrix *M;
  int rows;
  int cols;
  double k;
  Matrix *R;
  int start_row; // fila de inicio asignada a este hilo
  int end_row;   // fila final asignada a este hilo (no procesada)
} ScalarMatrixArgs;

void *scalar_matrix_thread(void *arg)
{
  ScalarMatrixArgs *args = (ScalarMatrixArgs *)arg;

  scalar_matrix(args->M, args->end_row - args->start_row, args->cols, args->k, args->R);

  return NULL;
}

// Aquí va el proceso concurrente
void scalar_matrix_concurrent(const Matrix *M, int rows, int cols, double k, Matrix *R)
{
  const int NUM_THREADS = 4; // Numero de hilos
  pthread_t threads[NUM_THREADS];
  ScalarMatrixArgs args[NUM_THREADS];

  int rows_per_thread = rows / NUM_THREADS;
  int last_thread_rows = rows_per_thread + rows % NUM_THREADS;

  for (int t = 0; t < NUM_THREADS; ++t)
  {
    args[t].M = M;
    args[t].rows = rows;
    args[t].cols = cols;
    args[t].k = k;
    args[t].R = R;
    args[t].start_row = t * rows_per_thread;
    args[t].end_row = args[t].start_row + rows_per_thread;

    pthread_create(&threads[t], NULL, scalar_matrix_thread, &args[t]);
  }

  for (int t = 0; t < NUM_THREADS; ++t)
  {
    pthread_join(threads[t], NULL);
  }
}