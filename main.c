#include <stdio.h>
#include "matrix.h"
#include "matrix.c"
#include <pthread.h>
#include <time.h>



int main()
{

	clock_t start, end;
	double cpu_time_used_secuencial, cpu_time_used_parallel;

	// Define el número de filas y columnas de una matriz
	int rows = 4;
	int cols = 4;

	// Crea una matriz con la cantidad de filas y columnas especificadas
	Matrix *M = create_matrix(rows, cols);

	// Inicializa la matriz con valores aleatorios
	init_matrix_rand(M);

	// Imprime la matriz en la consola
	printf("\nMatriz inicial:\n");
	print_matrix(M);

	// 1. Calcula la media de cada columna de la matriz y almacena los resultados en el vector creado, se muestra el tiempo de forma secuencial y paralela.
	//Comienza ejecución secuencial
	printf("\n1. Media de cada columna de una matriz:\n");
	start = clock();
	Vector *col_means = matrix_col_mean(M);
	end = clock();
	cpu_time_used_parallel = ((double)(end - start)) / CLOCKS_PER_SEC;
	//Imprime el tiempo de ejecución secuencial en segundos
	printf("%f\n", cpu_time_used_parallel);

	//Comienza ejecución paralela
	start = clock();
	col_means = matrix_col_mean_parallel(M);
	end = clock();
	cpu_time_used_parallel = ((double)(end - start)) / CLOCKS_PER_SEC;
	//Imprime el tiempo de ejecución paralela en segundos
	printf("%f\n", cpu_time_used_parallel);
	// Imprime el vector en la consola
	print_vector(col_means);

	// 2. Calcula la varianza de cada columna de la matriz y almacena los resultados en el vector creado, se muestra el tiempo de forma secuencial y paralela.
	printf("\n2. Varianza de cada columna de una matriz:\n");
	//Comienza ejecución secuencial
	start = clock();
	Vector *variances = matrix_col_vrz(M);
	end = clock();
	cpu_time_used_parallel = ((double)(end - start)) / CLOCKS_PER_SEC;
	//Imprime el tiempo de ejecución secuencial en segundos
	printf("%f\n", cpu_time_used_parallel);

	//Comienza ejecución paralela
	start = clock();
	variances = matrix_col_vrz_parallel(M);
	end = clock();
	cpu_time_used_parallel = ((double)(end - start)) / CLOCKS_PER_SEC;
	//Imprime el tiempo de ejecución paralela en segundos
	printf("%f\n", cpu_time_used_parallel);
	print_vector(variances);

	//3. Calcular la desviacion estandar de cada columna de una matriz
	printf("\n3. Desviacion estandar de cada columna de una matriz:\n");
	//Comienza ejecución secuencial
	start = clock();
	Vector *std = matrix_col_std(M);
	end = clock();
	cpu_time_used_parallel = ((double)(end - start)) / CLOCKS_PER_SEC;
	//Imprime el tiempo de ejecución secuencial en segundos
	printf("%f\n", cpu_time_used_parallel);

	//Comienza ejecución paralela
	start = clock();
	std = matrix_col_std_parallel(M);
	end = clock();
	cpu_time_used_parallel = ((double)(end - start)) / CLOCKS_PER_SEC;
	//Imprime el tiempo de ejecución paralela en segundos
	printf("%f\n", cpu_time_used_parallel);
	print_vector(std);

	//4. Calcular el valor mínimo y el valor máximo de cada columna de una matriz
	printf("\n4. Calcular el valor mínimo y el valor máximo de cada columna de una matriz\n");

	printf("\nEjecución secuencial:\n");
	start = clock();
	min_max(M);
	end = clock();
	cpu_time_used_parallel = ((double)(end - start)) / CLOCKS_PER_SEC;
	printf("\nTiempo secuencial:\n");
	printf("%f\n", cpu_time_used_parallel);

	printf("\nEjecución paralela:\n");
	start = clock();
	min_max_parallel(M, 4);
	end = clock();
	cpu_time_used_parallel = ((double)(end - start)) / CLOCKS_PER_SEC;
	printf("\nTiempo paralelo:\n");
	printf("%f\n", cpu_time_used_parallel);

	//5. Calcular la suma de dos matrices
	printf("\n5. Suma de dos matrices:\n");
	Matrix *Matrix1 = create_matrix(4, 4);
	init_matrix_rand(Matrix1);

	start = clock();
	Matrix *sum = add_matrix(M, Matrix1);
	end = clock();

	cpu_time_used_parallel = ((double)(end - start)) / CLOCKS_PER_SEC;
	//Imprime el tiempo de ejecución secuencial en segundos
	printf("%f\n", cpu_time_used_parallel);
	print_matrix(sum);

	start = clock();
	Matrix *sum_parallel = add_matrix_parallel(M, Matrix1, 4);
	end = clock();

	cpu_time_used_parallel = ((double)(end - start)) / CLOCKS_PER_SEC;
	//Imprime el tiempo de ejecución paralela en segundos
	printf("%f\n", cpu_time_used_parallel);
	print_matrix(sum_parallel);


	//6. Calcular el producto punto
	printf("\n6. Producto punto entre matrices:\n");
	
    //int n = 1000, i, j;
    //double cpu_time_used_secuencial, cpu_time_used_parallel;

    // Create matrices
    //matrix* mat1 = create_matrix2(n, n);
    //matrix* mat2 = create_matrix2(n, n);

    // Initialize matrices
    //matrix_init(&mat1, size, size);
    //matrix_init(&mat2, size, size);

    // Calculate product sequentially
    //matrix* result_seq = matrix_product_seq(&mat1, &mat2, &result_seq);
    //cpu_time_used_secuencial = ((double) (clock())) / CLOCKS_PER_SEC;

    // Calculate product in parallel
	
    //matrix* result_par = matrix_product_par(&mat1, &mat2, &result_par);
    //cpu_time_used_parallel = ((double) (clock())) / CLOCKS_PER_SEC;

    // Print matrices and result
    //printf("Matrix 1:\n");
    //for (i = 0; i < n; i++) {
    //    for (j = 0; j < n; j++) {
    //        printf("%lf ", mat1->data[i][j]);
    //    }
    //    printf("\n");
    //}

	//printf("\nMatrix 2:\n");
    //for (i = 0; i < n; i++) {
    //    for (j = 0; j < n; j++) {
    //        printf("%lf ", mat2->data[i][j]);
    //    }
    //    printf("\n");
    //}

    // Verify result
    //for (i = 0; i < n; i++) {
    //    for (j = 0; j < n; j++) {
    //        if (fabs(result_seq->data[i][j] - result_par->data[i][j]) > 0.001) {
    //            printf("Error: matrices are not equal\n");
    //            return -1;
    //        }
    //    }
    //}

	//printf("\nResultado del producto secuencial:\n");
    //for (i = 0; i < n; i++) {
    //    for (j = 0; j < n; j++) {
    //        printf("%lf ", result_seq->data[i][j]);
    //    }
    //    printf("\n");
    //}

    //printf("\nResultado del producto en paralelo:\n");
    //for (i = 0; i < n; i++) {
    //    for (j = 0; j < n; j++) {
    //        printf("%lf ", result_par->data[i][j]);
    //    }
    //    printf("\n");
    //}

    // Print results
    //printf("Time taken for sequential execution: %lf sec.\n", cpu_time_used_secuencial);
    //printf("Time taken for parallel execution: %lf sec.\n", cpu_time_used_parallel);


    // Free memory
	//free_matrix2(&mat1);
    //free_matrix2(&mat2);
    //free_matrix2(&result_seq);
    //free_matrix2(&result_par);
	////aaaaaaaaaq

	


	//7. Calcular la multiplicación de una matriz por un escalar
	printf("\n7. Multiplicar de una matriz por un escalar:\n");
	Matrix *Matrix2 = create_matrix(12, 4);
	init_matrix_rand(Matrix2);

	start = clock();
	scalar_matrix(Matrix2, 2.0);
	end = clock();

	cpu_time_used_parallel = ((double)(end - start)) / CLOCKS_PER_SEC;
	//Imprime el tiempo de ejecución secuencial en segundos
	printf("%f\n", cpu_time_used_parallel);
	print_matrix(Matrix2);

	start = clock();
	scalar_matrix_parallel(Matrix2, 2.0);
	end = clock();

	cpu_time_used_parallel = ((double)(end - start)) / CLOCKS_PER_SEC;
	//Imprime el tiempo de ejecución paralela en segundos
	printf("%f\n", cpu_time_used_parallel);
	print_matrix(Matrix2);

	// 8. Normalizar una matriz columna por columna x´ = (x - x(min) / (x(max) - x(min)))
	printf("\n8. Normalizar una matriz con los valores max y min:\n");
	Vector *max_numbers = matrix_col_max(M);
	Vector *min_numbers = matrix_col_min(M);

	printf("\nEjecución secuencial:\n");
	Matrix *normalize_max_min = M;
	start = clock();
	normalize_matrix(normalize_max_min, min_numbers, max_numbers);
	end = clock();
	cpu_time_used_parallel = ((double)(end - start)) / CLOCKS_PER_SEC;
	printf("\nTiempo secuencial:\n");
	printf("%f\n", cpu_time_used_parallel);

	printf("\nEjecución paralela:\n");
	start = clock();
	Matrix *normalize_max_min2 = M;
	normalize_matrix_parallel(normalize_max_min2, min_numbers, max_numbers, 4);
	end = clock();
	cpu_time_used_parallel = ((double)(end - start)) / CLOCKS_PER_SEC;
	printf("\nTiempo paralelo:\n");
	printf("%f\n", cpu_time_used_parallel);

	//9. Normalizar una matriz columna por columna de acuerdo con la siguiente formula: x'=(x-u)/r, donde x’ es el nuevo valor que tomara cada elemento de la matriz, u es la media de cada columna y r es la desviacion estandar de cada columna.
	printf("\n9. Normalizar una matriz columna por columna:\n");

	Matrix* M2 = create_matrix(rows, cols);
	// Inicializa la matriz con valores aleatorios
	init_matrix_rand(M2);
	//Comienza ejecución secuencial
	start = clock();
	Matrix* normal = normalize_matrix_2(M2);
	print_matrix(normal);
	end = clock();
	cpu_time_used_parallel = ((double)(end - start)) / CLOCKS_PER_SEC;
	//Imprime el tiempo de ejecución secuencial en segundos
	printf("%f\n", cpu_time_used_parallel);

	//Comienza ejecución paralela
	start = clock();
	normal = normalize_matrix_parallel_2(M2);
	end = clock();
	cpu_time_used_parallel = ((double)(end - start)) / CLOCKS_PER_SEC;
	//Imprime el tiempo de ejecución paralela en segundos
	printf("%f\n", cpu_time_used_parallel);
	print_matrix(normal);
}
