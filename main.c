#include <stdio.h>
#include "matrix.h"
#include "matrix.c"
#include <pthread.h>
#include <time.h>

int main(){

	clock_t start, end;
	double cpu_time_used_secuencial, cpu_time_used_parallel;

	// Define el número de filas y columnas de una matriz
	int rows = 3;
	int cols = 4;

	// Crea una matriz con la cantidad de filas y columnas especificadas
	Matrix* M = create_matrix(rows, cols);

	// Inicializa la matriz con valores aleatorios
	init_matrix_rand(M);

	// Imprime la matriz en la consola
	printf("\nMatriz inicial:\n");
	print_matrix(M);

	// 1. Calcula la media de cada columna de la matriz y almacena los resultados en el vector creado, se muestra el tiempo de forma secuencial y paralela.
	//Comienza ejecución secuencial
	printf("\nMedia de cada columna de una matriz:\n");
	start = clock();
	Vector* col_means = matrix_col_mean(M);
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
	printf("\nVarianza de cada columna de una matriz:\n");
	//Comienza ejecución secuencial
	start = clock();
	Vector* variances = matrix_col_vrz(M);
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
}
