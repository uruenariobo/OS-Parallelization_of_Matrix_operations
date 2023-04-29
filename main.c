#include <stdio.h>
#include "matrix.h"
#include "matrix.c"
#include <pthread.h>
#include <time.h>

int main(){

	clock_t start, end;
	double cpu_time_used_secuencial, cpu_time_used_parallel;

	// Define el número de filas y columnas de una matriz
	int rows = 4;
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

	//3. Calcular la desviacion estandar de cada columna de una matriz
	printf("\nDesviacion estandar de cada columna de una matriz:\n");
	//Comienza ejecución secuencial
	start = clock();
	Vector* std = matrix_col_std(M);
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
	printf("\nCalcular el valor mínimo y el valor máximo de cada columna de una matriz\n");

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

	// 8. Normalizar una matriz columna por columna x´ = (x - x(min) / (x(max) - x(min)))
	printf("\nNormalizar una matriz con los valores max y min:\n");
    Vector* max_numbers = matrix_col_max(M);
    Vector* min_numbers = matrix_col_min(M);

	printf("\nEjecución secuencial:\n");
	start = clock();
	Matrix* normalize_max_min = M;
    normalize_matrix(normalize_max_min, min_numbers, max_numbers);
	end = clock();
	cpu_time_used_parallel = ((double)(end - start)) / CLOCKS_PER_SEC;
	printf("\nTiempo secuencial:\n");
	printf("%f\n", cpu_time_used_parallel);

	printf("\nEjecución paralela:\n");
	start = clock();
	Matrix* normalize_max_min2 = M;
    normalize_matrix_parallel(normalize_max_min2, min_numbers ,max_numbers, 4);
	end = clock();
	cpu_time_used_parallel = ((double)(end - start)) / CLOCKS_PER_SEC;
	printf("\nTiempo paralelo:\n");
	printf("%f\n", cpu_time_used_parallel);


	//9. Normalizar una matriz columna por columna de acuerdo con la siguiente formula: x'=(x-u)/r, donde x’ es el nuevo valor que tomara cada elemento de la matriz, u es la media de cada columna y r es la desviacion estandar de cada columna.
	printf("\nNormalizar una matriz columna por columna:\n");
	//Comienza ejecución secuencial
	start = clock();
	Matrix* normal = normalize_matrix_2(M);
	end = clock();
	cpu_time_used_parallel = ((double)(end - start)) / CLOCKS_PER_SEC;
	//Imprime el tiempo de ejecución secuencial en segundos
	printf("%f\n", cpu_time_used_parallel);
	
	//Comienza ejecución paralela
	start = clock();
	normal = normalize_matrix_parallel_2(M);
	end = clock();
	cpu_time_used_parallel = ((double)(end - start)) / CLOCKS_PER_SEC;
	//Imprime el tiempo de ejecución paralela en segundos
	printf("%f\n", cpu_time_used_parallel);
	print_matrix(normal);
}
