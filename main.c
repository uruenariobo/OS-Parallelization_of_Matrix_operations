#include <stdio.h>
#include "matrix.h"
#include "matrix.c"

int main(){

	printf("Hello world!\n");

	// Define el número de filas y columnas de una matriz
	int rows = 3;
	int cols = 4;

	// Crea una matriz con la cantidad de filas y columnas especificadas
	Matrix* M = create_matrix(rows, cols);

	// Inicializa la matriz con valores aleatorios
	init_matrix_rand(M);

	// Imprime la matriz en la consola
	print_matrix(M);

	// Crea un vector con la misma cantidad de columnas que la matriz
	Vector* mean = create_vector(cols);

	// Calcula la media de cada columna de la matriz y almacena los resultados en el vector creado
	mean_columns(M, rows, cols, mean);
	

	// Imprime el vector en la consola
	print_vector(mean);

	print_vector(calculate_std_deviation(M));
}
