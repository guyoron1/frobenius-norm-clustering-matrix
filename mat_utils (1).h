#ifndef MAT_UTILS_H
#define MAT_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DELIMITER ','

/**
 * @brief Calculate the element-wise difference between two matrices.
 *
 * @param output Output matrix to store the difference.
 * @param mat1 First input matrix.
 * @param mat2 Second input matrix.
 * @param rows Number of rows in the matrices.
 * @param cols Number of columns in the matrices.
 */
void calc_mat_difference(double **output, double** mat1, double** mat2, const int rows, const int cols);

/**
 * @brief Multiply two matrices.
 *
 * @param mat1 First input matrix.
 * @param mat2 Second input matrix.
 * @param rows1 Number of rows in the first matrix.
 * @param cols1 Number of columns in the first matrix.
 * @param cols2 Number of columns in the second matrix.
 * @return Resulting matrix after multiplication.
 */
double** multiply_matrixes(double** mat1, double** mat2, const int rows1, const int cols1, const int cols2);

/**
 * @brief Copy the contents of one matrix to another.
 *
 * @param dest Destination matrix.
 * @param src Source matrix.
 * @param rows Number of rows in the matrices.
 * @param cols Number of columns in the matrices.
 */
void copy_matrix(double **dest, double **src, const int rows, const int cols);

/**
 * @brief Calculate the transpose of a matrix.
 *
 * @param mat Input matrix.
 * @param rows Number of rows in the input matrix.
 * @param cols Number of columns in the input matrix.
 * @return Transposed matrix.
 */
double** calc_transpose(double** mat, const int rows, const int cols);

/**
 * @brief Given D, diagonal matrix, it calculates D^(-0.5).
 *
 * @param D Diagonal matrix.
 * @param n Dimension of the matrix.
 * @return Resulting matrix after power operation.
 */
double** calc_inverse_sqrt_diagonal(double **D, int n);

/**
 * @brief Calculate the Frobenius squared norm of a matrix.
 *
 * @param mat Input matrix.
 * @param rows Number of rows in the matrix.
 * @param cols Number of columns in the matrix.
 * @return Frobenius squared norm.
 */
double calc_frobenius_squared_norm(double **mat, const int rows, const int cols);

/**
 * @brief Calculate the squared Euclidean distance between two vectors.
 *
 * @param x First input vector.
 * @param y Second input vector.
 * @param d Dimension of the vectors.
 * @return Squared Euclidean distance.
 */
double calc_squared_euclidean_distance(double *x, double *y, int d);

/**
 * @brief Print a matrix to the standard output.
 *
 * @param matrix Input matrix.
 * @param rows Number of rows in the matrix.
 * @param cols Number of columns in the matrix.
 */
void print_matrix(double **matrix, int rows, int cols);

/**
 * @brief Allocate memory for a matrix.
 *
 * @param rows Number of rows in the matrix.
 * @param cols Number of columns in the matrix.
 * @return Allocated matrix.
 */
double **allocate_matrix(int rows, int cols);

/**
 * @brief Free the allocated memory for a matrix.
 *
 * @param matrix Input matrix.
 * @param rows Number of rows in the matrix.
 */
void free_matrix(double **matrix, int rows);

/**
 * @brief Counts the number of rows and columns in a filename that cointais the data.
 *
 * @param filename A string representing the name of the file to read.
 * @param rows A pointer to an integer where the function will store the number of rows.
 * @param cols A pointer to an integer where the function will store the number of columns.
 *
 * @return An integer indicating the success of the function.
 *         Returns 0 if the file was successfully read and dimensions were counted.
 *         Returns 1 if there was an error opening the file.
 */
int count_dimensions(const char *filename, int *num_rows, int *num_cols);

/**
 * @brief Read data from a file into a matrix.
 *
 * @param filename Name of the file to read data from.
 * @param n Pointer to store the number of data points (rows).
 * @param d Pointer to store the dimension of each data point (columns).
 * @return Data matrix read from the file.
 */
double **read_data(const char *filename, int *n, int *d);

#endif
