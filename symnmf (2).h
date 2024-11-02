#ifndef SYMNMF_H
#define SYMNMF_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mat_utils.h"

#define ERROR_MESSAGE "An Error Has Occurred"


/**
 * @brief Calculate the symmetric similarity matrix.
 *
 * @param X Input data matrix.
 * @param n Number of data points.
 * @param d Dimension of each data point.
 * @return Symmetric similarity matrix.
 */
double **sym(double **X, const int n, const int d);


/**
 * @brief Calculate the diagonal degree matrix from input data.
 *
 * @param X Input data matrix.
 * @param n Number of data points.
 * @param d Dimension of each data point.
 * @return Diagonal degree matrix.
 */
double **ddg(double **X, const int n, const int d);

/**
 * @brief Calculate the normalized symmetric matrix from input data.
 *
 * @param X Input data matrix.
 * @param n Number of data points.
 * @param d Dimension of each data point.
 * @return Normalized symmetric matrix.
 */
double **norm(double **X, const int n, const int d);

/**
 * @brief Perform the Symmetric Non-negative Matrix Factorization (SymNMF).
 *
 * @param X Input data matrix.
 * @param H Initial matrix H.
 * @param W Normalized symmetric matrix.
 * @param n Number of data points.
 * @param d Dimension of each data point.
 * @param k Number of clusters.
 * @return Factorized matrix H.
 */
double** symnmf(double **H, double **W, const int n, const int k);

#endif