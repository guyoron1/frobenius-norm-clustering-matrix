#include "mat_utils.h"

void calc_mat_difference(double **output, double** mat1, double** mat2, const int rows, const int cols){
    int i, j;
    for (i = 0; i < rows; i++){
        for (j = 0; j < cols; j++){
             output[i][j] = mat1[i][j] - mat2[i][j];
        }
    }
}

double** multiply_matrixes(double** mat1, double** mat2, const int rows1, const int cols1, const int cols2){
     int i, j, k;
    double **result = allocate_matrix(rows1, cols2);
    if (result == NULL) {
        return NULL;
    }
    for (i = 0; i < rows1; i++) {
        for (j = 0; j < cols2; j++) {
            result[i][j] = 0.0;
            for (k = 0; k < cols1; k++) {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    return result;
}

void copy_matrix(double **dest, double **src, const int rows, const int cols){
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            dest[i][j] = src[i][j];
        }
    }
}

double** calc_transpose(double** mat, const int rows, const int cols){
    int i, j;
    double **transposed = allocate_matrix(cols, rows);
    if (transposed == NULL){
        return NULL;
    }
    for (i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
            transposed[j][i] = mat[i][j];
        }
    }
    return transposed;
}

double** calc_inverse_sqrt_diagonal(double **D, int n){
    int i, j;
    double** Q = allocate_matrix(n, n);
    if (Q == NULL){
        return NULL;
    }

    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (i == j) {
                if (D[i][j] == 0) {
                    printf("Error: Diagonal element is zero, cannot take inverse square root.\n");
                    exit(1);
                }
                Q[i][j] = 1.0 / sqrt(D[i][j]);
            } else {
                Q[i][j] = 0.0;
            }
        }
    }
    return Q;
}

double calc_frobenius_squared_norm(double **mat, const int rows, const int cols){
    int i, j;
    double norm = 0.0;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            norm += mat[i][j] * mat[i][j];
        }
    }
    return norm;
}

double calc_squared_euclidean_distance(double *x, double *y, int d){
    int k = 0;
    double sum = 0.0;
    for (k = 0; k < d; k++) {
        sum += pow(x[k] - y[k], 2);
    }
    return sum;
}

void print_matrix(double **matrix, int rows, int cols) {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            printf("%.4f", matrix[i][j]);
            if (j < cols - 1) {
                printf(",");
            }
        }
        printf("\n");
    }
}

double **allocate_matrix(int rows, int cols) {
    int i;
    double **mat = (double **)malloc(rows * sizeof(double *));
    if (mat == NULL) {
        return NULL;
    }
    for (i = 0; i < rows; i++) {
        mat[i] = (double *)malloc(cols * sizeof(double));
        if (mat[i] == NULL) {
            free_matrix(mat, rows);
            return NULL;
        }
    }
    return mat;
}

void free_matrix(double **mat, int rows) {
    int i;
    for (i = 0; i < rows; i++) {
        free(mat[i]);
    }
    free(mat);
}

int count_dimensions(const char *filename, int *rows, int *cols) {
    FILE *file;
    int c;
    double dummy;

    (*rows) = 0;
    (*cols) = 0;
    
    file = fopen(filename, "r");
    if (file == NULL) {
        return 1;
    }

    
    while (fscanf(file, "%lf", &dummy) != EOF){
        c = fgetc(file);

        if ((*rows) == 0)
            (*cols)++;

        if (c == '\n' || c == EOF)
            (*rows)++;
    }
    fclose(file);
    return 0;
}

double **read_data(const char *filename, int *n, int *d) {
    FILE *file;
    int c, row, col;
    double** X;

    if (count_dimensions(filename, n, d) != 0){
        return NULL;
    }
    
    X = allocate_matrix(*n, *d);
    if (X == NULL){
        return NULL;
    }
    
    file = fopen(filename, "r");
    if (file == NULL) {
        return NULL;
    }

    row = 0;
    col = 0;
    while (fscanf(file, "%lf", &(X[row][col++])) != EOF)
    {
        c = fgetc(file);
        if (c == '\n'){
            row++;
            col = 0;
        }
        else if (c != DELIMITER && c != EOF)
        {
            free_matrix(X, *n);
            return NULL;
        }
    }
    fclose(file);
    return X;
}

