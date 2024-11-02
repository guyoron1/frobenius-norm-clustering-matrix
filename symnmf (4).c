#include "symnmf.h"
#define EPS 1e-4
#define MAX_ITER 300
#define BETA 0.5

double** sym(double **X, const int n, const int d){
    int i,j;
    double **A = allocate_matrix(n, n);
    double dist;

    if (A == NULL) {
        return NULL;
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i != j) {
                dist = calc_squared_euclidean_distance(X[i], X[j], d);
                A[i][j] = exp((-dist)/2.0);
            } else {
                A[i][j] = 0.0;
            }
        }
    }
    return A;
}

/**
 * @brief Calculate the diagonal degree matrix.
 *
 * @param A Symmetric similarity matrix.
 * @param n Number of data points.
 * @return Diagonal degree matrix.
 */
double** calc_diagonal_degree_mat(double**A, const int n){
    double sum;
    int i, j;
    double **D = allocate_matrix(n, n);
    if (D == NULL){
        return NULL;
    }
    
    for (i = 0; i < n; i++) {
        sum = 0.0;
        for (j = 0; j < n; j++) {
            sum += A[i][j];
            D[i][j]= 0.0;
        }
        D[i][i] = sum;
    }
    return D;
}

/**
 * @brief Calculate the normalized symmetric matrix.
 *
 * @param A Symmetric similarity matrix.
 * @param D Diagonal degree matrix.
 * @param n Number of data points.
 * @return Normalized symmetric matrix.
 */
double** calc_normalized_sym(double **A, double**D, const int n){
    double **Q, **temp, **W;
    Q = calc_inverse_sqrt_diagonal(D, n);
    if (Q == NULL){
        return NULL;
    }
    temp = multiply_matrixes(Q, A, n, n, n);
    if (temp == NULL){
        free_matrix(Q, n);
        return NULL;
    }
    W = multiply_matrixes(temp, Q, n, n, n);
    if (W == NULL){
        free_matrix(Q, n);
        free_matrix(temp, n);
        return NULL;
    }

    free_matrix(Q, n);
    free_matrix(temp, n);
    return W;
}

double** ddg(double **X, const int n, const int d){
    double **D;
    double **A = sym(X, n, d);
    if (A==NULL){
        return NULL;
    }
    D = calc_diagonal_degree_mat(A, n);
    free_matrix(A, n);
    return D;
}

double** norm(double **X, const int n, const int d){
    double **A, **D, **W;
    A = sym(X, n, d);
    if (A==NULL){
        return NULL;
    }

    D = calc_diagonal_degree_mat(A, n);
    if (D==NULL){
        free_matrix(A, n);
        return NULL;
    }

    W = calc_normalized_sym(A, D, n);
    free_matrix(A, n);
    free_matrix(D, n);
    return W;
}

/**
 * @brief Calculate WH and HHtH matrices for the update rule.
 *
 * @param H Matrix H.
 * @param W Normalized symmetric matrix.
 * @param n Number of data points.
 * @param k Number of clusters.
 * @param WH_out Output WH matrix.
 * @param HHtH_out Output HHtH matrix.
 * @return 0 on success, 1 on failure.
 */
int calc_WH_HHth(double **H, double **W, const int n, const int k, double ***WH_out, double ***HHtH_out){
    double **WH, **HHtH, **Ht, **HHt;

    Ht = calc_transpose(H, n, k);
    if (Ht == NULL){
        return 1;
    }

    WH = multiply_matrixes(W, H, n, n, k);
    if (WH == NULL){
        free_matrix(Ht, k);
        return 1;
    }

    HHt = multiply_matrixes(H, Ht, n, k, n);
    if (HHt == NULL){
        free_matrix(Ht, k);
        free_matrix(WH, n);
        return 1;
    }

    HHtH = multiply_matrixes(HHt, H, n, n, k);
    if (HHtH == NULL){
        free_matrix(Ht, k);
        free_matrix(WH, n);
        free_matrix(HHt, n);
        return 1;
    }
    
    free_matrix(HHt, n);
    free_matrix(Ht, k);

    *WH_out = WH;
    *HHtH_out = HHtH;
    return 0;
}


/**
 * @brief Update the matrix H using the SymNMF update rule.
 *
 * @param H Matrix H.
 * @param W Normalized symmetric matrix.
 * @param n Number of data points.
 * @param k Number of clusters.
 * @return 0 on success, 1 on failure.
 */
int update_H(double **H, double **W, const int n, const int k){
    double **WH, **HHtH;
    int i, j;

    if (calc_WH_HHth(H, W, n, k, &WH, &HHtH) != 0) {
        return 1;
    }

    for (i = 0; i < n; i++){
        for(j = 0; j < k; j++){
            /*TODO: check if we need to handle the case where HHtH[i][j]=0 */
            H[i][j] = H[i][j] * (1-BETA + BETA * (WH[i][j] / HHtH[i][j]));
        }
    }

    free_matrix(WH, n);
    free_matrix(HHtH, n);
    return 0; 
}


double** symnmf(double **H, double **W, const int n, const int k){
    int iter;
    double f_norm_diff;
    double **H_old, **H_diff;

    H_old = allocate_matrix(n, k);
    if (H_old == NULL){
        return NULL;
    }

    H_diff = allocate_matrix(n, k); 
    if (H_diff == NULL){
        free_matrix(H_old, n);
        return NULL;
    }

    for (iter = 0; iter < MAX_ITER; iter++){
        printf("ITER = %d\n", iter);
        copy_matrix(H_old, H, n, k);
        if (update_H(H, W, n, k) != 0){
            free_matrix(H_old, n);
            free_matrix(H_diff, n);
            return NULL;
        }
        calc_mat_difference(H_diff, H, H_old, n, k);
        f_norm_diff = calc_frobenius_squared_norm(H_diff, n, k);
        if (f_norm_diff < EPS){
            break;
        } 
    }
    
    free_matrix(H_old, n);
    free_matrix(H_diff, n);
    return H;
}


int main(int argc, char *argv[]){
    int n, d;
    double **X;
    double **res = NULL;
    char *goal, *file_name;

    if(argc < 3){
        /** This will not happen because based on the instructions
         *  we can assume all the arguments are valid.
         *  we needed to use argc otherwise the code would not be able to compile.
         */
        printf("%s\n", ERROR_MESSAGE);
        return 1;
    }
    goal = argv[1];
    file_name = argv[2];
    X = read_data(file_name, &n, &d);

    if (X == NULL){
        printf("%s\n", ERROR_MESSAGE);
        return 1;
    }

    if (strcmp(goal, "sym") == 0) {
        res = sym(X, n, d);
    } else if (strcmp(goal, "ddg") == 0) {
        res = ddg(X, n, d);
    } else if (strcmp(goal, "norm") == 0) {
        res = norm(X, n, d);
    }
    
    free_matrix(X, n);
    if (res == NULL){
        printf("%s\n", ERROR_MESSAGE);
        return 1;
    }
    print_matrix(res, n, n);
    free_matrix(res, n);
    return 0;
}
