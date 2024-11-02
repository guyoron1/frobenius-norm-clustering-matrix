#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "symnmf.h"

enum Action{
    SYM,
    DDG,
    NORM,
};

int get_matrix_rows(PyObject *py_mat){
    if (!PyList_Check(py_mat)) {
        return 0;
    }

    PyObject *row0 = PyList_GetItem(py_mat, 0);
    if (PyList_Check(row0)) {
        return (int)PyList_Size(py_mat);
    } else {
        return 1;
    }
}

int get_matrix_cols(PyObject *py_mat){
     if (!PyList_Check(py_mat)) {
        return 0;
    }

    if (PyList_Size(py_mat) == 0) {
        return 0;
    }

    PyObject *row0 = PyList_GetItem(py_mat, 0);
    if (PyList_Check(row0)) {
        return (int)PyList_Size(row0);
    } else {
        return (int)PyList_Size(py_mat);
    }
}

double** PyObject_to_double_mat(PyObject *py_mat, int N, int d){
    int i,j;

    if (N == 0 || d == 0){
        return NULL;
    }
    
    double** mat = allocate_matrix(d, N);
    if (mat == NULL){
        return NULL;
    }
    for (i = 0; i < N; i++){
        PyObject* vector = PyList_GetItem(py_mat, i);
        if (!PyList_Check(vector)) {
            free_matrix(mat, N);
            PyErr_SetString(PyExc_TypeError, ERROR_MESSAGE);
            return NULL;
        }

        for (j = 0; j < d; j++) {
            PyObject* item = PyList_GetItem(vector, j);
            if (!PyFloat_Check(item)) {
                free_matrix(mat, N);
                PyErr_SetString(PyExc_TypeError, ERROR_MESSAGE);
                return NULL;
            }
            mat[i][j] = PyFloat_AsDouble(item);
        }
    }
    return mat;
}

PyObject* double_mat_to_PyObject(double** mat, int N, int d){
    int i, j;
    PyObject* pyList = PyList_New(N);
    if (!pyList) {
        return NULL;
    }
    for (i = 0; i < N; i++){
    
        PyObject* pyPoint = PyList_New(d);
        if (!pyPoint) {
            Py_DECREF(pyList);
            return NULL;
        }
        for (j = 0; j < d; j++){
            PyObject* pyElement = PyFloat_FromDouble(mat[i][j]);
            if (!pyElement) {
                Py_DECREF(pyPoint);
                Py_DECREF(pyList);
                return NULL;
            }
            PyList_SetItem(pyPoint, j, pyElement);
        }
        PyList_SetItem(pyList, i, pyPoint);
    }
    return pyList;
}



PyDoc_STRVAR(sym_doc,
"sym(arg1)\n"
"It returns the similarity matrix of the provided argument\n"
"\n"
"Parameters:\n"
"    arg1 (float[][]): X - data points.\n"
"\n"
"Returns:\n"
"    float[][]: similarity matrix A\n"
"\n"
"Preconditions:\n"
"    All the given data point are different \n");

static PyObject *py_sym(PyObject *self, PyObject *args){
    return execute_partial_symnmf(args, SYM);
}


PyDoc_STRVAR(ddg_doc,
"ddg(arg1)\n"
"It returns the diagonal degree matrix of the provided argument\n"
"\n"
"Parameters:\n"
"    arg1 (float[][]): X - data points.\n"
"\n"
"Returns:\n"
"    float[][]: diagonal degree matrix D\n"
"\n"
"Preconditions:\n"
"    All the given data point are different \n");

static PyObject *py_ddg(PyObject *self, PyObject *args){
    return execute_partial_symnmf(args, DDG);
}


PyDoc_STRVAR(norm_doc,
"norm(arg1)\n"
"It returns the normalized similarity matrix of the provided argument\n"
"\n"
"Parameters:\n"
"    arg1 (float[][]): X - data points.\n"
"\n"
"Returns:\n"
"    float[][]: normalized similarity matrix W\n"
"\n"
"Preconditions:\n"
"    All the given data point are different \n");

static PyObject *py_norm(PyObject *self, PyObject *args){
    return execute_partial_symnmf(args, NORM);
}

PyObject *execute_partial_symnmf(PyObject *args, enum Action action){
    PyObject *py_X, *py_res;
    double **X, **res_mat;
    int N, d;

    if (!PyArg_ParseTuple(args, "O", &py_X)) {
        return NULL;
    }

    if (!PyList_Check(py_X)) {
        return NULL;
    }
    
    N = get_matrix_rows(py_X);
    d = get_matrix_cols(py_X);
    X = PyObject_to_double_mat(py_X, N, d);
    
    if (X == NULL) {
        return NULL;
    }

    switch (action) {
        case SYM:
            res_mat = sym(X, N, d);
            break;
        case DDG:
            res_mat = ddg(X, N, d);
            break;
        case NORM:
            res_mat = norm(X, N, d);
            break;
    }

    if (res_mat == NULL) {
        free_matrix(X, N);
        PyErr_SetString(PyExc_RuntimeError, ERROR_MESSAGE);
        return NULL;
    }

    py_res = double_mat_to_PyObject(res_mat, N, N);
    free_matrix(X, N);
    free_matrix(res_mat, N);
    return py_res;
}

PyDoc_STRVAR(symnmf_doc,
"symnmf(arg1, arg2, arg3, arg4)\n"
"It solves the symNMF algorithm on the provided data\n"
"\n"
"Parameters:\n"
"    arg1 (float[][]): H - initial H.\n"
"    arg2 (float[][]): W - normalized similarity matrix.\n"
"    arg3 (float[][]): N - number of rows in the original data.\n"
"    arg4 (float[][]): k - number of required cluesters.\n"
"\n"
"Returns:\n"
"    float[][]: similarity matrix A\n"
"\n"
"Preconditions:\n"
"    All the given data point are different \n"
"    W is the normalized similarity matrix \n"
"    H is randomly initialized with the values from the interval [0, 2*sqrt(m/k)], where m is the average of all entrie of W \n"
"    k < N \n");

PyObject *py_symnmf(PyObject *self, PyObject *args){
    PyObject *py_H,*py_W, *py_res;
    double **H, **W, **updated_H;
    int N, k;

    if (!PyArg_ParseTuple(args, "OOii", &py_H, &py_W, &N, &k)) {
        return NULL;
    }

    if (!PyList_Check(py_H) || !PyList_Check(py_W)) {
        return NULL;
    }

    H = PyObject_to_double_mat(py_H, N, k);
    W = PyObject_to_double_mat(py_W, N, N);
    
    if (H == NULL || W == NULL) {
        free_matrix(H, N);
        free_matrix(W, N);
        return NULL;
    }

    updated_H = symnmf(H, W, N, k);

    if (updated_H == NULL) {
        free_matrix(H, N);
        free_matrix(W, N);
        PyErr_SetString(PyExc_RuntimeError, ERROR_MESSAGE);
        return NULL;
    }

    py_res = double_mat_to_PyObject(updated_H, N, k);
    free_matrix(H, N);
    free_matrix(W, N);
    free_matrix(updated_H, N);
    return py_res;
}


static PyMethodDef symnmfMethods[] = {
    {"sym", py_sym, METH_VARARGS, sym_doc},
    {"ddg", py_ddg, METH_VARARGS, ddg_doc},
    {"norm", py_norm, METH_VARARGS, norm_doc},
    {"symnmf", py_symnmf, METH_VARARGS,symnmf_doc},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef symnmfmoduledef = {
    PyModuleDef_HEAD_INIT,
    "symnmfmodule",
    NULL,
    -1,
    symnmfMethods
};

PyMODINIT_FUNC PyInit_symnmfmodule(void) {
    PyObject *m;
    m = PyModule_Create(&symnmfmoduledef);
    if (!m) {
        return NULL;
    }
    return m;
}