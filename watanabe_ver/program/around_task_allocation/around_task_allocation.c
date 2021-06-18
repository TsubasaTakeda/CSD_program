#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <direct.h>

#define LOGIT_DRIVER 5.0
#define CONVERGENCE 0.00001
#define ETA_ACCEL_PARAM 1.1
#define FIRST_LIPS 0.01
#define DISPLAY_ITERATION 100

// ベクトル構造体
typedef struct
{
    int num_elements;
    double *vector;
} Vector;
typedef struct
{
    int num_elements;
    int *vector;
} Vector_int;

// 行列構造体
typedef struct
{
    int num_row; // 行数
    int num_col; // 列数
    double **matrix;
} Matrix;
typedef struct
{
    int num_row;
    int num_col;
    int **matrix;
} Matrix_int;

// pのためだけの3次元行列構造体
typedef struct
{
    int num_row_1;
    int num_row_2;
    int num_row_3;
    double ***matrix;
} Matrix_3_dim;

// 加速勾配法の返り値構造体
typedef struct
{
    Vector sol;
    Matrix sol_matrix;
    double CPU_time_sec;
    int iteration;
    int num_calls_obj;
    int num_calls_nabla;
} Optimization;

// 問題情報の構造体
typedef struct
{
    Matrix cost_between_nodes;
    Matrix_int num_shippers;
    Matrix_int num_drivers;
    Vector_int depots_index;
    double LOGIT_param_driver;
    int num_nodes;
    int num_depots;
} Problem_struct;