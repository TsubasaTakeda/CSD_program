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

/*----------------------------------------------------------------------------------------------------------------------
行列関係
----------------------------------------------------------------------------------------------------------------------*/

/*
3次元行列メモリ領域確保の関数
第一引数：対象の3次元行列構造体
*/
Matrix get_matrix_memory_area(Matrix matrix)
{

    matrix.matrix = malloc(sizeof(double *) * matrix.num_row);
    for (int dim_1 = 0; dim_1 < matrix.num_row; dim_1++)
    {
        matrix.matrix[dim_1] = malloc(sizeof(double) * matrix.num_col);
    }

    return matrix;
}

/*
行列のメモリ領域を解放する関数
*/
int free_Matrix(Matrix A)
{
    for (int row = 0; row < A.num_row; row++)
    {
        free(A.matrix[row]);
    }
    free(A.matrix);

    return 0;
}

/*
行列のメモリ領域を解放する関数
*/
int free_Matrix_int(Matrix_int A)
{
    for (int row = 0; row < A.num_row; row++)
    {
        free(A.matrix[row]);
    }
    free(A.matrix);

    return 0;
}

/*
ベクトル(int型)の総和を計算する関数
*/
int sum_int_vector(const int* vector, const int num_elements)
{
    int sum = 0;
    for(int i = 0; i < num_elements; i++)
    {
        sum += vector[i];
    }
    return sum;
}

/*
ベクトル(double型)の最大要素indexを返す関数
*/
int max_index_vector(const double* vector, const int num_elements)
{
    double max = INT_MIN;
    int max_index = -1;
    for(int i = 0; i < num_elements; i++)
    {
        if(vector[i] > max){
            max = vector[i];
            max_index = i;
        }
    }
    if(max_index == -1){
        printf("Cannot find max_element");
    }
    return max_index;
}

/*----------------------------------------------------------------------------------------------------------------------
メイン関数
----------------------------------------------------------------------------------------------------------------------*/

int main(){

    Problem_struct data;

    // 実数値のタスク配分行列の読み込み
    char* filename = ".\\input\\task_allocation.csv";
    Matrix allocation;
    read_matrix_csv(filename, &allocation.matrix, &allocation.num_row, &allocation.num_col);
    // printf_matrix_double(allocation.matrix, allocation.num_row, allocation.num_col);

    // ドライバー行列の読み込み
    filename = ".\\input\\num_drivers.csv";
    read_int_matrix_csv(filename, &data.num_drivers.matrix, &data.num_drivers.num_row, &data.num_drivers.num_col);
    data.num_nodes = data.num_drivers.num_row;

    // 荷主行列の読み込み
    filename = ".\\input\\num_shippers.csv";
    read_int_matrix_csv(filename, &data.num_shippers.matrix, &data.num_shippers.num_row, &data.num_shippers.num_col);
    data.num_depots = data.num_shippers.num_row;


    // 整数値のタスク配分行列のメモリ領域確保 & 実数値の切り捨て数値を代入
    Matrix_int int_allocation;
    int_allocation.num_row = allocation.num_row;
    int_allocation.num_col = allocation.num_col;
    int_allocation.matrix = malloc(sizeof(int*) * int_allocation.num_row);
    for(int od = 0; od < int_allocation.num_row; od++)
    {
        int_allocation.matrix[od] = malloc(sizeof(int) * int_allocation.num_col);
        for(int rs = 0; rs < int_allocation.num_col; rs++)
        {
            int_allocation.matrix[od][rs] = (int)allocation.matrix[od][rs];
            // printf("int_allocation[%d][%d] = %d\n", od, rs, int_allocation.matrix[od][rs]);
        }
    }


    // 実数値-整数値の差分行列
    Matrix variation;
    variation.num_row = allocation.num_row;
    variation.num_col = allocation.num_col;
    variation = get_matrix_memory_area(variation);
    for(int od = 0; od < variation.num_row; od++)
    {
        for(int rs = 0; rs < variation.num_col; rs++)
        {
            // printf("(double)int_allocation[%d][%d] = %lf\n", od, rs, (double)int_allocation.matrix[od][rs]);
            // printf("allocation[%d][%d] = %lf\n", od, rs, allocation.matrix[od][rs]);
            variation.matrix[od][rs] = allocation.matrix[od][rs] - (double)int_allocation.matrix[od][rs];
        }
    }

    // printf("variation_first:");
    // printf_matrix_double(variation.matrix, variation.num_row, variation.num_col);


    // ドライバーODペア足りない人数を計算
    Vector_int shortage_driver;
    shortage_driver.num_elements = allocation.num_row;
    shortage_driver.vector = malloc(sizeof(int) * shortage_driver.num_elements);
    for(int o = 0; o < data.num_nodes; o++)
    {
        for(int d = 0; d < data.num_nodes; d++)
        {
            int od = o*data.num_nodes + d;
            shortage_driver.vector[od] = data.num_drivers.matrix[o][d] - sum_int_vector(int_allocation.matrix[od], int_allocation.num_col);
            // printf("shortage_driver[%d] = %d\n", od, shortage_driver.vector[od]);
        }
    }
    
    // odドライバーの人数を調整完了
    for(int od = 0; od < allocation.num_row; od++)
    {
        while(shortage_driver.vector[od] > 0)
        {
            int max_rs = max_index_vector(variation.matrix[od], variation.num_col);
            variation.matrix[od][max_rs] -= 1.0;
            int_allocation.matrix[od][max_rs]++;
            shortage_driver.vector[od]--;
        }
        // printf_vector_double(variation.matrix[od], variation.num_col);
    }

    free(shortage_driver.vector);
    printf("variation_second:");
    printf_matrix_double(variation.matrix, variation.num_row, variation.num_col);


    free_Matrix(allocation);
    free_Matrix_int(int_allocation);
    free_Matrix(variation);

    free_Matrix_int(data.num_drivers);
    free_Matrix_int(data.num_shippers);

    return 0;
}