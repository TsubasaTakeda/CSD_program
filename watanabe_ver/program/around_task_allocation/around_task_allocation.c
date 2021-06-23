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

/*
行列(double型)のcol列要素の最大要素indexを返す関数
*/
int max_index_matrix_dim_2(const Matrix matrix, const int col)
{
    int max_index = -1;
    double max = DBL_MIN;
    for(int dim_1 = 0; dim_1 < matrix.num_row; dim_1++)
    {
        if(matrix.matrix[dim_1][col] > max)
        {
            max = matrix.matrix[dim_1][col];
            max_index = dim_1;
        }
        if(max_index == -1)
        {
            printf("Cannot max_element.");
        }
    }
    return max_index;
}

/*
行列(int型)の列要素の総和をあらわすベクトルを返す関数
*/
Vector_int sum_int_matrix_dim_2(Matrix_int matrix)
{
    Vector_int sum_int;
    sum_int.num_elements = matrix.num_col;
    sum_int.vector = malloc(sizeof(int) * sum_int.num_elements);


    for(int dim_2 = 0; dim_2 < matrix.num_col; dim_2++)
    {
        int sum = 0;
        for(int dim_1 = 0; dim_1 < matrix.num_row; dim_1++)
        {
            sum += matrix.matrix[dim_1][dim_2];
        }
        sum_int.vector[dim_2] = sum;
    }

    // printf("complete_2\n");

    return sum_int;
}

/*----------------------------------------------------------------------------------------------------------------------
odドライバー数に合わせた配分調整
----------------------------------------------------------------------------------------------------------------------*/
/*
odドライバーの足りない人数を計算する関数
*/
Vector_int calc_shortage_driver(Problem_struct data, Matrix_int int_allocation)
{
    // odドライバーの足りない人数を計算(real_driver - allocation)
    Vector_int shortage_driver;
    shortage_driver.num_elements = int_allocation.num_row;
    shortage_driver.vector = malloc(sizeof(int) * shortage_driver.num_elements);
    for (int o = 0; o < data.num_nodes; o++)
    {
        for (int d = 0; d < data.num_nodes; d++)
        {
            int od = o * data.num_nodes + d;
            shortage_driver.vector[od] = data.num_drivers.matrix[o][d] - sum_int_vector(int_allocation.matrix[od], int_allocation.num_col);
            // printf("shortage_driver[%d] = %d\n", od, shortage_driver.vector[od]);
        }
    }
    return shortage_driver;
}

/*
配分をodドライバー数に合わせる
*/
int Adjust_allocation_to_od_driver(Matrix_int int_allocation, Vector_int shortage_driver, Matrix variation)
{
    // odドライバーの人数を調整
    for (int od = 0; od < int_allocation.num_row; od++)
    {
        while (shortage_driver.vector[od] > 0)
        {
            int max_rs = max_index_vector(variation.matrix[od], variation.num_col);
            variation.matrix[od][max_rs] -= 1.0;
            int_allocation.matrix[od][max_rs]++;
            shortage_driver.vector[od]--;
        }
        // printf_vector_double(variation.matrix[od], variation.num_col);
    }
    return 0;
}

/*----------------------------------------------------------------------------------------------------------------------
rsタスク数に合わせた配分調整
----------------------------------------------------------------------------------------------------------------------*/

// rs荷主の足りない人数を計算(real_shipper - allocation)
Vector_int calc_shortage_shipper(Problem_struct data, Matrix_int int_allocation)
{
    // 現在の配分を列要素について足し合わせる
    Vector_int sum_allocation_dim_2;
    sum_allocation_dim_2 = sum_int_matrix_dim_2(int_allocation);


    Vector_int shortage_shipper;
    shortage_shipper.num_elements = int_allocation.num_col;
    shortage_shipper.vector = malloc(sizeof(int) * shortage_shipper.num_elements);
    for (int r = 0; r < data.num_depots; r++)
    {
        for (int s = 0; s < data.num_nodes; s++)
        {
            int rs = r * data.num_nodes + s;
            shortage_shipper.vector[rs] = data.num_shippers.matrix[r][s] - sum_allocation_dim_2.vector[rs];
        }
    }

    free(sum_allocation_dim_2.vector);

    return shortage_shipper;
}

/* 
excess → shortage へと配分を移す od を選択する関数
variation(増やしたい度合)の差(shortage - excess)が最も大きい od とする
*/
int select_od(Matrix variation, Matrix_int int_allocation, int shortage_rs, int excess_rs)
{
    double max = DBL_MIN;
    int select_od = -1;
    // printf("variation:");
    // printf_matrix_double(variation.matrix, variation.num_row, variation.num_col);
    // printf("shortage_rs = %d\texcess_rs = %d\n", shortage_rs, excess_rs);
    for(int od = 0; od < variation.num_row; od++)
    {
        if(int_allocation.matrix[od][excess_rs] > 0 && variation.matrix[od][shortage_rs] - variation.matrix[od][excess_rs] > max){
            max = variation.matrix[od][shortage_rs] - variation.matrix[od][excess_rs];
            select_od = od;
        }
    }
    if(select_od == -1){
        printf("Cannot find Appropreate od(allocation change excess to shortage)");
    }
    return select_od;
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


    // odドライバーの足りない人数を計算(real_driver - allocation)
    Vector_int shortage_driver = calc_shortage_driver(data, int_allocation);
    // for(int i = 0; i < shortage_driver.num_elements; i++)
    // {
    //     printf("%d\t", shortage_driver.vector[i]);
    // }
    
    // odドライバーの人数を調整
    Adjust_allocation_to_od_driver(int_allocation, shortage_driver, variation);

    free(shortage_driver.vector);
    // printf("variation_second:");
    // printf_matrix_double(variation.matrix, variation.num_row, variation.num_col);


    // rs荷主の足りない人数を計算(real_shipper - allocation)
    Vector_int shortage_shipper = calc_shortage_shipper(data, int_allocation);


    /*
    配分を増やす必要のあるindex(shortage_index)，
    配分を減らす必要のあるindex(excess_index)を計算
    */
    Vector_int shortage_index;
    Vector_int excess_index;
    shortage_index.num_elements = 0;
    excess_index.num_elements = 0;
    for(int rs = 0; rs < int_allocation.num_col; rs++)
    {
        if(shortage_shipper.vector[rs] > 0){
            shortage_index.num_elements++;
        } else if(shortage_shipper.vector[rs] < 0){
            excess_index.num_elements++;
        }
    }
    if(shortage_index.num_elements != 0)
    {
        shortage_index.vector = malloc(sizeof(int) * shortage_index.num_elements);
        excess_index.vector = malloc(sizeof(int) * excess_index.num_elements);
        int num_shortage = 0;
        int num_excess = 0;
        for(int rs = 0; rs < int_allocation.num_col; rs++)
        {
            if(shortage_shipper.vector[rs] > 0){
                shortage_index.vector[num_shortage] = rs;
                num_shortage++;
            }else if(shortage_shipper.vector[rs] < 0){
                excess_index.vector[num_excess] = rs;
                num_excess++;
            }
        }

        // 実行するrsタスク数を調整
        int search_shortage = 0;
        int search_excess = 0;
        while(1)
        {
            int search_rs_s = shortage_index.vector[search_shortage];
            int search_rs_e = excess_index.vector[search_excess];
            int od = select_od(variation, int_allocation, search_rs_s, search_rs_e);
            int_allocation.matrix[od][search_rs_s]++;
            int_allocation.matrix[od][search_rs_e]--;
            shortage_shipper.vector[search_rs_s]--;
            shortage_shipper.vector[search_rs_e]++;
            variation.matrix[od][search_rs_s] -= 1.0;
            variation.matrix[od][search_rs_e] -= 1.0;
            // printf("complete\n");
            if(shortage_shipper.vector[shortage_index.vector[search_shortage]] == 0)
            {
                // printf("end\n");
                if(search_shortage == shortage_index.num_elements-1)break;
                search_shortage++;
            }
            if(shortage_shipper.vector[excess_index.vector[search_excess]] == 0)
            {
                search_excess++;
            }
        }
        free(shortage_index.vector);
        free(excess_index.vector);
    }

    free(shortage_shipper.vector);

    mkdir("output");
    filename = "output\\int_allocation.csv";
    write_int_matrix_csv(filename, int_allocation.matrix, int_allocation.num_row, int_allocation.num_col);

    printf("complete\n");

    free_Matrix(allocation);
    free_Matrix_int(int_allocation);
    free_Matrix(variation);

    free_Matrix_int(data.num_drivers);
    free_Matrix_int(data.num_shippers);

    return 0;
}