#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <time.h>

// ベクトル構造体
typedef struct {
    int num_elements;
    double* vector;
} Vector;
typedef struct {
    int num_elements;
    int* vector;
} Vector_int;

// 行列構造体
typedef struct {
    int num_row; // 行数
    int num_col; // 列数
    double** matrix;
} Matrix;
typedef struct {
    int num_row;
    int num_col;
    int** matrix;
} Matrix_int;

// 加速勾配法の返り値構造体
typedef struct {
    Vector sol;
    double CPU_time_sec;
    int iteration;
    int num_calls_obj;
    int num_calls_nabla;
} Optimization;

// 問題情報の構造体
typedef struct {
    Matrix cost_between_nodes;
    Matrix_int num_shippers;
    Matrix_int num_drivers;
    Vector_int depots_index;
    double LOGIT_param_driver;
    int num_nodes;
    int num_depots;
} Problem_struct;




/*----------------------------------------------------------------------------------------------------
加速勾配法
----------------------------------------------------------------------------------------------------*/

/*
backtracking (AlgorithmファイルのFISTAのアルゴリズムを参照)
第一引数：現在のリプシッツ定数，第二引数：リプシッツ定数更新用の定数，
第三引数：目的関数，第四引数：勾配関数，第五引数：利用する暫定解，
第六引数：目的関数の呼び出し回数ポインタ，第七引数：勾配関数の呼び出し回数ポインタ
返り値：iota(本来は自然数だが，便宜上double型で返す)
*/
double backtracking(const double L, const double eta, const double (*function)(const Vector), const Vector (*nabla_function)(const Vector), const Vector x, int *num_calls_obj, int *num_calls_nabla)
{

    double F, Q;

    double iota = 0.0;
    Vector y;
    y.vector = malloc(sizeof(double) * x.num_elements);
    y.num_elements = x.num_elements;

    double obj = function(x);
    num_calls_obj[0]++;

    Vector nabla = nabla_function(x);
    num_calls_nabla[0]++;

    double nolm = nolm_2(nabla.vector, x.num_elements);

    // 確認用
    // printf_vector_double(nabla, x_elements);
    // printf("obj = %f\nnolm = %f\n", obj, nolm);

    while (1)
    {

        Q = obj - nolm / (2 * pow(eta, iota) * L);

        // Fの移動先のベクトルを取得
        copy_vector_double(y.vector, nabla.vector, x.num_elements);
        k_times_double(y.vector, -1.0 / (pow(eta, iota) * L), x.num_elements);
        sum_b_to_a_double(y.vector, x.vector, x.num_elements);

        F = function(y);
        num_calls_obj[0]++;

        if (F <= Q)
        {
            break;
        }

        iota += 1.0;
    }

    // printf("F = %f\nQ = %f\n", F, Q);

    free(nabla.vector);
    free(y.vector);

    return iota;
}

/*
FISTA with リスタート
第一引数：収束判定値(勾配の最大要素)，第二引数：リプシッツ定数更新用の定数(1.0より大きい)，第三引数：初期リプシッツ定数(ステップサイズの逆数，0.0より大きい)，
第四引数：目的関数，第五引数：勾配関数，第六引数：初期解
返り値：0
*/
Optimization FISTA_with_restart(const double epsilon, const double eta, double L, const double (*function)(const Vector), const Vector (*nabla_function)(const Vector), const Vector x_first)
{

    Optimization solution;

    int x_elements = x_first.num_elements;

    // 暫定解等の領域確保
    Vector x, y;
    x.vector = malloc(sizeof(double) * x_elements);
    x.num_elements = x_elements;
    y.vector = malloc(sizeof(double) * x_elements);
    y.num_elements = x_elements;
    Vector delta_x;
    delta_x.vector = malloc(sizeof(double) * x_elements);
    delta_x.num_elements = x_elements;
    Vector x_nabla;
    x_nabla.num_elements = x_elements;
    Vector y_nabla;
    y_nabla.vector = malloc(sizeof(double) * x_elements);
    y_nabla.num_elements = x_elements;

    // 初期解の代入
    copy_vector_double(x.vector, x_first.vector, x_elements);
    copy_vector_double(y.vector, x_first.vector, x_elements);
    zeros_vector_double(delta_x.vector, x_elements);

    // iota(勾配方向のステップサイズを決定する変数)の保管場所
    double iota;

    // 慣性方向のステップサイズを決定するための変数の置き場所
    double t_0 = 1.0;
    double t_1 = 1.0;

    // 反復回数，目的関数呼び出し回数，勾配関数呼び出し回数
    int iteration = 0;
    int num_calls_obj = 0;
    int num_calls_nabla = 0;

    long cpu_start_time, cpu_end_time;
    double sec;

    cpu_start_time = clock();

    while (1)
    {

        // backtracking によるリプシッツ定数の取得
        iota = backtracking(L, eta, function, nabla_function, y, &num_calls_obj, &num_calls_nabla);
        L = pow(eta, iota) * L;

        // リプシッツ定数(ステップサイズの逆数)の確認
        // printf("Lipsitz = %f\n", L);

        // 暫定解の更新
        free(y_nabla.vector);
        y_nabla = nabla_function(y);
        num_calls_nabla++;
        k_times_double(y_nabla.vector, -1 / L, x_elements);
        copy_vector_double(x.vector, y.vector, x_elements);
        sum_b_to_a_double(x.vector, y_nabla.vector, x_elements);
        sum_b_to_a_double(delta_x.vector, y_nabla.vector, x_elements);

        iteration++;

        // 勾配の計算
        x_nabla = nabla_function(x);
        num_calls_nabla++;

        // 勾配を確認したいときは on に
        // printf_vector_double(x_nabla, x_elements);
        // printf("max = %f\nmin = %f\n", max_vector_double(x_nabla, x_elements), min_vector_double(x_nabla, x_elements));

        // 収束判定
        double max_vector = max_vector_double(x_nabla.vector, x_elements);
        double min_vector = min_vector_double(x_nabla.vector, x_elements);
        double conv;
        if(max_vector > - min_vector){
            conv = max_vector;
        } else {
            conv = - min_vector;
        }
        if (conv < epsilon)
        {
            break;
        }

        // 目的関数値を計算したいときは on
        if(iteration%100 == 0){
            printf("Iteration:%d\n", iteration);
            printf("Conv = %f \tObj = %f \tLips = %f\n\n", conv, function(x), L);
        }
        // printf("objective = %f", function(x));

        // 次ステップのための慣性方向の移動
        if (inner_product(x_nabla.vector, delta_x.vector, x_elements) > 0.0)
        {
            t_1 = 1.0;
        }
        t_0 = t_1;
        t_1 = (1.0 + sqrt(1.0 + 4.0 * t_0 * t_0)) / 2.0;
        k_times_double(delta_x.vector, (t_0 - 1) / t_1, x_elements);
        copy_vector_double(y.vector, x.vector, x_elements);
        sum_b_to_a_double(y.vector, delta_x.vector, x_elements);

        free(x_nabla.vector);

        // 暫定解を書き出ししたいときは on
        // printf_vector_double(x, x_elements);

    }

    cpu_end_time = clock();
    sec = (double)(cpu_end_time - cpu_start_time) / CLOCKS_PER_SEC;

    // 解の取得
    solution.sol = x;
    // printf_vector_double(x.vector, x_elements);
    solution.iteration = iteration;
    solution.num_calls_nabla = num_calls_nabla;
    solution.num_calls_obj = num_calls_obj;
    solution.CPU_time_sec = sec;



    free(x_nabla.vector);
    free(y_nabla.vector);
    free(y.vector);

    return solution;
}















/*----------------------------------------------------------------------------------------------------------------------
目的関数関係
----------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------
行列のメモリ領域を解放する関数
----------------------------------------------------------------------------------------------------------------------*/
int free_Matrix(Matrix A)
{
    for(int row = 0; row < A.num_row; row++)
    {
        free(A.matrix[row]);
    }
    free(A.matrix);

    return 0;
}


/*----------------------------------------------------------------------------------------------------------------------
r_index(0スタート)，s_index(0スタート)から
rs_index(0スタート)をえる関数
----------------------------------------------------------------------------------------------------------------------*/
int calc_rs_index_from_r_s(const int r, const int s, const int num_s)
{
    int rs = 0;

    rs += r * num_s + ( s - 1 );

    return rs;
}

/*----------------------------------------------------------------------------------------------------------------------
mu_os計算関数
第一引数：暫定解，第二引数：問題構造
返り値：mu_os行列(行がo，列がsで0スタート)
----------------------------------------------------------------------------------------------------------------------*/
Matrix calc_mu_os(const Vector now_sol, const Problem_struct data)
{
    Matrix mu_os;
    mu_os.num_row = data.num_nodes;
    mu_os.num_col = data.num_nodes;
    mu_os.matrix = malloc(sizeof(double *) * data.num_nodes);
    for (int o = 0; o < data.num_nodes; o++)
    {
        mu_os.matrix[o] = malloc(sizeof(double) * data.num_nodes);
        for (int s = 0; s < data.num_nodes; s++)
        {
            double sum = 0.0;

            for (int r = 0; r < data.num_depots; r++)
            {
                int rs = calc_rs_index_from_r_s(r, s, data.num_nodes);
                sum += exp(- data.LOGIT_param_driver * (data.cost_between_nodes.matrix[o][data.depots_index.vector[r]] - now_sol.vector[rs] + data.cost_between_nodes.matrix[data.depots_index.vector[r]][s]));
            }

            mu_os.matrix[o][s] = - log(sum) / data.LOGIT_param_driver;
        }
    }

    return mu_os;
}

/*----------------------------------------------------------------------------------------------------------------------
mu_od計算関数
第一引数：mu_os行列，第二引数：問題構造
返り値：mu_od行列(行がo，列がdで0スタート)
----------------------------------------------------------------------------------------------------------------------*/
Matrix calc_mu_od(const Matrix mu_os, const Problem_struct data)
{
    Matrix mu_od;
    mu_od.num_row = data.num_nodes;
    mu_od.num_col = data.num_nodes;
    mu_od.matrix = malloc(sizeof(double *) * data.num_nodes);
    for (int o = 0; o < data.num_nodes; o++)
    {
        mu_od.matrix[o] = malloc(sizeof(double) * data.num_nodes);
        for (int d = 0; d < data.num_nodes; d++)
        {
            double sum = 0.0;

            for (int s = 0; s < data.num_nodes; s++)
            {
                sum += exp(-data.LOGIT_param_driver * (mu_os.matrix[o][s] + data.cost_between_nodes.matrix[s][d]));
            }

            mu_od.matrix[o][d] = -log(sum) / data.LOGIT_param_driver;
        }
    }

    return mu_od;
}


/*----------------------------------------------------------------------------------------------------------------------
目的関数
第一引数：暫定解，第二引数：問題のデータ
返り値：目的関数値
----------------------------------------------------------------------------------------------------------------------*/
double obj_function(const Vector now_sol, const Problem_struct data)
{

    double obj = 0.0;

    Matrix mu_os;
    mu_os = calc_mu_os(now_sol, data);

    Matrix mu_od;
    mu_od = calc_mu_od(mu_os, data);


    for(int r = 0; r < data.num_depots; r++)
    {
        for (int s = 0; s < data.num_nodes; s++)
        {
            int rs = calc_rs_index_from_r_s(r, s, data.num_nodes);
            obj -= data.num_shippers.matrix[r][s] * now_sol.vector[rs];
        }
        
    }

    for(int o = 0; o < data.num_nodes; o++)
    {
        for(int d = 0; d < data.num_nodes; d++)
        {
            obj -= data.num_drivers.matrix[o][d] * mu_od.matrix[o][d];
        }
    }

    free_Matrix(mu_os);
    free_Matrix(mu_od);


    return obj;
}















/*----------------------------------------------------------------------------------------------------------------------
勾配関数
----------------------------------------------------------------------------------------------------------------------*/
// pのためだけの3次元行列構造体
typedef struct {
    int num_row_1;
    int num_row_2;
    int num_row_3;
    double*** matrix;
} Matrix_3_dim;

/*----------------------------------------------------------------------------------------------------------------------
3次元行列メモリ領域確保の関数
第一引数：対象の3次元行列構造体
----------------------------------------------------------------------------------------------------------------------*/
Matrix_3_dim keep_3_dim_memory_area(Matrix_3_dim matrix)
{

    matrix.matrix = malloc(sizeof(double**) * matrix.num_row_1);
    for(int row_1 = 0; row_1 < matrix.num_row_1; row_1++){
        matrix.matrix[row_1] = malloc(sizeof(double*) * matrix.num_row_2);
        for(int row_2 = 0; row_2 < matrix.num_row_2; row_2++){
            matrix.matrix[row_1][row_2] = malloc(sizeof(double) * matrix.num_row_3);
        }
    }

    return matrix;
}

/*----------------------------------------------------------------------------------------------------------------------
3次元行列のメモリ領域を解放する関数
----------------------------------------------------------------------------------------------------------------------*/
int free_Matrix_3_dim(Matrix_3_dim A)
{
    for (int dim_1 = 0; dim_1 < A.num_row_1; dim_1++)
    {
        for(int dim_2 = 0; dim_2 < A.num_row_2; dim_2++)
        {
            free(A.matrix[dim_1][dim_2]);
        }
        free(A.matrix[dim_1]);
    }
    free(A.matrix);

    return 0;
}

/*----------------------------------------------------------------------------------------------------------------------
p_osd計算関数
第一引数：問題のデータ，第二引数：mu_os，第三引数：mu_od
返り値：p_osd(3次元行列 0スタート)
----------------------------------------------------------------------------------------------------------------------*/
Matrix_3_dim calc_p_osd(const Problem_struct data, const Matrix mu_os, const Matrix mu_od)
{

    Matrix_3_dim p_osd;
    p_osd.num_row_1 = data.num_nodes;
    p_osd.num_row_2 = data.num_nodes;
    p_osd.num_row_3 = data.num_nodes;
    p_osd = keep_3_dim_memory_area(p_osd);

    for(int o = 0; o < p_osd.num_row_1; o++){
        for(int s = 0; s < p_osd.num_row_2; s++){
            for(int d = 0; d < p_osd.num_row_3; d++){
                p_osd.matrix[o][s][d] = exp(-data.LOGIT_param_driver * (mu_os.matrix[o][s] + data.cost_between_nodes.matrix[s][d] - mu_od.matrix[o][d]));
            }
        }
    }

    return p_osd;
}

/*----------------------------------------------------------------------------------------------------------------------
p_ors計算関数
第一引数：問題のデータ，第二引数：暫定解，第三引数：mu_os
返り値：p_ors(3次元行列 0スタート)
----------------------------------------------------------------------------------------------------------------------*/
Matrix_3_dim calc_p_ors(const Problem_struct data, const Vector now_sol, const Matrix mu_os)
{
    Matrix_3_dim p_ors;
    p_ors.num_row_1 = data.num_nodes;
    p_ors.num_row_2 = data.num_depots;
    p_ors.num_row_3 = data.num_nodes;
    p_ors = keep_3_dim_memory_area(p_ors);

    for(int o = 0; o < p_ors.num_row_1; o++){
        for(int r = 0; r < p_ors.num_row_2; r++){
            for(int s = 0; s < p_ors.num_row_3; s++){
                int rs = calc_rs_index_from_r_s(r, s, data.cost_between_nodes.num_row);
                p_ors.matrix[o][r][s] = exp(- data.LOGIT_param_driver * (data.cost_between_nodes.matrix[o][r] + data.cost_between_nodes.matrix[r][s] - now_sol.vector[rs] - mu_os.matrix[o][s]));
            }
        }
    }

    return p_ors;
}

/*----------------------------------------------------------------------------------------------------------------------
X_os計算関数
第一引数：問題のデータ，第二引数：p_osd
返り値：X_os行列
----------------------------------------------------------------------------------------------------------------------*/
Matrix calc_x_os(const Problem_struct data, const Matrix_3_dim p_osd){

    Matrix x_os;
    x_os.num_col = data.num_nodes;
    x_os.num_row = data.num_nodes;
    x_os.matrix = malloc(sizeof(double*) * x_os.num_row);
    for(int o = 0; o < x_os.num_row; o++){
        x_os.matrix[o] = malloc(sizeof(double) * x_os.num_col);
        for(int s = 0; s < x_os.num_col; s++){
            x_os.matrix[o][s] = 0.0;
            for(int d = 0; d < data.num_nodes; d++){
                x_os.matrix[o][s] += data.num_drivers.matrix[o][d] * p_osd.matrix[o][s][d]; 
            }
        }
    }

    return x_os;

}


/*----------------------------------------------------------------------------------------------------------------------
勾配関数
第一引数：暫定解，第二引数：問題のデータ
返り値：勾配ベクトル
----------------------------------------------------------------------------------------------------------------------*/
Vector nabla_function(const Vector now_sol, const Problem_struct data)
{

    Vector nabla;
    nabla.num_elements = now_sol.num_elements;
    nabla.vector = malloc(sizeof(double) * (now_sol.num_elements));

    Matrix mu_os = calc_mu_os(now_sol, data);
    Matrix mu_od = calc_mu_od(mu_os, data);

    Matrix_3_dim p_osd = calc_p_osd(data, mu_os, mu_od);
    Matrix_3_dim p_ors = calc_p_ors(data, now_sol, mu_os);

    Matrix x_os = calc_x_os(data, p_osd);

    for(int r = 0; r < data.num_depots; r++){
        for(int s = 0; s < data.num_nodes; s++){
            int rs = calc_rs_index_from_r_s(r, s, data.num_nodes);
            nabla.vector[rs] = - data.num_shippers.matrix[r][s];
            for(int o = 0; o < data.num_nodes; o++){
                nabla.vector[rs] += p_ors.matrix[o][r][s] * x_os.matrix[o][s];
            }
        }
    }

    free_Matrix(mu_os);
    free_Matrix(mu_od);
    free_Matrix(x_os);

    free_Matrix_3_dim(p_osd);
    free_Matrix_3_dim(p_ors);

    return nabla;
}