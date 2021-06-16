#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

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
        if (max_vector_double(x_nabla.vector, x_elements) < epsilon && min_vector_double(x_nabla.vector, x_elements) > -epsilon)
        {
            break;
        }

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

        // 目的関数値を計算したいときは on
        // printf("objective = %f", function(x));
    }

    // 解の取得
    solution.sol = x;
    // printf_vector_double(x.vector, x_elements);
    solution.iteration = iteration;
    solution.num_calls_nabla = num_calls_nabla;
    solution.num_calls_obj = num_calls_obj;



    free(x_nabla.vector);
    free(y_nabla.vector);
    free(y.vector);

    return solution;
}















/*----------------------------------------------------------------------------------------------------------------------
目的関数関係
----------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------
r_index(0スタート)，s_index(0スタート)から
rs_index(0スタート)をえる関数
----------------------------------------------------------------------------------------------------------------------*/
void free_Matrix(Matrix A)
{
    for(int row = 0; row < A.num_row; row++)
    {
        free(A.matrix[row]);
    }
    free(A.matrix);
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
    mu_os.num_row = data.cost_between_nodes.num_row;
    mu_os.num_col = data.cost_between_nodes.num_row;
    mu_os.matrix = malloc(sizeof(double *) * data.cost_between_nodes.num_row);
    for (int o = 0; o < data.cost_between_nodes.num_row; o++)
    {
        mu_os.matrix[o] = malloc(sizeof(double) * data.cost_between_nodes.num_row);
        for (int s = 0; s < data.cost_between_nodes.num_row; s++)
        {
            double sum = 0.0;

            for (int r = 0; r < data.depots_index.num_elements; r++)
            {
                int rs = calc_rs_index_from_r_s(r, s, data.cost_between_nodes.num_row);
                sum += exp(-data.LOGIT_param_driver * (data.cost_between_nodes.matrix[o][data.depots_index.vector[r]] - now_sol.vector[rs] + data.cost_between_nodes.matrix[data.depots_index.vector[r]][s]));
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
    mu_od.num_row = data.cost_between_nodes.num_row;
    mu_od.num_col = data.cost_between_nodes.num_row;
    mu_od.matrix = malloc(sizeof(double *) * data.cost_between_nodes.num_row);
    for (int o = 0; o < data.cost_between_nodes.num_row; o++)
    {
        mu_od.matrix[o] = malloc(sizeof(double) * data.cost_between_nodes.num_row);
        for (int d = 0; d < data.cost_between_nodes.num_row; d++)
        {
            double sum = 0.0;

            for (int s = 0; s < data.cost_between_nodes.num_row; s++)
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
第一引数：暫定解，第二引数：タスク起点のインデックス行列，
第三引数：ノード間のコスト行列，第四引数；ドライバーLOGITパラメータ
----------------------------------------------------------------------------------------------------------------------*/
double obj_function(const Vector now_sol, const Problem_struct data)
{

    double obj = 0.0;

    Matrix mu_os;
    mu_os = calc_mu_os(now_sol, data);

    Matrix mu_od;
    mu_od = calc_mu_od(mu_os, data);

    for(int r = 0; r < data.depots_index.num_elements; r++)
    {
        for (int s = 0; s < data.cost_between_nodes.num_row; s++)
        {
            int rs = calc_rs_index_from_r_s(r, s, data.cost_between_nodes.num_row);
            obj -= data.num_shippers.matrix[r][s] * now_sol.vector[rs];
        }
        
    }

    for(int o = 0; o < data.cost_between_nodes.num_row; o++)
    {
        for(int d = 0; d < data.cost_between_nodes.num_row; d++)
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


Vector nabla_function(const Vector x)
{

    Vector nabla;
    nabla.vector = malloc(sizeof(double) * x.num_elements);
    nabla.num_elements = x.num_elements;

    zeros_vector_double(nabla.vector, nabla.num_elements);

    for (int i = 0; i < 10; i++)
    {
        // nabla[i] = 2*x[i];
        nabla.vector[i] = 2 * x.vector[i] + 1 + exp(x.vector[i]);
    }

    return nabla;
}