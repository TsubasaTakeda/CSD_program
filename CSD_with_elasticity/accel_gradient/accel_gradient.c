#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <direct.h>

#define LOGIT_DRIVER 5.0
#define LOGIT_SHIPPER 5.0
#define CONVERGENCE 1.0
#define ETA_ACCEL_PARAM 1.1
#define FIRST_LIPS 0.01
#define DISPLAY_ITERATION 100


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

// pのためだけの3次元行列構造体
typedef struct
{
    int num_row_1;
    int num_row_2;
    int num_row_3;
    double ***matrix;
} Matrix_3_dim;

// 加速勾配法の返り値構造体
typedef struct {
    Vector sol;
    Matrix sol_matrix;
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
    double LOGIT_param_shipper;
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
        matrix.matrix[dim_1] = malloc(sizeof(double ) * matrix.num_col);
    }

    return matrix;
}

/*
ベクトルのメモリ領域を解放する関数
*/
int free_Vector(Vector A)
{
    free(A.vector);
    return 0;
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
3次元行列メモリ領域確保の関数
第一引数：対象の3次元行列構造体
*/
Matrix_3_dim get_3_dim_memory_area(Matrix_3_dim matrix)
{

    matrix.matrix = malloc(sizeof(double **) * matrix.num_row_1);
    for (int row_1 = 0; row_1 < matrix.num_row_1; row_1++)
    {
        matrix.matrix[row_1] = malloc(sizeof(double *) * matrix.num_row_2);
        for (int row_2 = 0; row_2 < matrix.num_row_2; row_2++)
        {
            matrix.matrix[row_1][row_2] = malloc(sizeof(double) * matrix.num_row_3);
        }
    }

    return matrix;
}

/*
3次元行列のメモリ領域を解放する関数
*/
int free_Matrix_3_dim(Matrix_3_dim A)
{
    for (int dim_1 = 0; dim_1 < A.num_row_1; dim_1++)
    {
        for (int dim_2 = 0; dim_2 < A.num_row_2; dim_2++)
        {
            free(A.matrix[dim_1][dim_2]);
        }
        free(A.matrix[dim_1]);
    }
    free(A.matrix);

    return 0;
}

char **read_filename_csv(const char *filename, int *numfile)
{

    char chData;

    int num_row, num_col;

    // 行数と列数を読み取る
    if (csv_CountRowCol(filename, &num_row, &num_col))
    {
        // return -1;
    }
    // printf("%d, %d\n", num_row, num_col);
    *numfile = num_row * num_col;

    FILE *fp = fopen(filename, "r");

    // 領域確保
    char **filelist = malloc(sizeof(char *) * num_row * num_col);
    for (int i = 0; i < num_row * num_col; i++)
    {
        filelist[i] = malloc(sizeof(char) * 256);
    }
    // printf("%s\n", filelist[0]);
    // printf("%d", num_row);

    // char s[num_row*num_col][100];

    for (int i = 0; i < num_row * num_col; i++)
    {
        // printf("%d", i);

        fscanf(fp, "%s", filelist[i]);

        // printf("%s", s[i]);

        // filelist[i] = &s[i][0];

        // printf(filelist[i]);
        // printf("\n");
    }
    // printf("end");

    // printf_matrix_double(matrix[0], *num_row, *num_col);

    fclose(fp);
    return filelist;
}

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
double backtracking(const double L, const double eta, const double (*function)(const Vector, const Problem_struct), const Vector (*nabla_function)(const Vector, const Problem_struct), const Vector x, int *num_calls_obj, int *num_calls_nabla, const Problem_struct data)
{

    double F, Q;

    double iota = 0.0;
    Vector y;
    y.vector = malloc(sizeof(double) * x.num_elements);
    y.num_elements = x.num_elements;

    // printf("now_sol:");
    // printf_vector_double(x.vector, x.num_elements);

    double obj = function(x, data);
    // printf("obj = %lf\n", obj);
    num_calls_obj[0]++;

    Vector nabla = nabla_function(x, data);
    // printf("now_nabla:");
    // printf_vector_double(nabla.vector, nabla.num_elements);
    num_calls_nabla[0]++;

    double nolm = nolm_2(nabla.vector, x.num_elements);
    // printf("nolm_2 = %lf\n", nolm);

    // 確認用
    // printf_vector_double(nabla, x_elements);
    // printf("obj = %f\nnolm = %f\n", obj, nolm);

    while (1)
    {

        Q = obj - nolm / (2 * pow(eta, iota) * L);
        // printf("Q = %lf\n", Q);

        // Fの移動先のベクトルを取得
        copy_vector_double(y.vector, nabla.vector, x.num_elements);
        k_times_double(y.vector, -1.0 / (pow(eta, iota) * L), x.num_elements);
        sum_b_to_a_double(y.vector, x.vector, x.num_elements);

        F = function(y, data);
        // printf("F = %lf\n", F);
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
Optimization FISTA_with_restart(const double epsilon, const double eta, double L, const double (*function)(const Vector, const Problem_struct), const Vector (*nabla_function)(const Vector, const Problem_struct), const Vector x_first, const Problem_struct data)
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
        iota = backtracking(L, eta, function, nabla_function, y, &num_calls_obj, &num_calls_nabla, data);
        L = pow(eta, iota) * L;

        // リプシッツ定数(ステップサイズの逆数)の確認
        // printf("Lipsitz = %f\n", L);

        // 暫定解の更新
        free(y_nabla.vector);
        y_nabla = nabla_function(y, data);
        num_calls_nabla++;
        k_times_double(y_nabla.vector, -1.0 / L, x_elements);
        copy_vector_double(x.vector, y.vector, x_elements);
        sum_b_to_a_double(x.vector, y_nabla.vector, x_elements);
        sum_b_to_a_double(delta_x.vector, y_nabla.vector, x_elements);

        iteration++;

        // 勾配の計算
        x_nabla = nabla_function(x, data);
        // printf("now_sol:");
        // printf_vector_double(x.vector, x.num_elements);
        num_calls_nabla++;

            // printf("complete_3\n");

        // 勾配を確認したいときは on に
        // printf_vector_double(x_nabla.vector, x_elements);
        // printf("max = %f\nmin = %f\n", max_vector_double(x_nabla.vector, x_elements), min_vector_double(x_nabla.vector, x_elements));

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
        if(iteration%(DISPLAY_ITERATION) == 0){
            printf("Iteration:%d\n", iteration);
            printf("Conv = %f \tObj = %f \tLips = %f\n", conv, function(x, data), L);
            printf("now_nabla:");
            printf_vector_double(x_nabla.vector, x_nabla.num_elements);
            printf("\n");
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

/*
r_index(0スタート)，s_index(0スタート)から
rs_index(0スタート)をえる関数
*/
int calc_rs_index_from_r_s(const int r, const int s, const int num_s)
{
    int rs = 0;

    rs += r * num_s + s;

    return rs;
}

/*
mu_os計算関数
第一引数：暫定解，第二引数：問題構造
返り値：mu_os行列(行がo，列がsで0スタート)
*/
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
            // printf("mu_os.matrix[%d][%d] = %lf\n", o, s, mu_os.matrix[o][s]);
        }
    }

    // printf_matrix_double(mu_os.matrix, mu_os.num_row, mu_os.num_col);

    return mu_os;
}

/*
mu_od計算関数
第一引数：mu_os行列，第二引数：問題構造
返り値：mu_od行列(行がo，列がdで0スタート)
*/
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
                // printf("%lf\n", mu_os.matrix[o][s] + data.cost_between_nodes.matrix[s][d]);
                // printf("mu_os = \n");
                // printf_matrix_double(mu_os.matrix, mu_os.num_row, mu_os.num_col);
                // printf("cost = \n");
                // printf_matrix_double(data.cost_between_nodes.matrix, data.cost_between_nodes.num_row, data.cost_between_nodes.num_col);
                // printf("sum = %lf\n", sum);
            }

            sum += exp(-data.LOGIT_param_driver * data.cost_between_nodes.matrix[o][d]);

            mu_od.matrix[o][d] = -log(sum) / data.LOGIT_param_driver;
            // printf("%lf\t%lf\n", mu_od.matrix[o][d], mu_os.matrix[o][d]);
        }
    }
    // printf_matrix_double(mu_od.matrix, mu_od.num_row, mu_od.num_col);

    return mu_od;
}


Matrix calc_V_rs(const Vector now_sol, const Problem_struct data)
{
    Matrix V_rs;
    V_rs.num_row = data.num_depots;
    V_rs.num_col = data.num_nodes;
    V_rs.matrix = malloc(sizeof(double * ) * V_rs.num_row);

    for(int r = 0; r < data.num_depots; r++)
    {
        V_rs.matrix[r] = malloc(sizeof(double) * V_rs.num_col);

        for(int s = 0; s < data.num_nodes; s++)
        {

            double sum = 0.0;

            sum += exp(-data.LOGIT_param_shipper * ( data.cost_between_nodes.matrix[data.depots_index.vector[r]][s] ));
            int rs = calc_rs_index_from_r_s(r, s, data.num_nodes);
            sum += exp(-data.LOGIT_param_shipper * ( now_sol.vector[rs] ));

            V_rs.matrix[r][s] = -log(sum) / data.LOGIT_param_shipper;

        }

    }

    return V_rs;

}


/*
目的関数
第一引数：暫定解，第二引数：問題のデータ
返り値：目的関数値
*/
double obj_function(const Vector now_sol, const Problem_struct data)
{

    double obj = 0.0;

    Matrix mu_os;
    mu_os = calc_mu_os(now_sol, data);

    Matrix mu_od;
    mu_od = calc_mu_od(mu_os, data);

    Matrix V_rs;
    V_rs = calc_V_rs(now_sol, data);


    for(int r = 0; r < data.num_depots; r++)
    {
        for (int s = 0; s < data.num_nodes; s++)
        {
            obj -= data.num_shippers.matrix[r][s] * V_rs.matrix[r][s];
        }
        
    }

    for(int o = 0; o < data.num_nodes; o++)
    {
        for(int d = 0; d < data.num_nodes; d++)
        {
            obj -= data.num_drivers.matrix[o][d] * mu_od.matrix[o][d];
        }
    }
    // printf_matrix_double(mu_od.matrix, mu_od.num_row, mu_od.num_col);

    free_Matrix(mu_os);
    free_Matrix(mu_od);
    free_Matrix(V_rs);


    return obj;
}















/*----------------------------------------------------------------------------------------------------------------------
勾配関数
----------------------------------------------------------------------------------------------------------------------*/

/*
p_osd計算関数
第一引数：問題のデータ，第二引数：mu_os，第三引数：mu_od
返り値：p_osd(3次元行列 0スタート)
*/
Matrix_3_dim calc_p_osd(const Problem_struct data, const Matrix mu_os, const Matrix mu_od)
{

    Matrix_3_dim p_osd;
    p_osd.num_row_1 = data.num_nodes;
    p_osd.num_row_2 = data.num_nodes;
    p_osd.num_row_3 = data.num_nodes;
    p_osd = get_3_dim_memory_area(p_osd);

    for(int o = 0; o < p_osd.num_row_1; o++){
        for(int s = 0; s < p_osd.num_row_2; s++){
            for(int d = 0; d < p_osd.num_row_3; d++){
                p_osd.matrix[o][s][d] = exp(-data.LOGIT_param_driver * (mu_os.matrix[o][s] + data.cost_between_nodes.matrix[s][d] - mu_od.matrix[o][d]));
            }
        }
    }

    return p_osd;
}

/*
p_ors計算関数
第一引数：問題のデータ，第二引数：暫定解，第三引数：mu_os
返り値：p_ors(3次元行列 0スタート)
*/
Matrix_3_dim calc_p_ors(const Problem_struct data, const Vector now_sol, const Matrix mu_os)
{
    Matrix_3_dim p_ors;
    p_ors.num_row_1 = data.num_nodes;
    p_ors.num_row_2 = data.num_depots;
    p_ors.num_row_3 = data.num_nodes;
    p_ors = get_3_dim_memory_area(p_ors);
    // printf_vector_double(now_sol.vector, now_sol.num_elements);

    for(int o = 0; o < p_ors.num_row_1; o++){
        for(int r = 0; r < p_ors.num_row_2; r++){
            // printf("start\n");
            for(int s = 0; s < p_ors.num_row_3; s++){
                int rs = calc_rs_index_from_r_s(r, s, data.num_nodes);
                // printf("num_nodes = %d\n", data.num_nodes);
                // printf("num_depots = %d\n", data.num_depots);
                // printf("rs = %d\n", rs);
                // printf("p_ors[%d][%d][%d] = %lf\n", o, r, s, data.cost_between_nodes.matrix[o][data.depots_index.vector[r]]);
                p_ors.matrix[o][r][s] = exp( - data.LOGIT_param_driver * (data.cost_between_nodes.matrix[o][data.depots_index.vector[r]] + data.cost_between_nodes.matrix[data.depots_index.vector[r]][s] - now_sol.vector[rs] - mu_os.matrix[o][s]));
            }
            // printf("complete\n");
        }
    }

    return p_ors;
}

/*
X_os計算関数
第一引数：問題のデータ，第二引数：p_osd
返り値：X_os行列
*/
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

/*
荷主のLOGIT配分
*/
Matrix calc_z_rs(const Vector now_sol, const Problem_struct data){

    Matrix z_rs;
    z_rs.num_row = data.num_depots;
    z_rs.num_col = data.num_nodes;
    z_rs.matrix = malloc(sizeof(double*) * z_rs.num_row );
    for(int r = 0; r < data.num_depots; r++){
        z_rs.matrix[r] = malloc(sizeof(double) * z_rs.num_col);
        for(int s = 0; s < data.num_nodes; s++){
            int rs = calc_rs_index_from_r_s(r, s, data.num_nodes);
            z_rs.matrix[r][s] = data.num_shippers.matrix[r][s] / (1.0 + exp( - data.LOGIT_param_shipper * (data.cost_between_nodes.matrix[data.depots_index.vector[r]][s] - now_sol.vector[rs])));
        }
    }

    return z_rs;

}


/*
勾配関数
第一引数：暫定解，第二引数：問題のデータ
返り値：勾配ベクトル
*/
Vector nabla_function(const Vector now_sol, const Problem_struct data)
{

    Vector nabla;
    nabla.num_elements = now_sol.num_elements;
    nabla.vector = malloc(sizeof(double) * (now_sol.num_elements));

    Matrix mu_os = calc_mu_os(now_sol, data);
    // printf_matrix_double(mu_os.matrix, mu_os.num_row, mu_os.num_col);
    Matrix mu_od = calc_mu_od(mu_os, data);
    // printf_matrix_double(mu_od.matrix, mu_od.num_row, mu_od.num_col);

    Matrix_3_dim p_osd = calc_p_osd(data, mu_os, mu_od);
    Matrix_3_dim p_ors = calc_p_ors(data, now_sol, mu_os);

    Matrix x_os = calc_x_os(data, p_osd);

    

    Matrix z_rs = calc_z_rs(now_sol, data);

    

    for(int r = 0; r < data.num_depots; r++)
    {
        for(int s = 0; s < data.num_nodes; s++){
            int rs = calc_rs_index_from_r_s(r, s, data.num_nodes);
            nabla.vector[rs] = - z_rs.matrix[r][s];
            for(int o = 0; o < data.num_nodes; o++){
                nabla.vector[rs] +=  p_ors.matrix[o][r][s] * x_os.matrix[o][s];
            }
        }
        // printf("complete_mid\n");
    }


    free_Matrix(mu_os);
    free_Matrix(mu_od);
    free_Matrix(x_os);
    free_Matrix(z_rs);

    free_Matrix_3_dim(p_osd);
    free_Matrix_3_dim(p_ors);

    // printf_vector_double(nabla.vector, nabla.num_elements);


    return nabla;
}







/*----------------------------------------------------------------------------------------------------------------------
ラスト配分関数
----------------------------------------------------------------------------------------------------------------------*/

/*
ラスト配分関数(ドライバー)
第一引数：問題データ，第二引数：p_osd，第三引数：p_ors
返り値：タスク配分行列(od(num_nodes*num_nodes) * rs(num_depots*num_nodes))
*/
Matrix allocate_task_to_driver(Problem_struct data, Vector sol)
{

    Matrix allocate;
    allocate.num_row = data.num_nodes*data.num_nodes;
    allocate.num_col = data.num_depots*data.num_nodes + 1;
    allocate = get_matrix_memory_area(allocate);

    Matrix mu_os = calc_mu_os(sol, data);
    Matrix mu_od = calc_mu_od(mu_os, data);

    Matrix_3_dim p_osd = calc_p_osd(data, mu_os, mu_od);
    Matrix_3_dim p_ors = calc_p_ors(data, sol, mu_os);


    for(int o = 0; o < data.num_nodes; o++)
    {
        for(int d = 0; d < data.num_nodes; d++)
        {
            int od = calc_rs_index_from_r_s(o, d, data.num_nodes);

            allocate.matrix[od][data.num_depots * data.num_nodes] = data.num_drivers.matrix[o][d];

            // printf("od = %d\n", od);
            for(int r = 0; r < data.num_depots; r++)
            {
                for(int s = 0; s < data.num_nodes; s++)
                {
                    int rs = calc_rs_index_from_r_s(r, s, data.num_nodes);
                    // printf("s = %d", s);
                    // printf("rs = %d\n", rs);
                    // printf("num_driver = %d,", data.num_drivers.matrix[o][d]);
                    // printf("\tp_osd = %lf,", p_osd.matrix[o][s][0]);
                    // printf("\tp_ors = %lf\n", p_ors.matrix[o][r][s]);
                    // printf_matrix_double(p_osd.matrix[o], p_osd.num_row_2, p_osd.num_row_3);
                    allocate.matrix[od][rs] = (double)(data.num_drivers.matrix[o][d]) * p_osd.matrix[o][s][d] * p_ors.matrix[o][r][s];
                    allocate.matrix[od][data.num_depots*data.num_nodes] -= allocate.matrix[od][rs];
                    // printf_matrix_double(p_osd.matrix[0], p_osd.num_row_2, p_osd.num_row_3);
                    // printf("allocate[%d][%d][%d][%d] = %lf\n", o, d, r, s, allocate.matrix[od][rs]); 
                }
            }


        }
    }

    free_Matrix(mu_os);
    free_Matrix(mu_od);
    free_Matrix_3_dim(p_osd);
    free_Matrix_3_dim(p_ors);

    return allocate;

}

/*
ラスト配分関数(荷主)
第一引数：問題データ，第二引数：p_osd，第三引数：p_ors
返り値：タスク配分行列(od(num_nodes*num_nodes) * rs(num_depots*num_nodes))
*/
Matrix allocate_task(Problem_struct data, Vector sol)
{

    Matrix allocate_task;
    allocate_task.num_row = data.num_depots * data.num_nodes;
    allocate_task.num_col = 2;

    allocate_task = get_matrix_memory_area(allocate_task);


    Matrix num_task = calc_z_rs(sol, data);

    for(int r = 0; r < data.num_depots; r++){

        for(int s = 0; s < data.num_nodes; s++){

            int rs = calc_rs_index_from_r_s(r, s, data.num_nodes);

            allocate_task.matrix[rs][1] = num_task.matrix[r][s];
            allocate_task.matrix[rs][0] = (double)data.num_shippers.matrix[r][s] - allocate_task.matrix[rs][1];

        }

    }

    free_Matrix(num_task);
    
    return allocate_task;
}











/*----------------------------------------------------------------------------------------------------------------------
メイン関数
----------------------------------------------------------------------------------------------------------------------*/
int main(){

    char *filename = malloc(sizeof(char) * 256);
    filename = "..//filelist.csv";
    int numfile;
    char **filelist;
    filelist = read_filename_csv(filename, &numfile);

    for (int i = 0; i < numfile; i++){

        for (int exp = 0; exp < 10; exp++){

            Problem_struct data;

            char buff[256];
            char file[256];

            strcpy(buff, filelist[i]);

            strcat(buff, "exp_");
            // // printf("end");

            char exp_str[10];
            itoa(exp, exp_str, 10);
            // printf(exp_str);
            strcat(buff, exp_str);



            // コスト行列読み込み (num_node * num_node)
            strcpy(file, buff);
            strcat(file, "\\input\\cost.csv");
            read_matrix_csv(file, &data.cost_between_nodes.matrix, &data.cost_between_nodes.num_row, &data.cost_between_nodes.num_col);
            // printf_matrix_double(data.cost_between_nodes.matrix, data.cost_between_nodes.num_row, data.cost_between_nodes.num_col);

            // 荷主数行列読み込み(num_node * num_node)
            strcpy(file, buff);
            strcat(file, "\\input\\num_shippers.csv");
            read_int_matrix_csv(file, &data.num_shippers.matrix, &data.num_shippers.num_row, &data.num_shippers.num_col);

            // ドライバー数行列読み込み(num_node * num_node)
            strcpy(file, buff);
            strcat(file, "\\input\\num_drivers.csv");
            read_int_matrix_csv(file, &data.num_drivers.matrix, &data.num_drivers.num_row, &data.num_drivers.num_col);

            // タスク起点インデックス読み込み(num_depots)
            strcpy(file, buff);
            strcat(file, "\\input\\depots_index.csv");
            read_int_vector_csv(file, &data.depots_index.vector, &data.depots_index.num_elements);

            printf("complete read files.\n\n");

            int sum_driver = 0;
            int sum_shipper = 0;


            // タスク起点インデックスを(0スタート)から(1スタート)に変換
            for (int r = 0; r < data.depots_index.num_elements; r++)
            {
                data.depots_index.vector[r] += -1;
            }

            // LOGITパラメータを設定
            data.LOGIT_param_driver = LOGIT_DRIVER;
            data.LOGIT_param_shipper = LOGIT_SHIPPER;
            data.num_nodes = data.cost_between_nodes.num_row;
            data.num_depots = data.depots_index.num_elements;

            double (*obj)(const Vector, const Problem_struct);
            Vector (*nabla)(const Vector, const Problem_struct);

            obj = obj_function;
            nabla = nabla_function;

            Vector sol_first;
            sol_first.num_elements = data.num_depots * data.num_nodes;
            sol_first.vector = malloc(sizeof(double) * sol_first.num_elements);
            zeros_vector_double(sol_first.vector, sol_first.num_elements);

            // printf("complete_2\n");

            // 勾配関数の動作確認 -------------------------------------------------------------------------------------
            // Vector nabla_first = nabla_function(sol_first, data);
            // printf_vector_double(nabla_first.vector, nabla_first.num_elements);
            // -------------------------------------------------------------------------------------------------------

            // 目的関数の動作確認 -------------------------------------------------------------------------------------
            // double obj_first = obj_function(sol_first, data);
            // printf("obj = %lf\n", obj_first);
            // -------------------------------------------------------------------------------------------------------

            // printf("complete_mid\n");

            Optimization solution;

            solution = FISTA_with_restart(CONVERGENCE, ETA_ACCEL_PARAM, FIRST_LIPS, obj, nabla, sol_first, data);
            printf("complete FISTA\n");

            printf("Iteration = %d\n\nSolution:", solution.iteration);
            // printf_vector_double(solution.sol.vector, solution.sol.num_elements);

            solution.sol_matrix.num_row = data.num_depots;
            solution.sol_matrix.num_col = data.num_nodes;
            solution.sol_matrix.matrix = malloc(sizeof(double *) * solution.sol_matrix.num_row);
            for (int dim_1 = 0; dim_1 < solution.sol_matrix.num_row; dim_1++)
            {
                solution.sol_matrix.matrix[dim_1] = malloc(sizeof(double) * solution.sol_matrix.num_col);
                for (int dim_2 = 0; dim_2 < solution.sol_matrix.num_col; dim_2++)
                {
                    solution.sol_matrix.matrix[dim_1][dim_2] = solution.sol.vector[dim_1 * solution.sol_matrix.num_col + dim_2];
                }
            }
            printf_matrix_double(solution.sol_matrix.matrix, solution.sol_matrix.num_row, solution.sol_matrix.num_col);

            printf("time = %lf [s]\n", solution.CPU_time_sec);
            printf("num of calls obj_function = %d\n", solution.num_calls_obj);
            printf("num of calls nabla_function = %d\n", solution.num_calls_nabla);

            strcpy(file, buff);
            strcat(file, "\\output");
            mkdir(file);

            strcpy(file, buff);
            strcat(file, "\\output\\sol_price.csv");
            write_matrix_csv(file, solution.sol_matrix.matrix, solution.sol_matrix.num_row, solution.sol_matrix.num_col);

            strcpy(file, buff);
            strcat(file, "\\time");
            mkdir(file);

            strcpy(file, buff);
            strcat(file, "\\time\\accel_CPU_time_sec.csv");
            write_vector_csv(file, &solution.CPU_time_sec, 1);

            strcpy(file, buff);
            strcat(file, "\\time\\accel_num_call_obj.csv");
            write_int_vector_csv(file, &solution.num_calls_obj, 1);

            strcpy(file, buff);
            strcat(file, "\\time\\accel_num_call_nabla.csv");
            write_int_vector_csv(file, &solution.num_calls_nabla, 1);

            strcpy(file, buff);
            strcat(file, "\\time\\accel_iteration.csv");
            write_int_vector_csv(file, &solution.iteration, 1);

            printf("\n\n");

            Matrix task_allocation = allocate_task_to_driver(data, solution.sol);
            Matrix num_task = allocate_task(data, solution.sol);
            // printf("task_allocation:");
            // printf_matrix_double(task_allocation.matrix, task_allocation.num_row, task_allocation.num_col);
            strcpy(file, buff);
            strcat(file, "\\middle");
            mkdir(file);

            strcpy(file, buff);
            strcat(file, "\\middle\\task_allocation_to_driver.csv");
            write_matrix_csv(file, task_allocation.matrix, task_allocation.num_row, task_allocation.num_col);
            
            strcpy(file, buff);
            strcat(file, "\\middle\\num_tasks.csv");
            write_matrix_csv(file, num_task.matrix, num_task.num_row, num_task.num_col);

            free(sol_first.vector);
            free(solution.sol.vector);
            free_Matrix(solution.sol_matrix);
            free_Matrix(task_allocation);
            free_Matrix(num_task);

            printf("\n\ncomplete!!!!!!!!!!!!!!!!!!!!!!!!!\n");
            printf("%s___exp_%d", filelist[i], exp);
            printf("\n\n\n\n");
        }
    }
    free(filename);

    for (int i = 0; i < numfile; i++){
        printf("%s\n", filelist[i]);
        free(filelist[i]);
    }
    free(filelist);

    return 0;

}