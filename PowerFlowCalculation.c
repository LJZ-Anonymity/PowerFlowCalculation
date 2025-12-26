// #include <stdbool.h> 减少调用标准库，用int代替bool
#include <windows.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

 /*====================== 基本数据结构定义 ======================*/

const double PI = 3.14159265358979323846;

typedef enum
{
    NODE_PQ,
    NODE_PV,
    NODE_SLACK
} NodeType; // 节点类型

typedef struct
{
    int from_bus; // 起始节点
    int to_bus;   // 终止节点
    double series_R; // 电阻
    double series_X; // 电抗
    double shunt_B; // 并联电纳
} LineArg; // 线路参数

typedef struct
{
    int node_num;
    NodeType type;
    double P;   /* 节点净注入有功（标幺），PG-PL */
    double Q;   /* 节点净注入无功（标幺），QG-QL */
    double Pg;  /* 节点发电有功（标幺） */
    double Qg;  /* 节点发电无功（标幺） */
    double Pl;  /* 节点负荷有功（标幺） */
    double Ql;  /* 节点负荷无功（标幺） */
    double Gs;  /* 节点并联电导（标幺） */
    double Bs;  /* 节点并联电纳（标幺） */
    double V;   /* 给定电压幅值（PV/平衡节点用） */
    double Theta;
} NodeArg; // 节点参数

typedef struct
{
    int node_num; /* 节点编号 */
    double e;     /* 节点电压幅值 */
    double f;     /* 节点电压相角 */
} InitVal; // 初始值

typedef struct
{
    LineArg* lines; /* 线路参数 */
    size_t line_count; /* 线路参数数量 */

    NodeArg* nodes; /* 节点参数 */
    size_t node_count; /* 节点参数数量 */

    InitVal* init_values; /* 初始值 */
    size_t init_count; /* 初始值数量 */

    int order;      /* 节点总数 */

    double* G;      /* 节点导纳矩阵实部 */
    double* B;      /* 节点导纳矩阵虚部 */

    double* e;      /* 节点电压实部 */
    double* f;      /* 节点电压虚部 */
} NetworkInfo; // 网络信息

/* 稀疏矩阵结构（CSR格式：Compressed Sparse Row） */
typedef struct
{
    int nrows;          /* 矩阵行数 */
    int ncols;          /* 矩阵列数 */
    int nnz;            /* 非零元素个数 */
    int* row_ptr;       /* 行指针数组，长度为 nrows+1 */
    int* col_idx;       /* 列索引数组，长度为 nnz */
    double* values;     /* 非零元素值数组，长度为 nnz */
} SparseMatrix;

typedef struct
{
    NetworkInfo* info; /* 网络信息 */

    int slack_index; /* 平衡节点索引 */
    int pq_count;    /* PQ节点个数 */
    int pv_count;    /* PV节点个数 */
    int non_slack;  /* 非平衡节点个数 */

    int eq_count;   /* 方程个数 = 2*non_slack（PQ: ΔP,ΔQ；PV: ΔP,Δ|U|） */

    int* var_indices;   /* 非平衡节点在网络中的索引顺序：先 PQ 再 PV */
    double* delta;      /* 常数项向量 */
    SparseMatrix J_sparse; /* 稀疏雅可比矩阵 */

    double tolerance; /* 收敛容差 */
    int max_iter;    /* 最大迭代次数 */
} NLIteration; // 牛顿-拉夫逊迭代

/*====================== 工具函数 ======================*/
/// <summary>
/// 确保分配内存成功
/// </summary>
/// <param name="ptr">要分配内存的指针</param>
/// <param name="label">内存标签</param>
static void ensure_alloc(void* ptr, const char* label)
{
    if (!ptr)
    {
        fprintf(stderr, "%s allocation failed\n", label); /* 分配失败 */
        exit(EXIT_FAILURE);
    }
}

/*====================== 稀疏矩阵操作 ======================*/
/// <summary>
/// 初始化稀疏矩阵
/// </summary>
static void sparse_matrix_init(SparseMatrix* mat, int nrows, int ncols)
{
    if (!mat) return;
    memset(mat, 0, sizeof(*mat));
    mat->nrows = nrows;
    mat->ncols = ncols;
    mat->nnz = 0;
    mat->row_ptr = (int*)calloc((size_t)(nrows + 1), sizeof(int));
    ensure_alloc(mat->row_ptr, "sparse row_ptr");
}

/// <summary>
/// 释放稀疏矩阵
/// </summary>
static void sparse_matrix_free(SparseMatrix* mat)
{
    if (!mat) return;
    if (mat->row_ptr) free(mat->row_ptr);
    if (mat->col_idx) free(mat->col_idx);
    if (mat->values) free(mat->values);
    memset(mat, 0, sizeof(*mat));
}

/// <summary>
/// 从三元组格式构建CSR稀疏矩阵（先收集所有非零元素，再转换为CSR）
/// </summary>
/// <param name="mat">输出稀疏矩阵</param>
/// <param name="triplets">三元组数组，每个元素为 (row, col, value)</param>
/// <param name="nnz">非零元素个数</param>
static void sparse_matrix_from_triplets(SparseMatrix* mat, 
    const int* row_idx, const int* col_idx, const double* values, int nnz)
{
    int i, j;
    int* row_counts = NULL;

    if (!mat || !row_idx || !col_idx || !values || nnz <= 0) return;

    /* 统计每行的非零元素个数 */
    row_counts = (int*)calloc((size_t)mat->nrows, sizeof(int));
    ensure_alloc(row_counts, "row_counts");
    
    for (i = 0; i < nnz; ++i)
    {
        if (row_idx[i] >= 0 && row_idx[i] < mat->nrows &&
            col_idx[i] >= 0 && col_idx[i] < mat->ncols)
        {
            row_counts[row_idx[i]]++;
        }
    }

    /* 构建row_ptr */
    mat->row_ptr[0] = 0;
    for (i = 0; i < mat->nrows; ++i)
    {
        mat->row_ptr[i + 1] = mat->row_ptr[i] + row_counts[i];
    }
    mat->nnz = mat->row_ptr[mat->nrows];

    /* 分配col_idx和values */
    mat->col_idx = (int*)calloc((size_t)mat->nnz, sizeof(int));
    mat->values = (double*)calloc((size_t)mat->nnz, sizeof(double));
    ensure_alloc(mat->col_idx, "sparse col_idx");
    ensure_alloc(mat->values, "sparse values");

    /* 重置row_counts作为每行的当前位置计数器 */
    memset(row_counts, 0, (size_t)mat->nrows * sizeof(int));

    /* 填充col_idx和values */
    for (i = 0; i < nnz; ++i)
    {
        int r = row_idx[i];
        int c = col_idx[i];
        if (r >= 0 && r < mat->nrows && c >= 0 && c < mat->ncols)
        {
            int pos = mat->row_ptr[r] + row_counts[r];
            mat->col_idx[pos] = c;
            mat->values[pos] = values[i];
            row_counts[r]++;
        }
    }

    /* 对每行的列索引进行排序（可选，但有助于数值稳定性） */
    if (mat->col_idx && mat->values && mat->nnz > 0)
    {
        /* 保存nnz值到局部变量，帮助静态分析工具理解边界 */
        const int nnz = mat->nnz;
        for (i = 0; i < mat->nrows; ++i)
        {
            int start = mat->row_ptr[i];
            int end = mat->row_ptr[i + 1];
            /* 严格的边界保护，防止静态分析认为可能越界 */
            if (start < 0 || start >= nnz)
                continue; /* 跳过无效起始位置 */
            if (end < start || end > nnz)
                continue; /* 跳过无效结束位置 */
            if (start >= end)
                continue; /* 跳过空行 */

            /* 简单的冒泡排序（每行元素通常很少） */
            /* 此时 start 在 [0, nnz) 范围内，end 在 (start, nnz] 范围内 */
            /* 在循环条件中明确限制 j 和 k 的范围，确保不超过 nnz */
            for (j = start; j < end - 1 && j < nnz - 1; ++j)
            {
                int k;
                /* 在循环条件中明确限制 k 的范围，确保不超过 nnz */
                for (k = j + 1; k < end && k < nnz; ++k)
                {
                    /* 此时 j 和 k 都在 [0, nnz) 范围内，可以安全访问数组 */
                    if (mat->col_idx[j] > mat->col_idx[k])
                    {
                        int tmp_col = mat->col_idx[j];
                        double tmp_val = mat->values[j];
                        mat->col_idx[j] = mat->col_idx[k];
                        mat->values[j] = mat->values[k];
                        mat->col_idx[k] = tmp_col;
                        mat->values[k] = tmp_val;
                    }
                }
            }
        }
    }
    if (row_counts) free(row_counts);
}

/// <summary>
/// 获取稀疏矩阵指定位置的元素值
/// </summary>
static double sparse_matrix_get(const SparseMatrix* mat, int row, int col)
{
    int i;
    if (!mat || row < 0 || row >= mat->nrows || col < 0 || col >= mat->ncols)
        return 0.0;
    
    for (i = mat->row_ptr[row]; i < mat->row_ptr[row + 1]; ++i)
    {
        if (mat->col_idx[i] == col)
            return mat->values[i];
    }
    return 0.0;
}

/// <summary>
/// 获取矩阵中指定位置的元素值
/// </summary>
/// <param name="mat">矩阵</param>
/// <param name="order">矩阵阶数</param>
/// <param name="row">行索引</param>
/// <param name="col">列索引</param>
/// <returns>矩阵中指定位置的元素值</returns>
static inline double mat_get(const double* mat, int order, int row, int col)
{
    return mat[row * order + col];
}

/// <summary>
/// 设置矩阵中指定位置的元素值
/// </summary>
/// <param name="mat">矩阵</param>
/// <param name="order">矩阵阶数</param>
/// <param name="row">行索引</param>
/// <param name="col">列索引</param>
/// <param name="value">要设置的值</returns>
static inline void mat_set(double* mat, int order, int row, int col, double value)
{
    mat[row * order + col] = value;
}

/*====================== 文本输入读取（caseX.txt） ======================*/

/// <summary>
/// 根据节点编号查找其在数组中的索引
/// </summary>
/// <param name="nodes">节点数组</param>
/// <param name="count">节点数量</param>
/// <param name="node_num">节点编号</param>
/// <returns>对应索引，未找到返回 -1</returns>
static int find_node_index_by_num(const NodeArg* nodes, size_t count, int node_num)
{
    size_t idx;

    if (!nodes) return -1;
    for (idx = 0; idx < count; ++idx)
    {
        if (nodes[idx].node_num == node_num)
            return (int)idx;
    }
    return -1;
}

/// <summary>
/// 从 m 文件读取网络数据，格式参考 Input/case4.m
/// </summary>
/// <param name="filename">文件名</param>
/// <param name="out_lines">输出线路参数</param>
/// <param name="out_line_count">输出线路参数数量</param>
/// <param name="out_nodes">输出节点参数</param>
/// <param name="out_node_count">输出节点参数数量</param>
/// <param name="out_inits">输出初始值</param>
/// <returns>0: 成功, -1: 失败</returns>
static int read_case_m(const char* filename,
    LineArg** out_lines, size_t* out_line_count,
    NodeArg** out_nodes, size_t* out_node_count,
    InitVal** out_inits, size_t* out_init_count,
    double* baseMVA)
{
    FILE* fp;
    char buf[512];
    int nBus = 0, slackBus = 1;
    double Us = 1.0, eps = 1e-4;
    int i = 0;

    if (!filename || !out_lines || !out_line_count || !out_nodes || !out_node_count || !out_inits || !out_init_count)
        return -1;

    errno_t err = fopen_s(&fp, filename, "r");
    if (err != 0 || !fp)
    {
        fprintf(stderr, "Failed to open input file: %s\n", filename);
        return -1;
    }

    // 读取 baseMVA
    *baseMVA = 100.0;  // 默认值
    while (fgets(buf, sizeof(buf), fp))
    {
        if (strstr(buf, "mpc.baseMVA") != NULL)
        {
            if (sscanf_s(buf, "mpc.baseMVA = %lf", baseMVA) != 1)
            {
                *baseMVA = 100.0;
            }
            break;
        }
    }

    // 读取 bus 数据
    while (fgets(buf, sizeof(buf), fp))
    {
        if (strstr(buf, "mpc.bus = [") != NULL)
        {
            // 读取节点数据
            int bus_count = 0;
            while (fgets(buf, sizeof(buf), fp))
            {
                if (strstr(buf, "];") != NULL)
                    break;
                int bus_i = 0;
                int type = 0;
                double Pd = 0.0;
                double Qd = 0.0;
                if (sscanf_s(buf, "%d %d %lf %lf", &bus_i, &type, &Pd, &Qd) == 4)
                    bus_count++;
            }
            nBus = bus_count;
            break;
        }
    }

    // 读取 branch 数据
    size_t line_count = 0;
    while (fgets(buf, sizeof(buf), fp))
    {
        if (strstr(buf, "mpc.branch = [") != NULL)
        {
            // 读取线路数据
            while (fgets(buf, sizeof(buf), fp))
            {
                if (strstr(buf, "];") != NULL)
                    break;
                int fbus = 0;
                int tbus = 0;
                double r = 0.0;
                double x = 0.0;
                double b_half = 0.0; // 线路充电电纳 b/2
                if (sscanf_s(buf, "%d %d %lf %lf %lf", &fbus, &tbus, &r, &x, &b_half) >= 5)
                    line_count++;
            }
            break;
        }
    }

    // 分配内存
    *out_lines = (LineArg*)calloc(line_count, sizeof(LineArg));
    *out_nodes = (NodeArg*)calloc(nBus, sizeof(NodeArg));
    *out_inits = (InitVal*)calloc(nBus, sizeof(InitVal));
    ensure_alloc(*out_lines, "lines");
    ensure_alloc(*out_nodes, "nodes");
    ensure_alloc(*out_inits, "inits");

    rewind(fp); // 重新读取文件，填充具体数据

    // 读取节点数据
    while (fgets(buf, sizeof(buf), fp))
    {
        if (strstr(buf, "mpc.bus = [") != NULL)
        {
            int idx = 0;
            while (fgets(buf, sizeof(buf), fp))
            {
                if (strstr(buf, "];") != NULL)
                    break;
                int bus_i = 0, type = 0;
                double Pd = 0.0, Qd = 0.0;
                double Gs = 0.0, Bs = 0.0, area = 0.0;
                double Vm = 1.0, Va = 0.0;
                double baseKV = 0.0, zone = 0.0, Vmax = 1.1, Vmin = 0.9;
                int matched = sscanf_s(buf,
                    "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                    &bus_i, &type, &Pd, &Qd, &Gs, &Bs, &area,
                    &Vm, &Va, &baseKV, &zone, &Vmax, &Vmin);
                if (matched >= 9)
                {
                    double theta_rad = Va * PI / 180.0;
                    double e0 = Vm * cos(theta_rad);
                    double f0 = Vm * sin(theta_rad);

                    (*out_nodes)[idx].node_num = bus_i;
                    (*out_nodes)[idx].type = (type == 3) ? NODE_SLACK :
                        (type == 2) ? NODE_PV : NODE_PQ;
                    (*out_nodes)[idx].Pl = Pd / *baseMVA;
                    (*out_nodes)[idx].Ql = Qd / *baseMVA;
                    (*out_nodes)[idx].P = (*out_nodes)[idx].Pg - (*out_nodes)[idx].Pl;
                    (*out_nodes)[idx].Q = (*out_nodes)[idx].Qg - (*out_nodes)[idx].Ql;

                    (*out_nodes)[idx].Gs = Gs / *baseMVA;
                    (*out_nodes)[idx].Bs = Bs / *baseMVA;

                    (*out_nodes)[idx].V = Vm;
                    (*out_nodes)[idx].Theta = theta_rad;
                    (*out_inits)[idx].node_num = bus_i;
                    (*out_inits)[idx].e = e0;
                    (*out_inits)[idx].f = f0;
                    idx++;
                }
            }
            break;
        }
    }

    // 读取发电机数据，累加至对应节点
    while (fgets(buf, sizeof(buf), fp))
    {
        if (strstr(buf, "mpc.gen = [") != NULL)
        {
            while (fgets(buf, sizeof(buf), fp))
            {
                if (strstr(buf, "];") != NULL)
                    break;
                int gen_bus = 0;
                double Pg = 0.0, Qg = 0.0;
                double Qmax = 0.0, Qmin = 0.0, Vg = 0.0;
                double mBase = 0.0;
                int status = 0;
                if (sscanf_s(buf, "%d %lf %lf %lf %lf %lf %lf %d",
                    &gen_bus, &Pg, &Qg, &Qmax, &Qmin, &Vg, &mBase, &status) >= 3)
                {
                    int node_idx = find_node_index_by_num(*out_nodes, nBus, gen_bus);
                    if (node_idx >= 0)
                    {
                        (*out_nodes)[node_idx].Pg += Pg / *baseMVA;
                        (*out_nodes)[node_idx].Qg += Qg / *baseMVA;
                        (*out_nodes)[node_idx].P = (*out_nodes)[node_idx].Pg - (*out_nodes)[node_idx].Pl;
                        (*out_nodes)[node_idx].Q = (*out_nodes)[node_idx].Qg - (*out_nodes)[node_idx].Ql;
                    }
                }
            }
            break;
        }
    }

    // 读取线路数据
    while (fgets(buf, sizeof(buf), fp))
    {
        if (strstr(buf, "mpc.branch = [") != NULL)
        {
            int idx = 0;
            while (fgets(buf, sizeof(buf), fp))
            {
                if (strstr(buf, "];") != NULL)
                    break;
                int fbus = 0, tbus = 0;
                double r = 0.0, x = 0.0, b_half = 0.0; // 线路充电电纳 b/2
                if (*out_lines && idx < line_count && line_count > 0)
                {
                    int matched = sscanf_s(buf, "%d %d %lf %lf %lf", &fbus, &tbus, &r, &x, &b_half);
                    if (matched >= 5)
                    {
                        // 只有成功读取数据后才写入
                        (*out_lines)[idx].from_bus = fbus;
                        (*out_lines)[idx].to_bus = tbus;
                        (*out_lines)[idx].series_R = r;
                        (*out_lines)[idx].series_X = x;
                        (*out_lines)[idx].shunt_B = b_half; // 存储线路充电电纳 B/2
                        idx++;
                    }
                }
            }
            break;
        }
    }

    fclose(fp);
    *out_line_count = line_count;
    *out_node_count = nBus;
    *out_init_count = nBus;

    return 0;
}

/*====================== 网络初始化：形成 Ybus ======================*/

/// <summary>
/// 将阻抗转换为导纳
/// </summary>
/// <param name="R">电阻</param>
/// <param name="X">电抗</param>
/// <param name="G">导纳实部</param>
/// <param name="B">导纳虚部</param>
static void impedance_to_admittance(double R, double X, double* G, double* B)
{
    double mod = R * R + X * X; /* 阻抗模平方，不开方减少额外运算，提高精度与效率 */
    if (mod == 0.0)
    {
        *G = 0.0;
        *B = 0.0;
        return;
    }
    *G = R / mod;
    *B = -X / mod;
}

/// <summary>
/// 复制并转换线路参数
/// </summary>
/// <param name="info">网络信息</param>
/// <param name="lines">线路参数</param>
/// <param name="count">线路参数数量</param>
static void copy_and_convert_lines(NetworkInfo* info, const LineArg* lines, size_t count)
{
    size_t i;

    if (!info || (!lines && count > 0)) return;

    info->lines = (LineArg*)calloc(count, sizeof(LineArg));
    ensure_alloc(info->lines, "LineArg");
    /* 仅为安抚静态分析器，逻辑上 ensure_alloc 已经保证非 NULL */
    if (!info->lines) return;
    info->line_count = count;

    for (i = 0; i < count; ++i)
    {
        double G = 0.0, B = 0.0;
        info->lines[i] = lines[i];
        impedance_to_admittance(lines[i].series_R, lines[i].series_X, &G, &B);
        /* 在结构里直接存成导纳的实部/虚部 */
        info->lines[i].series_R = G;
        info->lines[i].series_X = B;
    }
}

/// <summary>
/// 发现节点总数
/// </summary>
/// <param name="nodes">节点参数</param>
/// <param name="node_count">节点参数数量</param>
/// <returns>节点总数</returns>
static int discover_order(const NodeArg* nodes, size_t node_count)
{
    int max_id = 0;
    size_t i;
    for (i = 0; i < node_count; ++i)
    {
        if (nodes[i].node_num > max_id)
            max_id = nodes[i].node_num;
    }
    return max_id;
}

/// <summary>
/// 重新排序节点参数
/// </summary>
/// <param name="info">网络信息</param>
/// <param name="nodes">节点参数</param>
/// <param name="count">节点参数数量</param>
static void reorder_nodes(NetworkInfo* info, const NodeArg* nodes, size_t count)
{
    int idx;

    if (!info || (!nodes && count > 0)) return;

    info->nodes = (NodeArg*)calloc(info->order, sizeof(NodeArg));
    ensure_alloc(info->nodes, "NodeArg");
    if (!info->nodes) return;

    for (idx = 0; idx < info->order; ++idx)
    {
        int target = idx + 1;
        int found = 0;
        size_t j;
        for (j = 0; j < count; ++j)
        {
            if (nodes[j].node_num == target)
            {
                info->nodes[idx] = nodes[j];
                found = 1;
                break;
            }
        }
        if (!found)
        {
            info->nodes[idx].node_num = target;
            info->nodes[idx].type = NODE_PQ;
            info->nodes[idx].P = 0.0;
            info->nodes[idx].Q = 0.0;
            info->nodes[idx].V = 1.0;
            info->nodes[idx].Theta = 0.0;
        }
    }
    info->node_count = info->order;
}

/// <summary>
/// 重新排序初始值
/// </summary>
/// <param name="info">网络信息</param>
/// <param name="init_vals">初始值</param>
/// <param name="count">初始值数量</param>
static void reorder_init_values(NetworkInfo* info, const InitVal* init_vals, size_t count)
{
    int idx;

    if (!info || (!init_vals && count > 0)) return;

    info->init_values = (InitVal*)calloc(info->order, sizeof(InitVal));
    ensure_alloc(info->init_values, "InitVal");
    if (!info->init_values) return;

    for (idx = 0; idx < info->order; ++idx)
    {
        int target = idx + 1;
        int found = 0;
        size_t j;
        for (j = 0; j < count; ++j)
        {
            if (init_vals[j].node_num == target)
            {
                info->init_values[idx] = init_vals[j];
                found = 1;
                break;
            }
        }
        if (!found)
        {
            info->init_values[idx].node_num = target;
            info->init_values[idx].e = 1.0;
            info->init_values[idx].f = 0.0;
        }
    }
    info->init_count = info->order;
}

/// <summary>
/// 分配状态变量
/// </summary>
/// <param name="info">网络信息</param>
static void allocate_state(NetworkInfo* info)
{
    int i;

    if (!info) return;
    info->e = (double*)calloc(info->order, sizeof(double));
    info->f = (double*)calloc(info->order, sizeof(double));
    ensure_alloc(info->e, "state e");
    ensure_alloc(info->f, "state f");
    if (!info->e || !info->f) return;
    for (i = 0; i < info->order; ++i) {
        info->e[i] = info->init_values[i].e;
        info->f[i] = info->init_values[i].f;
    }
}

/// <summary>
/// 分配导纳矩阵
/// </summary>
/// <param name="info">网络信息</param>
static void allocate_admittance(NetworkInfo* info)
{
    if (!info) return;
    size_t mat_size = (size_t)info->order * (size_t)info->order;
    info->G = (double*)calloc(mat_size, sizeof(double));
    info->B = (double*)calloc(mat_size, sizeof(double));
    ensure_alloc(info->G, "G matrix");
    ensure_alloc(info->B, "B matrix");
}

/// <summary>
/// 填充导纳矩阵
/// </summary>
/// <param name="info">网络信息</param>
static void fill_admittance(NetworkInfo* info)
{
    int i, j;
    if (!info || !info->G || !info->B) return;

    for (i = 0; i < info->order; ++i) // 行 (From Bus Index - 1)
    {
        for (j = 0; j < info->order; ++j) // 列 (To Bus Index - 1)
        {
            double val_real = 0.0;
            double val_imag = 0.0;
            if (i == j)
            {
                size_t k;
                for (k = 0; k < info->line_count; ++k)
                {
                    const LineArg* line = &info->lines[k];
                    // 检查线路 k 是否连接到当前节点 i
                    if (line->from_bus == (i + 1) || line->to_bus == (i + 1))
                    {
                        // 累加线路串联导纳的实部和虚部 (G + jB)
                        val_real += line->series_R;
                        val_imag += line->series_X;

                        // 累加线路并联电纳 B/2
                        val_imag += line->shunt_B;
                    }
                }

                // 累加节点自身的并联导纳 Gs + jBs
                const NodeArg* node = &info->nodes[i];
                val_real += node->Gs; // 节点并联电导 Gs (实部)
                val_imag += node->Bs; // 节点并联电纳 Bs (虚部)
            }
            else
            {
                size_t k;
                // 寻找连接节点 i 和 j 的支路
                for (k = 0; k < info->line_count; ++k)
                {
                    const LineArg* line = &info->lines[k];
                    int bus_i = i + 1;
                    int bus_j = j + 1;

                    // 检查是否是 i-j 支路
                    if ((line->from_bus == bus_i && line->to_bus == bus_j) ||
                        (line->from_bus == bus_j && line->to_bus == bus_i))
                    {
                        // 非对角元素 Yij = -Yseries
                        val_real = -line->series_R;
                        val_imag = -line->series_X;
                        break; // 找到支路后退出循环
                    }
                }
            }
            mat_set(info->G, info->order, i, j, val_real);
            mat_set(info->B, info->order, i, j, val_imag);
        }
    }
}

/// <summary>
/// 初始化网络
/// </summary>
/// <param name="info">网络信息</param>
/// <param name="lines">线路参数</param>
/// <param name="line_count">线路参数数量</param>
/// <param name="nodes">节点参数</param>
/// <param name="node_count">节点参数数量</param>
/// <param name="init_vals">初始值</param>
/// <param name="init_count">初始值数量</param>
static void init_network(NetworkInfo* info,
    const LineArg* lines,
    size_t line_count,
    const NodeArg* nodes,
    size_t node_count,
    const InitVal* init_vals,
    size_t init_count)
{
    if (!info) return;

    memset(info, 0, sizeof(*info));
    copy_and_convert_lines(info, lines, line_count);
    info->order = discover_order(nodes, node_count);
    reorder_nodes(info, nodes, node_count);
    reorder_init_values(info, init_vals, node_count);
    allocate_state(info);
    allocate_admittance(info);
    fill_admittance(info);
}

/// <summary>
/// 释放网络
/// </summary>
/// <param name="info">网络信息</param>
static void free_network(NetworkInfo* info)
{
    if (!info) return;

    if (info->lines)       free(info->lines);
    if (info->nodes)       free(info->nodes);
    if (info->init_values) free(info->init_values);
    if (info->G)           free(info->G);
    if (info->B)           free(info->B);
    if (info->e)           free(info->e);
    if (info->f)           free(info->f);
    memset(info, 0, sizeof(*info));
}

/*====================== 牛顿-拉夫逊迭代相关 ======================*/

/// <summary>
/// 定位平衡节点
/// </summary>
/// <param name="info">网络信息</param>
/// <returns>平衡节点索引</returns>
static int locate_slack(const NetworkInfo* info)
{
    int i;
    for (i = 0; i < info->order; ++i)
    {
        if (info->nodes[i].type == NODE_SLACK) return i;
    }
    return -1;
}

/// <summary>
/// 构建索引集合
/// </summary>
/// <param name="ctx">牛顿-拉夫逊迭代上下文</param>
static void build_index_sets(NLIteration* ctx)
{
    int i;

    if (!ctx || !ctx->info) return;
    ctx->pq_count = 0;
    ctx->pv_count = 0;

    for (i = 0; i < ctx->info->order; ++i)
    {
        if (ctx->info->nodes[i].type == NODE_SLACK)
            continue;
        if (ctx->info->nodes[i].type == NODE_PQ)
            ctx->pq_count++;
        else if (ctx->info->nodes[i].type == NODE_PV)
            ctx->pv_count++;
    }

    ctx->non_slack = ctx->pq_count + ctx->pv_count;
    ctx->eq_count = ctx->non_slack * 2;

    /* 为了避免静态分析器认为会越界，这里按全部节点数分配，实际只用前 non_slack 个 */
    ctx->var_indices = (int*)calloc((size_t)ctx->info->order, sizeof(int));
    ensure_alloc(ctx->var_indices, "var indices");
    if (!ctx->var_indices) return;

    /* 序号：先 PQ 再 PV */
    {
        int cursor = 0;
        for (i = 0; i < ctx->info->order; ++i)
        {
            if (ctx->info->nodes[i].type == NODE_PQ)
                ctx->var_indices[cursor++] = i;
        }
        for (i = 0; i < ctx->info->order; ++i)
        {
            if (ctx->info->nodes[i].type == NODE_PV)
                ctx->var_indices[cursor++] = i;
        }
    }
}

/// <summary>
/// 初始化牛顿-拉夫逊迭代
/// </summary>
/// <param name="ctx">牛顿-拉夫逊迭代上下文</param>
/// <param name="info">网络信息</param>
/// <returns>0: 成功, -1: 失败</returns>
static int init_iteration(NLIteration* ctx, NetworkInfo* info)
{
    if (!ctx || !info) return -1;

    memset(ctx, 0, sizeof(*ctx));
    ctx->info = info;
    ctx->slack_index = locate_slack(info);
    if (ctx->slack_index < 0)
    {
        fprintf(stderr, "Slack bus not found\n\n");
        return -1;
    }

    build_index_sets(ctx);

    ctx->delta = (double*)calloc((size_t)ctx->eq_count, sizeof(double));
    ensure_alloc(ctx->delta, "delta");
    
    /* 初始化稀疏矩阵 */
    sparse_matrix_init(&ctx->J_sparse, ctx->eq_count, ctx->eq_count);
    
    return 0;
}

/// <summary>
/// 释放牛顿-拉夫逊迭代
/// </summary>
/// <param name="ctx">牛顿-拉夫逊迭代上下文</param>
static void free_iteration(NLIteration* ctx)
{
    if (!ctx) return;

    if (ctx->var_indices) free(ctx->var_indices);
    if (ctx->delta)       free(ctx->delta);
    sparse_matrix_free(&ctx->J_sparse);
    memset(ctx, 0, sizeof(*ctx));
}

/// <summary>
/// 计算功率平衡
/// </summary>
/// <param name="info">网络信息</param>
/// <param name="node_idx">节点索引</param>
/// <param name="P_calc">计算有功功率</param>
/// <param name="Q_calc">计算无功功率</param>
static void compute_power_balance(const NetworkInfo* info,
    int node_idx,
    double* P_calc,
    double* Q_calc)
{
    int j;
    double e_i;
    double f_i;

    if (!info || !P_calc || !Q_calc || !info->e || !info->f) return;

    e_i = info->e[node_idx];
    f_i = info->f[node_idx];
    *P_calc = 0.0;
    *Q_calc = 0.0;

    for (j = 0; j < info->order; ++j)
    {
        double e_j = info->e[j];
        double f_j = info->f[j];
        double Gij = mat_get(info->G, info->order, node_idx, j);
        double Bij = mat_get(info->B, info->order, node_idx, j);
        /* 优化计算顺序，减少舍入误差 */
        double Gij_ej = Gij * e_j;
        double Gij_fj = Gij * f_j;
        double Bij_fj = Bij * f_j;
        double Bij_ej = Bij * e_j;
        *P_calc += (e_i * (Gij_ej - Bij_fj)) + (f_i * (Gij_fj + Bij_ej));
        *Q_calc += (f_i * (Gij_ej - Bij_fj)) - (e_i * (Gij_fj + Bij_ej));
    }
}

/// <summary>
/// 批量计算所有节点的功率
/// </summary>
/// <param name="info">网络信息</param>
/// <param name="P_calc">输出数组，计算的有功功率</param>
/// <param name="Q_calc">输出数组，计算的无功功率</param>
static void compute_all_power_balance(const NetworkInfo* info, double* P_calc, double* Q_calc)
{
    int i, j;

    if (!info || !P_calc || !Q_calc) return;

    /* 初始化 */
    for (i = 0; i < info->order; ++i)
    {
        P_calc[i] = 0.0;
        Q_calc[i] = 0.0;
    }

    /* 批量计算：类似 MATLAB 的 V .* conj(Ybus * V) */
    for (i = 0; i < info->order; ++i)
    {
        double e_i = info->e[i];
        double f_i = info->f[i];

        for (j = 0; j < info->order; ++j)
        {
            double e_j = info->e[j];
            double f_j = info->f[j];
            double Gij = mat_get(info->G, info->order, i, j);
            double Bij = mat_get(info->B, info->order, i, j);

            /* 优化计算顺序，减少舍入误差 */
            double Gij_ej = Gij * e_j;
            double Gij_fj = Gij * f_j;
            double Bij_fj = Bij * f_j;
            double Bij_ej = Bij * e_j;
            /* P = e_i * (Gij * e_j - Bij * f_j) + f_i * (Gij * f_j + Bij * e_j) */
            P_calc[i] += (e_i * (Gij_ej - Bij_fj)) + (f_i * (Gij_fj + Bij_ej));
            /* Q = f_i * (Gij * e_j - Bij * f_j) - e_i * (Gij * f_j + Bij * e_j) */
            Q_calc[i] += (f_i * (Gij_ej - Bij_fj)) - (e_i * (Gij_fj + Bij_ej));
        }
    }
}

/// <summary>
/// 计算增量
/// </summary>
/// <param name="ctx">牛顿-拉夫逊迭代上下文</param>
static void calc_delta(NLIteration* ctx)
{
    int eq_pos = 0;
    int i;

    /* 批量计算所有节点的功率（避免重复计算） */
    double* P_calc = (double*)calloc((size_t)ctx->info->order, sizeof(double));
    double* Q_calc = (double*)calloc((size_t)ctx->info->order, sizeof(double));

    if (P_calc && Q_calc)
    {
        compute_all_power_balance(ctx->info, P_calc, Q_calc);

        /* PQ 节点：ΔP, ΔQ */
        for (i = 0; i < ctx->pq_count; ++i)
        {
            int node_idx = ctx->var_indices[i];
            ctx->delta[eq_pos++] = ctx->info->nodes[node_idx].P - P_calc[node_idx];
            ctx->delta[eq_pos++] = ctx->info->nodes[node_idx].Q - Q_calc[node_idx];
        }

        /* PV 节点：ΔP, Δ(|U|^2) */
        for (i = 0; i < ctx->pv_count; ++i)
        {
            int node_idx = ctx->var_indices[ctx->pq_count + i];
            double e_i = ctx->info->e[node_idx];
            double f_i = ctx->info->f[node_idx];

            ctx->delta[eq_pos++] = ctx->info->nodes[node_idx].P - P_calc[node_idx];
            ctx->delta[eq_pos++] = ctx->info->nodes[node_idx].V * ctx->info->nodes[node_idx].V
                - (e_i * e_i + f_i * f_i);
        }

        free(P_calc);
        free(Q_calc);
    }
    else
    {
        /* 如果内存分配失败，回退到原始方法 */
        /* PQ 节点：ΔP, ΔQ */
        for (i = 0; i < ctx->pq_count; ++i)
        {
            int node_idx = ctx->var_indices[i];
            double P_calc_single = 0.0, Q_calc_single = 0.0;
            compute_power_balance(ctx->info, node_idx, &P_calc_single, &Q_calc_single);
            ctx->delta[eq_pos++] = ctx->info->nodes[node_idx].P - P_calc_single;
            ctx->delta[eq_pos++] = ctx->info->nodes[node_idx].Q - Q_calc_single;
        }

        /* PV 节点：ΔP, Δ(|U|^2) */
        for (i = 0; i < ctx->pv_count; ++i)
        {
            int node_idx = ctx->var_indices[ctx->pq_count + i];
            double P_calc_single = 0.0, Q_calc_single = 0.0;
            double e_i, f_i;

            compute_power_balance(ctx->info, node_idx, &P_calc_single, &Q_calc_single);
            e_i = ctx->info->e[node_idx];
            f_i = ctx->info->f[node_idx];

            ctx->delta[eq_pos++] = ctx->info->nodes[node_idx].P - P_calc_single;
            ctx->delta[eq_pos++] = ctx->info->nodes[node_idx].V * ctx->info->nodes[node_idx].V
                - (e_i * e_i + f_i * f_i);
        }
    }
}

/// <summary>
/// 预先计算所有节点的A_i和B_i（用于雅可比矩阵对角元素）
/// </summary>
/// <param name="ctx">牛顿-拉夫逊迭代上下文</param>
/// <param name="A">输出数组，A[i] = sum_k(Gik*fk + Bik*ek)</param>
/// <param name="B">输出数组，B[i] = sum_k(Gik*ek - Bik*fk)</param>
static void precompute_AB(NLIteration* ctx, double* A, double* B)
{
    const NetworkInfo* info = ctx->info;
    int i, k;

    if (!A || !B) return;

    for (i = 0; i < info->order; ++i)
    {
        A[i] = 0.0;
        B[i] = 0.0;
        double e_i = info->e[i];
        double f_i = info->f[i];

        for (k = 0; k < info->order; ++k)
        {
            double e_k = info->e[k];
            double f_k = info->f[k];
            double Gik = mat_get(info->G, info->order, i, k);
            double Bik = mat_get(info->B, info->order, i, k);

            // B_i 对应实部 Gik*ek - Bik*fk
            B[i] += Gik * e_k - Bik * f_k;
            // A_i 对应虚部 Gik*fk + Bik*ek
            A[i] += Gik * f_k + Bik * e_k;
        }
    }
}

/// <summary>
/// 生成稀疏雅可比矩阵
/// </summary>
/// <param name="ctx">牛顿-拉夫逊迭代上下文</param>
static void gen_sparse_jacobian(NLIteration* ctx)
{
    int row, col;
    int* row_idx = NULL;
    int* col_idx = NULL;
    double* values = NULL;
    int nnz_count = 0;
    int max_nnz_estimate;

    if (!ctx) return;

    /* 估算最大非零元素个数：每个节点最多与所有其他节点相连 */
    /* 对于雅可比矩阵，每个节点对最多有4个非零元素（2x2块） */
    max_nnz_estimate = ctx->eq_count * ctx->eq_count / 4; /* 保守估计 */
    if (max_nnz_estimate < ctx->eq_count * 10) /* 至少每行10个非零元素 */
        max_nnz_estimate = ctx->eq_count * 10;

    row_idx = (int*)calloc((size_t)max_nnz_estimate, sizeof(int));
    col_idx = (int*)calloc((size_t)max_nnz_estimate, sizeof(int));
    values = (double*)calloc((size_t)max_nnz_estimate, sizeof(double));
    ensure_alloc(row_idx, "sparse row_idx");
    ensure_alloc(col_idx, "sparse col_idx");
    ensure_alloc(values, "sparse values");

    /* 预先计算所有节点的A_i和B_i */
    double* A = (double*)calloc((size_t)ctx->info->order, sizeof(double));
    double* B = (double*)calloc((size_t)ctx->info->order, sizeof(double));
    if (!A || !B)
    {
        free(row_idx);
        free(col_idx);
        free(values);
        return;
    }
    precompute_AB(ctx, A, B);

    /* 收集所有非零元素 */
    for (row = 0; row < ctx->non_slack; ++row)
    {
        int i = ctx->var_indices[row]; /* 当前节点在网络中的索引 */
        const NetworkInfo* info = ctx->info;
        double e_i = info->e[i];
        double f_i = info->f[i];
        int eq_row1 = 2 * row;     /* ΔP 行 */
        int eq_row2 = 2 * row + 1;  /* ΔQ 或 Δ|U|^2 行 */

        for (col = 0; col < ctx->non_slack; ++col)
        {
            int j = ctx->var_indices[col]; /* 变量对应的节点索引 */
            double Gij = mat_get(info->G, info->order, i, j);
            double Bij = mat_get(info->B, info->order, i, j);
            double dP_dfj, dP_dej, dQ_dfj, dQ_dej;

            if (i == j)
            {
                /* 对角元素 */
                double A_i = A[i];
                double B_i = B[i];

                dP_dfj = A_i + f_i * Gij - e_i * Bij;
                dP_dej = B_i + e_i * Gij + f_i * Bij;

                if (info->nodes[i].type == NODE_PV)
                {
                    dQ_dfj = 2.0 * f_i;
                    dQ_dej = 2.0 * e_i;
                }
                else
                {
                    dQ_dfj = B_i - f_i * Bij - e_i * Gij;
                    dQ_dej = -A_i - e_i * Bij + f_i * Gij;
                }
            }
            else
            {
                /* 非对角元素 */
                dP_dfj = -e_i * Bij + f_i * Gij;
                dP_dej = e_i * Gij + f_i * Bij;

                if (info->nodes[i].type == NODE_PV)
                {
                    dQ_dfj = 0.0;
                    dQ_dej = 0.0;
                }
                else
                {
                    dQ_dfj = -e_i * Gij - f_i * Bij;
                    dQ_dej = -e_i * Bij + f_i * Gij;
                }
            }

            int eq_col1 = 2 * col;     /* f_j 列 */
            int eq_col2 = 2 * col + 1;  /* e_j 列 */

            /* 只存储非零元素 */
            if (fabs(dP_dfj) > 1e-15)
            {
                if (nnz_count >= max_nnz_estimate) break;
                row_idx[nnz_count] = eq_row1;
                col_idx[nnz_count] = eq_col1;
                values[nnz_count] = dP_dfj;
                nnz_count++;
            }
            if (fabs(dP_dej) > 1e-15)
            {
                if (nnz_count >= max_nnz_estimate) break;
                row_idx[nnz_count] = eq_row1;
                col_idx[nnz_count] = eq_col2;
                values[nnz_count] = dP_dej;
                nnz_count++;
            }
            if (fabs(dQ_dfj) > 1e-15)
            {
                if (nnz_count >= max_nnz_estimate) break;
                row_idx[nnz_count] = eq_row2;
                col_idx[nnz_count] = eq_col1;
                values[nnz_count] = dQ_dfj;
                nnz_count++;
            }
            if (fabs(dQ_dej) > 1e-15)
            {
                if (nnz_count >= max_nnz_estimate) break;
                row_idx[nnz_count] = eq_row2;
                col_idx[nnz_count] = eq_col2;
                values[nnz_count] = dQ_dej;
                nnz_count++;
            }
        }
    }

    free(A);
    free(B);

    /* 释放旧的稀疏矩阵 */
    sparse_matrix_free(&ctx->J_sparse);

    /* 构建新的稀疏矩阵 */
    sparse_matrix_init(&ctx->J_sparse, ctx->eq_count, ctx->eq_count);
    sparse_matrix_from_triplets(&ctx->J_sparse, row_idx, col_idx, values, nnz_count);

    free(row_idx);
    free(col_idx);
    free(values);
}

/// <summary>
/// 生成雅可比矩阵（保留密集版本作为备用）
/// </summary>
/// <param name="ctx">牛顿-拉夫逊迭代上下文</param>
static void gen_jacobian(NLIteration* ctx)
{
    /* 使用稀疏版本 */
    gen_sparse_jacobian(ctx);
}

/// <summary>
/// 稀疏LU分解求解器（使用部分主元选择）
/// </summary>
/// <param name="A">稀疏系数矩阵</param>
/// <param name="b">常数项向量</param>
/// <param name="x">解向量</param>
static void solve_sparse_linear_system(const SparseMatrix* A, const double* b, double* x)
{
    int i, j, k, n;
    int* perm = NULL;  /* 行置换数组 */
    double* y = NULL;

    if (!A || !b || !x) return;
    n = A->nrows;
    
    /* 分配密集矩阵用于LU分解（初始化为 0） */
    double* dense = (double*)calloc((size_t)n * (size_t)n, sizeof(double));
    ensure_alloc(dense, "dense matrix for LU");
    if (!dense) return;

    /* 
     * 将稀疏矩阵高效转换为密集格式：
     * 只遍历非零元，复杂度 O(nnz)，避免反复调用 sparse_matrix_get 的额外开销
     */
    for (i = 0; i < n; ++i)
    {
        int row_start = A->row_ptr[i];
        int row_end   = A->row_ptr[i + 1];
        for (j = row_start; j < row_end; ++j)
        {
            int col = A->col_idx[j];
            if (col >= 0 && col < n)
            {
                dense[i * n + col] = A->values[j];
            }
        }
    }

    /* 分配行置换数组 */
    perm = (int*)calloc((size_t)n, sizeof(int));
    ensure_alloc(perm, "permutation");
    if (!perm) { free(dense); return; }
    for (i = 0; i < n; ++i) perm[i] = i;

    /* 执行LU分解（带部分主元选择） */
    for (k = 0; k < n - 1; ++k)
    {
        /* 寻找主元 */
        int pivot = k;
        double max_val = fabs(dense[k * n + k]);
        for (i = k + 1; i < n; ++i)
        {
            double val = fabs(dense[i * n + k]);
            if (val > max_val)
            {
                max_val = val;
                pivot = i;
            }
        }

        if (fabs(max_val) < 1e-12)
        {
            fprintf(stderr, "Jacobian matrix is singular at column %d\n", k);
            free(dense);
            free(perm);
            exit(EXIT_FAILURE);
        }

        /* 交换行 */
        if (pivot != k)
        {
            int tmp_perm = perm[k];
            perm[k] = perm[pivot];
            perm[pivot] = tmp_perm;
            for (j = 0; j < n; ++j)
            {
                double tmp = dense[k * n + j];
                dense[k * n + j] = dense[pivot * n + j];
                dense[pivot * n + j] = tmp;
            }
        }

        /* 消元 */
        double pivot_val = dense[k * n + k];
        for (i = k + 1; i < n; ++i)
        {
            double factor = dense[i * n + k] / pivot_val;
            dense[i * n + k] = factor; /* 存储L矩阵的元素 */
            for (j = k + 1; j < n; ++j)
            {
                dense[i * n + j] -= factor * dense[k * n + j];
            }
        }
    }

    /* 前向替换：Ly = Pb */
    y = (double*)calloc((size_t)n, sizeof(double));
    ensure_alloc(y, "y vector");
    if (!y) { free(dense); free(perm); return; }

    for (i = 0; i < n; ++i)
    {
        double sum = b[perm[i]];
        for (j = 0; j < i; ++j)
        {
            sum -= dense[i * n + j] * y[j];
        }
        y[i] = sum;
    }

    /* 后向替换：Ux = y */
    for (i = n - 1; i >= 0; --i)
    {
        double sum = y[i];
        for (j = i + 1; j < n; ++j)
        {
            sum -= dense[i * n + j] * x[j];
        }
        if (fabs(dense[i * n + i]) < 1e-12)
        {
            fprintf(stderr, "U matrix is singular at row %d\n", i);
            free(dense);
            free(perm);
            free(y);
            exit(EXIT_FAILURE);
        }
        x[i] = sum / dense[i * n + i];
    }

    free(dense);
    free(perm);
    free(y);
}

/// <summary>
/// 应用校正
/// </summary>
/// <param name="ctx">牛顿-拉夫逊迭代上下文</param>
/// <param name="correction">校正向量</param>
/// <returns>0: 应用成功, 1: 无需应用（已收敛）</returns>
static int apply_correction(NLIteration* ctx, const double* correction)
{
    int i;

    if (!ctx || !correction) return 1;

    for (i = 0; i < ctx->non_slack; ++i)
    {
        int node_idx = ctx->var_indices[i];
        ctx->info->f[node_idx] += correction[2 * i];
        ctx->info->e[node_idx] += correction[2 * i + 1];
    }
    return 0;
}

/// <summary>
/// 开始迭代
/// </summary>
/// <param name="ctx">牛顿-拉夫逊迭代上下文</param>
/// <param name="tolerance">容差，用于判断迭代是否收敛</param>
/// <param name="max_iter">最大迭代次数，防止无限迭代</param>
/// <param name="converged_out">迭代是否收敛</param>
/// <param name="iter_out">完成迭代的次数</param>
/// <param name="runtime_out">迭代所用的运行时间（秒）</param>
static void start_iteration(NLIteration* ctx,
    double tolerance,
    int max_iter,
    int* converged_out,
    int* iter_out,
    double* runtime_out)
{
    double* correction;
    int converged = 0;
    int iteration = 0;
    double runtime = 0.0;

    if (!ctx) return;

    ctx->tolerance = tolerance;
    ctx->max_iter = max_iter;

    /* 为校正向量分配内存，大小为eq_count个double类型的空间 */
    correction = (double*)calloc((size_t)ctx->eq_count, sizeof(double));
    ensure_alloc(correction, "correction vector");

    /* 高精度计时 */
    LARGE_INTEGER freq, t0, t1;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&t0);

    /* 初始检查：计算初始功率不匹配 */
    calc_delta(ctx);
    {
        int i;
        double max_mismatch = 0.0;
        for (i = 0; i < ctx->eq_count; ++i)
        {
            double val = fabs(ctx->delta[i]);
            if (val > max_mismatch)
                max_mismatch = val;
        }
        if (max_mismatch < tolerance)
        {
            converged = 1;
            if (converged_out) *converged_out = converged;
            if (iter_out) *iter_out = iteration;
            QueryPerformanceCounter(&t1);
            runtime = (double)(t1.QuadPart - t0.QuadPart) / (double)freq.QuadPart;
            if (runtime_out) *runtime_out = runtime;
            free(correction);
            return;
        }
    }

    while (!converged && iteration < ctx->max_iter)
    {
        int i;
        gen_jacobian(ctx);
        solve_sparse_linear_system(&ctx->J_sparse, ctx->delta, correction);
        apply_correction(ctx, correction);
        
        /* 计算更新后的功率不匹配 */
        calc_delta(ctx);
        
        /* 检查收敛：使用功率不匹配的无穷范数（与MATLAB版本一致） */
        {
            double max_mismatch = 0.0;
            for (i = 0; i < ctx->eq_count; ++i)
            {
                double val = fabs(ctx->delta[i]);
                if (val > max_mismatch)
                    max_mismatch = val;
            }
            if (max_mismatch < tolerance)
                converged = 1;
        }
        
        iteration++;
    }

    QueryPerformanceCounter(&t1);
    runtime = (double)(t1.QuadPart - t0.QuadPart) / (double)freq.QuadPart;

    if (converged)
        printf("Converged : YES\n");
    else
        printf("Converged : NO  (max iterations = %d)\n", ctx->max_iter);
    printf("Iterations: %d\n", iteration);
    printf("Runtime   : %.3f ms\n\n", runtime * 1000.0);

    if (converged_out) *converged_out = converged;
    if (iter_out) *iter_out = iteration;
    if (runtime_out) *runtime_out = runtime;

    free(correction);
}

/*====================== 结果输出（仿样例） ======================*/

/// <summary>
/// 更新平衡节点和PV节点的发电量（根据功率平衡方程计算）
/// </summary>
/// <param name="info">网络信息</param>
static void update_generation(NetworkInfo* info)
{
    int i;
    double P_calc, Q_calc;

    if (!info) return;

    /* 对于平衡节点和PV节点，逐个计算功率并更新发电量 */
    for (i = 0; i < info->order; ++i)
    {
        if (info->nodes[i].type == NODE_SLACK || info->nodes[i].type == NODE_PV)
        {
            /* 计算节点的注入功率 */
            compute_power_balance(info, i, &P_calc, &Q_calc);

            /* 发电量 = 计算出的注入功率 + 负荷 */
            /* P_injection = P_gen - P_load, 所以 P_gen = P_injection + P_load */
            info->nodes[i].Pg = P_calc + info->nodes[i].Pl;
            info->nodes[i].Qg = Q_calc + info->nodes[i].Ql;

            /* 更新净注入功率 */
            info->nodes[i].P = info->nodes[i].Pg - info->nodes[i].Pl;
            info->nodes[i].Q = info->nodes[i].Qg - info->nodes[i].Ql;
        }
    }
}

/// <summary>
/// 输出潮流计算结果，输出到控制台
/// </summary>
/// <param name="info">网络信息</param>
/// <param name="converged">是否收敛</param>
/// <param name="runtime">运行时间</param>
static void print_report(const NetworkInfo* info,
    int converged,
    double runtime,
    double baseMVA)
{
    int i;
    double totalPg = 0.0, totalQg = 0.0;
    double totalPl = 0.0, totalQl = 0.0;
    double sumPloss = 0.0, sumQloss = 0.0; // 用于存储累加的支路串联损耗

    if (!info) return;

    /* 在输出前，更新平衡节点和PV节点的发电量（根据功率平衡计算） */
    update_generation((NetworkInfo*)info);

    for (i = 0; i < info->order; ++i)
    {
        totalPg += info->nodes[i].Pg;
        totalQg += info->nodes[i].Qg;
        totalPl += info->nodes[i].Pl;
        totalQl += info->nodes[i].Ql;
    }

    printf("%% ===============================================================================\n");
    printf("%% |     System Summary                                                           |\n");
    printf("%% ===============================================================================\n");
    printf("%% How many?                How much?              P (MW)            Q (MVAr)\n");
    printf("%% ---------------------    -------------------  -------------  -----------------\n");
    printf("%% Buses              %2d     Total Gen Capacity       -                 -\n",
        info->order);

    int genBus = 0, loadBus = 0;
    for (i = 0; i < info->order; ++i)
    {
        if (info->nodes[i].Pg > 1e-6) genBus++;
        if (info->nodes[i].Pl > 1e-6) loadBus++;
    }
    printf("%% Generators         %2d     Generation (actual)  %11.2f       %11.2f\n",
        genBus, totalPg * baseMVA, totalQg * baseMVA);
    printf("%% Loads              %2d     Load                %11.2f       %11.2f\n",
        loadBus, totalPl * baseMVA, totalQl * baseMVA);

    printf("%% Shunts             0     Shunt (inj)              0.00              0.00\n");

    /* 计算并累加总损耗 (仅使用串联部分损耗) */
    if (info->lines && info->line_count > 0)
    {
        size_t k;
        for (k = 0; k < info->line_count; ++k)
        {
            int fb = info->lines[k].from_bus - 1;
            int tb = info->lines[k].to_bus - 1;
            double G = info->lines[k].series_R;
            double B_series = info->lines[k].series_X;

            double Vi_re = info->e[fb], Vi_im = info->f[fb];
            double Vj_re = info->e[tb], Vj_im = info->f[tb];

            /* 从 i 到 j 的串联电流 */
            double dV_re_ij = Vi_re - Vj_re;
            double dV_im_ij = Vi_im - Vj_im;
            double Iser_ij_re = G * dV_re_ij - B_series * dV_im_ij;
            double Iser_ij_im = G * dV_im_ij + B_series * dV_re_ij;

            /* 从 j 到 i 的串联电流 */
            double dV_re_ji = Vj_re - Vi_re;
            double dV_im_ji = Vj_im - Vi_im;
            double Iser_ji_re = G * dV_re_ji - B_series * dV_im_ji;
            double Iser_ji_im = G * dV_im_ji + B_series * dV_re_ji;

            /* 从 i 到 j 的串联功率 Sser_ij = Vi * conj(Iser_ij) */
            /* 使用更稳定的计算顺序 */
            double Pser_ij = (Vi_re * Iser_ij_re) + (Vi_im * Iser_ij_im);
            double Qser_ij = (Vi_im * Iser_ij_re) - (Vi_re * Iser_ij_im);

            /* 从 j 到 i 的串联功率 Sser_ji = Vj * conj(Iser_ji) */
            double Pser_ji = (Vj_re * Iser_ji_re) + (Vj_im * Iser_ji_im);
            double Qser_ji = (Vj_im * Iser_ji_re) - (Vj_re * Iser_ji_im);

            // 损耗 (I^2 * Z) = Sser_ij + Sser_ji
            sumPloss += (Pser_ij + Pser_ji) * baseMVA;
            sumQloss += (Qser_ij + Qser_ji) * baseMVA;
        }
    }

    printf("%% Branches        %4zu     Losses (I^2 * Z)      %11.2f       %11.2f\n\n",
        info->line_count, sumPloss, sumQloss);


    printf("%% ===============================================================================\n");
    printf("%% |     Bus Data                                                                 |\n");
    printf("%% ===============================================================================\n");
    printf("%%  Bus      Voltage          Generation             Load        \n");
    printf("%%   #   Mag(pu) Ang(deg)   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)\n");
    printf("%% ----- ------- --------  --------  --------  --------  --------\n");

    for (i = 0; i < info->order; ++i)
    {
        double e = info->e[i];
        double f = info->f[i];
        double mag = sqrt(e * e + f * f);
        double ang = atan2(f, e) * 180.0 / PI;
        printf("%% %5d %7.3f %8.3f  %8.2f  %8.2f  %8.2f  %8.2f\n",
            i + 1,
            mag,
            ang,
            info->nodes[i].Pg * baseMVA,
            info->nodes[i].Qg * baseMVA,
            info->nodes[i].Pl * baseMVA,
            info->nodes[i].Ql * baseMVA);
    }

    printf("%%                         --------  --------  --------  --------\n");
    printf("%%                Total:  %8.2f  %8.2f  %8.2f  %8.2f\n",
        totalPg * baseMVA,
        totalQg * baseMVA,
        totalPl * baseMVA,
        totalQl * baseMVA);

    printf("\n");

    printf("%% ===============================================================================\n");
    printf("%% |     Branch Data                                                              |\n");
    printf("%% ===============================================================================\n");
    printf("%% Brnch   From   To    From Bus Injection   To Bus Injection     Loss (I^2 * Z)  \n");
    printf("%%   #     Bus    Bus    P (MW)   Q (MVAr)   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)\n");
    printf("%% -----  -----  -----  --------  --------  --------  --------  --------  --------\n");

    /* 计算并打印 Branch Data */
    if (info->lines && info->line_count > 0)
    {
        size_t k;
        // 重新设置损耗累加器
        sumPloss = 0.0;
        sumQloss = 0.0;

        for (k = 0; k < info->line_count; ++k)
        {
            int fb = info->lines[k].from_bus - 1;
            int tb = info->lines[k].to_bus - 1;
            /* 注意：series_R 和 series_X 在 copy_and_convert_lines 中已被转换为导纳的实部和虚部 */
            double Yseries_G = info->lines[k].series_R;  // 支路串联导纳的实部
            double Yseries_B = info->lines[k].series_X;  // 支路串联导纳的虚部
            double b2 = info->lines[k].shunt_B;  // 并联电纳 B/2

            double Vi_re = info->e[fb], Vi_im = info->f[fb];
            double Vj_re = info->e[tb], Vj_im = info->f[tb];

            /* 计算支路功率，直接使用支路参数计算（更准确） */

            /* 计算电压差 */
            double dV_re = Vi_re - Vj_re;
            double dV_im = Vi_im - Vj_im;

            /* 从节点i流入支路的串联电流：Iser_ij = Y_series * (Vi - Vj) */
            double Iser_ij_re = Yseries_G * dV_re - Yseries_B * dV_im;
            double Iser_ij_im = Yseries_G * dV_im + Yseries_B * dV_re;

            /* 节点i侧的并联电流：Ishunt_i = j * (b/2) * Vi */
            double Ish_ij_re = -b2 * Vi_im;
            double Ish_ij_im = b2 * Vi_re;

            /* 从节点i流入支路的总电流：Iij = Iser_ij + Ishunt_i */
            double Iij_re = Iser_ij_re + Ish_ij_re;
            double Iij_im = Iser_ij_im + Ish_ij_im;

            /* 从节点i流入支路的功率：Sij = Vi * conj(Iij) */
            /* 使用更稳定的计算顺序以减少舍入误差 */
            /* Sij = (Vi_re + j*Vi_im) * (Iij_re - j*Iij_im) */
            /*     = Vi_re*Iij_re + Vi_im*Iij_im + j*(Vi_im*Iij_re - Vi_re*Iij_im) */
            /* 优化：先计算乘积项，再求和，减少舍入误差 */
            double Pij = (Vi_re * Iij_re) + (Vi_im * Iij_im);
            double Qij = (Vi_im * Iij_re) - (Vi_re * Iij_im);

            /* 从节点j流入支路的串联电流：Iser_ji = Y_series * (Vj - Vi) = -Y_series * (Vi - Vj) */
            double Iser_ji_re = -Iser_ij_re;
            double Iser_ji_im = -Iser_ij_im;

            /* 节点j侧的并联电流：Ishunt_j = j * (b/2) * Vj */
            double Ish_ji_re = -b2 * Vj_im;
            double Ish_ji_im = b2 * Vj_re;

            /* 从节点j流入支路的总电流：Iji = Iser_ji + Ishunt_j */
            double Iji_re = Iser_ji_re + Ish_ji_re;
            double Iji_im = Iser_ji_im + Ish_ji_im;

            /* 从节点j流入支路的功率：Sji = Vj * conj(Iji) */
            /* 使用更稳定的计算顺序以减少舍入误差 */
            /* Sji = (Vj_re + j*Vj_im) * (Iji_re - j*Iji_im) */
            /*     = Vj_re*Iji_re + Vj_im*Iji_im + j*(Vj_im*Iji_re - Vj_re*Iji_im) */
            double Pji = (Vj_re * Iji_re) + (Vj_im * Iji_im);
            double Qji = (Vj_im * Iji_re) - (Vj_re * Iji_im);

            /* 计算串联损耗 */
            /* 从 i 到 j 的串联功率 Sser_ij = Vi * conj(Iser_ij) */
            /* 使用更稳定的计算顺序 */
            double Pser_ij = (Vi_re * Iser_ij_re) + (Vi_im * Iser_ij_im);
            double Qser_ij = (Vi_im * Iser_ij_re) - (Vi_re * Iser_ij_im);
            /* 从 j 到 i 的串联功率 Sser_ji = Vj * conj(Iser_ji) */
            double Pser_ji = (Vj_re * Iser_ji_re) + (Vj_im * Iser_ji_im);
            double Qser_ji = (Vj_im * Iser_ji_re) - (Vj_re * Iser_ji_im);

            /* 计算最终打印的 Loss */
            double Ploss = (Pser_ij + Pser_ji) * baseMVA;
            double Qloss = (Qser_ij + Qser_ji) * baseMVA;

            sumPloss += Ploss;
            sumQloss += Qloss;

            printf("%% %5zu  %5d  %5d  %8.2f  %8.2f  %8.2f  %8.2f  %8.3f  %8.3f\n",
                k + 1,
                info->lines[k].from_bus,
                info->lines[k].to_bus,
                Pij * baseMVA, // Total Injection
                Qij * baseMVA, // Total Injection
                Pji * baseMVA, // Total Injection
                Qji * baseMVA, // Total Injection
                Ploss,         // Series Loss
                Qloss);        // Series Loss
        }

        printf("%%                                                              --------  --------\n");
        printf("%%                                                     Total:  %8.3f  %8.3f\n",
            sumPloss,
            sumQloss);
    }

    printf("\n");
}

/*====================== main：程序入口 ======================*/
int main(void)
{
    const char* input_path = "Input/case4.m";
    LineArg* line_args = NULL;
    NodeArg* node_args = NULL;
    InitVal* init_vals = NULL;
    size_t line_count = 0, node_count = 0, init_count = 0;
    double baseMVA = 100.0;

    NetworkInfo info;
    NLIteration iter;
    int converged = 0;
    int iterations = 0;
    double runtime = 0.0;

    if (read_case_m(input_path, &line_args, &line_count, &node_args, &node_count, &init_vals, &init_count, &baseMVA) != 0)
        return 1;

    init_network(&info, line_args, line_count, node_args, node_count, init_vals, init_count);

    if (init_iteration(&iter, &info) != 0)
    {
        free(line_args);
        free(node_args);
        free(init_vals);
        free_network(&info);
        return 1;
    }

    /* 提高收敛精度以减少数值误差，MATPOWER默认使用1e-8 */
    start_iteration(&iter, 1e-8, 1000, &converged, &iterations, &runtime);

    /* 按课程给的 printpf 样式输出结果（输出到控制台） */
    print_report(&info, converged, runtime, baseMVA);

    free_iteration(&iter);
    free_network(&info);
    free(line_args);
    free(node_args);
    free(init_vals);

    return 0;
}