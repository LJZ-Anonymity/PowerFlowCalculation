// #include <stdbool.h> 减少调用标准库，converged用int代替bool
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
    double shunt_G; // 并联电导
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
    double* B;      /* 节点导纳矩阵虚部（电纳） */

    double* e;      /* 节点电压实部 */
    double* f;      /* 节点电压虚部 */
} NetworkInfo; // 网络信息

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
    double* J;          /* 雅可比矩阵 */

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
                if (sscanf_s(buf, "%d %d %lf %lf", &fbus, &tbus, &r, &x) == 4)
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
                    (*out_nodes)[idx].P = -(*out_nodes)[idx].Pl;
                    (*out_nodes)[idx].Q = -(*out_nodes)[idx].Ql;
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
                if (sscanf_s(buf, "%d %lf %lf", &gen_bus, &Pg, &Qg) >= 3)
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
                double r = 0.0, x = 0.0;
                if (*out_lines && idx < line_count && line_count > 0)
                {
                    // 先尝试读取数据
                    if (sscanf_s(buf, "%d %d %lf %lf", &fbus, &tbus, &r, &x) == 4)
                    {
                        // 只有成功读取数据后才写入
                        (*out_lines)[idx].from_bus = fbus;
                        (*out_lines)[idx].to_bus = tbus;
                        (*out_lines)[idx].series_R = r;
                        (*out_lines)[idx].series_X = x;
                        (*out_lines)[idx].shunt_G = 0.0;
                        (*out_lines)[idx].shunt_B = 0.0;
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
    *B = X / mod;
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
    for (i = 0; i < info->order; ++i)
    {
        for (j = 0; j < info->order; ++j)
        {
            double val_real = 0.0;
            double val_imag = 0.0;
            if (i == j)
            {
                size_t k;
                for (k = 0; k < info->line_count; ++k)
                {
                    const LineArg* line = &info->lines[k];
                    if (line->from_bus == (j + 1) || line->to_bus == (j + 1))
                    {
                        val_real += line->series_R + line->shunt_G;
                        val_imag += line->series_X + line->shunt_B;
                    }
                }
            }
            else
            {
                size_t k;
                for (k = 0; k < info->line_count; ++k)
                {
                    const LineArg* line = &info->lines[k];
                    if ((line->from_bus == (i + 1) && line->to_bus == (j + 1)) ||
                        (line->to_bus == (i + 1) && line->from_bus == (j + 1)))
                    {
                        val_real += -line->series_R;
                        val_imag += -line->series_X;
                    }
                }
            }
            mat_set(info->G, info->order, i, j, val_real);
            /* 注意这里 B 存的是 -Im(Y) */
            mat_set(info->B, info->order, i, j, -val_imag);
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
    reorder_init_values(info, init_vals, init_count);
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
    ctx->J = (double*)calloc((size_t)ctx->eq_count * (size_t)ctx->eq_count,
        sizeof(double));
    ensure_alloc(ctx->delta, "delta");
    ensure_alloc(ctx->J, "Jacobian");
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
    if (ctx->J)           free(ctx->J);
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
        *P_calc += e_i * (Gij * e_j - Bij * f_j) + f_i * (Gij * f_j + Bij * e_j);
        *Q_calc += f_i * (Gij * e_j - Bij * f_j) - e_i * (Gij * f_j + Bij * e_j);
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

    /* PQ 节点：ΔP, ΔQ */
    for (i = 0; i < ctx->pq_count; ++i)
    {
        int node_idx = ctx->var_indices[i];
        double P_calc = 0.0, Q_calc = 0.0;
        compute_power_balance(ctx->info, node_idx, &P_calc, &Q_calc);
        ctx->delta[eq_pos++] = ctx->info->nodes[node_idx].P - P_calc;
        ctx->delta[eq_pos++] = ctx->info->nodes[node_idx].Q - Q_calc;
    }

    /* PV 节点：ΔP, Δ(|U|^2) */
    for (i = 0; i < ctx->pv_count; ++i)
    {
        int node_idx = ctx->var_indices[ctx->pq_count + i];
        double P_calc = 0.0, Q_calc = 0.0;
        double e_i, f_i;

        compute_power_balance(ctx->info, node_idx, &P_calc, &Q_calc);
        e_i = ctx->info->e[node_idx];
        f_i = ctx->info->f[node_idx];

        ctx->delta[eq_pos++] = ctx->info->nodes[node_idx].P - P_calc;
        ctx->delta[eq_pos++] = ctx->info->nodes[node_idx].V * ctx->info->nodes[node_idx].V
            - (e_i * e_i + f_i * f_i);
    }
}

/// <summary>
/// 填充行对
/// </summary>
/// <param name="ctx">牛顿-拉夫逊迭代上下文</param>
/// <param name="row_idx">行索引</param>
/// <param name="row1">第一行</param>
/// <param name="row2">第二行</param>
static void fill_row_pair(NLIteration* ctx, int row_idx, double* row1, double* row2)
{
    int node_idx = ctx->var_indices[row_idx];
    const NetworkInfo* info = ctx->info;
    double e_i = info->e[node_idx];
    double f_i = info->f[node_idx];
    double H_ii = 0.0, N_ii = 0.0;
    int k;
    double Gii, Bii, J_ii, L_ii;
    int col;

    for (k = 0; k < info->order; ++k)
    {
        double e_k = info->e[k];
        double f_k = info->f[k];
        double Gik = mat_get(info->G, info->order, node_idx, k);
        double Bik = mat_get(info->B, info->order, node_idx, k);
        H_ii += Gik * f_k + Bik * e_k;
        N_ii += Gik * e_k - Bik * f_k;
    }

    Gii = mat_get(info->G, info->order, node_idx, node_idx);
    Bii = mat_get(info->B, info->order, node_idx, node_idx);
    J_ii = N_ii - Bii * f_i - Gii * e_i;
    L_ii = -H_ii + Gii * f_i - Bii * e_i;

    H_ii += -Bii * e_i + Gii * f_i;
    N_ii += Gii * e_i + Bii * f_i;

    for (col = 0; col < ctx->non_slack; ++col)
    {
        int node_j = ctx->var_indices[col];
        if (node_idx == node_j)
        {
            row1[col * 2] = H_ii;
            row1[col * 2 + 1] = N_ii;
            if (info->nodes[node_idx].type == NODE_PV)
            {
                row2[col * 2] = 2.0 * f_i;
                row2[col * 2 + 1] = 2.0 * e_i;
            }
            else
            {
                row2[col * 2] = J_ii;
                row2[col * 2 + 1] = L_ii;
            }
        }
        else
        {
            double Gij = mat_get(info->G, info->order, node_idx, node_j);
            double Bij = mat_get(info->B, info->order, node_idx, node_j);
            double H_ij = -Bij * e_i + Gij * f_i;
            double N_ij = Gij * e_i + Bij * f_i;
            double J_ij = -N_ij;
            double L_ij = H_ij;

            row1[col * 2] = H_ij;
            row1[col * 2 + 1] = N_ij;
            if (info->nodes[node_idx].type == NODE_PV)
            {
                row2[col * 2] = 0.0;
                row2[col * 2 + 1] = 0.0;
            }
            else
            {
                row2[col * 2] = J_ij;
                row2[col * 2 + 1] = L_ij;
            }
        }
    }
}

/// <summary>
/// 生成雅可比矩阵
/// </summary>
/// <param name="ctx">牛顿-拉夫逊迭代上下文</param>
static void gen_jacobian(NLIteration* ctx)
{
    int row;

    if (!ctx || !ctx->J) return;
    memset(ctx->J, 0, (size_t)ctx->eq_count * (size_t)ctx->eq_count
        * sizeof(double));
    for (row = 0; row < ctx->non_slack; ++row)
    {
        double* row1 = ctx->J + (size_t)(2 * row) * (size_t)ctx->eq_count;
        double* row2 = ctx->J + (size_t)(2 * row + 1) * (size_t)ctx->eq_count;
        fill_row_pair(ctx, row, row1, row2);
    }
}

/// <summary>
/// 求解线性系统
/// </summary>
/// <param name="A">系数矩阵</param>
/// <param name="b">常数项向量</param>
/// <param name="x">解向量</param>
/// <param name="n">方程个数</param>
static void solve_linear_system(const double* A, const double* b, double* x, int n)
{
    int i, j, col;

    if (!A || !b || !x || n <= 0) return;
    double* aug = (double*)calloc((size_t)n * (size_t)(n + 1), sizeof(double));
    ensure_alloc(aug, "augmented matrix");
    if (!aug) return;

    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
            aug[i * (n + 1) + j] = A[i * n + j];
        aug[i * (n + 1) + n] = b[i];
    }

    for (col = 0; col < n; ++col)
    {
        int pivot = col;
        double max_val = fabs(aug[col * (n + 1) + col]);
        int row;

        for (row = col + 1; row < n; ++row)
        {
            double val = fabs(aug[row * (n + 1) + col]);
            if (val > max_val)
            {
                max_val = val;
                pivot = row;
            }
        }
        if (fabs(max_val) < 1e-12)
        {
            fprintf(stderr, "Jacobian matrix is singular\n");
            free(aug);
            exit(EXIT_FAILURE);
        }

        if (pivot != col)
        {
            int k;
            for (k = col; k <= n; ++k)
            {
                double tmp = aug[col * (n + 1) + k];
                aug[col * (n + 1) + k] = aug[pivot * (n + 1) + k];
                aug[pivot * (n + 1) + k] = tmp;
            }
        }

        double pivot_val = aug[col * (n + 1) + col];
        int k;
        for (k = col; k <= n; ++k)
            aug[col * (n + 1) + k] /= pivot_val;
        for (row = 0; row < n; ++row)
        {
            if (row == col) continue;
            double factor = aug[row * (n + 1) + col];
            for (k = col; k <= n; ++k)
                aug[row * (n + 1) + k] -= factor * aug[col * (n + 1) + k];
        }
    }

    for (i = 0; i < n; ++i)
        x[i] = aug[i * (n + 1) + n];
    if (aug) free(aug);
}

/// <summary>
/// 应用校正
/// </summary>
/// <param name="ctx">牛顿-拉夫逊迭代上下文</param>
/// <param name="correction">校正向量</param>
/// <returns>0: 成功, -1: 失败</returns>
static int apply_correction(NLIteration* ctx, const double* correction)
{
    int i;
    double max_abs = 0.0;

    if (!ctx || !correction) return 1;
    for (i = 0; i < ctx->eq_count; ++i)
    {
        double val = fabs(correction[i]);
        if (val > max_abs)
            max_abs = val;
    }
    if (max_abs <= ctx->tolerance) return 1;

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

    while (!converged && iteration < ctx->max_iter)
    {
        calc_delta(ctx);
        gen_jacobian(ctx);
        solve_linear_system(ctx->J, ctx->delta, correction, ctx->eq_count);
        converged = apply_correction(ctx, correction);
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
/// 输出潮流计算结果，输出到控制台
/// </summary>
/// <param name="info">网络信息</param>
/// <param name="lines">线路参数</param>
/// <param name="line_count">线路参数数量</param>
/// <param name="converged">是否收敛</param>
/// <param name="runtime">运行时间</param>
static void print_report(const NetworkInfo* info,
    const LineArg* lines,
    size_t line_count,
    int converged,
    double runtime)
{
    const double baseMVA = 100.0; /* 可根据需要修改或从文件读取 */
    int i;
    double totalPg = 0.0, totalQg = 0.0;
    double totalPl = 0.0, totalQl = 0.0;

    if (!info) return;

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

    /* 这里简单用“有发电”的节点个数当作机组数 */
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

    double pLoss = (totalPg - totalPl) * baseMVA;
    double qLoss = (totalQg - totalQl) * baseMVA;
    printf("%% Shunts             0     Shunt (inj)              0.00              0.00\n");
    printf("%% Branches        %4zu     Losses (I^2 * Z)      %11.2f       %11.2f\n\n",
        line_count, pLoss, qLoss);

    printf("%% ===============================================================================\n");
    printf("%% |     Bus Data                                                                 |\n");
    printf("%% ===============================================================================\n");
    printf("%%  Bus      Voltage          Generation             Load        \n");
    printf("%%   #   Mag(pu) Ang(deg)   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)\n");
    printf("%% ----- ------- --------  --------  --------  --------  --------\n");

    for (i = 0; i < info->order; ++i)
    {
        double e = info->e[i]; /* 电压为虚部 */
        double f = info->f[i]; /* 功率为实部 */
        double mag = sqrt(e * e + f * f); /* 功率转电压 */
        double ang = atan2(f, e) * 180.0 / PI; /* 弧度转角度 */
        printf("%% %5d %7.3f %8.3f  %8.2f  %8.2f  %8.2f  %8.2f\n",
            i + 1,
            mag,
            ang,
            info->nodes[i].Pg * baseMVA,
            info->nodes[i].Qg * baseMVA,
            info->nodes[i].Pl * baseMVA,
            info->nodes[i].Ql * baseMVA);
    }

    /* Bus 部分 Total 行 */
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

    if (lines && line_count > 0)
    {
        size_t k;
        double sumPloss = 0.0, sumQloss = 0.0;
        for (k = 0; k < line_count; ++k)
        {
            int fb = lines[k].from_bus - 1;
            int tb = lines[k].to_bus - 1;
            double G = lines[k].series_R;   /* 此时已是导纳实部 */
            double B = -lines[k].series_X;  /* 导纳虚部为 -series_X */
            double b2 = lines[k].shunt_B;   /* 对端 shunt b/2，标幺 */

            double Vi_re = info->e[fb], Vi_im = info->f[fb];
            double Vj_re = info->e[tb], Vj_im = info->f[tb];

            double dV_re = Vi_re - Vj_re;
            double dV_im = Vi_im - Vj_im;

            /* 从 i 节点到 j 节点的串联电流：(G + jB)*(Vi - Vj) */
            double Iser_ij_re = G * dV_re - B * dV_im;
            double Iser_ij_im = G * dV_im + B * dV_re;
            /* i 节点的并联支路电流：j*b2*Vi */
            double Ish_ij_re = -b2 * Vi_im;
            double Ish_ij_im = b2 * Vi_re;

            double Iij_re = Iser_ij_re + Ish_ij_re;
            double Iij_im = Iser_ij_im + Ish_ij_im;

            /* i 节点总电流 = 串联电流 + 并联电流 */
            double Sij_re = Vi_re * Iij_re + Vi_im * Iij_im;
            double Sij_im = Vi_im * Iij_re - Vi_re * Iij_im;

            /* 从 j 节点到 i 节点的计算 */
            dV_re = Vj_re - Vi_re;
            dV_im = Vj_im - Vi_im;
            double Iser_ji_re = G * dV_re - B * dV_im;
            double Iser_ji_im = G * dV_im + B * dV_re;
            double Ish_ji_re = -b2 * Vj_im;
            double Ish_ji_im = b2 * Vj_re;
            double Iji_re = Iser_ji_re + Ish_ji_re;
            double Iji_im = Iser_ji_im + Ish_ji_im;
            double Sji_re = Vj_re * Iji_re + Vj_im * Iji_im;
            double Sji_im = Vj_im * Iji_re - Vj_re * Iji_im;

            double Ploss = (Sij_re + Sji_re) * baseMVA;
            double Qloss = (Sij_im + Sji_im) * baseMVA;

            sumPloss += Ploss;
            sumQloss += Qloss;

            printf("%% %5zu  %5d  %5d  %8.2f  %8.2f  %8.2f  %8.2f  %8.3f  %8.3f\n",
                k + 1,
                lines[k].from_bus,
                lines[k].to_bus,
                Sij_re * baseMVA,
                Sij_im * baseMVA,
                Sji_re * baseMVA,
                Sji_im * baseMVA,
                Ploss,
                Qloss);
        }

        /* Branch 部分 Total 行 */
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

    start_iteration(&iter, 1e-4, 1000, &converged, &iterations, &runtime);

    /* 按课程给的 printpf 样式输出结果（输出到控制台） */
    print_report(&info, line_args, line_count, converged, runtime);

    free_iteration(&iter);
    free_network(&info);
    free(line_args);
    free(node_args);
    free(init_vals);

    return 0;
}