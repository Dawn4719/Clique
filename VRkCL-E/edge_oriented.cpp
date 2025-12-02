#include <set>
#include <algorithm>
#include <unordered_map>
#include <omp.h>
#include <cassert>
#include <chrono>
#include <fstream>
#include <map>
#include "set_operation.h"
#include "edge_oriented.h"
#include "truss/util/graph/graph.h"
#include "truss/util/log/log.h"
#include "truss/util/timer.h"
#include "truss/decompose/parallel_all_edge_cnc.h"
#include "truss/util/reordering/reorder_utils.h"
#include "truss/decompose/iter_helper.h"
#include <future>
#include <thread>
#include <atomic>
#include <queue>
#include <mutex>
#include <sys/time.h>
extern const int K, L;
extern unsigned long long N;
std::map<int, std::map<int, int>> mp;
EBBkC_Graph_t::EBBkC_Graph_t() = default;

EBBkC_Graph_t::~EBBkC_Graph_t() {
    int i, k = K, node_size = v_size, link_size = e_size;

    if (is_sub) {
        k = K - 2;
        node_size = truss_num;
        link_size = truss_num * (truss_num - 1) / 2;
    }

    // delete [] edges;

    if (T) {
        for (i = 0; i < link_size; i++) delete [] T[i];
        delete [] T;
    }

    if (C) {
        for (i = 0; i < link_size; i++) delete [] C[i];
        delete [] C;
    }

    delete [] T_size;

    delete [] C_size;

    if (sub_v) {
        for (i = 0; i <= k; i++) delete [] sub_v[i];
        delete [] sub_v;
    }

    if (sub_e) {
        for (i = 0; i <= k; i++) delete [] sub_e[i];
        delete [] sub_e;
    }

    delete [] sub_v_size;

    delete [] sub_e_size;

    delete [] lab;

    if (DAG_deg) {
        for (i = 0; i <= k; i++) delete [] DAG_deg[i];
        delete [] DAG_deg;

    }

    if (G_deg) {
        for (i = 0; i <= k; i++) delete [] G_deg[i];
        delete [] G_deg;
    }


    delete [] col;

    if (DAG_adj) {
        for (i = 0; i < node_size; i++) delete [] DAG_adj[i];
        delete [] DAG_adj;
    }

    if (G_adj) {
        for (i = 0; i < node_size; i++) delete [] G_adj[i];
        delete [] G_adj;
    }

    if (used) {
        for (i = 0; i <= k; i++) delete [] used[i];
        delete [] used;
    }


    // delete [] v_lab;

    // if (out_v_size) {
    //     for (i = 0; i <= k; i++) delete [] out_v_size[i];
    //     delete [] out_v_size;
    // }
    //
    // if (out_e_size) {
    //     for (i = 0; i <= k; i++) delete [] out_e_size[i];
    //     delete [] out_e_size;
    // }
    //
    // delete [] F;
    //
    // delete [] P;
    //
    // delete [] lack_size;
    //
    // if (lack) {
    //     for (i = 0; i < node_size; i++) delete [] lack[i];
    //     delete [] lack;
    // }
    //
    // delete [] lev;
    //
    // delete [] loc;
}

void EBBkC_Graph_t::build(bool sub) {
    int i, k = K, node_size = v_size, link_size = e_size;

    is_sub = sub;

    if (sub) {
        k = K - 2;
        node_size = truss_num;
        link_size = (truss_num) * (truss_num - 1) / 2;
    }

    sub_v = new int* [k + 1];

    sub_e = new int* [k + 1];

    canAdd = new int* [k + 1];

    candidate = new int* [k + 1];

    // alreadyChoose = new int* [k + 1];

    sub_e_size = new int [k + 1];

    sub_v_size = new int [k + 1];

    canAddSize = new int [k + 1];

    candidateSize = new int [k + 1];

    // alreadyChooseSize = new int [k + 1];

    for (i = 0; i < k; i++) sub_v[i] = new int [truss_num + 1];
    sub_v[k] = new int [node_size];

    for (i = 0; i < k; i++) sub_e[i] = new int [truss_num * (truss_num - 1) / 2];
    sub_e[k] = new int [link_size];

    for(i = 0; i < k; i++) canAdd[i] = new int[truss_num + 1];
    canAdd[k] = new int [node_size];

    for(i = 0; i < k; i++) candidate[i] = new int[truss_num + 1];
    candidate[k] = new int [node_size];

    // for(i = 0; i < k; i++) alreadyChoose[i] = new int[truss_num + 1];
    // alreadyChoose[k] = new int[node_size];

    sub_v_size[k] = 0;

    sub_e_size[k] = 0;

    lab = new int [node_size];
    for (i = 0; i < node_size; i++) lab[i] = k;

    DAG_deg = new int* [k + 1];
    for (i = 0; i <= k; i++) DAG_deg[i] = new int [node_size];

    G_deg = new int* [k + 1];
    for (i = 0; i <= k; i++) G_deg[i] = new int [node_size];

    col = new int [node_size];

    DAG_adj = new int* [node_size];
    for (i = 0; i < node_size; i++) DAG_adj[i] = new int [truss_num + 1];

    G_adj = new int* [node_size];
    for (i = 0; i < node_size; i++) G_adj[i] = new int  [truss_num + 1];

    used = new bool* [k + 1];
    for (i = 0; i <= k; i++) used[i] = new bool [node_size + 1]();

    // v_lab = new int [node_size];
    // for (i = 0; i < node_size; i++) v_lab[i] = k;

    // out_v_size = new int* [k + 1];
    // for (i = 0; i <= k; i++) out_v_size[i] = new int [link_size];
    //
    // out_e_size = new int* [k + 1];
    // for (i = 0; i <= k; i++) out_e_size[i] = new int [link_size];
    //
    // F = new int [truss_num + 1];
    //
    // P = new int [truss_num + 1];
    //
    // lack_size = new int [node_size];
    //
    // lack = new int* [node_size];
    // for (i = 0; i < node_size; i++) lack[i] = new int [L + 1];
    //
    // lev = new int [node_size]();
    //
    // loc = new int [node_size];
}

void EBBkC_Graph_t::branch(int e, EBBkC_Graph_t* g) {
    // if (C_size[e] < K - 2) {
    //     g->v_size = 0;
    //     return;
    // }
    int c, i, j, k, p, e_, u, v, w, s, t, end, dist, l = K;
    int *old2new = new int[v_size];

    g->v_size = 0;
    g->e_size = 0;
    g->edges = new Edge_t[C_size[e] + 1];
    g->new2old = vector<int>(T_size[e]);

    for (i = 0; i < T_size[e]; i++) {
        u = T[e][i];
        old2new[u] = g->v_size;
        g->new2old[g->v_size++] = u;
    }

    for (i = 0; i < C_size[e]; i++) {
        e_ = C[e][i];
        s = edges[e_].s;
        t = edges[e_].t;
        g->edges[g->e_size].s = old2new[s];
        g->edges[g->e_size++].t = old2new[t];
    }

    delete [] old2new;

    if (l == 3 || l == 4) {
        g->sub_v_size[l - 2] = g->v_size;
        g->sub_e_size[l - 2] = g->e_size;
        return;
    }

    if (g->v_size < l - 2 || g->e_size < (l - 2) * (l - 3) / 2) {
        g->sub_v_size[l - 2] = g->v_size;
        g->sub_e_size[l - 2] = g->e_size;
        return;
    }

    for (j = 0; j < g->v_size; j++) {
        g->col[j] = 0;
        g->DAG_deg[l - 2][j] = 0;
        g->G_deg[l - 2][j] = 0;
    }

    for (j = 0; j < g->e_size; j++) {
        s = g->edges[j].s;
        t = g->edges[j].t;
        g->G_adj[s][g->G_deg[l - 2][s]++] = t;
        g->G_adj[t][g->G_deg[l - 2][t]++] = s;
    }

    auto *list = new KeyVal_t[truss_num + 1];
    for (j = 0; j < g->v_size; j++) {
        list[j].key = j;
        list[j].val = g->G_deg[l - 2][j];
    }
    sort(list, list + g->v_size);

    for (j = 0; j < g->v_size; j++) {
        u = list[j].key;
        for (k = 0; k < g->G_deg[l - 2][u]; k++) {
            v = g->G_adj[u][k];
            g->used[l - 2][g->col[v]] = true;
        }
        for (c = 1; g->used[l - 2][c]; c++);
        g->col[u] = c;
        if (c > g->Col) g->Col = c;
        for (k = 0; k < g->G_deg[l - 2][u]; k++) {
            v = g->G_adj[u][k];
            g->used[l - 2][g->col[v]] = false;
        }
    }
    delete[] list;

    g->sub_v_size[l - 2] = 0;

    for (j = 0; j < g->v_size; j++) {
        g->sub_v[l - 2][g->sub_v_size[l - 2]++] = j;
    }

    g->sub_e_size[l - 2] = 0;
    for (j = 0; j < g->e_size; j++) {
        g->sub_e[l - 2][g->sub_e_size[l - 2]++] = j;
        s = g->edges[j].s;
        t = g->edges[j].t;
        g->edges[j].s = (g->col[s] > g->col[t]) ? s : t;
        g->edges[j].t = (g->col[s] > g->col[t]) ? t : s;
        s = g->edges[j].s;
        t = g->edges[j].t;

        g->DAG_adj[s][g->DAG_deg[l - 2][s]++] = t;
    }

    return;

}

timeval time_start;
timeval time_end;
map<string, size_t> get_index_mem() {
    FILE *fp = fopen("/proc/self/status", "r");
    char line[128];
    map<string, size_t> res;
    while (fgets(line, 128, fp) != NULL) {
        //        if (strncmp(line, "VmPeak", 2) == 0)
        //        {
        //            cout << line << endl;
        ////            printf("当前进程占用虚拟内存大小为：%d KB\n", atoi(line + 6));
        //        }
        if (strncmp(line, "VmRSS:", 6) == 0) {
            string p = line;
            res["now"] = size_t(stoull(p.substr(6)));
            // cout << line;
        }
        if (strncmp(line, "VmPeak:", 7) == 0) {
            string p = line;
            res["pk"] = size_t(stoull(p.substr(7)));
            // cout << line;
        }
    }
    fclose(fp);
    return res;
}
#define Get_T() std::chrono::high_resolution_clock::now();

// auto ti = Get_T();
size_t origin_graphsize;
void EBBkC_Graph_t::truss_decompose(const char *dir)  {
    graph_t g;

    //load the graph from file
    auto pk1 = get_index_mem()["pk"];
    Graph G(dir, K);
    cout << g.m << " " << K << endl;
    // origin_graphsize = get_index_mem()["now"] - pk1;
    // cout << "--------------------------------------" << get_index_mem()["now"] - pk1 << " " << get_index_mem()["pk"] - pk1 << endl;
    gettimeofday(&time_start, NULL);
    // ti = Get_T();
    g.adj = G.edge_dst;
    g.num_edges = G.node_off;
    g.n = G.nodemax;
    g.m = G.edgemax;

    std::cout << "|V| " << g.n << " |E| " << g.m << std::endl;

    string reorder_method("core");

    vector <int32_t> new_vid_dict;
    vector <int32_t> old_vid_dict;
    ReorderWrapper(g, dir, reorder_method, new_vid_dict, old_vid_dict);

    //edge list array
    Timer get_eid_timer;

    Edge *edgeIdToEdge = (Edge *) malloc((g.m / 2) * sizeof(Edge));
    assert(edgeIdToEdge != nullptr);
    log_info("Malloc Time: %.9lf s", get_eid_timer.elapsed());
    get_eid_timer.reset();

    //Populate the edge list array
    getEidAndEdgeList(&g, edgeIdToEdge);
    log_info("Init Eid Time: %.9lf s", get_eid_timer.elapsed());
    get_eid_timer.reset();

    int *EdgeSupport = (int *) malloc(g.m / 2 * sizeof(int));
    assert(EdgeSupport != nullptr);

    auto max_omp_threads = 1;
    // omp_get_max_threads();
    omp_set_num_threads(max_omp_threads);
    log_info("Max Threads: %d", max_omp_threads);
#pragma omp parallel for
    for (auto i = 0; i < max_omp_threads; i++) {
        auto avg = g.m / 2 / max_omp_threads;
        auto iter_beg = avg * i;
        auto iter_end = (i == max_omp_threads - 1) ? g.m / 2 : avg * (i + 1);
        memset(EdgeSupport + iter_beg, 0, (iter_end - iter_beg) * sizeof(int));
    }
    log_info("Init EdgeSupport Time: %.9lf s", get_eid_timer.elapsed());
    get_eid_timer.reset();

    Timer global_timer;
    truss_num = PKT_intersection(&g, EdgeSupport, edgeIdToEdge);

#pragma omp single
    {
        int u, v, w;
        int *old2new = new int[N_NODES];
        for (int i = 0; i < N_NODES; i++) old2new[i] = -1;

        e_size = 0;
        edges = new Edge_t[g.m / 2];

        T = new int *[g.m / 2];
        T_size = new int[g.m / 2];

        for (int i = 0; i < g.m / 2; i++) {
            if (g.edge_truss[i] <= K) continue;

            Edge e = edgeIdToEdge[i];
            u = e.u;
            v = e.v;

            if (old2new[u] == -1) {
                old2new[u] = (int) new2old.size();
                new2old.push_back(u);
            }

            if (old2new[v] == -1) {
                old2new[v] = (int) new2old.size();
                new2old.push_back(v);
            }

            edges[e_size] = Edge_t(old2new[u], old2new[v], false);
            edge2id.insert(edges[e_size], e_size);
            rank.push_back(g.edge_rank[i]);

            int sz = g.v_set[i].size();
            T_size[e_size] = sz;

            T[e_size] = new int[truss_num + 1];
            for (int j = 0; j < sz; j++) {
                int w = g.v_set[i][j];

                if (old2new[w] == -1) {
                    old2new[w] = (int) new2old.size();
                    new2old.push_back(w);
                }

                T[e_size][j] = old2new[w];
            }

            e_size++;
        }

        C = new int *[e_size];
        C_size = new int[e_size];

        for (int i = 0; i < e_size; i++) {

            C_size[i] = 0;
            int sz = T_size[i];
            C[i] = new int[sz * (sz - 1) / 2];

            for (int j = 0; j < T_size[i]; j++) {
                for (int k = j + 1; k < T_size[i]; k++) {

                    Edge_t e = Edge_t(T[i][j], T[i][k], false);
                    int idx = edge2id.exist(e);

                    if (idx != -1 && rank[idx] > rank[i]) {
                        C[i][C_size[i]++] = idx;
                    }
                }
            }
        }
    };

    v_size = new2old.size();

    printf("|V| = %d, |E| = %d\n", v_size, e_size);
    printf("Truss number = %d\n", truss_num - 2);

    //Free memory
    free_graph(&g);
    free(edgeIdToEdge);
    free(EdgeSupport);
}

int dfs_count;
vector<int> comblist;

inline void EBBkC_Comb_list2(int list_size, int start, int picked, int k, int m, unsigned long long* cliques) {
    // return;
    int l1 = list_size - (k - picked), l2, l3, i, j, q, w;
    if (picked == k - 2) {
        l2 = l1 + 1;
            for (i = start; i <= l1; i++) {
                for (j = i + 1; j <= l2; j++)
                    (*cliques) += m;
            }
        return;
    }

    for (i = start; i <= l1; i++) {
        EBBkC_Comb_list2(list_size, i + 1, picked + 1, k, m, cliques);
    }
}

inline void EBBkC_Comb_list(int *list, int list_size, int start, int picked, int k, unsigned long long* cliques) {
    int l1 = list_size - (k - picked), l2, l3, i, j, q, w;
    if (picked == k) {
        (*cliques)++;
        return;
    }
    else
    if (picked == k - 1) {
            for (i = start; i <= l1; i++)
                (*cliques)++;
        return;
    }
    else
    if (picked == k - 2) {
        l2 = l1 + 1;
        for (i = start; i <= l1; i++) {
            for (j = i + 1; j <= l2; j++)
                (*cliques)++;
        }
        return;
    }

    for (i = start; i <= l1; i++) {
        EBBkC_Comb_list(list, list_size, i + 1, picked + 1, k, cliques);
    }
}

using ull = unsigned long long;

vector<vector<int>> his;
vector<int> IDX;

void PT() {
    fstream fs("../tt.txt", ios::app);
    for (auto i : his) {
        if (i[0] == -2) {
            fs << "(";
            for (int j = 1; j < i.size(); j++) {
                fs << i[j];
                if (j + 1 < i.size())
                    fs << " ";
            }
            fs << ") ";
        }
        else if (i[0] == -1) {
            fs << "<";
            for (int j = 1; j < i.size(); j++) {
                fs << i[j];
                if (j + 1 < i.size())
                    fs << " ";
            }
            fs << "> ";
        }
        else {
            fs << "(2222)";
        }
    }
    fs << endl;
    fs.close();
}

map<int, int> DEPTH;
int depp;

struct Node {
    int id=-1;
    int fa=-1;
    vector<int> val{};
    vector<int> sons{};

    Node(int _id, int _fa, vector<int> _val): id(_id), fa(_fa) {
        val = _val;
    }

    bool operator==(const Node& other) const {
        return id == other.id &&
               fa == other.fa &&
               val == other.val &&
               sons == other.sons;
    }
};

struct searchTree{
    int cnt = 0;
    uint64_t hashVal = 0;
    // vector<vector<vector<int>>> node;
    vector<Node> tree;
    vector<string> nodesets;

    inline void init() {
        cnt = 0;
        hashVal = 0;
        if (!tree.empty()) tree.clear();
        nodesets.clear();
    }

    void increase(vector<string>& v) {
        cnt++;
        // for (int i = 0; i < v.size(); i++) {
        //     nodesets[i] += " " + v[i];
        // }
    }

    void storeV() {
        // for (auto i : tree) {
        //     for (int j = 1; j < i.val.size(); j++) {
        //         nodesets.emplace_back(to_string(i.val[j]));
        //     }
        // }
    }

    bool operator==(const searchTree& o) const {
        bool same = true;
        for (int i = 0; i < tree.size(); i++) {
            if (tree[i].val[0] == o.tree[i].val[0] && tree[i].val.size() == o.tree[i].val.size()) continue;
            same = false;
        }
        return same;
    }

    bool operator!=(const searchTree& o) const {
        for (int i = 0; i < tree.size(); i++) {
            if (tree[i].val[0] != o.tree[i].val[0] || tree[i].val.size() != o.tree[i].val.size()) return true;
        }
        return false;
    }
};

searchTree st;
int treeIdx;

bool l2;

// bool skip;
int maxdepth;
bool hasAdd;

typedef unsigned long long ULL;

const int P = 131;
vector<ULL> h, p;

ULL hash(string& s)
{
    p[0] = 1;
    for (int i = 1; i <= s.size(); i ++ )
    {
        h[i] = h[i - 1] * P + s[i - 1];
        p[i] = p[i - 1] * P;
    }
    return h[s.size()];
}

inline uint64_t fnv1a_64(const std::string& s) {
    uint64_t hash = 1469598103934665603ULL;      // offset basis
    for (unsigned char c : s) {
        hash ^= c;
        hash *= 1099511628211ULL;                // FNV prime
    }
    return hash;
}

inline void cal(int l, int lz, int depp) {
    // int l = K - 2, lz = 0;
    if (depp == 0) {
        l = K - 2, lz = 0;
    }
    for (int i = depp; i < his.size(); i++) {
        if (l == 2) {
            if (l <= lz)
                if (l == 0 || l == lz)
                    N++;
                else if (l == 1 || l + 1 == lz)
                    for (int q = 0; q < lz; q++)
                        N++;
                else
                    EBBkC_Comb_list2(lz, 0, 0,  min(l, lz - l), 1, &N);
            break;
        }
        if (his[i][0] == -1) {
            lz += his[i].size() - 1;
        }
        else {
            l--;
        }
        if (l == 0) {
            N++;
            break;
        }

        if (l <= lz) {
            if (l == 0 || l == lz)
                N++;
            else if (l == 1 || l + 1 == lz)
                for (int q = 0; q < lz; q++)
                    N++;
            else
                EBBkC_Comb_list2(lz, 0, 0,  min(l, lz - l), 1, &N);
        }

        // if (his[i][0] == -1 && i + 1 == his.size()) {
        //     break;
        // }
    }
}
// 
string code;
string encode(int idx) {
    string tmp;
    vector<pair<int, int>> v;
    for (auto i : st.tree[idx].sons) {
        if (st.tree[i].val[0] == -1) v.emplace_back(st.tree[i].val.size() - 1, i);
        else v.emplace_back(0, i);
    }
    sort(v.begin(), v.end());

    for (auto i : v) {
        tmp += to_string(i.first);
        if (st.tree[i.second].sons.empty()) {
            continue;
        }
        tmp += '(';
        tmp += encode(i.second);
        tmp += ')';
    }

    return tmp;
}

int treeBranchCnt;
string vertexLabel;

map<uint64_t, int> vertexHash;
vector<searchTree> ST;
int Maxdepth;

inline void search(const int cnt, const vector<Node>& tree, int idx, int l, int lz, unsigned long long * cliques) {
    // if (l == 2) {
    //     if (lz >= 2) {
    //         // auto ti = Get_T();
    //
    //         if (l <= lz)
    //             if (l == 0 || l == lz)
    //                 N += cnt;
    //             else if (l == 1 || l + 1 == lz)
    //                 for (int j = 0; j < lz; j++)
    //                     N += cnt;
    //             else
    //                 EBBkC_Comb_list2(lz, 0, 0,  min(l, lz - l), cnt, &N);
    //
    //         // time2 += Duration(ti);
    //     }
    //     return;
    // }

    if (tree[idx].val[0] == -1) {
        lz += tree[idx].val.size() - 1;
    }

    if (l == 0) {
        N += cnt;
        if (tree[idx].val[0] == -1) lz -= tree[idx].val.size() - 1;
        // cout << "++clique " << N << endl;
        return;
    }

    if (tree[idx].val[0] == -1 && tree[idx].sons.empty()) {
        // auto ti = Get_T();
        if (l <= lz)
            if (l == 0 || l == lz)
                N += cnt;
            else if (l == 1 || l + 1 == lz)
                for (int j = 0; j < lz; j++) N += cnt;
            else
                EBBkC_Comb_list2(lz, 0, 0,  min(l, lz - l), cnt, &N);

        // cout << "++comblist candidateSize[l] == 0 " << lz << " " << min(l, lz - l) << " " << N << endl;
        // time2 += Duration(ti);
        if (tree[idx].val[0] == -1) lz -= tree[idx].val.size() - 1;
        // cout << "candidate == 0 " << N << endl;
        return;
    }

    if (l <= lz) {
        // auto ti = Get_T();
            if (l == 0 || l == lz)
                N += cnt;
            else if (l == 1 || l + 1 == lz)
                for (int j = 0; j < lz; j++) N += cnt;
            else
                EBBkC_Comb_list2(lz, 0, 0,  min(l, lz - l), cnt, &N);
        // time2 += Duration(ti);
    }

    for (auto node : tree[idx].sons) {
        // if (tree[node].val[0] == -2) {

        if (idx == 0 && tree[node].val[0] == -1) {
            search(cnt, tree, node, l, lz, cliques);
        }
        else if (tree[node].sons.size() == 1 && tree[tree[node].sons[0]].val[0] == -1) {
            // numberOfVertex += cnt;
            // numberOfNode += cnt;
            search(cnt, tree, tree[node].sons[0], l - 1, lz, cliques);
        }
        else {
            search(cnt, tree, node, l - 1, lz, cliques);
        }
    }
    if (tree[idx].val[0] == -1) lz -= tree[idx].val.size() - 1;
}

int lz_his;
int minlz;

inline void add(int l, int lz) {
    // PT();
    hasAdd = true;

    // for (const auto& ii : his) {
    //     for (auto jj : ii)
    //         cout << jj << " ";
    //     cout << " - ";
    // }
    // cout << endl;
    if (his.size() == 1) return;
    if (depp == 0 && !st.tree.empty()) {
        // assert(maxdepth <= 4);
        lz_his = lz;
        if (maxdepth > Maxdepth || lz < minlz) {
            // if (st.tree[0].val[0] == -1)
            search(1, st.tree, 0, K - 2 - (st.tree[0].val[0] == -2), 0, &N);
        }
        else {
            string str;
            if (st.tree[0].val[0] == -1) {
                str = to_string(st.tree[0].val.size() - 1) + "(";
            }
            else {
                str = "0(";
            }
            str += encode(0) + ")";
            auto hashVal = fnv1a_64(str);

            // if (st.tree[0].val[0] == -2)
            // search(1, st.tree, 0, K - 2 - (st.tree[0].val[0] == -2), 0);
            // else {
            st.hashVal = hashVal;

            if (!vertexHash.count(hashVal)) {
                vertexHash[hashVal] = ST.size();
                // st.storeV();
                // st.cnt = 1;
                ST.emplace_back(st);
            }
            else {
                // st.storeV();
                ST[vertexHash[hashVal]].increase(st.nodesets);
            }
        }
        st.init();
        maxdepth = 0;
        lz_his = 0;
        // cout << N << endl;
    }

    // assert(his.back()[0] == -1);
    if (!st.tree.empty()) {
        for (int i = depp; i < his.size(); i++) {
            assert(!his[i].empty());
            st.tree[IDX.back()].sons.emplace_back(st.tree.size());
            st.tree.emplace_back(st.tree.size(), IDX.back(), his[i]);
            IDX.emplace_back(st.tree.size() - 1);
        }
        // IDX.pop_back();
    }
    else {
        st.tree.emplace_back(0, -1, his[0]);
        IDX.emplace_back(0);

        for (int i = 1; i < his.size(); i++) {
            st.tree[i - 1].sons.emplace_back(i);
            IDX.emplace_back(i);
            st.tree.emplace_back(st.tree.size(), st.tree.size() - 1, his[i]);
        }
    }

    if (his.size() > maxdepth) maxdepth = his.size();

    depp = his.size();
    DEPTH[depp]++;
}

void EBBkC_Graph_t::EBBkC_plus_plus(int l, unsigned long long* cliques, int& lz) {
    // cout << "l = " << l << endl;
    int c, cc, i, j, k, p, e, e_, u, v, w, s, t, end, dist;
    /*sub_v_size : g->ns[]    sub_e_size: g->es[]*/

    // if (l-lz < 0 || sub_v_size[l] < l - lz || sub_e_size[l] < (l-lz) * (l - lz - 1) / 2) return;
    if (sub_v_size[l] < max(0, l - lz) || sub_e_size[l] < max(0, l - lz) * max(0, l - lz - 1) / 2) return;
    assert(depp >= 0);
    dfs_count++;

    if (K == 3) {
        // cout << epoch << " K = 3" <<endl;
        (*cliques) += sub_v_size[l];
        return;
    }
    if (K == 4) {
        // cout << epoch << " K = 4" <<endl;
        (*cliques) += sub_e_size[l];
        return;
    }
    if (l == 2) {
        // cout << "l == 2" << endl;

        if (lz >= 2) {
            // choose 2 AC
            // if(!skip)
            // if (l < lz - l)
            //     EBBkC_Comb_list(comblist.data(), lz, 0, 0, l, cliques);
            // else
            //     EBBkC_Comb_list(comblist.data(), lz, 0, 0, lz - l, cliques);

            if (l <= lz) {
                if (lz < minlz) {
                    if (l < lz - l) {
                        EBBkC_Comb_list(comblist.data(), lz, 0, 0, l, cliques);
                    }
                    else {
                        EBBkC_Comb_list(comblist.data(), lz, 0, 0, lz - l, cliques);
                    }
                }
            }

            // choose 1 AC
            for (i = 0; i < sub_v_size[l]; i++) {
                u = sub_v[l][i];
                // if (u >= b.size()) b.resize(u + 1);
                // b[u] = '1';
                // if(!skip)
                for (c = 0; c < lz; c++) {
                    (*cliques) ++;
                }
            }
            // bv.emplace_back(b);

            // choose 0 AC
            for (i = 0; i < sub_v_size[l]; i++) {
                u = sub_v[l][i];
                // b.clear();
                for (j = 0; j < DAG_deg[l][u]; j++) {
                    v = DAG_adj[u][j];
                    (*cliques) ++;
                    // if (v >= b.size()) b.resize(v + 1);
                    // b[v] = '1';
                }
                // bv.emplace_back(b);
            }
        }
        else if (lz == 1) {
            // if (sub_v_size[l] == 0) {
            //     // st.node.emplace_back(his);
            //     add(l, lz);
            //     // IDX.pop_back();
            //     // depp--;
            //     return;
            // }

            // choose 1 AC
            for (i = 0; i < sub_v_size[l]; i++) {
                u = sub_v[l][i];
                // if (u >= b.size()) b.resize(u + 1);
                // b[u] = '1';

                for (c = 0; c < lz; c++) {
                    (*cliques) ++;
                }
            }
            // bv.emplace_back(b);

            // choose 0 AC
            for (i = 0; i < sub_v_size[l]; i++) {
                u = sub_v[l][i];
                // b.clear();
                for (j = 0; j < DAG_deg[l][u]; j++) {
                    v = DAG_adj[u][j];
                    (*cliques) ++;
                    // if (v >= b.size()) b.resize(v + 1);
                    // b[v] = '1';
                }
                // bv.emplace_back(b);
            }
        }
        else if (lz == 0) {
            // if (sub_v_size[l] == 0) {
            //     // st.node.emplace_back(his);
            //     add(l, lz);
            //     // IDX.pop_back();
            //     // depp--;
            //     return;
            // }

            // choose 0 AC
            // b.clear();
            // for (i = 0; i < sub_v_size[l]; i++) {
            //     u = sub_v[l][i];
            //     if (u >= b.size()) b.resize(u + 1);
            //     b[u] = '1';
            // }
            // b.resize(b.size() + 1);
            // b.back() = 'x';
            // bv.emplace_back(b);

            for (i = 0; i < sub_v_size[l]; i++) {
                u = sub_v[l][i];
                // b.clear();
                for (j = 0; j < DAG_deg[l][u]; j++) {
                    v = DAG_adj[u][j];
                    (*cliques) ++;
                    // if (v >= b.size()) b.resize(v + 1);
                    // b[v] = '1';
                }
                // bv.emplace_back(b);
            }
        }

        // if (l <= lz) {
        //     if (lz < minlz) {
        //         if (l < lz - l) {
        //             EBBkC_Comb_list(comblist.data(), lz, 0, 0, l, cliques);
        //         }
        //         else {
        //             EBBkC_Comb_list(comblist.data(), lz, 0, 0, lz - l, cliques);
        //         }
        //     }
        //     else {
        if (his.size() > maxdepth) {
            maxdepth = his.size();
        }
        if (lz >= minlz)
            add(l, lz);
    // }
        // }

        // st.bv.emplace_back(bv);
        // assert((int)st.bv.size() > 0);
        // his.emplace_back(vector<int>{(int)st.bv.size() - 1});
        // st.node.emplace_back(his);


        // add(l, lz);
        // his.pop_back();
        // IDX.pop_back();
        // depp--;
        // his.emplace_back(vector<int>{-2, 2222});
        // PT();
        // his.pop_back();
        // DEPTH[depp--]++;
        return;
    }

    canAddSize[l] = candidateSize[l] = 0;
    // assert(l <= K);

    for(i = 0; i < sub_v_size[l]; i++) {
        u = sub_v[l][i];
        if(G_deg[l][u] == sub_v_size[l] - 1) {
            canAdd[l][canAddSize[l]++] = u;
        }
        else {
            candidate[l][candidateSize[l]++] = u;
        }
    }
    int curTreeId;
    // if (l == K - 2 && canAddSize[l] == 0 && !l2) {
    //     his.resize(his.size() + 1);
    //     his.back().emplace_back(-1);
    // }
    if (canAddSize[l]) {
        // if (l == K - 2) {
        //     st.tree.resize(treeIdx + 1);
        //     st.tree.back().val.reserve(canAddSize[l]);
        //     st.tree.back().id = treeIdx++;
        //
        //     for (i = 0; i < canAddSize[l]; i++) st.tree.back().val.emplace_back(canAdd[l][i]);
            // his.emplace_back(vector<int>(canAdd[l], canAdd[l] + canAddSize[l]));
        // }
        // else {
        his.resize(his.size() + 1);
        his.back().reserve(canAddSize[l] + 1);
        his.back().emplace_back(-1);
        for (i = 0; i < canAddSize[l]; i++) his.back().emplace_back(canAdd[l][i]);
        assert(his.back().size() > 1);
        // }
    }

#ifdef DBG
    cout << "AC: ";
    for (i = 0; i < canAddSize[l]; i++) cout << canAdd[l][i] << " ";
    cout << endl;
    cout << "Can: ";
    for (i = 0; i < candidateSize[l]; i++) cout << candidate[l][i] << " ";
    cout << endl;
#endif

    lz += canAddSize[l];
    lz_his = lz;
    if (l == 0) {
        // if (l2) (*cliques)++;
        lz -= canAddSize[l];
        // cout << "++clique" << endl;
#ifdef DBG
        for (i = 0; i < canAddSize[l]; i++)
            comblist.pop_back();
#endif
        // PT();

        if (lz < minlz) {
            (*cliques)++;
        }
        else {
            add(l, lz);
        }

            if (canAddSize[l] && !his.empty()) {
                assert(!his.empty());
                assert(!IDX.empty());
                // if (his.size() == IDX.size()) {
                //     IDX.pop_back();
                    // depp--;
                // }
                his.pop_back();
                if (his.size() < depp) {
                    depp = his.size();
                    IDX.resize(depp);
                }
            }
        // }
        // DEPTH[depp--]++;
        return;
    }

    if (candidateSize[l] == 0) {
        assert(l <= lz);
        // if (l2)
        // if (l < lz - l)
        //     EBBkC_Comb_list(comblist.data(), lz, 0, 0, l, cliques);
        //     // else DFS_work_stealing(lz, l, cliques);
        // else
        //     EBBkC_Comb_list(comblist.data(), lz, 0, 0, lz - l, cliques);

        if (lz < minlz) {
            if (l < lz - l)
                EBBkC_Comb_list(comblist.data(), lz, 0, 0, l, cliques);
            else {
                EBBkC_Comb_list(comblist.data(), lz, 0, 0, lz - l, cliques);
            }
        }
        else {
            add(l, lz);
        }

        // cout << "++comblist candidateSize[l] == 0" << lz << " " << min(l, lz - l) << endl;
        lz -= canAddSize[l];
#ifdef DBG
        for (i = 0; i < canAddSize[l]; i++)
            comblist.pop_back();
#endif
        // PT();
        // if (!l2) {
            // st.node.resize(st.node.size() + 1);
            // st.node.back() = his;

            if (canAddSize[l] && !his.empty()) {
                // if (his.size() == IDX.size()) {
                //     IDX.pop_back();
                    // depp--;
                his.pop_back();
                if (his.size() < depp) {
                    depp = his.size();
                    IDX.resize(depp);
                }
                // }

            }
        // }
        // DEPTH[depp--]++;
        return;
    }
    // if (l2)
    // if (l <= lz) {
    //     if (l < lz - l)
    //         EBBkC_Comb_list(comblist.data(), lz, 0, 0, l, cliques);
    //     else
    //         EBBkC_Comb_list(comblist.data(), lz, 0, 0, lz - l, cliques);
    //     // cout << "++comblist " << lz << " " << min(l, lz - l) << endl;
    // }
    for(i=0; i<candidateSize[l]; i++) {
        u = candidate[l][i];
        // cout << "can: " << u << endl;
        if (col[u] < l - lz) continue;
        sub_v_size[l - 1] = 0;
        dist = 0;

        for (j = 0; j < DAG_deg[l][u]; j++) {
            // assert(u < v_size);
            // assert(j < truss_num);
            v = DAG_adj[u][j];
            // if(G_deg[l][v] == sub_v_size[l] - 1 || G_deg[l][v] == sub_v_size[l] - 2) continue;
            if(G_deg[l][v] == sub_v_size[l] - 1) continue;
            lab[v] = l - 1;

            sub_v[l - 1][sub_v_size[l - 1]++] = v;
            DAG_deg[l - 1][v] = 0;
            G_deg[l - 1][v] = 0;

            if (!used[l][col[v]]) {
                used[l][col[v]] = true;
                dist++;
            }
        }
        if (dist >= l - lz - 1) {
            sub_e_size[l - 1] = 0;
            /*sub_v_size: g->sub[]*/
            for (int jj = 0; jj < sub_v_size[l - 1]; jj++) {
                v = sub_v[l - 1][jj];

                end = DAG_deg[l][v];
                for (k = 0; k < end; k++) {
                    w = DAG_adj[v][k];
                    if (lab[w] == l - 1) {
                        // assert(G_deg[l][w] != sub_v_size[l] - 1);
                        DAG_deg[l - 1][v]++;
                        sub_e_size[l - 1]++;

                        // just for early-termination
                        /*G_deg[l-1][j] 第j个节点的度数，入度+出度*/
                        G_deg[l - 1][v]++;
                        G_deg[l - 1][w]++;

                    } else {
                        DAG_adj[v][k--] = DAG_adj[v][--end];
                        DAG_adj[v][end] = w;
                    }
                }
            }
#ifdef DBG
            cout << "choose: " << u << " " << l << endl;
            res.emplace_back(u);
#endif
            // if (!l2) {
            his.emplace_back(vector<int>{-2, u});
            // }
            // hasAdd = false;

            EBBkC_plus_plus(l - 1, cliques, lz);
            // if (!l2) {
            his.pop_back();
            if (his.size() < depp) {
                depp = his.size();
                IDX.resize(depp);
            // }
            }
#ifdef DBG
            res.pop_back();
            cout << "del: " << u << " " << l << endl;
#endif
        }
        for (j = 0; j < sub_v_size[l - 1]; j++) {
            v = sub_v[l - 1][j];
            lab[v] = l;
            used[l][col[v]] = false;
        }
    }

    if (canAddSize[l] && !his.empty() && !l2) {
        // if (his.size() == IDX.size()) {
        //     IDX.pop_back();
            // depp--;
        // }
        his.pop_back();
        if (his.size() < depp) {
            depp = his.size();
            IDX.resize(depp);
        }
    }

    lz -= canAddSize[l];
#ifdef DBG
    for (i = 0; i < canAddSize[l]; i++)
        comblist.pop_back();
#endif
}

double time2 = 0;

int dbgcnt;
// void getres(const map<int, map<int, map<int, vector<searchTree>>>>& stree) {
// void getres(searchTree& stree) {
//     for (const auto& acCount : stree) {
//         for (const auto& V : acCount.second) {
//             for (const auto& E : V.second) {
//                 for (const auto& _st : E.second) {
//                     search(_st.cnt + 1, _st.tree, 0, K - 2, 0);
//                     // search(st.cnt, st.bv, st.tree, 0, K - 2, 0);
//                 }
//             }
//         }
//     }
// }

// 获取 std::string 的堆内存大小（考虑 SSO）
inline size_t string_heap_memory(const std::string& s) {
    // 典型 SSO 阈值（GCC/Clang 为 15，MSVC 为 16）
    constexpr size_t SSO_THRESHOLD = 15;
    return (s.size() > SSO_THRESHOLD) ? s.capacity() : 0;
}

// 计算单个 searchTree 的内存占用
size_t calculate_searchTree_memory(const searchTree& st) {
    size_t total = sizeof(searchTree); // 对象本身大小（含成员变量）

    // 计算 tree 的内存
    total += st.tree.capacity() * sizeof(Node); // Node 数组
    for (size_t i = 0; i < st.tree.size(); ++i) {
        const Node& node = st.tree[i];
        total += node.val.capacity() * sizeof(int);   // val 的堆内存
        total += node.sons.capacity() * sizeof(int);  // sons 的堆内存
    }

    // 计算 nodesets 的内存
    total += st.nodesets.capacity() * sizeof(std::string); // string 数组
    for (size_t i = 0; i < st.nodesets.size(); ++i) {
        total += string_heap_memory(st.nodesets[i]); // 每个 string 的堆内存
    }

    return total;
}

// 计算 vector<searchTree> 的总内存占用
size_t calculate_vector_searchTree_memory(const std::vector<searchTree>& vec) {
    size_t total = sizeof(vec); // vector 自身开销（24 字节 on 64-bit）
    total += vec.capacity() * sizeof(searchTree); // searchTree 数组

    // 逐个计算每个 searchTree 的动态内存
    for (size_t i = 0; i < vec.size(); ++i) {
        total += calculate_searchTree_memory(vec[i]) - sizeof(searchTree);
    }
    return total;
}
double EBBkC_t::list_k_clique(const char *file_name, long long& dfs_count_) {
    EBBkC_Graph_t G, g;
    struct rusage end;

    printf("Reading edges from %s ...\n", file_name);
    timeval truss_decompose_time_end;

    G.truss_decompose(file_name);
    // auto s2 = std::chrono::high_resolution_clock::now();
    // double truss_decompose_time = std::chrono::duration_cast<std::chrono::microseconds>(s2 - ti).count() / (double)1000;
    gettimeofday(&truss_decompose_time_end, NULL);
    double truss_decompose_tim1e = (truss_decompose_time_end.tv_sec - time_start.tv_sec) * 1000.0 + (truss_decompose_time_end.tv_usec - time_start.tv_usec) / 1000.0;

    printf("Building necessary data structure ...\n");
    G.build(false);

    printf("Iterate over all cliques\n");

    g.truss_num = G.truss_num;

    g.build(true);

    vector<pair<int, int>> C;
    for (int i = 0; i < G.e_size; i++) {
        if (G.C_size[i] * 2 >= (K - 2) * (K - 3)) {
            C.emplace_back(i, G.C_size[i]);
        }
    }

    sort(C.begin(), C.end(), [](const pair<int, int>& a, const pair<int, int>& b) {return a.second > b.second;});


    unsigned long long NN = 0;
    int lz;
    vector<vector<int>> egg;
    map<int, int> mSet;
    // map<int, map<int, map<int, map<int, vector<searchTree>>>>> ST;

    ull N2 = 0;
    cout << C.size() << endl;
    // 1007211
    map<int, int> DEP;

    Maxdepth = 6;
    minlz = 3;
    int allBranch = 0;
    int calBranch = 0;
    map<int, int> CCCNT;
    int itr = 0;

    for (auto & i : C) {
    // for (int i = 0; i < G.e_size; i++) {
        lz = 0;
        G.branch(i.first, &g);
        // G.branch(i, &g);
        if (g.v_size < K - 2) continue;
        if (K == 3) {
            N += g.sub_v_size[K - 2];
            continue;
        }
        if (K == 4) {
            N += g.sub_e_size[K - 2];
            continue;
        }
        // assert(g.e_size == C[i].second);
        if (g.sub_e_size[K-2] * 2 == g.sub_v_size[K-2] * (g.sub_v_size[K-2] - 1)) {
            // continue;
            if (K - 2 == g.sub_v_size[K-2] || K - 2 == 0) {
                NN++;
                continue;
            }
            if (K - 2 == 1 || K - 2 + 1 == g.sub_v_size[K-2]) {
                for (int j = 0; j < g.sub_v_size[K-2]; j++) NN++;
                continue;
            }
            assert(g.new2old.size() == g.v_size);
            int minN = min(K - 2, g.sub_v_size[K-2] - (K - 2));
            // EBBkC_Comb_list(nullptr, g.sub_v_size[K-2], 0, 0, minN, &NN);
            // EBBkC_Comb_list2(g.sub_v_size[K-2], 0, 0, minN, 1, &NN);
            if (!egg.empty() && egg[0].size() != g.sub_v_size[K-2]) {
                for (auto nval : mSet) {
                    auto mn = min((int)egg[0].size() - nval.first, nval.first);
                        EBBkC_Comb_list2(egg[0].size(), 0, 0, mn, nval.second, &NN);
                }
                mSet.clear();
                egg.clear();
            }
            // egg.emplace_back(g.new2old);
            mSet[minN]++;
            continue;
        }
        // continue;

        {
            // cout << g.v_size << " " << g.e_size << endl;
            // maxdepth = 0;
            g.EBBkC_plus_plus(K - 2, &N, lz);
            // if (maxdepth) DEP[maxdepth]++;

            if (!st.tree.empty()) {
                if (maxdepth > Maxdepth || lz_his < minlz) {
                    search(1, st.tree, 0, K - 2 - (st.tree[0].val[0] == -2), 0, &N);
                }
                else {

                    string str;
                    if (st.tree[0].val[0] == -1) {
                        str = to_string(st.tree[0].val.size() - 1) + "(";
                    }
                    else {
                        str = "0(";
                    }
                    str += encode(0) + ")";
                    auto hashVal = fnv1a_64(str);

                    // search(1, st.tree, 0, K - 2 - (st.tree[0].val[0] == -2), 0);

                    st.hashVal = hashVal;
                    if (!vertexHash.count(hashVal)) {
                        vertexHash[hashVal] = ST.size();
                        // st.storeV();

                        ST.emplace_back(st);
                    }
                    else {
                        // st.storeV();
                        ST[vertexHash[hashVal]].increase(st.nodesets);
                    }
                }
                st.init();
            }
            // CCCNT[dfs_count]++;
            // dfs_count = 0;
        }
    }

    // getres2();

    int all = 0;
    size_t eggsz = 0;
    // for (const auto& i : egg) {
    //     eggsz += i.size() * sizeof(int);
    // }
    // auto rm = calculate_vector_searchTree_memory(ST);
    // std::cout << "ST memory usage: " << calculate_vector_searchTree_memory(ST) + eggsz << std::endl;
    if (!ST.empty()) {
        for (const auto& i : ST) {
            search(i.cnt + 1, i.tree, 0, K - 2 - (i.tree[0].val[0] == -2), 0, &N);
            allBranch += i.cnt + 1;
            calBranch++;
        }
        ST.clear();
    }
    // cout << "Same branch: " <<  calBranch << " " << allBranch - calBranch << " " << allBranch << endl;
    if (!egg.empty()) {
        // cal
        // auto ti = Get_T();

        for (auto nval : mSet) {
            EBBkC_Comb_list2(egg[0].size(), 0, 0, min((int)egg[0].size() - nval.first, nval.first), nval.second, &NN);

            allBranch += nval.second;
            calBranch++;
        }
        // time2 += Duration(ti);
        mSet.clear();
        egg.clear();
    }
    // cout << "DEP" << endl;
    // for (auto i : DEP) {
    //     cout << i.first << " " << i.second << endl;
    // }
    N += NN;
    // N = N2;
    gettimeofday(&time_end, NULL);
    double runtime = (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;

    // auto edd = std::chrono::high_resolution_clock::now();
    // double Alltime = std::chrono::duration_cast<std::chrono::microseconds>(edd - ti).count() / (double)1000;
    time_t time = chrono::system_clock::to_time_t(chrono::system_clock::now());
    string timestr = string(ctime(&time));
    cout << timestr << endl;
    timestr.pop_back();

    cout << "This is our" << endl;
    cout << "truss_decompose_time: " << truss_decompose_tim1e << endl;
    cout << runtime << "ms" << endl;
    return runtime;
}

double EBBkC_t::list_k_clique_parallel(const char *file_name, int threads) {
    EBBkC_Graph_t G, g;
    struct rusage end;

    printf("Reading edges from %s ...\n", file_name);
    timeval truss_decompose_time_end;

    G.truss_decompose(file_name);
    // auto s2 = std::chrono::high_resolution_clock::now();
    // double truss_decompose_time = std::chrono::duration_cast<std::chrono::microseconds>(s2 - ti).count() / (double)1000;
    gettimeofday(&truss_decompose_time_end, NULL);
    double truss_decompose_tim1e = (truss_decompose_time_end.tv_sec - time_start.tv_sec) * 1000.0 + (truss_decompose_time_end.tv_usec - time_start.tv_usec) / 1000.0;

    printf("Building necessary data structure ...\n");
    G.build(false);

    printf("Iterate over all cliques\n");

    g.truss_num = G.truss_num;

    g.build(true);

    vector<pair<int, int>> C;
    for (int i = 0; i < G.e_size; i++) {
        if (G.C_size[i] * 2 >= (K - 2) * (K - 3)) {
            C.emplace_back(i, G.C_size[i]);
        }
    }

    sort(C.begin(), C.end(), [](const pair<int, int>& a, const pair<int, int>& b) {return a.second > b.second;});


    unsigned long long NN = 0;
    int lz;
    vector<vector<int>> egg;
    map<int, int> mSet;

    ull N2 = 0;
    cout << C.size() << endl;
    // 1007211
    map<int, int> DEP;

    Maxdepth = 10;
    minlz = 3;
    int allBranch = 0;
    int calBranch = 0;
    auto tt0 = omp_get_wtime();
    for (auto & i : C) {
    // for (int i = 0; i < G.e_size; i++) {
        lz = 0;
        G.branch(i.first, &g);
        // G.branch(i, &g);
        if (g.v_size < K - 2) continue;
        if (K == 3) {
            N += g.sub_v_size[K - 2];
            continue;
        }
        if (K == 4) {
            N += g.sub_e_size[K - 2];
            continue;
        }
        // assert(g.e_size == C[i].second);
        if (g.sub_e_size[K-2] * 2 == g.sub_v_size[K-2] * (g.sub_v_size[K-2] - 1)) {
            // continue;
            if (K - 2 == g.sub_v_size[K-2] || K - 2 == 0) {
                N++;
                continue;
            }
            if (K - 2 == 1 || K - 2 + 1 == g.sub_v_size[K-2]) {
                for (int j = 0; j < g.sub_v_size[K-2]; j++) N++;
                continue;
            }
            // assert(g.new2old.size() == g.v_size);
            // int minN = min(K - 2, g.sub_v_size[K-2] - (K - 2));
            // // EBBkC_Comb_list(nullptr, g.sub_v_size[K-2], 0, 0, minN, &NN);
            // EBBkC_Comb_list2(g.sub_v_size[K-2], 0, 0, minN, 1, &NN);
            // if (!egg.empty() && egg[0].size() != g.sub_v_size[K-2]) {
            //     // cal
            //     // auto ti = Get_T();
            //     for (auto nval : mSet) {
            //         auto mn = min((int)egg[0].size() - nval.first, nval.first);
            //             //C(a, b)         a                             b     w
            //             EBBkC_Comb_list2(egg[0].size(), 0, 0, mn, nval.second, &NN);
            //         // allBranch += nval.second;
            //         // calBranch++;
            //         //
            //         // numberOfLF += nval.second;
            //         // numberOfNode += nval.second;
            //         // numberOfFCV += egg[0].size() * nval.second;
            //         // numberOfVertex += egg[0].size() * nval.second;
            //     }
            //     // time2 += Duration(ti);
            //     // cout << C[i].first << " " << NN << endl;
            //     mSet.clear();
            //     egg.clear();
            // }
            // egg.emplace_back(g.new2old);
            // mSet[minN]++;

            std::vector<int> vv (g.sub_v_size[K-2] + 1);
            vv[0] = -1;
            // Node a(0, -1, vv);
            st.tree.emplace_back(0, -1, vv);
            std::string str = std::to_string(g.sub_v_size[K-2] + 1);
            auto hashVal = fnv1a_64(str);

            // search(1, st.tree, 0, K - 2, 0, &N);

            st.hashVal = hashVal;
            if (!vertexHash.count(hashVal)) {
                vertexHash[hashVal] = ST.size();
                // st.storeV();
                ST.emplace_back(st);
            }
            else {
                // st.storeV();
                ST[vertexHash[hashVal]].increase(st.nodesets);
            }
            st.init();

            continue;
        }
        // continue;

        {
            // cout << g.v_size << " " << g.e_size << endl;
            // maxdepth = 0;
            g.EBBkC_plus_plus(K - 2, &N, lz);
            // if (maxdepth) DEP[maxdepth]++;

            if (!st.tree.empty()) {
                if (maxdepth > Maxdepth || lz_his < minlz) {
                    search(1, st.tree, 0, K - 2 - (st.tree[0].val[0] == -2), 0, &N);
                }
                else {
                    string str;
                    if (st.tree[0].val[0] == -1) {
                        str = to_string(st.tree[0].val.size() - 1) + "(";
                    }
                    else {
                        str = "0(";
                    }
                    str += encode(0) + ")";
                    auto hashVal = fnv1a_64(str);

                    // search(1, st.tree, 0, K - 2 - (st.tree[0].val[0] == -2), 0);

                    st.hashVal = hashVal;
                    if (!vertexHash.count(hashVal)) {
                        vertexHash[hashVal] = ST.size();
                        // st.storeV();

                        ST.emplace_back(st);
                    }
                    else {
                        // st.storeV();
                        ST[vertexHash[hashVal]].increase(st.nodesets);
                    }
                }
                st.init();
            }
        }
    }

    // getres2();

    int all = 0;
    size_t eggsz = 0;
    // for (const auto& i : egg) {
    //     eggsz += i.size() * sizeof(int);
    // }
    auto rm = calculate_vector_searchTree_memory(ST);
    omp_set_num_threads(threads);
    int n_edges;
    cout << ST.size() << endl;

    auto tt1 = omp_get_wtime();
#pragma omp parallel private(n_edges) reduction(+:N)
    {
        n_edges = 0;
#pragma omp for schedule(dynamic, 1) nowait
        for (auto i : ST) {
            search(i.cnt + 1, i.tree, 0, K - 2 - (i.tree[0].val[0] == -2), 0, &N);
            n_edges++;
        }
    // cout << N3 << endl;
    printf("Thread: %d, handled %d edges\n", omp_get_thread_num(), n_edges);
    }
    auto pt1 = (omp_get_wtime() - tt0) * 1e3;
    auto pt2 = (omp_get_wtime() - tt1) * 1e3;

    N += NN;
    // N = N2;
    gettimeofday(&time_end, NULL);
    double runtime = (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;

    // auto edd = std::chrono::high_resolution_clock::now();
    // double Alltime = std::chrono::duration_cast<std::chrono::microseconds>(edd - ti).count() / (double)1000;
    time_t time = chrono::system_clock::to_time_t(chrono::system_clock::now());
    string timestr = string(ctime(&time));
    cout << timestr << endl;
    timestr.pop_back();

    cout << "This is our-P" << endl;
    cout << "truss_decompose_time: " << truss_decompose_tim1e << endl;
    cout << pt1 << "ms" << " " << pt2 << "ms" << endl;
    return runtime;
}
