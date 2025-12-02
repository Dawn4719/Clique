#ifndef DESCOL_EBBKC_H
#define DESCOL_EBBKC_H
#include "def.h"
#include <vector>
#include <fstream>

#include "truss/util/graph/graph.h"
using namespace std;
struct CoreGraph {
    uint32_t nodemax;
    uint32_t edgemax;
    uint32_t *node_off;
    int *sub;
    int *edge_dst;
    int *lab;
};
class EBBkC_Graph_t {

public:
    long long ffff = 0;
    // for both EBBkC & EBBkC+ framework
    int v_size = 0;
    int e_size = 0;
    int truss_num = 0;
    bool is_sub = false;
    HashMap_t edge2id;
    Edge_t* edges = nullptr;

    int* sub_e_size = nullptr;    // number of edges in the sub-branch
    int* sub_v_size = nullptr;    // number of nodes in the sub-branch
    int *canAddSize = nullptr;
    int *canAddFSize = nullptr;
    int *candidateSize = nullptr;
    // int *alreadyChooseSize = nullptr;
    int** sub_v = nullptr;        // nodes in the sub-branch
    int** sub_e = nullptr;        // edges in the sub-branch
    int** canAdd = nullptr;
    int** canAddF = nullptr;
    int** candidate = nullptr;
    // int** alreadyChoose = nullptr;

    // int* nowSize = nullptr;

    int** T = nullptr;            // T[e]: a set of nodes
    int* T_size = nullptr;        // T_size[e]: the size of T[e]
    int** C = nullptr;            // C[e]: a set of edges
    int* C_size = nullptr;        // C_size[e]: the size of C[e]

    vector<int> new2old;
    vector<int> rank;

    // for EBBkC framework only
    int* v_lab = nullptr;
    int* e_lab = nullptr;
    int** out_v_size = nullptr;
    int** out_e_size = nullptr;
    // for EBBkC+ framework only
    // int* col = nullptr;           // col[u]: the color of vertex 'u'
    int* col = nullptr;           // col[u]: the color of vertex 'u'
    int* lab = nullptr;           // lab[u]: the level of vertex 'u'
    int** G_deg = nullptr;        // G_deg[l][u]: the number of neighbors of vertex 'u' at level 'l'
    int** G_adj = nullptr;        // G_adj[l][u]: the set of neighbors of vertex 'u' at level 'l'
    int** DAG_deg = nullptr;
    int** DAG_adj = nullptr;      // DAG_adj[l][u]: the set of out-neighbors of vertex 'u' at level 'l'

    int* GBit = nullptr;
    unsigned** NeiBit = nullptr;
    int BitSize;
    int** DAG_out_adj = nullptr;
    bool** used = nullptr;

    int** vSame = nullptr; // stores same node
    int* vSameSize = nullptr;
    int Col = 0;

    // for EBBkC++ (early-termination)
    // int P_size = 0;
    // int P_act = 0;
    // int F_size = 0;
    // int *P = nullptr;
    // int *F = nullptr;
    // int *lack_size = nullptr;
    // int **lack = nullptr;
    // int *lev = nullptr;
    // int *loc = nullptr;

    // std::ifstream fin;
    // std::ofstream fout;
    // FILE *file;


    EBBkC_Graph_t();
    ~EBBkC_Graph_t();

//    void read_edges_from_file(const char* file_name);
    void truss_decompose(const char* w_file_name);
    void truss_decompose(const CoreGraph* G);
//    void read_ordered_edges_from_file(const char* file_name);
    void build(bool sub);
    void build(int n, int m, int truss);

//    void EBBkC(int l, unsigned long long *cliques);
//    void EBBkC_plus(int, unsigned long long *cliques);
//    void EBBkC_plus_plus(int l, unsigned long long *cliques);

    void branch(int e, EBBkC_Graph_t* g);
    void EBBkC_plus_plus(int l, unsigned long long* cliques, int& lz);
    void EBBkC_plus_plus_AC(int l, unsigned long long *cliques);
    void EBBKC_Comb_list_dfs(int l, unsigned long long *cliques);
    void EBBKC_Comb_list_dfs(int l, int *list, int list_size, int start, int picked, int kk, unsigned long long *cliques);
    void dfs(int l, int now,unsigned long long *cliques);
    void calcr(int l, unsigned long long *cliques);
    // bool can_terminate(int l, unsigned long long* cliques);
    // void list_in_plex(int start, int p, int q, unsigned long long* cliques);
    void EBBKC_My_Comb_list(int ii, int start, int picked, unsigned long long *cliques);
    void EnumSameList1and0AC(int picked, unsigned long long* cliques);
    void EnumSameList0AC(int picked, unsigned long long* cliques);


    // void EBBkC_Comb_list(int* list, int list_size, int start, int picked, int k, unsigned long long* cliques);
    // void EBBkC_Comb_listWithSame(int* list, int list_size, int start, int picked, int k, unsigned long long* cliques);
    // void EnumSameList2(int picked, unsigned long long* cliques);
    // void EnumSameList(int* list, int list_size, int picked, int k, unsigned long long* cliques);
    // void EBBkC_Comb_listWithOutSame(int* list, int list_size, int start, int picked, int k, unsigned long long* cliques);
    // void EBBkC_Comb_listWithOutSamePlus(int* list, int list_size, int start, int picked, int k, unsigned long long* cliques);
    // void EBBKC_My_Comb_list(int ii, int start, int picked, int need1, int need2, unsigned long long *cliques);
    // void EBBKC_My_Comb_list2(int ii, int start, int picked, int need1, int need2,unsigned long long *cliques);
};


class EBBkC_t {
public:
//    static double truss_order(const char* r_file_name, const char* w_file_name);
    static double list_k_clique(const char* file_name, long long& dfs_count);
    static double list_k_clique(const char* file_name, std::string s, long long& dfs_count);
    static double list_k_clique_parallel(const char* file_name, int threads);
};


#endif //DESCOL_EBBKC_H
