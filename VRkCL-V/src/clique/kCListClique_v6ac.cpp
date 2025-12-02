#include "kCListClique.h"
#include "../tools/listLinearHeap.hpp"
#include <algorithm>
#include <cassert>
#include <random>
#include <map>
#include <unordered_map>


// #define DEBUG
// #define DDEBUG

// #define PRINT_SG

#ifdef PRINT_SG
#define PRINT_AVG_SG
double sRaw = 0.0, sCCr = 0.0, cntNodes__ = 0;
#endif

// #define PRINT_ARBORICITY

#ifdef PRINT_ARBORICITY
std::vector<std::pair<ui, ui>> edgesRandom;
std::vector<ui> r;
std::vector<ui> Pava, reordered;

#define PRING_AVG_ARBORICITY
double sRaw = 0.0, sCCr = 0.0, cntNodes__ = 0;
#endif

// #define PRINT_RATIO_OF_KCLIQUES
#ifdef PRINT_RATIO_OF_KCLIQUES
double numberOfKcliquesInCoreClique = 0;
#endif
int lz_his;
std::vector<std::vector<int>> his;
std::vector<int> IDX;
std::map<int, int> DEPTH;
int depp;
struct Node {
    int id=-1;
    int fa=-1;
    std::vector<int> val{};
    std::vector<int> sons{};

    Node(int _id, int _fa, std::vector<int> _val): id(_id), fa(_fa) {
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
    std::vector<Node> tree;
    std::vector<std::string> nodesets;

    inline void init() {
        cnt = 0;
        hashVal = 0;
        if (!tree.empty()) tree.clear();
        nodesets.clear();
    }

    void increase(std::vector<std::string>& v) {
        cnt++;
        for (int i = 0; i < v.size(); i++) {
            nodesets[i] += " " + v[i];
        }
    }

    void storeV() {
        for (auto i : tree) {
            for (int j = 1; j < i.val.size(); j++) {
                nodesets.emplace_back(std::to_string(i.val[j]));
            }
        }
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
std::vector<std::vector<unsigned>> egg;
std::map<int, int> mSet;
double time2;

inline uint64_t fnv1a_64(const std::string& s) {
    uint64_t hash = 1469598103934665603ULL;      // offset basis
    for (unsigned char c : s) {
        hash ^= c;
        hash *= 1099511628211ULL;                // FNV prime
    }
    return hash;
}


inline void Comb_list(int list_size, int start, int picked, int k, unsigned long long* cliques) {
    // return;
    if (picked == k) {
        (*cliques)++;
        return;
    }

    for (int i = start; i <= list_size - (k - picked); i++) {
        Comb_list(list_size, i + 1, picked + 1, k, cliques);
    }
}

unsigned long long kclistClique::run_v6() {
    printf("kClistClique_v6.cpp\n");

    g.initHash();
    printf("initHash\n");

    edAdj.resize(k);
    for (ui i = 0; i < k; i++)
        edAdj[i].resize(g.n);
    deg.resize(k);
    for (ui i = 0; i < k; i++)
        deg[i].resize(g.n);

    candidateSize.resize(k);
    candidate.resize(k);
    for (ui i = 0; i < k; i++)
        candidate[i].resize(g.coreNumber + 1);
    canAddSize.resize(k);
    canAdd.resize(k);
    for (ui i = 0; i < k; i++)
        canAdd[i].resize(g.coreNumber + 1);

    adj.resize(g.n);
    for (ui i = 0; i < g.n; i++) {
        adj[i].resize(g.pIdx[i + 1] - g.pIdx2[i]);
        for (ui j = 0; j < g.pIdx[i + 1] - g.pIdx2[i]; j++) {
            adj[i][j] = g.pEdge[g.pIdx2[i] + j];
        }
    }

    level.resize(g.n);
    vis.resize(g.n);

    std::vector<ui> newCliqueNei;

    std::vector<ui> keys(g.coreNumber);
    std::vector<ui> heads(g.coreNumber);
    std::vector<ui> nxts(g.coreNumber);

    std::vector<ui> id_s(g.coreNumber + 1);
    std::vector<ui> degree(g.coreNumber);
    std::vector<bool> removed(g.coreNumber);

    std::vector<ui> pIdx(g.coreNumber + 1);
    std::vector<ui> pEdge(g.coreNumber * g.coreNumber);

    for (ui u = 0; u < g.n; u++) if (g.pIdx[u + 1] - g.pIdx2[u] >= k - 1) {
        // std::cout << "Induce: " << u << std::endl;

        for (ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            level[g.pEdge[i]] = i - g.pIdx2[u] + 1;
        }
        ui n = g.pIdx[u + 1] - g.pIdx2[u];
        for (ui i = 0; i <= n; i++)
            pIdx[i] = 0;
        // std::cout << "Cout Origin subgraph" << std::endl;
        for (ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            ui v = g.pEdge[i];

            for (ui j = g.pIdx2[v]; j < g.pIdx[v + 1]; j++) {
                ui w = g.pEdge[j];
                if (level[w] > 0) {
                    // std::cout << v << " " << w << std::endl;
                    pIdx[i - g.pIdx2[u] + 1]++;
                    pIdx[level[w]]++;
                }
            }
        }

        for (ui i = 1; i <= n; i++)
            pIdx[i] += pIdx[i - 1];
        for (ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            ui v = g.pEdge[i];
            for (ui j = g.pIdx2[v]; j < g.pIdx[v + 1]; j++) {
                ui w = g.pEdge[j];
                if (level[w]) {
                    pEdge[pIdx[i - g.pIdx2[u]]++] = level[w] - 1;
                    pEdge[pIdx[level[w] - 1]++] = i - g.pIdx2[u];
                }
            }
        }
        for (ui i = n; i > 0; i--)
            pIdx[i] = pIdx[i - 1];
        pIdx[0] = 0;

        clique.clear();
        std::vector<ui> &nxtC = nodes[k - 1];
        nxtC.clear();

        for (ui i = g.pIdx2[u]; i < g.pIdx[u + 1]; i++) {
            deg[k - 1][g.pEdge[i]] = 0;
            nxtC.pb(g.pEdge[i]);
            level[g.pEdge[i]] = 0;
        }

        for (auto v : nxtC)
            level[v] = k-1;

        // 重建图
        int edges = 0;
        for (auto v : nxtC) {
            ui ed = adj[v].size();
            for (ui i = 0; i < ed;) {
                if (level[adj[v][i]] != k-1)
                    std::swap(adj[v][i], adj[v][--ed]);
                else {
                    deg[k - 1][v]++;
                    deg[k - 1][adj[v][i]]++;
                    i++, edges++;
                }
            }
            edAdj[k - 1][v] = ed;
        }
        int lz = 0;
        if (edges * 2 == nxtC.size() * (nxtC.size() - 1)) {
            Comb_list(nxtC.size(), 0, 0, std::min(k - 1, (int)nxtC.size() - (k - 1)), &answer);
            for (auto u2 : clique)
                vis[u2] = 0;
            for (ui v : nxtC)
                level[v] = 0;
            continue;
        }
        listing(k - 1, clique.size(), lz);

        for (auto u2 : clique)
            vis[u2] = 0;
        for (ui v : nxtC)
            level[v] = 0;
        // std::cout << u << " " << answer << std::endl;
        // answer = 0;
    }

    return answer;
}



void kclistClique::listing(ui deep, ui edClique, int& lz) {
#ifdef DDEBUG
printf("    deep %u\n", deep);
printf("C:");
for(auto v: nodes[deep]) printf("%u ", v);printf("\n");
printf("clique:");
for(ui i = 0; i < edClique; i++) printf("%u ", clique[i]);printf("\n");
#endif

    std::vector<ui> & C = nodes[deep];
    if(deep == 2) {
        /*for(auto v:C)
            answer += edAdj[deep][v];*/
        if (lz >= 2) {
            if (deep <= lz) {
                if (deep < lz - deep)
                    Comb_list(lz, 0, 0, deep, &answer);
                else
                    Comb_list(lz, 0, 0, lz - deep, &answer);
            }

            // choose 1 AC
            for (auto i : C) {
                for (int c = 0; c < lz; c++) {
                    answer ++;
                }
            }

            // choose 0 AC
            for (auto u : C) {
                for(ui j = 0; j < edAdj[deep][u]; j++) {
                    ui v = adj[u][j];
                    answer ++;
                }
            }
        }
        else if (lz == 1) {
            // choose 1 AC
            for (auto i : C) {
                for (int c = 0; c < lz; c++) {
                    answer ++;
                }
            }

            // choose 0 AC
            for (auto u : C) {
                for(ui j = 0; j < edAdj[deep][u]; j++) {
                    ui v = adj[u][j];
                    answer ++;
                }
            }
        }
        else if (lz == 0) {
            // choose 0 AC
            for (auto u : C) {
                for(ui j = 0; j < edAdj[deep][u]; j++) {
                    ui v = adj[u][j];
                    answer ++;
                }
            }
        }

        // if (his.size() > maxdepth) {
        // 	maxdepth = his.size();
        // }
        // if (lz >= minlz)
        // 	add(l, lz);
        depp = his.size();
        // for (const auto& ii : his) {
        // 	for (auto jj : ii)
        // 		cout << jj << " ";
        // 	cout << " - ";
        // }
        // cout << endl;
        // cout << "l=2 " << l << " " << lz << " " << *n << endl;
        return;
    }

    if (deep > C.size() + lz)
        return;

    canAddSize[deep] = 0;
    candidateSize[deep] = 0;
    for (auto u : C) {
        if(deg[deep][u] == C.size() - 1)
            canAdd[deep][canAddSize[deep]++] = u;
        else
            candidate[deep][candidateSize[deep]++] = u;
    }
    assert(candidateSize[deep] < g.coreNumber);
    lz += canAddSize[deep];
    lz_his = lz;
    if (deep == 0) {
        lz -= canAddSize[deep];
        // if (lz < minlz) {
        answer++;
        // }
        // else {
        // add(l, lz);
        // }
        depp = his.size();
        // for (const auto& ii : his) {
        // 	for (auto jj : ii)
        // 		cout << jj << " ";
        // 	cout << " - ";
        // }
        // cout << endl;

        if (canAddSize[deep] && !his.empty()) {
            his.pop_back();
            if (his.size() < depp) {
                depp = his.size();
                IDX.resize(depp);
            }
        }

        // cout << "l=0 " << l << " " << lz << " " << *n << endl;
        return;
    }

    if (candidateSize[deep] == 0) {
        // if (lz < minlz) {
        if (deep <= lz) {
            if (deep < lz - deep)
                Comb_list(lz, 0, 0, deep, &answer);
            else
                Comb_list(lz, 0, 0, lz - deep, &answer);
        }
        // }
        // else {
        // add(deep, lz);
        // }
        lz -= canAddSize[deep];
        depp = his.size();

        // for (const auto& ii : his) {
        // 	for (auto jj : ii)
        // 		cout << jj << " ";
        // 	cout << " - ";
        // }

        if (canAddSize[deep] && !his.empty()) {
            // if (his.size() == IDX.size()) {
            //     IDX.pop_back();
            // depp--;
            his.pop_back();
            if (his.size() < depp) {
                depp = his.size();
                IDX.resize(depp);
            }
        }

        // cout << endl;
        // cout << "candidate=0 " << l << " " << lz << " " << *n << endl;
        return;
    }

    if (deep <= lz)
        if (deep < lz - deep)
            Comb_list(lz, 0, 0, deep, &answer);
        else
            Comb_list(lz, 0, 0, lz - deep, &answer);

    std::vector<ui> & nxtC = nodes[deep - 1];

    for(int i = 0; i < candidateSize[deep]; i++) {
        int u = candidate[deep][i];
        if(edAdj[deep][u] + 1 + lz < deep) continue;

        nxtC.clear();
        for(ui j = 0; j < edAdj[deep][u]; j++) {
            ui v = adj[u][j];
            if (deg[deep][v] == C.size() - 1) continue;
            nxtC.push_back(v);
            level[v] = deep - 1;
            deg[deep - 1][v] = 0;
        }

        for(ui j = 0; j < nxtC.size(); j++) {
            ui v = nxtC[j];
            ui & ed = edAdj[deep - 1][v];
            ed = edAdj[deep][v];

            for(ui l = 0; l < ed; ) {
                ui w = adj[v][l];
                if(level[w] == deep - 1) {
                    l++;
                    deg[deep - 1][v]++;
                    deg[deep - 1][w]++;
                }
                else std::swap(adj[v][l], adj[v][--ed]);
            }
        }
        his.emplace_back(std::vector<int>{-2, (int)u});
        listing(deep - 1, 0, lz);
        his.pop_back();
        if (his.size() < depp) {
            depp = his.size();
            IDX.resize(depp);
        }

        for(auto v : nxtC) level[v] = deep;
    }

    if (canAddSize[deep] && !his.empty()) {
        his.pop_back();
        if (his.size() < depp) {
            depp = his.size();
            IDX.resize(depp);
        }
    }
    lz -= canAddSize[deep];
}
