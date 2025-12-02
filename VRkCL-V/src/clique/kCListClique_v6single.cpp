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
int K;
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

inline void EBBkC_Comb_list2(int list_size, int start, int picked, int k, int m, unsigned long long* cliques) {
    return;
    int l1 = list_size - (k - picked), l2, l3, i, j, q, w;
    // if (picked == k) {
    // 	for (w = 0; w < m; w++)
    // 		(*cliques)++;
    // 	return;
    // }
    // else if (picked == k - 1) {
    // 	for (w = 0; w < m; w++)
    // 		for (i = start; i <= l1; i++)
    // 			(*cliques)++;
    // 	return;
    // }
    // else
    if (picked == k - 2) {
        l2 = l1 + 1;
        for (w = 0; w < m; w++)
            for (i = start; i <= l1; i++) {
                for (j = i + 1; j <= l2; j++)

                    (*cliques)++;
            }
        return;
    }

    for (i = start; i <= l1; i++) {
        EBBkC_Comb_list2(list_size, i + 1, picked + 1, k, m, cliques);
    }
}

inline void cal(int l, int lz, int depp, unsigned long long* cliques) {
	if (depp == 0) {
		l = K - 1, lz = 0;
	}
	for (int i = depp; i < his.size(); i++) {
		if (l == 2) {
			if (l <= lz)
				if (l == 0 || l == lz)
					(*cliques)++;
				else if (l == 1 || l + 1 == lz)
					for (int q = 0; q < lz; q++)
						(*cliques)++;
				else
					EBBkC_Comb_list2(lz, 0, 0,  std::min(l, lz - l), 1, cliques);
			break;
		}
		if (his[i][0] == -1) {
			lz += his[i].size() - 1;
		}
		else {
			l--;
		}
		if (l == 0) {
			(*cliques)++;
			break;
		}

		if (l <= lz) {
			if (l == 0 || l == lz)
				(*cliques)++;
			else if (l == 1 || l + 1 == lz)
				for (int q = 0; q < lz; q++)
					(*cliques)++;
			else
				EBBkC_Comb_list2(lz, 0, 0,  std::min(l, lz - l), 1, cliques);
		}

		// if (his[i][0] == -1 && i + 1 == his.size()) {
		//     break;
		// }
	}
}
//
std::string code;
std::string encode(int idx) {
	std::string tmp;
	std::vector<std::pair<int, int>> v;
	for (auto i : st.tree[idx].sons) {
		if (st.tree[i].val[0] == -1) v.emplace_back(st.tree[i].val.size() - 1, i);
		else v.emplace_back(0, i);
	}
	sort(v.begin(), v.end());

	for (auto i : v) {
		tmp += std::to_string(i.first);
		if (st.tree[i.second].sons.empty()) {
			continue;
		}
		tmp += '(';
		tmp += encode(i.second);
		tmp += ')';
	}
	return tmp;
}

std::unordered_map<uint64_t, int> vertexHash;
std::vector<searchTree> ST;
int Maxdepth;

inline void search(const int cnt, const std::vector<Node>& tree, int idx, int l, int lz, unsigned long long* cliques) {
    if (tree[idx].val[0] == -1) {
        lz += tree[idx].val.size() - 1;
        // idx = tree[idx].sons[0];
    }

    if (l == 0) {
        (*cliques) += cnt;
        if (tree[idx].val[0] == -1) lz -= tree[idx].val.size() - 1;
        // cout << "++clique " << N << endl;
        return;
    }

    if (tree[idx].val[0] == -1 && tree[idx].sons.empty()) {
        // auto ti = Get_T();
        if (l <= lz)
            if (l == 0 || l == lz)
                (*cliques) += cnt;
            else if (l == 1 || l + 1 == lz)
                for (int j = 0; j < lz; j++) (*cliques) += cnt;
            else
                EBBkC_Comb_list2(lz, 0, 0,  std::min(l, lz - l), cnt, cliques);

        // cout << "++comblist candidateSize[l] == 0 " << lz << " " << min(l, lz - l) << " " << N << endl;
        // time2 += Duration(ti);
        if (tree[idx].val[0] == -1) lz -= tree[idx].val.size() - 1;
        // cout << "candidate == 0 " << N << endl;
        return;
    }

    if (l <= lz) {
        // auto ti = Get_T();
            if (l == 0 || l == lz)
                (*cliques) += cnt;
            else if (l == 1 || l + 1 == lz)
                for (int j = 0; j < lz; j++) (*cliques) += cnt;
            else
                EBBkC_Comb_list2(lz, 0, 0,  std::min(l, lz - l), cnt, cliques);
        // time2 += Duration(ti);
    }

    for (auto node : tree[idx].sons) {
        // if (tree[node].val[0] == -2) {

        if (idx == 0 && tree[node].val[0] == -1) {
            search(cnt, tree, node, l, lz, cliques);
        }
        else if (tree[node].sons.size() == 1 && tree[tree[node].sons[0]].val[0] == -1) {
            search(cnt, tree, tree[node].sons[0], l - 1, lz, cliques);
        }
        else {
            search(cnt, tree, node, l - 1, lz, cliques);
        }
    }
    if (tree[idx].val[0] == -1) lz -= tree[idx].val.size() - 1;
}
int maxdepth;
bool hasAdd;
int minlz;

inline void add(int l, int lz, unsigned long long* cliques) {
    // PT();
    hasAdd = true;

    // for (const auto& ii : his) {
    //     for (auto jj : ii)
    //         std::cout << jj << " ";
    //     std::cout << " - ";
    // }
    // std::cout << std::endl;
    if (his.size() == 1) return;
    if (depp == 0 && !st.tree.empty()) {
        // assert(maxdepth <= 4);
        lz_his = lz;
        if (maxdepth > Maxdepth || lz < minlz) {
            // if (st.tree[0].val[0] == -1)
			search(1, st.tree, 0, K - 1 - (st.tree[0].val[0] == -2), 0, cliques);
        }
        else {
            std::string str;
            if (st.tree[0].val[0] == -1) {
                str = std::to_string(st.tree[0].val.size() - 1) + "(";
            }
            else {
                str = "0(";
            }
            str += encode(0) + ")";
            auto hashVal = fnv1a_64(str);

            // if (st.tree[0].val[0] == -2)
            // search(1, st.tree, 0, K - 1 - (st.tree[0].val[0] == -2), 0);
            // else {
            st.hashVal = hashVal;

            if (!vertexHash.count(hashVal)) {
                vertexHash[hashVal] = ST.size();
                st.storeV();
                // st.cnt = 1;
                ST.emplace_back(st);
            }
            else {
                st.storeV();
                ST[vertexHash[hashVal]].increase(st.nodesets);
            }
        }
        st.init();
        maxdepth = 0;
        lz_his = 0;
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

// std::size_t memory_usage(const Node& n) {
//     std::size_t m = sizeof(Node); // Node 对象本身（通常在 vector<Node> 那块大内存里）
//
//     // val 和 sons 是两个动态数组
//     m += n.val.capacity()  * sizeof(int);
//     m += n.sons.capacity() * sizeof(int);
//
//     return m;
// }
//
// std::size_t memory_usage(const searchTree& t) {
//     std::size_t m = sizeof(searchTree); // searchTree 对象本身
//
//     // 1) tree: vector<Node>
//     //   a) 底层 Node 数组
//     m += t.tree.capacity() * sizeof(Node);
//
//     //   b) 每个 Node 里面的 val/sons 的动态数组
//     for (const auto& node : t.tree) {
//         m += node.val.capacity()  * sizeof(int);
//         m += node.sons.capacity() * sizeof(int);
//     }
//
//     // 2) nodesets: vector<string>
//     //   a) 底层 string 数组
//     m += t.nodesets.capacity() * sizeof(std::string);
//
//     //   b) 每个 string 内部的字符缓冲区
//     for (const auto& s : t.nodesets) {
//         m += s.capacity(); // 每个 char 按 1 字节算
//     }
//
//     return m;
// }
//
// // 估算 vector<searchTree> 占用的内存
// std::size_t memory_usage(const std::vector<searchTree>& vec) {
//     std::size_t m = sizeof(vec); // vector 对象本身（通常在栈或别处）
//
//     // 1) vector<searchTree> 底层数组：capacity() * sizeof(searchTree)
//     m += vec.capacity() * sizeof(searchTree);
//
//     // 2) 每个 searchTree 里的内部动态内存
//     // 注意：memory_usage(t) 里面已经加了 sizeof(searchTree)，
//     // 如果不想重复，可以在这里减掉一份：
//     for (const auto& t : vec) {
//         m += memory_usage(t) - sizeof(searchTree);
//     }
//
//     return m;
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

std::map<std::string, size_t> get_index_mem2() {
    FILE *fp = fopen("/proc/self/status", "r");
    char line[128];
    std::map<std::string, size_t> res;
    while (fgets(line, 128, fp) != NULL) {
        //        if (strncmp(line, "VmPeak", 2) == 0)
        //        {
        //            cout << line << endl;
        ////            printf("当前进程占用虚拟内存大小为：%d KB\n", atoi(line + 6));
        //        }
        if (strncmp(line, "VmRSS:", 6) == 0) {
            std::string p = line;
            res["now"] = size_t(stoull(p.substr(6)));
            std::cout << line;
        }
        if (strncmp(line, "VmPeak:", 7) == 0) {
            std::string p = line;
            res["pk"] = size_t(stoull(p.substr(7)));
            std::cout << line;
        }
    }
    fclose(fp);
    return res;
}

unsigned long long kclistClique::run_v6() {
    printf("kClistClique_v6.cpp\n");

    // g.initHash();
    printf("initHash\n");
    K = k;
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

    Maxdepth = 10;
    minlz = 3;
    gv.reserve(g.n);
    for (ui u = 0; u < g.n; u++) {
        if (g.pIdx[u + 1] - g.pIdx2[u] >= k - 1)
        gv.emplace_back(g.pIdx[u + 1] - g.pIdx2[u], u);
    }
    sort(gv.begin(), gv.end(), [](std::pair<ui, ui> a, std::pair<ui, ui> b) {
        return a.first > b.first;
    });
    // for (ui u = 0; u < g.n; u++) if (g.pIdx[u + 1] - g.pIdx2[u] >= k - 1) {
        // std::cout << "Induce: " << u << std::endl;

    std::cout << "frst " << std::endl;
    auto m1 = get_index_mem2()["now"];
    for (auto [d, u] : gv) {
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
            if (k - 1 == nxtC.size() || k - 1 == 0) {
                answer++;
            }
            else if (k - 1 == 1 || k == nxtC.size()) {
                for (int j = 0; j < nxtC.size(); j++) answer++;
            }
            else {
                auto minN = std::min(k - 1, (int)nxtC.size() - (k - 1));
                assert(minN >= 0 && minN <= 150);
                // EBBkC_Comb_list2(nxtC.size(), 0, 0, minN, 1, &answer);
                if (!egg.empty() && egg[0].size() != nxtC.size()) {
                    for (auto nval : mSet) {
                        auto mn = std::min((int)egg[0].size() - nval.first, nval.first);
                        EBBkC_Comb_list2(egg[0].size(), 0, 0, mn, nval.second, &answer);
                    }
                    mSet.clear();
                    egg.clear();
                }
                egg.emplace_back(nxtC);
                mSet[minN]++;
            }

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

        if (!st.tree.empty()) {
            if (maxdepth > Maxdepth || lz_his < minlz) {
                search(1, st.tree, 0, K - 1 - (st.tree[0].val[0] == -2), 0, &answer);
            }
            else {
                std::string str;
                if (st.tree[0].val[0] == -1) {
                    str = std::to_string(st.tree[0].val.size() - 1) + "(";
                }
                else {
                    str = "0(";
                }
                str += encode(0) + ")";
                auto hashVal = fnv1a_64(str);

                // search(1, st.tree, 0, K - 1 - (st.tree[0].val[0] == -2), 0);

                st.hashVal = hashVal;
                if (!vertexHash.count(hashVal)) {
                    vertexHash[hashVal] = ST.size();
                    st.storeV();
                    ST.emplace_back(st);
                }
                else {
                    st.storeV();
                    ST[vertexHash[hashVal]].increase(st.nodesets);
                }
            }
            st.init();
        }

    }
    size_t eggsz = 0;
    // for (const auto& i : egg) {
    //     eggsz += i.size() * sizeof(int);
    // }
    // std::cout << "ST memory usage: " << (memory_usage(ST) + eggsz) / 1024 << std::endl;

    if (!egg.empty()) {
        for (auto nval : mSet) {
            EBBkC_Comb_list2(egg[0].size(), 0, 0, std::min((int)egg[0].size() - nval.first, nval.first), nval.second, &answer);
        }
        mSet.clear();
        egg.clear();
    }
    vertexHash.clear();
    auto m2 = get_index_mem2()["now"];
    answer = calculate_vector_searchTree_memory(ST);
    std::cout << "ST memory usage: " << calculate_vector_searchTree_memory(ST) << " " << (m2-m1) << std::endl;

    if (!ST.empty()) {
        for (const auto& i1 : ST) {
            search(i1.cnt + 1, i1.tree, 0, K - 1 - (i1.tree[0].val[0] == -2), 0, &answer);
        }
        ST.clear();
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
            if (lz < minlz && deep <= lz) {
            // if (deep <= lz) {
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

        if (his.size() > maxdepth) {
        	maxdepth = his.size();
        }
        if (lz >= minlz)
        	add(deep, lz, &answer);
        // for (const auto& ii : his) {
        //     for (auto jj : ii)
        //         std::cout << jj << " ";
        //     std::cout << " - ";
        // }
        // std::cout << std::endl;
        // std::cout << "l=2 " << deep << " " << lz << " " << answer << std::endl;
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
    // assert(candidateSize[deep] < g.coreNumber);
    if (canAddSize[deep]) {
        his.resize(his.size() + 1);
        his.back().reserve(canAddSize[deep] + 1);
        his.back().emplace_back(-1);
        for (int i = 0; i < canAddSize[deep]; i++) his.back().emplace_back(canAdd[deep][i]);
        assert(his.back().size() > 1);
    }
    lz += canAddSize[deep];
    lz_his = lz;
    if (deep == 0) {
        lz -= canAddSize[deep];
        if (lz < minlz) {
            answer++;
        }
        else {
            add(deep, lz, &answer);
        }
        // depp = his.size();

        // for (const auto& ii : his) {
        //     for (auto jj : ii)
        //         std::cout << jj << " ";
        //     std::cout << " - ";
        // }
        // std::cout << std::endl;

        if (canAddSize[deep] && !his.empty()) {
            his.pop_back();
            if (his.size() < depp) {
                depp = his.size();
                IDX.resize(depp);
            }
        }

        // std::cout << "l=0 " << deep << " " << lz << " " << answer << std::endl;
        return;
    }

    if (candidateSize[deep] == 0) {
        if (lz < minlz) {
            if (deep <= lz) {
                if (deep < lz - deep)
                    Comb_list(lz, 0, 0, deep, &answer);
                else
                    Comb_list(lz, 0, 0, lz - deep, &answer);
            }
        }
        else {
            add(deep, lz, &answer);
        }
        lz -= canAddSize[deep];
        // depp = his.size();

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

    // if (deep <= lz)
    //     if (deep < lz - deep)
    //         Comb_list(lz, 0, 0, deep, &answer);
    //     else
    //         Comb_list(lz, 0, 0, lz - deep, &answer);

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
