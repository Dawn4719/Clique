#include "graph.h"

#include <sys/time.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <chrono>
#include <fstream>
#include <cassert>

#include "../util.h"
#include "../log/log.h"
#include "../search/search_util.h"
#include "../../decompose/parallel_all_edge_cnc.h"

void free_graph(graph_t *g) {
    if (g->adj != nullptr)
        free(g->adj);

    if (g->num_edges != nullptr)
        free(g->num_edges);

    if (g->eid != nullptr)
        free(g->eid);

    if (g->edge_rank != nullptr)
        free(g->edge_rank);

    if (g->edge_truss != nullptr)
        free(g->edge_truss);
}

double timer() {
    struct timeval tp{};
    gettimeofday(&tp, nullptr);
    return ((double) (tp.tv_sec) + tp.tv_usec * 1e-6);
}

Graph::Graph(const char *dir_cstr) {
    // dir = string(dir_cstr);
    dir = string("/home/qsl/exp/EBBkC/dataset/") + string(dir_cstr);
    std::cout << dir << std::endl;
    ReadDegree();
    ReadAdjacencyList();
    CheckInputGraph();
}

Graph::Graph(const char *dir_cstr, int K_) {
    dir = string(dir_cstr);
    std::cout << dir << std::endl;
    string pa;
    std::fstream fs = std::fstream("/home/qsl/exp/EBBkC/dataset/" + dir + "/" + dir + ".clean", fstream::in);
    std::cout << "/home/qsl/exp/EBBkC/dataset/" + dir + "/" + dir + ".clean" << std::endl;
    // fs >> nodemax >> edgemax;
    // std::cout << nodemax << " " << edgemax << " " << K_ << std::endl;
    char c;
    int x;
    // for (int i = 0; i < nodemax; i++) {
    //     fs >> c >> x >> x;
    // }

    // int* vMaxClique = new int [nodemax];
    // int* eMaxClique = new int [edgemax];
    vector<int> vMaxClique(1);
    vector<int> eMaxClique(1);
    std::fstream fs2 = std::fstream("/home/qsl/exp/EBBkC/dataset/ifm/" + dir + "_V.txt", fstream::in);
    std::cout << "/home/qsl/exp/EBBkC/dataset/ifm/" + dir + "_V.txt" << std::endl;
    double t; fs2 >> t;
    while (fs2 >> vMaxClique.back()) {
        vMaxClique.resize(vMaxClique.size() + 1);
    }
    fs2.close();

    // fs2 = std::fstream("/home/qsl/exp/EBBkC/dataset/ifm/" + dir + "_E.txt", fstream::in);
    // std::cout << "/home/qsl/exp/EBBkC/dataset/ifm/" + dir + "_E.txt" << std::endl;
    //
    // while (fs2 >> eMaxClique.back()) {
    //     eMaxClique.resize(eMaxClique.size() + 1);
    // }
    // fs2.close();

    nodemax = vMaxClique.size() - 1;
    edgemax = eMaxClique.size() - 1;
    std::cout << nodemax << " -1- " << edgemax << std::endl;

    degree.resize(nodemax);
    int* ls = new int[nodemax + 1];
    memset(ls, -1, sizeof(int) * (nodemax + 1));

    int cnt_ = 0;
	int eMaxClique_cnt = 0;
	int vertex_cnt = 0;

    node_off = (uint32_t *) malloc(sizeof(uint32_t) * (nodemax + 1));


    // std::pair<int, int>* edges = new std::pair<int, int>[edgemax];
    // vector<pair<int, int>> ee;
    // for (int i = 0; i < edgemax; i++) {
    //     int u, v;
    //     fs >> u >> v;
    //     if (vMaxClique[u - 1] >= K_ && vMaxClique[v - 1] >= K_ && eMaxClique[i] >= K_) {
    //         if (ls[u] == -1) ls[u] = vertex_cnt++;
    //         if (ls[v] == -1) ls[v] = vertex_cnt++;
    //
    //         // edges[i].first = u;
    //         // edges[i].second = v;
    //         // degree[u]++;
    //         // degree[v]++;
    //
    //         edges[eMaxClique_cnt].first = ls[u];
    //         edges[eMaxClique_cnt].second = ls[v];
    //         degree[ls[u]]++;
    //         degree[ls[v]]++;
    //         eMaxClique_cnt++;
    //         // ee.emplace_back(ls[u], ls[v]);
    //     }
    // }

    int u, v;
    int i = 0;
    vector<std::pair<int, int>> edges;
    while (fs >> u >> v) {
        if (vMaxClique[u - 1] >= K_ && vMaxClique[v - 1] >= K_) {
            if (ls[u] == -1) ls[u] = vertex_cnt++;
            if (ls[v] == -1) ls[v] = vertex_cnt++;
            // edges[i].first = u;
            // edges[i].second = v;
            // degree[u]++;
            // degree[v]++;
            if (ls[u] >= degree.size()) degree.resize(ls[u] + 1);
            if (ls[v] >= degree.size()) degree.resize(ls[v] + 1);
            degree[ls[u]]++;
            degree[ls[v]]++;
            edges.emplace_back(ls[u], ls[v]);
            eMaxClique_cnt++;
            // ee.emplace_back(ls[u], ls[v]);
        }
    }
    // fstream ot("../" + dir + to_string(K_) + ".txt", fstream::out);
    // sort(ee.begin(), ee.end());
    // for (auto [u, v] : ee)
    //     ot << u << " " << v << endl;
    // ot.close();
    nodemax = vertex_cnt;
    edgemax = eMaxClique_cnt;

    cout << nodemax << " -- " << edgemax << std::endl;

    for (auto i = 0u; i < nodemax; i++) { node_off[i + 1] = node_off[i] + degree[i]; degree[i] = 0;}
    if (node_off[nodemax] == edgemax * 2) {
        edgemax = node_off[nodemax];
    }
    assert(node_off[nodemax] == edgemax);

    edge_dst = static_cast<int *>(malloc(sizeof(int) * static_cast<uint64_t>(edgemax + 16)));

    for (int i = 0; i < edgemax / 2; i++) {
        edge_dst[node_off[edges[i].first] + degree[edges[i].first]++] = edges[i].second;
        edge_dst[node_off[edges[i].second] + degree[edges[i].second]++] = edges[i].first;
    }

    for (int i = 0; i < nodemax; i++) {
        std::sort(edge_dst + node_off[i], edge_dst + node_off[i + 1]);
    }

    // delete[] vMaxClique;
    // delete[] eMaxClique;
    delete[] ls;
    // delete[] edges;

    // ReadDegree();
    // ReadAdjacencyList();
    CheckInputGraph();
    std::cout << "Finish" << std::endl;
}

using namespace std::chrono;

void Graph::ReadDegree() {
    auto start = high_resolution_clock::now();

    ifstream deg_file(dir + string("/b_degree.bin"), ios::binary);
    int int_size;
    deg_file.read(reinterpret_cast<char *>(&int_size), 4);

    deg_file.read(reinterpret_cast<char *>(&nodemax), 4);
    deg_file.read(reinterpret_cast<char *>(&edgemax), 4);
    log_info("int size: %d, n: %s, m: %s", int_size, FormatWithCommas(nodemax).c_str(),
             FormatWithCommas(edgemax).c_str());

    degree.resize(static_cast<unsigned long>(nodemax));
    deg_file.read(reinterpret_cast<char *>(&degree.front()), sizeof(int) * nodemax);

    auto end = high_resolution_clock::now();
    log_info("read degree file time: %.3lf s", duration_cast<milliseconds>(end - start).count() / 1000.0);
}

void Graph::ReadAdjacencyList() {
    auto start = high_resolution_clock::now();
    ifstream adj_file(dir + string("/b_adj.bin"), ios::binary);

    // csr representation
    node_off = (uint32_t *) malloc(sizeof(uint32_t) * (nodemax + 1));

    // prefix sum
    node_off[0] = 0;
    for (auto i = 0u; i < nodemax; i++) { node_off[i + 1] = node_off[i] + degree[i]; }
    if (node_off[nodemax] == edgemax * 2) {
        edgemax = node_off[nodemax];
    }
    assert(node_off[nodemax] == edgemax);

    edge_dst = static_cast<int *>(malloc(sizeof(int) * static_cast<uint64_t>(edgemax + 16)));

    string dst_v_file_name = dir + string("/b_adj.bin");
    auto dst_v_fd = open(dst_v_file_name.c_str(), O_RDONLY, S_IRUSR | S_IWUSR);
    int *buffer = (int *) mmap(0, static_cast<uint64_t >(edgemax) * 4u, PROT_READ, MAP_PRIVATE, dst_v_fd, 0);

    auto end = high_resolution_clock::now();
    log_info("malloc, and sequential-scan time: %.3lf s", duration_cast<milliseconds>(end - start).count() / 1000.0);
    // load dst vertices into the array
#pragma omp parallel for schedule(dynamic, 1000)
    for (auto i = 0u; i < nodemax; i++) {
        // copy to the high memory bandwidth mem
        for (uint64_t offset = node_off[i]; offset < node_off[i + 1]; offset++) {
            edge_dst[offset] = buffer[offset];
        }
        // inclusive
        degree[i]++;
    }
    munmap(buffer, static_cast<uint64_t >(edgemax) * 4u);

#ifdef VERIFY_INPUT
    // Verify.
#pragma omp parallel for schedule(dynamic, 1000)
    for (auto u = 0u; u < nodemax; u++) {
        for (size_t offset = node_off[u]; offset < node_off[u + 1]; offset++) {
            auto v = edge_dst[offset];
            if (BranchFreeBinarySearch(edge_dst, node_off[v], node_off[v + 1], (int) u) == node_off[v + 1]) {
                log_fatal("CSR not correct...");
                exit(-1);
            }
        }
    }
    log_info("CSR verify pass");
#endif

    auto end2 = high_resolution_clock::now();
    log_info("read adjacency list file time: %.3lf s", duration_cast<milliseconds>(end2 - end).count() / 1000.0);
}

void Graph::CheckInputGraph() {
    auto start = high_resolution_clock::now();

#pragma omp parallel for schedule(dynamic, 5000)
    for (auto i = 0u; i < nodemax; i++) {
        for (auto j = node_off[i]; j < node_off[i + 1]; j++) {
            if (edge_dst[j] == static_cast<int>(i)) {
                log_error("Self loop of v: %d", i);
                exit(1);
            }
            if (j > node_off[i] && edge_dst[j] <= edge_dst[j - 1]) {
                log_error("Edges not sorted in increasing id order!\nThe program may not run properly!");
                exit(1);
            }
        }
    }
    auto end = high_resolution_clock::now();
    log_info("check input graph file time: %.3lf s", duration_cast<milliseconds>(end - start).count() / 1000.0);
}
