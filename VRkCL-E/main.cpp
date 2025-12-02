#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <set>
#include <algorithm>
#include <string>
#include <omp.h>
#include <cassert>
#include "set_operation.h"
#include "def.h"
#include "edge_oriented.h"
#include <fstream>
#include <chrono>
#define Get_Time() std::chrono::high_resolution_clock::now()
#define Duration(start) std::chrono::duration_cast< \
std::chrono::microseconds>(Get_Time() - start).count() / (float)1000
#define Print_Time(str, start) std::cout << str << Duration(start) << " ms" << std::endl
using namespace std;
int  K, L;
unsigned long long N = 0;

int main(int argc, char** argv) {
    string src_filename;
    long long dfs_count = 0;
    double runtime;
    src_filename = argv[1];
    K = atoi(argv[2]);
    runtime = EBBkC_t::list_k_clique(src_filename.c_str(), dfs_count); // argv[2]
    printf("Number of %u-cliques: %llu\n", K, N);
    printf("EBBkC+ET (t = %d) runtime %.2lf ms\n\n", L, runtime);
    
    // // EBBkC+ET (parallel)
    // string src_filename(argv[2]);
    // K = atoi(argv[3]);
    // omp_set_num_threads(atoi(argv[4]));

    // runtime = EBBkC_t::list_k_clique_parallel(src_filename.c_str(), atoi(argv[4]));
    // // runtime = EBBkC_t::list_k_clique_parallel(src_filename.c_str(), 1);
    // printf("Number of %u-cliques: %llu\n", K, N);
    // printf("EBBkC+ET runtime %.2lf ms\n\n", runtime);
    return 0;
}