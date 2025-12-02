#include "../graph/graph.h"
#include "../tools/getArgs.hpp"
#include "../clique/kCListClique.h"
#include <iostream>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <map>
#include <sys/time.h>
#include <omp.h>

void printUsage() {
    std::cout << "bin/run " << std::endl;
    std::cout << "-f graph file directory(edge.bin & idx.bin)" << std::endl;
    std::cout << "-f_txt graph file text file, each edge exists one time" << std::endl;
    std::cout << "-f_txtD graph file text file, each edge exists two times" << std::endl;
    std::cout << "-k" << std::endl;
}

std::map<std::string, size_t> get_index_mem() {
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

int main(int argc, char * argv[])
{
    argsController ac(argc, argv);

    // char* argstr[100] = {"bin/kcc", "noUVM", "-f_txt", "ba", "-k", "6"};
    // argsController ac(6, argstr);
    if(!ac.exist("-f_txt") && !ac.exist("-f") && !ac.exist("-f_txtD")) {
        printUsage();
        return 0;
    }
    
    bool noVUM = false;
    if(ac.exist("noUVM")) noVUM = true;

    Graph g;
    g.k = std::stoi(ac["-k"]);
    auto pk1 = get_index_mem()["pk"];
    if(ac.exist("-f_txt")) {

        // std::string s = "/home/qsl/exp/coreCliqueRemoval/t.txt";
        // g.readFromText(ac["-f_txt"], noVUM);
        g.readFromTextDel(ac["-f_txt"], noVUM);
    }
    else if(ac.exist("-f_txtD")) g.readFromTextDoubleEdges(ac["-f_txtD"]);
    else g.readFromBin(ac["-f"]);

    auto graphsize = get_index_mem()["pk"];

    std::cout << "load graph: n " << g.n << " m " << g.m << " maxD " << g.maxD << std::endl;

    auto s1 = std::chrono::high_resolution_clock::now();
    timeval time_start;
    timeval time_end;
    gettimeofday(&time_start, NULL);
    g.changeToCoreOrder();

    auto s2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(s2 - s1);
    double coreTime = duration.count() / (double)1000;
    std::cout << "changeToCoreOrder:" << duration.count() << "ms" << std::endl;
    std::cout << "coreNumber:" << g.coreNumber << std::endl;

    int k = 5;
    if(ac.exist("-k")) k =  std::stoi(ac["-k"]);
    std::cout << "k:" << k << std::endl;

    int p=1;
    if(ac.exist("-p")) p =  std::stoi(ac["-p"]);
    // std::fstream f("/home/qsl/exp/coreCliqueRemoval-backup/CCCCR_" + ac["-f_txt"] + ".csv", std::ios::app);
    // f << "CCR," << k << "," << 0 << "," << coreTime << "," << 0 << "," << get_index_mem()["pk"] << "," << graphsize - pk1 << " " << g.m << std::endl;
    // f.close();
    // return 0;

    kclistClique pP(std::move(g), k);

    auto t1 = std::chrono::high_resolution_clock::now();

    double cnt = 0;
    omp_set_num_threads(p);
    auto p1 = omp_get_wtime();
    cnt = pP.run_v6();
    auto ptime = (omp_get_wtime() - p1) * 1e3;
    gettimeofday(&time_end, NULL);
    double runtime = (time_end.tv_sec - time_start.tv_sec) * 1000.0 + (time_end.tv_usec - time_start.tv_usec) / 1000.0;

    auto t2 = std::chrono::high_resolution_clock::now();
    auto durationt = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
    double listTime = duration.count() / (double)1000;
    std::cout << "time:" << listTime << "ms" << std::endl;

    std::cout << std::fixed << std::setprecision(0) <<
             k << "-clique:" << cnt << std::endl;

    durationt = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - s1);
    double ALLTIME = std::chrono::duration_cast<std::chrono::microseconds>(t2 - s1).count() / (double)1000;
    std::cout << "ALL time:" << ALLTIME << "ms " << runtime << " " << get_index_mem()["pk"] << std::endl;
    std::cout << "/home/qsl/exp/coreCliqueRemoval-backup/CCCCR-P" + std::to_string(p) + "_" + ac["-f_txt"] + ".csv" << std::endl;
    std::fstream f("/home/qsl/exp/coreCliqueRemoval-backup/CCCCR_P" + ac["-f_txt"] + ".csv", std::ios::app);
    // std::fstream f("/home/qsl/exp/coreCliqueRemoval-backup/CCCCRTreeSize" + std::to_string(p) + "_" + ac["-f_txt"] + ".csv", std::ios::app);
    // f << "CCR," << k << "," << cnt / 1024 << "," << get_index_mem()["pk"] << std::endl;
    // int mn = 0x3f3f3f3f, mx = 0, ct = 0;
    // unsigned long long sum = 0;
    // for (auto i : pP.CNT) {
    //     sum += (unsigned long long)i.first * i.second;
    //     if (i.first > mx) mx = i.first;
    //     if (i.first < mn) mn = i.first;
    //     ct += i.second;
    // }
    //
    // std::cout << mx << " " << mn << " " << sum << " " << ct << " " << 1.0 * sum / ct << std::endl;
    // f << "CCR," << ac["-f_txt"] << "," << k << "," << (unsigned long long)cnt << "," << mx << "," << mn << "," << sum / ct << std::endl;
    f << "CCR," << k << "," << (unsigned long long)cnt << "," << ALLTIME << "," << pP.pt1 << "," << pP.pt2 << "," << p << "," << get_index_mem()["pk"] << std::endl;
    // f << "CCR," << k << "," << (unsigned long long)cnt << "," << ALLTIME << "," << get_index_mem()["pk"] << std::endl;
    f.close();
    return 0;
}