#include <cstddef>
#include <cstring>
#include <iostream>
#include <omp.h>
#include <stdlib.h>
#include <string>
#include "src/index.h"
#include <chrono>
#include <unordered_map>
#include <vector>
#include <sstream>

int main(int argc, char** argv)
{
    auto s = omp_get_wtime();
    int dim = atoi(argv[1]);
    int nums = atoi(argv[2]);
    char* index_name = argv[3];
    char* data_type = argv[4];
    char* partition_name = argv[5];
    char* neigh_freq_file = argv[6];
    char* org_vamana_file = argv[7];
    char* rearr_vamana_file = argv[8];

    int thread_nums = atoi(argv[9]);
    int ivf_times = atoi(argv[10]);
    int descent_times = atoi(argv[11]);
    int strategy_num = atoi(argv[12]);
    char* neigh_freq_file = argv[13];
    int BS = 1;
    bool debug = false;
    if(atoi(argv[14]) == 1){
        debug = true;
    }    
    if(argc == 15){
        BS = atoi(argv[11]);
    }
    omp_set_num_threads(thread_nums);
    Index index(dim, nums, index_name, data_type, true, debug, false, BS);
    if (strategy_num != 0){
        index.rearrange_vamana(org_vamana_file, neigh_freq_file, rearr_vamana_file );   
    }
    double e = omp_get_wtime();
    // index.GreedPartition();
    index.graph_partition(partition_name, ivf_times, strategy_num, rearr_vamana_file);
    double ss = omp_get_wtime() -e;
    index.partition_statistic();
    std::cout << "TIME: " <<ss << std::endl;
    index.re_id2pid();
    index.save_partition(partition_name);
    return 0;
}
