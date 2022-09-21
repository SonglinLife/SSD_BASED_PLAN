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

    int thread_nums = atoi(argv[6]);
    int ivf_times = atoi(argv[7]);
    int descent_times = atoi(argv[8]);
    int BS = 1;
    bool debug = false;
    if(atoi(argv[9]) == 1){
        debug = true;
    }    
    if(argc == 11){
        BS = atoi(argv[10]);
    }
    omp_set_num_threads(thread_nums);
    Index index(dim, nums, index_name, data_type, true, debug, false, BS);
    double e = omp_get_wtime();
    // index.GreedPartition();
    index.graph_partition(partition_name, ivf_times);
    double ss = omp_get_wtime() -e;
    index.partition_statistic();
    std::cout << "TIME: " <<ss << std::endl;
    index.re_id2pid();
    index.save_partition(partition_name);
    return 0;
}
