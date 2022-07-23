#include <cstddef>
#include <cstring>
#include <iostream>
#include <omp.h>
#include <stdlib.h>
#include <string>
#include "src/index.h"
#include <chrono>
#include <metis.h>
#include <unordered_map>
#include <vector>
#include <sstream>
#include <kaHIP_interface.h>
void convert2CSR(std::vector<std::vector<unsigned>>& adj, unsigned* xadj, unsigned* adjncy){
    xadj[0] = 0;
    for(unsigned i=0; i< adj.size(); i++){
       xadj[i+1] = xadj[i] + adj[i].size();
       memcpy((char*)(adjncy + xadj[i]), (char*)adj[i].data(), sizeof(unsigned)* adj[i].size());
    }
}

void MetisRecursive(Index& index){
    idx_t nvtxs = index.nd;
    idx_t ncon = 1;
    idx_t* xadj = new idx_t[nvtxs+1];
    idx_t* adjncy = new idx_t[index.E];
    xadj[0] = 0;
    for(unsigned i=0; i<nvtxs;i++){
        xadj[i+1] = xadj[i] + index.undirect_graph[i].size();
        auto s = xadj[i];
        for(auto &n: index.undirect_graph[i]){
            adjncy[s++] = n;
        }
    }
    idx_t nparts = index._partition_number;
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 0;
    options[METIS_OPTION_UFACTOR] = 1;
    options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;
    options[METIS_OPTION_CTYPE] = METIS_CTYPE_RM;
    options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_GROW;
    options[METIS_OPTION_RTYPE] = METIS_RTYPE_GREEDY;
    options[METIS_OPTION_CONTIG] = 1;
    idx_t *part = new idx_t[index.nd];
    idx_t edge_cut = 0;
    idx_t npart = (unsigned)index._partition_number;
    std::cout << "test" << std::endl;
    METIS_PartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, NULL, NULL, NULL, &npart, NULL, NULL, options,&edge_cut, part);
    std::unordered_map<unsigned, std::vector<unsigned>> _part;
    for(unsigned i=0; i<nvtxs; i++){
        if(_part.count(part[i])){
            _part[part[i]].push_back(i);
        }else
        {
            _part[part[i]] = std::vector<unsigned>();
        }
    }
    index._partition.resize(index._partition_number);
    unsigned s = 0;
    for(auto &n:_part){
       index._partition[s++] = n.second; 
    }
    index.partition_statistic();
    index.save_partition("metis.pa");
}

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

    bool debug = false;
    if(atoi(argv[9]) == 1){
        debug = true;
    }    
    omp_set_num_threads(thread_nums);
    Index index(dim, nums, index_name, data_type, true, debug);
    double e = omp_get_wtime();
    // index.GreedPartition();
    index.graph_partition(partition_name, ivf_times);
    double ss = omp_get_wtime() -e;
    index.partition_statistic();
    std::cout << "METIS RECURSIVE TIME: " <<ss << std::endl;
    index.re_id2pid();
    index.save_partition(partition_name);
    // index.graph_partition(partition_name, ivf_times, descent_times);
   //  index._partition.resize(index._partition_number);
   //  std::cout << "test" << std::endl;
   //  std::cout<<index.E<<std::endl;

   //  int *xadj = new int[index._nd + 1];
   //  int *adjncy = new int[index.E];
   //  int *part = new int[index._nd*2];
   //  xadj[0] = 0;
   //  for(unsigned i=0; i< index.undirect_graph.size(); i++){
   //     xadj[i+1] = xadj[i] + index.undirect_graph[i].size();
   //     unsigned s = xadj[i];
   //     for(auto n: index.undirect_graph[i]){
   //      adjncy[s++] = n;
   //     }
   //  }
   //  std::cout << "start npart" << std::endl;
   //  int nd = index._nd;
   //  int npart = index._partition_number;
   //  int edge_cut;
   //  int ncon = 1;
   //  double imbalanced  = 0;
   //  // METIS_PartGraphKway(&nd,&ncon , xadj, adjncy, NULL, NULL, NULL, &npart,NULL, NULL, NULL, &edge_cut, part);
   //  kaffpa(&nd, NULL, xadj, NULL, adjncy, &npart, &imbalanced, false, 0, FAST, &edge_cut, part); 

   //  auto e = omp_get_wtime() - s;
   //  std::cout << "partition time: " <<e << std::endl;
   //  std::cout << "end partition." << std::endl;

   //  std::unordered_map<unsigned, std::vector<unsigned>> tmp;
   //  for(unsigned i=0; i< index._nd; i++){
   //     tmp[part[i]].push_back(i);
   //  }
   //  index._partition.clear();
   //  for(auto n: tmp){
   //     index._partition.push_back(n.second); 
   //  }
   //  index.partition_statistic();
    return 0;
}
