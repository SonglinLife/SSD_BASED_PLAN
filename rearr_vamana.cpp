#include "src/rearr_vamana.h"
#include <iostream>
#include <omp.h>
#include <string>
#include "src/rearr_vamana.h"
#include "src/index.h"
#include <chrono>
#include <sys/types.h>

int main(int argc, char ** argv){
    char *origin_vamana = argv[1];
    char *freq_file = argv[2];
    char *rearr_vamana = argv[3];
    
    int nd = atoi(argv[4]);
    int dim = atoi(argv[5]);

    std::priority_queue<std::pair<unsigned, unsigned>, std::vector<std::pair<unsigned, unsigned>>,
                    CompareByFirst> visited_freq;
    std::vector<std::priority_queue<std::pair<unsigned, unsigned>, std::vector<std::pair<unsigned, unsigned>>,
            CompareByFirst> > visited_neighbor_count;
    std::vector<std::vector<unsigned>> final_graph; // for generate rearranged_vamana, a neighbor list
    
   
    Index index(dim, nd);
    
    index.load_vamana(origin_vamana);
    final_graph = index.direct_graph;
    unsigned width = index._width;
    unsigned ep = index._ep;

    stat_neighbor_frequency(freq_file, nd, &visited_freq, &visited_neighbor_count);
    rearrangement_neighbor(nd, width, &visited_neighbor_count, &final_graph);
    save_vamana_index(rearr_vamana, &final_graph, nd, width, ep);

}