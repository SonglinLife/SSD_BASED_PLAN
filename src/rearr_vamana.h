#ifndef SSD_BASED_PLAN_REARR_VAMANA_H
#define SSD_BASED_PLAN_REARR_VAMANA_H

#include <atomic>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <numeric>
#include <ostream>
#define INF 0xffffffff
#include <shared_mutex>
#include <queue>
#include <memory>
#include <set>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <cstring>
#include <ctime>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <utility>
#include <omp.h>
#include <cmath>
#include <mutex>
#include <queue>
#include <random>
#include <oneapi/tbb/concurrent_queue.h>

typedef unsigned long int _u64;
typedef unsigned long int _u64;
using VecT = uint8_t;

struct CompareByFirst {
        constexpr bool operator()(std::pair<unsigned, unsigned> const &a,
                                    std::pair<unsigned, unsigned> const &b) const noexcept {
            return a.second < b.second;
        }
    }; 
       
void rearrangement_neighbor(int nd, int _width,
    std::vector<std::priority_queue<std::pair<unsigned, unsigned>, std::vector<std::pair<unsigned, unsigned>>, CompareByFirst>> *visited_neighbor_count,
    std::vector<std::vector<unsigned>> *_final_graph) {
    
    for (int i = 0; i < nd; i++) {
        if (!(*visited_neighbor_count[i]).empty()) {
            std::unordered_map <unsigned, bool> neighbor_tmp;
            for (auto& each_neigh: (*_final_graph)[i]) {
                neighbor_tmp[each_neigh] = true;
            }
            (*_final_graph)[i].clear();
            while (!(*visited_neighbor_count)[i].empty()) {
                (*_final_graph)[i].push_back((*visited_neighbor_count)[i].top().first);
                neighbor_tmp[(*visited_neighbor_count)[i].top().first] = false;
                (*visited_neighbor_count)[i].pop();
            }
            for (auto& tmp_res : neighbor_tmp) {
                if (tmp_res.second) {
                    (*_final_graph)[i].push_back(tmp_res.first);
                }
            }
        }
        assert(_final_graph[i].size() == _width);
    }

}
// This is a function for read frequency
void stat_neighbor_frequency(const char *freq_file, int nd,
            std::priority_queue<std::pair<unsigned, unsigned>, std::vector<std::pair<unsigned, unsigned>>, CompareByFirst> *visited_freq,
            std::vector<std::priority_queue<std::pair<unsigned, unsigned>, std::vector<std::pair<unsigned, unsigned>>, CompareByFirst> > *visited_neighbor_count) {

    std::ifstream reader(freq_file, std::ios::binary | std::ios::out);
    std::cout << "read visited neighbors information: " << freq_file << std::endl;
    unsigned num = 0;
    reader.read((char *)&num, 4);
    (*visited_neighbor_count).resize(num);
    unsigned tmp_count = 0;
    unsigned n_size = 0;
    assert(num == nd);
    for (size_t i = 0; i < num; i++){
        unsigned v_freq = 0;
        reader.read((char *)v_freq, sizeof(unsigned));
        (*visited_freq).emplace(std::make_pair(i, v_freq));
    }
    for (size_t i = 0; i < num; i++){
        reader.read((char *) &n_size, sizeof(unsigned));
        for  (size_t j = 0; j < n_size; j++){
            unsigned nbr_id, visit_nbr_id_freq;
            reader.read((char *)&nbr_id, sizeof(unsigned));
            reader.read((char *)&visit_nbr_id_freq, sizeof(unsigned));
            (*visited_neighbor_count)[i].emplace(std::make_pair(visit_nbr_id_freq, nbr_id));                
        }
    }
        
}
void save_vamana_index(const char *rearr_vamana, std::vector<std::vector<unsigned>> *_final_graph,
                            unsigned nd, unsigned _width, unsigned _ep) {
    long long     total_gr_edges = 0;
    size_t        index_size = 0;
    std::ofstream out;
    out.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    out.open(std::string(rearr_vamana), std::ios::binary | std::ios::out);
    out.write((char *) &index_size, sizeof(uint64_t));
    out.write((char *) &_width, sizeof(unsigned));
    out.write((char *) &_ep, sizeof(unsigned));
    for (unsigned i = 0; i < nd; i++) {
        unsigned GK = (unsigned) (*_final_graph)[i].size();
        out.write((char *) &GK, sizeof(unsigned));
        out.write((char *) (*_final_graph)[i].data(), GK * sizeof(unsigned));
        total_gr_edges += GK;
    }
    index_size = out.tellp();
    out.seekp(0, std::ios::beg);
    out.write((char *) &index_size, sizeof(uint64_t));
    out.close();

    std::cout << "Avg degree: "
                << ((float) total_gr_edges) /
                        ((float) (nd))
                << std::endl;
}
#endif

