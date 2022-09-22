//
// Created by Zilliz_wmz on 2022/6/6.
//

#ifndef SSD_BASED_PLAN_INDEX_H
#define SSD_BASED_PLAN_INDEX_H
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
using concurrent_queue= oneapi::tbb::concurrent_bounded_queue<unsigned>;
#define READ_U64(stream, val) stream.read((char *)&val, sizeof(_u64))
#define ROUND_UP(X, Y) \
    (((uint64_t)(X) / (Y)) + ((uint64_t)(X) % (Y) != 0)) * (Y)
#define SECTOR_LEN (_u64) 4096
#define DEBUG true

typedef unsigned long int _u64;
typedef unsigned long int _u64;
using VecT = uint8_t;

struct pnode
{
    float cnt;
    unsigned id;
    pnode(float cnt, unsigned id) : cnt(cnt), id(id) {}
    friend bool operator<(pnode a, pnode b)
    {
        return a.cnt < b.cnt;
    }
};
struct ListNode{
    unsigned id;
    unsigned num;
    ListNode* next;
    ListNode* pre;
    ListNode(unsigned i){
        id = i;
        next = nullptr;
        pre = nullptr;
        num = 1;
    }
};
void AddListNode(ListNode* A, ListNode* B){
    B->next = A->next;
    if(A->next)
        A->next->pre = B;
    A->next = B;
    B->pre = A;
}

ListNode* EraseListNode(ListNode* A){
    if(A->pre)
        A->pre->next = A->next;
    if(A->next)
        A->next->pre = A->pre;
    A->next=nullptr;
    A->pre =nullptr;
    return A;
}
inline size_t get_file_size(const std::string &fname)
{
    std::ifstream reader(fname, std::ios::binary | std::ios::ate);
    if (!reader.fail() && reader.is_open())
    {
        size_t end_pos = reader.tellg();
        reader.close();
        return end_pos;
    }
    else
    {
        std::cout << "Could not open file: " << fname << std::endl;
        return 0;
    }
}

typedef std::priority_queue<pnode, std::vector<pnode>> ppq;
class Index
{
public:
    // void debug()
    // {
    //     if (debug)
    //     {
    //         std::cout << "debug" << std::endl;
    //     }
    // }
    Index(const size_t dimension, const size_t n) : _dim(dimension), nd(n)
    {
    }
    Index(const size_t dimension, const size_t n, const char *indexName, const char *data_type = "uint8", bool load_disk = false, bool debug = false, bool sample = false, int BS=1)
    {

        this->debug = debug;
        std::srand(static_cast<unsigned int>(std::time(nullptr)));
        std::ifstream index_meta(indexName, std::ios::binary);
        size_t actual_size = get_file_size(indexName);
        size_t expected_size;
        index_meta.read((char *)&expected_size, sizeof(size_t));
        if (actual_size != expected_size)
        {
            std::cout << "index file not match!" << std::endl;
            exit(-1);
        }
        index_meta.close();
        _dim = dimension;
        nd = n;
        // load_vamana(indexName, sample);
        _rd = new std::random_device();
        _gen = new std::mt19937((*_rd)());
        _dis = new std::uniform_real_distribution<>(0, 1);
        if (load_disk)
        {
            if (std::string(data_type) == std::string("uint8"))
            {
                load_disk_index<uint8_t>(indexName, BS);
            }
            else if (std::string(data_type) == std::string("float"))
            {
                load_disk_index<float>(indexName, BS);
            }
            else
            {
                std::cout << "not support type" << std::endl;
                exit(-1);
            }
        }
        else
        {
            load_vamana(indexName,sample);
        }
        cursize = nd / 1000;
        for(unsigned i =0; i<nd; i++){
            free_node.insert(i);
        }
    }

    ~Index()
    {
        for (auto p : pmutex)
        {
            delete p;
        }
    };
    const bool Debug = true;
    void cout_step()
    {
        if (!debug)
        {
            return;
        }

#pragma omp atomic
        cur++;
        if ((cur + 0) % cursize == 0)
        {
            std::cout << (double)(cur + 0) / nd * 100 << "%    \r";
            std::cout.flush();
        }
    }
    /**
     * load vamana graph index from disk
     * @param filename
     */
    void load_vamana(const char *filename, bool sample = false)
    {
        std::cout << "Reading index file: " << filename << "... " << std::flush;
        std::ifstream in;
        in.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        try
        {
            in.open(filename, std::ios::binary);
            size_t expected_file_size;
            in.read((char *)&expected_file_size, sizeof(uint64_t));
            in.read((char *)&_width, sizeof(unsigned));
            in.read((char *)&_ep, sizeof(unsigned));
            std::cout << "Loading vamana index " << filename << "..."
                      << std::flush;

            size_t cc = 0;
            unsigned nodes = 0;
            while (in.peek() != EOF)
            {
                unsigned k;
                in.read((char *)&k, sizeof(unsigned));
                cc += k;
                ++nodes;
                std::vector<unsigned> tmp(k);
                in.read((char *)tmp.data(), k * sizeof(unsigned));
                direct_graph.emplace_back(tmp);
                if (nodes % 10000000 == 0)
                    std::cout << "." << std::flush;
            }
            if (direct_graph.size() != nd)
            {
                std::cout << "graph vertex size error!\n";
                exit(-1);
            }
            if (sample)
            {
                std::cout << "cut adj" << std::endl;
                for(unsigned i=0; i< nd; i++){
                    std::vector<unsigned> tmp;
                    tmp.reserve(10);
                    for(unsigned j=0; j<20 && j<direct_graph[i].size(); j++){
                        tmp.push_back(direct_graph[i][j]);
                    }
                    direct_graph[i].clear();
                    direct_graph[i].assign(tmp.begin(), tmp.end());
                }
            }
            C = 12;
            _partition_number = ROUND_UP(nd, C) / C;
            reverse_graph.resize(nd);
            std::vector<std::mutex> ms(nd);
#pragma omp parallel for shared(reverse_graph, direct_graph)
            for (unsigned i = 0; i < nd; i++)
            {
                for (unsigned j = 0; j < direct_graph[i].size(); j++)
                {
                    std::lock_guard<std::mutex> lock(ms[direct_graph[i][j]]);
                    reverse_graph[direct_graph[i][j]].emplace_back(i);
                    
                }
            }
            std::cout << "done. Index has " << nodes << " nodes and " << cc
                      << " out-edges" << std::endl;
            en = cc;
            for (int i = 0; i < _partition_number; i++)
            {
                pmutex.push_back(new std::mutex);
            }
        }
        catch (std::system_error &e)
        {
            exit(-1);
        }
    }
    template <typename T>
    void load_disk_index(const char *index_name, int BS=1)
    {
        std::cout << "loading disk index file: " << index_name << "... " << std::flush;
        std::ifstream in;
        in.exceptions(std::ifstream::failbit | std::ifstream::badbit);

        try
        {
            // get metadata
            in.open(index_name, std::ios::binary);
            std::unique_ptr<_u64[]> meta_data = std::make_unique<_u64[]>(12);
            memset(meta_data.get(), 0, 12 * sizeof(_u64));
            in.read((char *)meta_data.get(), 10 * sizeof(_u64));
            in.seekg(SECTOR_LEN, std::ios::beg);
            _u64 actual_file_size = get_file_size(index_name);
            if (actual_file_size != meta_data[0])
            {
                std::cout << "file size mismatch!" << std::endl;
                exit(-1);
            }
            if (nd != meta_data[1])
            {
                std::cout <<"nd: "<<nd<<"meta_data: "<<meta_data[1]<< " npts mismatch" << std::endl;
                exit(-1);
            }
            max_node_len = meta_data[3];
            C = meta_data[4];
            _partition_number = ROUND_UP(nd, C) / C;

            std::unique_ptr<char[]> mem_index = std::make_unique<char[]>(_partition_number * SECTOR_LEN);
            in.read(mem_index.get(), _partition_number * SECTOR_LEN);
            direct_graph.resize(nd);
            _u64 des = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+: des)
            for (unsigned i = 0; i < _partition_number; i++)
            {
                std::unique_ptr<char[]> sector_buf = std::make_unique<char[]>(SECTOR_LEN);
                memcpy(sector_buf.get(), mem_index.get() + i * SECTOR_LEN, SECTOR_LEN);
                for (int j = 0; j < C && i * C + j < nd; j++)
                {
                    std::unique_ptr<char[]> node_buf = std::make_unique<char[]>(max_node_len);
                    memcpy(node_buf.get(), sector_buf.get() + j * max_node_len, max_node_len);
                    unsigned &nnbr = *(unsigned *)(node_buf.get() + _dim * sizeof(T));
                    unsigned *nhood_buf = (unsigned *)(node_buf.get() + (_dim * sizeof(T)) + sizeof(unsigned));
                    std::vector<unsigned> tmp(nnbr);
                    des += nnbr;
                    memcpy((char *)tmp.data(), nhood_buf, nnbr * sizeof(unsigned));
                    direct_graph[i * C + j].assign(tmp.begin(), tmp.end());
                }
            }
            std::cout << "avg degree: " << (double)des / nd << std::endl;
            mem_index.reset();
            std::vector<std::mutex> ms(nd);
            reverse_graph.resize(nd);
#pragma omp parallel for shared(reverse_graph, direct_graph)
            for (unsigned i = 0; i < nd; i++)
            {
                for (unsigned j = 0; j < direct_graph[i].size(); j++)
                {
                    std::lock_guard<std::mutex> lock(ms[direct_graph[i][j]]);
                    reverse_graph[direct_graph[i][j]].emplace_back(i);
                }
            }
            undirect_graph.resize(nd);
            E = 0;
#pragma omp parallel for schedule(dynamic, 100)
            for(unsigned i=0; i< nd; i++){
                std::set<unsigned> ne;
                for(auto n: direct_graph[i]){
                    ne.insert(n);
                }
                for(auto n:reverse_graph[i]){
                    ne.insert(n);
                }
                for(auto n:ne){
                    undirect_graph[i].push_back(n);
                }
#pragma omp atomic
                E += undirect_graph[i].size();
            }
            for(unsigned i=0; i< nd; i++){
                
            }
            C = (SECTOR_LEN*BS) / max_node_len;
            _partition_number = ROUND_UP(nd, C) / C;
            std::cout << "nd: " << nd << " _dim:" << _dim << " C:" << C << " pn:" << _partition_number << std::endl;
            for (int i = 0; i < _partition_number; i++)
            {
                pmutex.push_back(new std::mutex);
            }
            std::cout << "load index over." << std::endl;
        }
        catch (std::system_error &e)
        {
            std::cout << "open file " << index_name << " error!" << std::endl;
            exit(-1);
        }
    }

    /**
     * load binary format data from disk
     * @tparam T
     * @param filename
     * @param data
     * @param num
     * @param dim
     */
    template <typename T>
    void load_bin_data(const char *filename, T *&data, unsigned &num, unsigned &dim)
    {
        std::ifstream in(filename, std::ios::binary);
        if (!in.is_open())
        {
            std::cerr << "open file error" << std::endl;
            exit(-1);
        }
        in.read((char *)&num, 4);
        in.read((char *)&dim, 4);
        data = new T[(size_t)num * (size_t)dim];

        // in.seekg(0, std::ios::beg);
        for (size_t i = 0; i < num; i++)
        {
            in.read((char *)(data + i * dim), dim * sizeof(T));
        }
        in.close();
    }

    /**
     * load the exact KNNG from disk
     * @param filename
     * @param set_dim
     */
    void load_knng(const char *filename, unsigned set_dim)
    {
        unsigned *data = nullptr;
        unsigned num = 0;
        unsigned dim = 0;
        load_bin_data(filename, data, num, dim);
        for (int i = 0; i < num; i++)
        {
            std::vector<unsigned> tmp;
            for (int j = 0; j < dim; j++)
            {
                tmp.push_back(*(data + i * dim + j));
                if (tmp.size() == set_dim)
                    break;
            }
            direct_graph.emplace_back(tmp);
        }
        delete[] data;
    }

    /**
     * evaluate the neighbor overlap ratio (graph quality) for graph indexes
     * built by an approximate algorithm and brute force
     * @param filename
     */
    void evaluate_graph_quality(const char *filename)
    {
        float graph_quality = 0;
        unsigned *data = nullptr;
        unsigned num = 0;
        unsigned dim = 0;
        load_bin_data(filename, data, num, dim);
        for (int i = 0; i < num; i++)
        {
            float tmp = 0;
            for (int j = 0; j < dim; j++)
            {
                auto result = std::find(direct_graph[i].begin(), direct_graph[i].end(), *(data + i * dim + j));
                if (result != direct_graph[i].end())
                {
                    tmp += 1;
                }
                if ((j + 1) == direct_graph[i].size())
                    break;
            }
            graph_quality += tmp / (float)direct_graph[i].size();
        }
        std::cout << "graph quality: " << graph_quality / (float)num << std::endl;
    }

    /**
     * the number of id overlaps for two sets (a and b)
     * @param a
     * @param b
     * @return
     */
    static unsigned overlap_number(const std::vector<unsigned> &a, const std::vector<unsigned> &b)
    {
        unsigned count = 0;
        for (unsigned int i : a)
        {
            if (std::find(b.begin(), b.end(), i) != b.end())
            {
                count++;
            }
        }
        return count;
    }

    /**
     * initialize the parameters for graph partitioning
     * @param vector_number_each_block
     */
    void init_partition(unsigned vector_number_each_block)
    {
        C = vector_number_each_block;
        _partition_number = nd / C + ((nd % C == 0) ? 0 : 1);
        _partition.resize(_partition_number);
        for (int i = 0; i < _partition_number; i++)
        {
            pmutex.push_back(new std::mutex);
        }
        cur = 0;
    }

    /**
     * partition the vertex v to a proper partition by LDG algorithm
     * @param v
     * @return
     */
    unsigned vertex_partition(unsigned v)
    {
        std::vector<float> objective_function(_partition_number);
        float cur_obj_score = -1;
        unsigned cur_obj_partition = 0;
        std::map<unsigned, std::pair<unsigned, float>> record_0_score;
        for (size_t i = 0; i < _partition_number; i++)
        {
            if (_partition[i].size() == C)
                continue;
            unsigned o = overlap_number(direct_graph[v], _partition[i]);
            float w = 1 - (float)_partition.size() / (float)C;
            objective_function[i] = (float)o * w;
            if (objective_function[i] > cur_obj_score)
            {
                cur_obj_score = objective_function[i];
                cur_obj_partition = i;
            }
            record_0_score[i] = std::make_pair(o, w);
        }
        if (cur_obj_score == 0)
        {
            for (auto &rs : record_0_score)
            {
                if (rs.second.first == 0 && rs.second.second != 0)
                {
                    cur_obj_partition = rs.first;
                }
            }
        }
        return cur_obj_partition;
    }

    /**
     * save the partition result
     * @tparam T
     * @param filename
     * @param partition
     */
    void save_partition(const char *filename)
    {

        // re_id2pid();
        std::ofstream writer(filename, std::ios::binary | std::ios::out);
        std::cout << "writing bin: " << filename << std::endl;
        writer.write((char *)&C, sizeof(_u64));
        writer.write((char *)&_partition_number, sizeof(_u64));
        writer.write((char *)&nd, sizeof(_u64));
        std::cout << "_partition_num: " << _partition_number << " C: " << C << " nd: " << nd << std::endl;
        for (unsigned i = 0; i < _partition_number; i++)
        {
            auto p = _partition[i];
            unsigned s = p.size();
            writer.write((char *)&s, sizeof(unsigned));
            writer.write((char *)p.data(), sizeof(unsigned) * s);
        }
        std::vector<unsigned> id2pidv(nd);
        for (auto n : id2pid)
        {
            id2pidv[n.first] = n.second;
        }
        writer.write((char *)id2pidv.data(), sizeof(unsigned) * nd);
    }

    /**
     * load partition from disk
     * @param filename
     */
    void load_partition(const char *filename)
    {
        std::ifstream reader(filename, std::ios::binary);
        reader.read((char *)&C, sizeof(_u64));
        reader.read((char *)&_partition_number, sizeof(_u64));
        reader.read((char*)&nd, sizeof(_u64));
        std::cout << "load partition _partition_num: " << _partition_number << ", C: " << C << std::endl;
        _partition.clear();
        auto tmp = new unsigned[C];
        for (unsigned i = 0; i < _partition_number; i++)
        {
            unsigned c;
            reader.read((char *)&c, sizeof(unsigned));
            reader.read((char *)tmp, c * sizeof(unsigned));
            std::vector<unsigned> tt;
            tt.reserve(C);
            for (int j = 0; j < c; j++)
            {
                tt.push_back(*(tmp + j));
            }
            _partition.push_back(tt);
        }
        delete []tmp;
        re_id2pid();
    }
    void re_id2pid()
    {

        id2pid.clear();
        for (int i = 0; i < _partition_number; i++)
        {
            for (int j = 0; j < _partition[i].size(); j++)
            {
                id2pid[_partition[i][j]] = i;
            }
        }
    }
    /**
     * graph partitioning by LDG algorithm
     * @param filename
     */
    void graph_partition_LDG(const char *filename)
    {
        int cur = 0;
        std::vector<std::mutex> pmutex(_partition_number);
#pragma omp parallel
        {
#pragma omp for schedule(dynamic, 100)
            for (size_t i = 0; i < nd; i++)
            {
                unsigned partition_id = vertex_partition(i);
                pmutex[partition_id].lock();
                while (_partition[partition_id].size() == C)
                {
                    pmutex[partition_id].unlock();
                    partition_id = vertex_partition(i);
                    pmutex[partition_id].lock();
                }
                _partition[partition_id].emplace_back(i);
                pmutex[partition_id].unlock();
                cur++;
                if ((cur + 1) % (100) == 0)
                {
                    std::cout << (double)(cur + 1) / nd * 100 << "%    \r";
                    std::cout.flush();
                }
            }
        }
        save_partition(filename);
    }

    /**
     * count the id overlap according to the graph partitioning
     */
    void partition_statistic()
    {
        std::vector<unsigned> overlap(nd, 0);
        std::vector<unsigned> blk_neighbor_overlap(_partition_number, 0);
        double overlap_ratio = 0;
        unsigned cluter_count = 0;

#pragma omp parallel for schedule(dynamic,100) reduction(+ : overlap_ratio)
        for (size_t i = 0; i < _partition_number; i++)
        {
            std::unordered_set<unsigned> neighbors;
            unsigned blk_neighbor_num = 0;
            for (size_t j = 0; j < _partition[i].size(); j++)
            {
                blk_neighbor_num += direct_graph[_partition[i][j]].size();
                std::unordered_set<unsigned> ne;
                for (unsigned &x : direct_graph[_partition[i][j]])
                {
                    neighbors.insert(x);
                    ne.insert(x);
                }
                blk_neighbor_overlap[i] = blk_neighbor_num - neighbors.size();
                for (size_t z = 0; z < _partition[i].size(); z++)
                {
                    if (_partition[i][j] == _partition[i][z])
                        continue;
                    if(ne.find(_partition[i][z])!=ne.end()){
                        overlap[_partition[i][j]]++;
                    }
                    // for (size_t x = 0; x < direct_graph[_partition[i][j]].size(); x++)
                    // {
                    //     if (_partition[i][z] == direct_graph[_partition[i][j]][x])
                    //     {
                    //         overlap[_partition[i][j]]++;
                    //         break;
                    //     }
                    // }
                }
                overlap_ratio += (_partition[i].size() == 1 ? 0 : (1.0 * overlap[_partition[i][j]] / (_partition[i].size() - 1)));
            }
        }
        unsigned max_overlaps = 0;
        unsigned min_overlaps = std::numeric_limits<unsigned>::max();
        double ave_overlap_ratio = 0;
        std::map<unsigned, unsigned> overlap_count;
        for (size_t i = 0; i < nd; i++)
        {
            if (overlap_count.count(overlap[i]))
            {
                overlap_count[overlap[i]]++;
            }
            else
            {
                overlap_count[overlap[i]] = 1;
            }
            if (overlap[i] > max_overlaps)
                max_overlaps = overlap[i];
            if (overlap[i] < min_overlaps)
                min_overlaps = overlap[i];
        }
        ave_overlap_ratio = overlap_ratio / (double)nd;
        for (auto &it : overlap_count)
        {
            std::cout << "each id, overlap number " << it.first << ", count: " << it.second << std::endl;
        }
        std::cout << "each id, max overlaps: " << max_overlaps << std::endl;
        std::cout << "each id, min overlaps: " << min_overlaps << std::endl;
        std::cout << "each id, average overlap ratio: " << ave_overlap_ratio << std::endl;

        // unsigned max_neighbor_overlaps = 0;
        // unsigned min_neighbor_overlaps = std::numeric_limits<unsigned>::max();
        // float ave_neighbor_overlaps = 0;
        // std::map<unsigned, unsigned> neighbor_overlap_count;
        // for (size_t i = 0; i < _partition_number; i++)
        // {
        //     ave_neighbor_overlaps += (float)blk_neighbor_overlap[i];
        //     if (neighbor_overlap_count.count(blk_neighbor_overlap[i]))
        //     {
        //         neighbor_overlap_count[blk_neighbor_overlap[i]]++;
        //     }
        //     else
        //     {
        //         neighbor_overlap_count[blk_neighbor_overlap[i]] = 1;
        //     }
        //     if (blk_neighbor_overlap[i] > max_neighbor_overlaps)
        //         max_neighbor_overlaps = blk_neighbor_overlap[i];
        //     if (blk_neighbor_overlap[i] < min_neighbor_overlaps)
        //         min_neighbor_overlaps = blk_neighbor_overlap[i];
        // }
        // for (auto &it : neighbor_overlap_count)
        // {
        //     std::cout << "each block, neighbor overlap number " << it.first << ", count: " << it.second << std::endl;
        // }
        // ave_neighbor_overlaps /= (float)_partition_number;
        // std::cout << "each block, max neighbor overlaps: " << max_neighbor_overlaps << std::endl;
        // std::cout << "each block, min neighbor overlaps: " << min_neighbor_overlaps << std::endl;
        // std::cout << "each block, average neighbor overlaps: " << ave_neighbor_overlaps << std::endl;
    }

    /**
     * load bbann cluster information from disk
     * @param filename
     * @param bucket_num
     * @param cur_k1
     */
    void load_cluster(const char *filename, size_t &bucket_num, unsigned cur_k1)
    {
        std::ifstream in(filename, std::ios::binary);
        if (!in.is_open())
        {
            std::cout << "open file error" << std::endl;
            exit(-1);
        }
        in.seekg(0, std::ios::end);
        std::ios::pos_type ss = in.tellg();
        size_t fsize = (size_t)ss;
        bucket_num = (size_t)(fsize / 4096);
        _cluster_size[cur_k1] = bucket_num;
        _total_cluster_size += bucket_num;
        _cluster_data[cur_k1] = new char[bucket_num * 4096];
        in.seekg(0, std::ios::beg);
        in.read((char *)_cluster_data[cur_k1], bucket_num * 4096);
        in.close();
    }

    /**
     * initialize the parameters of bbann cluster
     */
    void init_cluster_data()
    {
        _cluster_data.resize(_k1);
        _cluster_size.resize(_k1);
    }

    /**
     * bbann first level kmeans parameter k1
     * @param k1
     */
    void set_k1(unsigned k1)
    {
        _k1 = k1;
    }

    /**
     * count the id overlap according to the bbann cluster
     */
    void point_statistic()
    {
        std::vector<unsigned> overlap(nd, 0);
        std::vector<unsigned> blk_overlap(_total_cluster_size, 0);
        std::vector<unsigned> blk_neighbor_overlap(_total_cluster_size, 0);
        double overlap_ratio = 0;
        double blk_overlap_ratio = 0;
        unsigned vector_size = sizeof(VecT) * _dim;
        unsigned id_size = sizeof(unsigned);
        unsigned cluter_count = 0;
        for (size_t i = 0; i < _k1; i++)
        {
            char *cd = _cluster_data[i];
            for (size_t j = 0; j < _cluster_size[i]; j++)
            {
                char *cur_bucket = cd;
                cd += 4096;
                // memcpy(cur_bucket, cd, 4096);
                unsigned bucket_size = *(unsigned *)cur_bucket;
                std::vector<unsigned> bucket_element(bucket_size);
                cur_bucket += sizeof(unsigned);
                for (size_t k = 0; k < bucket_size; k++)
                {
                    cur_bucket += vector_size;
                    bucket_element[k] = *(unsigned *)cur_bucket;
                    cur_bucket += id_size;
                }
                std::unordered_set<unsigned> neighbors;
                unsigned blk_neighbor_num = 0;
                for (size_t k = 0; k < bucket_size; k++)
                {
                    blk_neighbor_num += direct_graph[bucket_element[k]].size();
                    for (unsigned int &x : direct_graph[bucket_element[k]])
                    {
                        neighbors.insert(x);
                    }
                }
                blk_neighbor_overlap[cluter_count] = blk_neighbor_num - neighbors.size();
                for (size_t k = 0; k < bucket_size; k++)
                {
                    if (neighbors.count(bucket_element[k]))
                        blk_overlap[cluter_count]++;
                }
                blk_overlap_ratio += (bucket_size == 0 ? 0 : (1.0 * blk_overlap[cluter_count] / bucket_size));
                cluter_count++;
                for (size_t k = 0; k < bucket_size; k++)
                {
                    for (size_t z = 0; z < bucket_size; z++)
                    {
                        if (k == z)
                            continue;
                        // if (bucket_element[k] == 10)
                        // std::cout << "bucket number:" << bucket_element[z] << " " << direct_graph[bucket_element[k]].size() << std::endl;
                        for (size_t x = 0; x < direct_graph[bucket_element[k]].size(); x++)
                        {
                            // if (bucket_element[k] == 1) {
                            //     std::cout << "neighbors: " << direct_graph[k][x] << " ";
                            // }
                            if (bucket_element[z] == direct_graph[bucket_element[k]][x])
                            {
                                overlap[bucket_element[k]]++;
                                break;
                            }
                        }
                        // std::cout << std::endl;
                    }
                    overlap_ratio += (bucket_size == 1 ? 0 : (1.0 * overlap[bucket_element[k]] / (bucket_size - 1)));
                }
                // delete cur_bucket;
            }
        }

        unsigned max_overlaps = 0;
        unsigned min_overlaps = std::numeric_limits<unsigned>::max();
        double ave_overlap_ratio = 0;
        std::map<unsigned, unsigned> overlap_count;
        for (size_t i = 0; i < nd; i++)
        {
            if (overlap_count.count(overlap[i]))
            {
                overlap_count[overlap[i]]++;
            }
            else
            {
                overlap_count[overlap[i]] = 1;
            }
            if (overlap[i] > max_overlaps)
                max_overlaps = overlap[i];
            if (overlap[i] < min_overlaps)
                min_overlaps = overlap[i];
        }
        ave_overlap_ratio = overlap_ratio / (double)nd;
        for (auto &it : overlap_count)
        {
            std::cout << "each id, overlap number " << it.first << ", count: " << it.second << std::endl;
        }
        std::cout << "each id, max overlaps: " << max_overlaps << std::endl;
        std::cout << "each id, min overlaps: " << min_overlaps << std::endl;
        std::cout << "each id, average overlap ratio: " << ave_overlap_ratio << std::endl;
        unsigned max_blk_overlaps = 0;
        unsigned min_blk_overlaps = std::numeric_limits<unsigned>::max();
        double ave_blk_overlap_ratio = 0;
        std::map<unsigned, unsigned> blk_overlap_count;
        for (size_t i = 0; i < _total_cluster_size; i++)
        {
            if (blk_overlap_count.count(blk_overlap[i]))
            {
                blk_overlap_count[blk_overlap[i]]++;
            }
            else
            {
                blk_overlap_count[blk_overlap[i]] = 1;
            }
            if (blk_overlap[i] > max_blk_overlaps)
                max_blk_overlaps = blk_overlap[i];
            if (blk_overlap[i] < min_blk_overlaps)
                min_blk_overlaps = blk_overlap[i];
        }
        ave_blk_overlap_ratio = blk_overlap_ratio / (double)_total_cluster_size;
        for (auto &it : blk_overlap_count)
        {
            std::cout << "each block, overlap number " << it.first << ", count: " << it.second << std::endl;
        }
        std::cout << "each block, max overlaps: " << max_blk_overlaps << std::endl;
        std::cout << "each block, min overlaps: " << min_blk_overlaps << std::endl;
        std::cout << "each block, average overlap ratio: " << ave_blk_overlap_ratio << std::endl;

        unsigned max_neighbor_overlaps = 0;
        unsigned min_neighbor_overlaps = std::numeric_limits<unsigned>::max();
        float ave_neighbor_overlaps = 0;
        std::map<unsigned, unsigned> neighbor_overlap_count;
        for (size_t i = 0; i < _total_cluster_size; i++)
        {
            ave_neighbor_overlaps += (float)blk_neighbor_overlap[i];
            if (neighbor_overlap_count.count(blk_neighbor_overlap[i]))
            {
                neighbor_overlap_count[blk_neighbor_overlap[i]]++;
            }
            else
            {
                neighbor_overlap_count[blk_neighbor_overlap[i]] = 1;
            }
            if (blk_neighbor_overlap[i] > max_neighbor_overlaps)
                max_neighbor_overlaps = blk_neighbor_overlap[i];
            if (blk_neighbor_overlap[i] < min_neighbor_overlaps)
                min_neighbor_overlaps = blk_neighbor_overlap[i];
        }
        for (auto &it : neighbor_overlap_count)
        {
            std::cout << "each block, neighbor overlap number " << it.first << ", count: " << it.second << std::endl;
        }
        ave_neighbor_overlaps /= (float)_total_cluster_size;
        std::cout << "each block, max neighbor overlaps: " << max_neighbor_overlaps << std::endl;
        std::cout << "each block, min neighbor overlaps: " << min_neighbor_overlaps << std::endl;
        std::cout << "each block, average neighbor overlaps: " << ave_neighbor_overlaps << std::endl;
    }
    unsigned select_partition(unsigned i)
    {

#pragma omp atomic
        select_nums++;

        float maxn = 0.0;
        unsigned res = INF;
        std::unordered_map<unsigned, unsigned> pcount;
        unsigned tpid = 0;
        for (auto n : direct_graph[i])
        {
            unsigned pid = id2pid[n];
            if (pid == INF)
                continue;
            pcount[pid] = pcount[pid] + 1;
            if (tpid < pid)
            {
                tpid = pid;
            }
        }
        for (auto n : reverse_graph[i])
        {
            unsigned pid = id2pid[n];
            if (pid == INF)
                continue;
            pcount[pid] = pcount[pid] + 1;
            if (tpid < pid)
            {
                tpid = pid;
            }
        }
        for (auto c : pcount)
        {
            unsigned pid = c.first;
            float cnt = c.second;
            std::lock_guard<std::mutex> lock(*pmutex[pid]);
            double s = _partition[pid].size();
            cnt *= (1 - s / C) ;
            if (cnt > maxn && _partition[pid].size() < C)
            {
                res = pid;
                maxn = cnt;
            }
            // std::lock_guard<std::mutex> lock(*pmutex[pid]);
            // if (_partition[pid].size() == C)
            // {
            //     continue;
            // }
            // std::vector<unsigned> tmp(_partition[pid].begin(), _partition[pid].end());
            // std::unordered_map<unsigned, unsigned> mpp;
            // unsigned s = _partition[pid].size();
            // float start =(float) partition_ratio(tmp, pid, mpp) /( s*(s-1));
            // tmp.push_back(i);
            // float end = (float)partition_ratio(tmp, pid, mpp) / (s *(s+1));
            // p.push({end -start, pid});
        }
        // if (p.empty())
        // {
        //     res = getUnfilled();
        // }
        // else
        // {
        //     res = p.top().id;
        // }
        pcount.clear();
        // if( res == INF-1){
        //     for(auto &n :direct_graph[i]){
        //         for(auto &s :direct_graph[n]){
        //             unsigned pid = id2pid[s];
        //             if(pid == INF){
        //                 continue;
        //             }
        //             pcount[pid] = pcount[pid] + 1;
        //         }
        //     }
        //     for(auto &n :reverse_graph[i]){
        //         for(auto &s :reverse_graph[n]){
        //             unsigned pid = id2pid[s];
        //             if(pid == INF){
        //                 continue;
        //             }
        //             pcount[pid] = pcount[pid] + 1;
        //         }
        //     }
        //     for (auto& c : pcount)
        //     {
        //         unsigned pid = c.first;
        //         float cnt = c.second;
        //         std::lock_guard<std::mutex> lock(*pmutex[pid]);
        //         double s = _partition[pid].size();
        //         cnt *= (1 - s / C);
        //         if (cnt > maxn && _partition[pid].size() < C)
        //         {
        //             res = pid;
        //             maxn = cnt;
        //         }
        //     }
        // }
        
        if (res == INF)
        {
#pragma omp atomic
            select_free++;
            res = getUnfilled();
        }
        return res;
    }
    unsigned getUnfilled(unsigned tpid)
    {
        std::lock_guard<std::mutex> lock(flock);
        for (unsigned i = 0; i < _partition_number; i++)
        {
            unsigned pid = (tpid + i) % _partition_number;
            if (_partition[pid].size() < C)
            {
                return pid;
            }
        }
        return 0;
    }
    unsigned getUnfilled()
    {
        // std::lock_guard<std::mutex> lock(plock);
        // unsigned res = unfill.top().id;
        // unsigned cnt = unfill.top().cnt;
        // unfill.pop();
        // unfill.push({res, cnt++})
#pragma omp atomic
        getUnfilled_nums++;

        unsigned res;

        // std::shared_lock<std::shared_mutex> lock(smutex);
        // unsigned tmp = (*_dis)(*_gen) * unfilled_partition.size();
        // if(tmp > unfilled_partition.size()){
        //     std::cout << "error tmp: "<< tmp << std::endl;
        // }
        // auto r = unfilled_partition.begin();
        // std::advance(r, tmp);
        // res = r->first;
        // for (int i = 0; i < _partition_number; i++)
        // {
        //     int pid = (tmp + i) % _partition_number;
        //     std::lock_guard<std::mutex> lock(*pmutex[pid]);
        //     if (_partition[pid].size() < C)
        //     {
        //         res = pid;
        //         break;
        //     }
        // }
        do
        {
           free_q.pop(res); 
        } while(_partition[res].size() == C);
        return res;
    }
    void bfs(unsigned start)
    {
        std::queue<unsigned> q;
        vlock.lock();
        if (!visited[start])
        {
            q.push(start);
            visited[start] = true;
        }
        vlock.unlock();
        while (!q.empty())
        {
            unsigned i = q.front();
            q.pop();
            unsigned pid = sync(i);
            vlock.lock();
            cout_step();
            for (auto n : direct_graph[i])
            {
                if (!visited[n])
                {
                    q.push(n);
                    visited[n] = true;
                }
            }
            vlock.unlock();
        }
    }
    void graph_partition_bfs()
    {

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0, 1);

        _partition.clear();
        _partition.resize(_partition_number);
        cur = 0;
        visited.clear();
        std::cout << "start" << std::endl;
        // init

        std::set<unsigned> eps;
        while (eps.size() < 64)
        {
            unsigned a = nd * dis(gen) - 1;
            eps.insert(a);
        }
        for (int i = 0; i < _partition_number; i++)
        {
            unfilled_partition[i] = true;
        }
        auto start = omp_get_wtime();
#pragma omp parallel for
        for (int i = 0; i < 64; i++)
        {
            auto s = eps.begin();
            std::advance(s, i);
            unsigned ep = *s;
            bfs(ep);
        }
        std::cout << "bfs" << std::endl;
#pragma omp parallel for
        for (int i = 0; i < nd; i++)
        {
            if (visited[i])
                continue;
            sync(i);
            // vlock.lock();
            // if (!visited[i])
            // {
            //     vlock.unlock();
            //     sync(i);
            // }
            // else
            // {
            //     vlock.unlock();
            // }
        }
        auto end = omp_get_wtime();
        std::cout << "Round " << round << "bfs time: " << end - start << " finished." << std::endl;
    }
    void graph_partition(const char *filename, int k, int dn = 0)
    {
        for (int i = 0; i < nd; i++)
        {
            id2pid[i] = INF;
        }
        _partition.resize(_partition_number);
        if (true)
        {
            std::unordered_set<unsigned> vis;
            for (unsigned i = 0; i < _partition_number; i++)
            {
                for (unsigned j = 0; j < C && i*C +j < nd; j++)
                {
                    id2pid[i * C + j] = i;
                    _partition[i].push_back(i*C+j);
                }
            }
            // auto ivf_file_name = std::string(filename) + std::string(".ivf0");
            // save_partition(ivf_file_name.c_str());
        }
        // load_partition("../bigann/DEEP_10M_R48_L128_B0.225_greed");
        for (int i = 0; i < k; i++)
        {
            select_free=0;
            graph_partition_ivf();
            std::cout << "select free: "<<(double)select_free/_partition_number << std::endl;
            partition_statistic();
            auto ivf_file_name = std::string(filename) + std::string(".ivf") + std::to_string(i + 1);
            std::cout << "total ivf time: " << ivf_time << std::endl;
            save_partition(ivf_file_name.c_str());
            round++;
        }
        std::cout << "select pid nums" << select_nums << " get unfilled partition nums: " << getUnfilled_nums << std::endl;
        std::cout << "total ivf time: " << ivf_time << std::endl;
        std::bitset<10000000>bs;
#pragma omp parallel for
        for (size_t i = 0; i < _partition.size(); ++i) {
            for (auto n : _partition[i]) {
                assert(!bs.test(n));
                bs.set(n);
            }
        }
        if(dn==0)return;
        // save_partition(filename);
        // load_partition(filename);
        // re_id2pid();
        // partition_statistic();
        round = 0;
        for (int i = 0; i < _partition_number; i++)
        {
            partition_ratio(i, true);
        }
        auto ivf_file_name = std::string(filename) + std::string(".ivf");
        save_partition(ivf_file_name.c_str());
        std::cout << "gp_descent" << std::endl;
        for (int i = 0; i < dn; i++)
        {
            cur = 0;
            gp_descent();
            std::cout << "over descent round " << i << std::endl;
            partition_statistic();
            // for (int i = 0; i < _partition_number; i++)
            // {
            //     partition_ratio(i, true);
            // }
            // re_id2pid();
        }
        std::cout << "descent total time: " << de_time << std::endl;
        auto de_file_name = std::string(filename) + std::string(".de");
        save_partition(de_file_name.c_str());
    }
    void graph_partition_ivf()
    {
        _partition.clear();
        _partition.resize(_partition_number);
        cur = 0;
        visited.clear();
        std::cout << "start" << std::endl;
        std::vector<unsigned> stream(nd);
        std::iota(stream.begin(), stream.end(), 0);
        auto rng = std::default_random_engine {};
        std::shuffle(std::begin(stream), std::end(stream), rng);
        for (int i = 0; i < _partition_number; i++)
        {
            // unfilled_partition[i] = true;
            free_q.push(i);
        }
        // for(int i=0; i<nd; i++){
        //     id2pid[i] = INF;
        // }
        auto start = omp_get_wtime();
        size_t offset = (*_dis)(*_gen) * nd;
#pragma omp parallel for schedule(dynamic)
        for (unsigned i = 0; i < nd; i++)
        {
            // debug();
            size_t n = stream[i];
            sync(n);
            // unsigned pid = select_partition(i);
            // pmutex[pid]->lock();
            // while (_partition[pid].size() == C)
            // {
            //     pmutex[pid]->lock();
            //     pid = select_partition(i);
            //     pmutex[pid]->lock();
            // }
            // _partition[pid].emplace_back(i);
            // pmutex[pid]->unlock();
            // std::lock_guard<std::mutex> lock(flock);
            // unfilled_partition.erase(pid);
            // cur++;
            // if ((cur + 1) % (100) == 0)
            // {
            //     std::cout << (double)(cur + 1) / nd * 100 << "%    \r";
            //     std::cout.flush();
            // }
            // debug();
            cout_step();
        }
        auto end = omp_get_wtime();
        std::cout << "ivf time: " << end - start << " round: " << round << std::endl;
        ivf_time += end - start;
        round++;
    }
    unsigned sync(unsigned i)
    {
        unsigned pid = select_partition(i);
        pmutex[pid]->lock();

        while (_partition[pid].size() == C)
        {
            pmutex[pid]->unlock();
            pid = select_partition(i);
            pmutex[pid]->lock();
        }
        _partition[pid].emplace_back(i);
        id2pid[i] = pid;
        unsigned s = _partition[pid].size();
        pmutex[pid]->unlock();

        if (s != C)
        {
            // std::lock_guard<std::mutex> lock(flock);
            // std::unique_lock<std::shared_mutex>lock(smutex);
            // unfilled_partition.erase(pid);
            free_q.push(pid);    
        }

        return pid;
    }
    void gp_descent()
    {
        unsigned cur = 0;
        swap_cnt = 0;
        swap_min = 0;
        double start = omp_get_wtime();
#pragma omp parallel for schedule(dynamic)
        for (unsigned i = 0; i < nd; i++)
        {
            for (unsigned j = 0; j < direct_graph[i].size(); j++)
            {
                for (unsigned k = j + 1; k < direct_graph[i].size(); k++)
                {
                    swap_partition(id2pid[direct_graph[i][j]], id2pid[direct_graph[i][k]]);
                }
            }
            cout_step();
            // #pragma omp atomic
            //             cur++;
            //             if ((cur + 1) % 100 == 0)
            //             {
            //                 std::cout << (double)(cur + 1) / nd * 100 << "%    \r";
            //                 std::cout.flush();
            //             }
        }
        cur = 0;
        // for (unsigned i = 0; i < nd; i++)
        // {
        //     for (unsigned j = 0; j < reverse_graph[i].size(); j++)
        //     {
        //         for (unsigned k = j + 1; k < reverse_graph[i].size(); k++)
        //         {
        //             swap_partition(id2pid[reverse_graph[i][j]], id2pid[reverse_graph[i][k]]);
        //         }
        //     }
        //     cout_step();
        //     // #pragma omp atomic
        //     //             cur++;
        //     //             if ((cur + 1) % 100 == 0)
        //     //             {
        //     //                 std::cout << (double)(cur + 1) / nd * 100 << "%    \r";
        //     //                 std::cout.flush();
        //     //             }
        // }
        std::cout << "swap_cnt: " << swap_cnt << std::endl;
        std::cout << "swap_min: " << swap_min << std::endl;
        auto end = omp_get_wtime();
        std::cout << "round " << round++ << "decent time: " << end - start << std::endl;
        de_time += end - start;
    }

    unsigned partition_ratio(unsigned pid, bool update = false)
    {
        // std::lock_guard<std::mutex> lock(*pmutex[pid]);
        unsigned oldc = 0, newc = 0;
        for (unsigned i = 0; i < _partition[pid].size(); i++)
        {
            unsigned cnt = 0;
            for (auto n : direct_graph[_partition[pid][i]])
            {
                if (id2pid[n] == pid)
                {
                    cnt++;
                }
            }
            oldc += id2ratio[_partition[pid][i]];
            if (update)
            {
                id2ratio[_partition[pid][i]] = cnt;
            }
            newc += cnt;
        }
        // if (oldc != newc)
        //     std::cout << "partition: " << pid << " oldc: " << oldc << " newc: " << newc << std::endl;
        return newc;
    }
    unsigned partition_ratio(std::vector<unsigned> &partition, unsigned pid, std::unordered_map<unsigned, unsigned> &id2r)
    {
        unsigned res = 0;
        for (unsigned i = 0; i < partition.size(); i++)
        {
            unsigned cnt = 0;
            for (auto n : direct_graph[partition[i]])
            {
                if (id2pid[n] == pid)
                {
                    cnt++;
                }
            }
            id2r[partition[i]] = cnt;
            res += cnt;
        }
        return res;
    }
    unsigned get_partition_min_ratio_id(unsigned pid)
    {
        unsigned minn = INF, res = INF;
        for (auto n : _partition[pid])
        {
            if (minn > id2ratio[n])
            {
                minn = id2ratio[n];
                res = n;
            }
        }
        return res;
    }
#pragma oninline
    unsigned __attribute__((noinline)) pid_ratio(unsigned pid)
    {
        std::unordered_map<unsigned, unsigned> tes;
        unsigned cnt = 0;
        for (auto n : _partition[pid])
        {
            tes[n] = id2ratio[n];
            std::cout << n << ": " << tes[n] << " ";
            cnt += tes[n];
        }
        std::cout << std::endl;
        return cnt;
    }
    void swap_partition(unsigned jpid, unsigned kpid)
    {
        if (jpid == kpid)
            return;

        std::scoped_lock lock(*pmutex[kpid], *pmutex[jpid]);
        if (_partition[kpid].size() < C || _partition[jpid].size() < C)
            return;

        unsigned kmin = get_partition_min_ratio_id(kpid);
        unsigned jmin = get_partition_min_ratio_id(jpid);
        std::vector<unsigned> jp, kp;
        jp.emplace_back(kmin);
        kp.emplace_back(jmin);
        unsigned oldr = 0;
        for (auto n : _partition[jpid])
        {
            oldr += id2ratio[n];
            if (n == jmin)
                continue;
            jp.emplace_back(n);
        }
        for (auto n : _partition[kpid])
        {
            oldr += id2ratio[n];
            if (n == kmin)
                continue;
            kp.emplace_back(n);
        }
        std::unordered_map<unsigned, unsigned> mapj, mapk;
        id2pid[kmin] = jpid;
        id2pid[jmin] = kpid;
        unsigned newr = partition_ratio(jp, jpid, mapj) + partition_ratio(kp, kpid, mapk);
        if (newr > oldr)
        {
#pragma omp atomic
            swap_cnt++;
            _partition[jpid] = jp;
            _partition[kpid] = kp;
            for (auto n : jp)
            {
                id2ratio[n] = mapj[n];
            }
            for (auto n : kp)
            {
                id2ratio[n] = mapk[n];
            }
            if (get_partition_min_ratio_id(kpid) == jmin && get_partition_min_ratio_id(jpid) == kmin)
            {
#pragma omp atomic
                swap_min++;
            }
        }
        else
        {
            id2pid[kmin] = kpid;
            id2pid[jmin] = jpid;
        }
    }
    void swap_LDG_partition(unsigned jpid, unsigned kpid)
    {

        if (jpid == kpid)
            return;

        std::scoped_lock lock(*pmutex[kpid], *pmutex[jpid]);
        if (_partition[kpid].size() < C || _partition[jpid].size() < C)
            return;
        std::vector<unsigned> all;
        auto jr = sort_ratio(jpid);
        auto kr = sort_ratio(kpid);
        std::set<unsigned> a, b;
        for (int i = 0; i < C; i++)
        {
            if (i < 6)
            {
                a.insert(jr[i].id);
                b.insert(kr[i].id);
                continue;
            }
            all.push_back(jr[i].id);
            all.push_back(kr[i].id);
        }
        for (auto n : all)
        {
            if (a.size() == C)
            {
                b.insert(n);
                continue;
            }
            else if (b.size() == C)
            {
                a.insert(n);
                continue;
            }
            if (selectOne(n, a) > selectOne(n, b))
            {
                a.insert(n);
            }
            else
            {
                b.insert(n);
            }
        }
        _partition[jpid] = std::vector<unsigned>(a.begin(), a.end());
        _partition[kpid] = std::vector<unsigned>(b.begin(), b.end());
        partition_ratio(jpid, true);
        partition_ratio(kpid, true);
    }
    std::vector<pnode> sort_ratio(unsigned pid)
    {
        std::vector<pnode> a;
        for (auto n : _partition[pid])
            a.push_back({(float)id2ratio[n], n});
        std::sort(a.begin(), a.end());
        return a;
    }
    float selectOne(unsigned id, std::set<unsigned> &a)
    {
        unsigned cnt = 0;
        for (auto n : direct_graph[id])
        {
            if (a.find(n) != a.end())
                cnt++;
        }
        for (auto n : reverse_graph[id])
        {
            if (a.find(n) != a.end())
                cnt++;
        }
        return (float)cnt * (1 - (float)a.size() / C);
    }

    unsigned get_free_node(){
        auto it = free_node.begin();
        if(free_node.end() == it){
            std::cout << "false" << std::endl;
            exit(-1);
        }
        unsigned res = *it;
        return res;
    }

    unsigned search_priority_list(){
        for(unsigned i=C-1; i>=1; i--){
            if(priority_list[i]->next == nullptr){
                continue;
            }
            auto s = priority_list[i]->next;
            while(s ->next != nullptr){
                s = s->next;
            }
            s = EraseListNode(s);
            auto id = s->id;
            delete s;
            return id;
        }
        return INF;         
    }
    void free_list(){
        for(int i=0; i<C; i++){
            auto s = priority_list[i];
            ListNode* p;
            while(s != nullptr){
                p = s->next;
                delete s;
                s = p;
            }
        }
    }
    unsigned print_list(){
        for(unsigned i=C-1; i>=1; i--){
            std::cout << "i: "<<i<< "    ";
            auto s = priority_list[i];
            while(s->next !=nullptr){
                std::cout << " "<< s->next->id;
                s = s->next;
            }
            std::cout << std::endl;
        }
        return INF;
    }
    void update_priority_list(unsigned id, unsigned pid = INF){
        for(auto &n:undirect_graph[id]){
            if(!free_node.count(n))continue;
            if(pid == 1412){
                std::cout << "n: "<< n  << std::endl;
                print_list();
                if(n == 8618452){
                    std::cout << "test" << std::endl;
                }
            }
            // std::cout << n << std::endl;
            if(id2ptr.count(n)){
                auto t = EraseListNode(id2ptr[n]);
                t->num++;
                // std::cout << t->num << " "<<n << std::endl;
                AddListNode(priority_list[t->num], t);
                // print_list();
                if(pid == 1412){
                    print_list();
                }
                
            }else{
                auto t = new ListNode(n);
                id2ptr[n] = t;
                AddListNode(priority_list[1], t);
                if(pid == 1412){
                    print_list();
                }
            }
        }
    }

    std::vector<unsigned> GreedInsert(unsigned pid){
        // if(pid == 1412){
        //     std::cout << "test" << std::endl;
        // }
        for(unsigned i=0; i<C; i++){
           auto tmp = new ListNode(-1);
           priority_list.push_back(tmp); 
        }
        id2ptr.clear();
        std::vector<unsigned> part;
        unsigned init_node = get_free_node();
        free_node.erase(init_node);
        fix_node.insert(init_node);
        std::queue<unsigned> q;
        // print_list();
        for(auto &n: undirect_graph[init_node]){
            if(free_node.count(n)){
                id2ptr[n] = new ListNode(n);
                AddListNode(priority_list[1], id2ptr[n]);
            }
        }    
        // if(pid == 1412){
        //     print_list();
        // }
        // print_list();
        // part.push_back(init_node);
        q.push(init_node);
        for(unsigned i=1; i<C && (pid*C + i)<nd; i++){
            unsigned id = search_priority_list();

        // if(pid == 1412){
        //     print_list();
        // }
            std::queue<unsigned> hop3;
            if(id == INF){
                unsigned last = q.back();
                bool flag = false;
                for(auto &n : undirect_graph[last]){
                    for(auto &s :undirect_graph[n]){
                        hop3.push(s);
                        if(free_node.count(s)){
                            id = s;
                            flag = true;
                            break;
                        }
                    }
                    if(flag){
                        break;
                    }
                }
            }
            // if (id == INF){
            //     while(!hop3.empty()){
            //         unsigned n = hop3.front();
            //         hop3.pop();
            //         bool flag = false;
            //         for(auto &s : undirect_graph[n]){
            //            if(free_node.count(s)){
            //             id = s;
            //             flag = true;
            //             break;
            //            } 
            //         }
            //         if(flag)
            //             break;
            //     }
            // }
            if(id == INF){
                id = get_free_node();
                select_free++;
            }
            free_node.erase(id);
            fix_node.insert(id);
            // part.push_back(id);
            q.push(id);
        if(pid == 1412){
            print_list();
            std::cout << "id "<< id<<"adj list ";
            for(auto &n : undirect_graph[id]){
                std::cout << " "<<n;
            }
            std::cout << "" << std::endl;
        }
            update_priority_list(id, pid);
            if(pid == 1412){
                std::cout << "update over!" << std::endl;
                print_list();
            }
            // print_list();
        }
        // std::cout << "queue:" ;
        while(q.size()>C){
            auto id = q.front();
            // std::cout << " "<<id;
            q.pop();
            free_node.insert(id);
        }
        while(!q.empty()){
            auto id = q.front();
            q.pop();
            // std::cout << " "<<id ;
            part.push_back(id);
        }
        // std::cout<<std::endl;
        // print_list();
        // for(auto &n :id2ptr){
        //     delete n.second;
        // }
        // for(unsigned i=1; i<C; i++){
        //     priority_list[i]->next = nullptr;
        // }
        free_list();
        priority_list.clear();
        // print_list();
        return part;
    }
    void GreedPartition(){
        _partition.clear();
        select_free = 0;
        for(unsigned i=0; i<_partition_number; i++){
            auto t = GreedInsert(i);
            _partition.push_back(t);
            // std::cout << "pid: "<<i << std::endl;
            // for(auto &n: t){
            //     std::cout << n;
            //     std::cout << " ";
            // }
            // std::cout << std::endl;
            if(i % 10000 == 0){
                
                std::cout <<(float) i / _partition_number  << std::endl;
            }
        }
        std::cout << "select free per partition "<<(double)select_free / _partition_number << std::endl;
    }
    

    size_t _dim; // vector dimension
    _u64 nd;    // vector number
    _u64 max_node_len;
    unsigned _width{};                               // max out-degree
    unsigned _ep{};                                  // seed vertex id
    std::vector<std::vector<unsigned>> direct_graph; // neighbor list
    unsigned select_free;
    unsigned _k1 = 10; // cluster number for 0 level clustering (file number)
    std::vector<char *> _cluster_data;
    std::vector<unsigned> _cluster_size;
    unsigned _total_cluster_size = 0;
    _u64 C;                                                 // partition size threshold
    _u64 _partition_number = 0;                             // the number of partitions
    std::vector<std::vector<unsigned>> _partition{1000000}; // each partition set
    std::vector<std::mutex *> pmutex;
    std::mutex vlock;
    std::mutex flock;
    int cur = 0;
    std::vector<std::vector<unsigned>> reverse_graph;
    std::vector<std::vector<unsigned>> undirect_graph;
    std::unordered_map<unsigned, unsigned> id2pid;
    std::unordered_map<unsigned, unsigned> id2ratio;
    std::unordered_map<unsigned, bool> unfilled_partition;
    std::unordered_map<unsigned, bool> visited;
    std::mutex plock;
    unsigned en;
    ppq unfill;
    int T = 2;
    int round = 0;
    unsigned swap_cnt = 0;
    unsigned swap_min = 0;
    double ivf_time = 0.0;
    bool debug = false;
    double de_time = 0.0;
    unsigned cursize = 10000;
    uint64_t select_nums = 0;
    uint64_t getUnfilled_nums = 0;
   _u64 E;
   std::uniform_real_distribution<> *_dis;
   std::mt19937 *_gen;
   std::random_device *_rd;
   std::shared_mutex smutex;
   std::vector<ListNode*> priority_list;
   std::unordered_set<unsigned> free_node;
   concurrent_queue free_q;
   std::unordered_set<unsigned> fix_node;
   std::unordered_map<unsigned, ListNode*> id2ptr;
};
#endif // SSD_BASED_PLAN_INDEX_H
