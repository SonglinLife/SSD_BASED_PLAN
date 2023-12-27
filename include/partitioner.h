//
// Created by Zilliz_wmz on 2022/6/6.
//
#pragma once

#include <omp.h>
#include <oneapi/tbb/concurrent_queue.h>
#include <oneapi/tbb/concurrent_set.h>
#include <algorithm>
#include <atomic>
#include <bitset>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <ostream>
#include <queue>
#include <random>
#include <set>
#include <shared_mutex>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include "filesystem"
#include "freq_relayout.h"

#ifndef INF
#define INF 0xffffffff
#endif  // INF
#ifndef READ_U64
#define READ_U64(stream, val) stream.read((char *)&val, sizeof(_u64))
#endif  // !READ_U64
#ifndef ROUND_UP
#define ROUND_UP(X, Y) (((uint64_t)(X) / (Y)) + ((uint64_t)(X) % (Y) != 0)) * (Y)
#endif  // !ROUND_UP
#ifndef SECTOR_LEN
#define SECTOR_LEN (_u64)4096
#endif  // !SECTOR_LEN

namespace GP {

using concurrent_queue = oneapi::tbb::concurrent_bounded_queue<unsigned>;
namespace fs = std::filesystem;
using _u64 = unsigned long int;
using _u32 = unsigned int;
using VecT = uint8_t;

inline size_t get_file_size(const std::string &fname) {
  std::ifstream reader(fname, std::ios::binary | std::ios::ate);
  if (!reader.fail() && reader.is_open()) {
    size_t end_pos = reader.tellg();
    reader.close();
    return end_pos;
  } else {
    std::cout << "Could not open file: " << fname << std::endl;
    return 0;
  }
}

// DiskANN changes how the meta is stored in the first sector of
// the _disk.index file after commit id 8bb74ff637cb2a77c99b71368ade68c62b7ca8e0
// (exclusive) It returns <is_new_version, vector of metas uint64_ts>
std::pair<bool, std::vector<_u64>> get_disk_index_meta(const std::string &path) {
  std::ifstream fin(path, std::ios::binary);

  int meta_n, meta_dim;
  const int expected_new_meta_n = 9;
  const int expected_new_meta_n_with_reorder_data = 12;
  const int old_meta_n = 11;
  bool is_new_version = true;
  std::vector<_u64> metas;

  fin.read((char *)(&meta_n), sizeof(int));
  fin.read((char *)(&meta_dim), sizeof(int));

  if (meta_n == expected_new_meta_n || meta_n == expected_new_meta_n_with_reorder_data) {
    metas.resize(meta_n);
    fin.read((char *)(metas.data()), sizeof(_u64) * meta_n);
  } else {
    is_new_version = false;
    metas.resize(old_meta_n);
    fin.seekg(0, std::ios::beg);
    fin.read((char *)(metas.data()), sizeof(_u64) * old_meta_n);
  }
  fin.close();
  return {is_new_version, metas};
}

class graph_partitioner {
 public:
  graph_partitioner(const char *indexName, const char *data_type = "uint8", unsigned fixed_ratio = 0,
                    bool load_disk = true, unsigned BS = 1, bool visual = false,
                    std::string freq_file = std::string(""), unsigned cut = INF) {
    _visual = visual;
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    // check file size
    size_t actual_size = get_file_size(indexName);
    size_t expected_size;
    auto meta_pair = get_disk_index_meta(indexName);
    if (meta_pair.first) {
      expected_size = meta_pair.second.back();
    } else {
      expected_size = meta_pair.second.front();
    }
    if (actual_size != expected_size) {
      std::cout << "index file not match!" << std::endl;
      exit(-1);
    }

    _rd = new std::random_device();
    _gen = new std::mt19937((*_rd)());
    _dis = new std::uniform_real_distribution<>(0, 1);
    if (load_disk) {
      if (std::string(data_type) == std::string("uint8")) {
        load_disk_index<uint8_t>(indexName, BS);
      } else if (std::string(data_type) == std::string("float")) {
        load_disk_index<float>(indexName, BS);
      } else {
        std::cout << "not support type" << std::endl;
        exit(-1);
      }
    } else {
      load_vamana(indexName);
    }
    cursize = _nd / 1000;

    if (!freq_file.empty()) {
      if (!fs::exists(freq_file)) {
        std::cout << "No such freq file!" << std::endl;
        exit(-1);
      }
      read_freq(_freq_list, _freq_nei_list, freq_file);
      if (_freq_list.size() != _nd) {
        std::cout << "freq not match, freq file has node " << _freq_list.size() << " but index file nodes is " << _nd
                  << std::endl;
        exit(-1);
      }
      relayout_adj(_freq_nei_list, full_graph);
    }
    // copy to direct_graph
    direct_graph.clear();
    direct_graph.resize(full_graph.size());
#pragma omp parallel for
    for (unsigned i = 0; i < _nd; i++) {
      direct_graph[i].assign(full_graph[i].begin(), full_graph[i].end());
    }
    // cut graph
    if (cut != INF) {
      std::cout << "direct graph will be cut, it degree become " << cut << std::endl;
    }
#pragma omp parallel for
    for (unsigned i = 0; i < _nd; i++) {
      if (cut < direct_graph[i].size()) {
        direct_graph[i].resize(cut);
      }
    }
// cout graph E num 
    size_t _ep = 0;
#pragma omp parallel for reduction(+ : _ep) 
    for (unsigned i = 0; i < _nd; i++) {
      _ep += direct_graph[i].size();
    }
    std::cout << "graph vertex size: " << _nd << " graph edge size: " << _ep << std::endl;
    // reverse graph
    std::vector<std::mutex> ms(_nd);
    reverse_graph.resize(_nd);
#pragma omp parallel for shared(reverse_graph, direct_graph)
    for (unsigned i = 0; i < _nd; i++) {
      for (unsigned j = 0; j < direct_graph[i].size(); j++) {
        std::lock_guard<std::mutex> lock(ms[direct_graph[i][j]]);
        reverse_graph[direct_graph[i][j]].emplace_back(i);
      }
    }
    std::cout << "reverse graph done." << std::endl;
    for (unsigned i = 0; i < _partition_number; i++) {
      pmutex.push_back(std::make_unique<std::mutex>());
    }
    if (fixed_ratio != 0) {
      // make a random vector from 0 to np, the vector size is _nd * fixed_ratio / 100
      _fixed_sz = _nd * fixed_ratio / 100;
    }
  }

  void cout_step() {
    if (!_visual) {
      return;
    }

#pragma omp atomic
    cur++;
    if ((cur + 0) % cursize == 0) {
      std::cout << (double)(cur + 0) / _nd * 100 << "%    \r";
      std::cout.flush();
    }
  }
  /**
   * load vamana graph index from disk
   * @param filename
   */
  void load_vamana(const char *filename, bool sample = false) {
    std::cout << "Reading index file: " << filename << "... " << std::flush;
    std::ifstream in;
    in.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
      in.open(filename, std::ios::binary);
      size_t expected_file_size;
      in.read((char *)&expected_file_size, sizeof(uint64_t));
      in.read((char *)&_width, sizeof(unsigned));
      in.read((char *)&_ep, sizeof(unsigned));
      std::cout << "Loading vamana index " << filename << "..." << std::flush;

      size_t cc = 0;
      unsigned nodes = 0;
      while (in.peek() != EOF) {
        unsigned k;
        in.read((char *)&k, sizeof(unsigned));
        cc += k;
        ++nodes;
        std::vector<unsigned> tmp(k);
        in.read((char *)tmp.data(), k * sizeof(unsigned));
        direct_graph.emplace_back(tmp);
        if (nodes % 10000000 == 0) std::cout << "." << std::flush;
      }
      if (direct_graph.size() != _nd) {
        std::cout << "graph vertex size error!\n";
        exit(-1);
      }
      if (sample) {
        std::cout << "cut adj" << std::endl;
        for (unsigned i = 0; i < _nd; i++) {
          std::vector<unsigned> tmp;
          tmp.reserve(10);
          for (unsigned j = 0; j < 20 && j < direct_graph[i].size(); j++) {
            tmp.push_back(direct_graph[i][j]);
          }
          direct_graph[i].clear();
          direct_graph[i].assign(tmp.begin(), tmp.end());
        }
      }
      C = 12;
      _partition_number = ROUND_UP(_nd, C) / C;
      reverse_graph.resize(_nd);
      std::vector<std::mutex> ms(_nd);
#pragma omp parallel for shared(reverse_graph, direct_graph)
      for (unsigned i = 0; i < _nd; i++) {
        for (unsigned j = 0; j < direct_graph[i].size(); j++) {
          std::lock_guard<std::mutex> lock(ms[direct_graph[i][j]]);
          reverse_graph[direct_graph[i][j]].emplace_back(i);
        }
      }
      std::cout << "done. Index has " << nodes << " nodes and " << cc << " out-edges" << std::endl;
      for (unsigned i = 0; i < _partition_number; i++) {
        pmutex.push_back(std::make_unique<std::mutex>());
      }
    } catch (std::system_error &e) {
      exit(-1);
    }
  }

  template <typename T>
  void load_disk_index(const char *index_name, int BS = 1) {
    std::cout << "loading disk index file: " << index_name << "... " << std::flush;
    std::ifstream in;
    in.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    try {
      _u64 expected_npts;
      auto meta_pair = get_disk_index_meta(index_name);

      if (meta_pair.first) {
        // new version
        expected_npts = meta_pair.second.front();
      } else {
        expected_npts = meta_pair.second[1];
      }
      _nd = expected_npts;
      _dim = meta_pair.second[1];

      _max_node_len = meta_pair.second[3];
      C = meta_pair.second[4];

      _partition_number = ROUND_UP(_nd, C) / C;

      std::unique_ptr<char[]> mem_index = std::make_unique<char[]>(_partition_number * SECTOR_LEN);
      in.open(index_name, std::ios::binary);
      in.seekg(SECTOR_LEN, std::ios::beg);
      in.read(mem_index.get(), _partition_number * SECTOR_LEN);
      in.close();
      full_graph.resize(_nd);
      _u64 des = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : des)
      for (unsigned i = 0; i < _partition_number; i++) {
        std::unique_ptr<char[]> sector_buf = std::make_unique<char[]>(SECTOR_LEN);
        memcpy(sector_buf.get(), mem_index.get() + i * SECTOR_LEN, SECTOR_LEN);
        for (unsigned j = 0; j < C && i * C + j < _nd; j++) {
          std::unique_ptr<char[]> node_buf = std::make_unique<char[]>(_max_node_len);
          memcpy(node_buf.get(), sector_buf.get() + j * _max_node_len, _max_node_len);
          unsigned &nnbr = *(unsigned *)(node_buf.get() + _dim * sizeof(T));
          unsigned *nhood_buf = (unsigned *)(node_buf.get() + (_dim * sizeof(T)) + sizeof(unsigned));
          std::vector<unsigned> tmp(nnbr);
          des += nnbr;
          memcpy((char *)tmp.data(), nhood_buf, nnbr * sizeof(unsigned));
          full_graph[i * C + j].assign(tmp.begin(), tmp.end());
        }
      }
      std::cout << "avg degree: " << (double)des / _nd << std::endl;
      mem_index.reset();
      C = (SECTOR_LEN * BS) / _max_node_len;
      _partition_number = ROUND_UP(_nd, C) / C;
      std::cout << "_nd: " << _nd << " _dim:" << _dim << " C:" << C << " pn:" << _partition_number << std::endl;
      std::cout << "load index over." << std::endl;
    } catch (std::system_error &e) {
      std::cout << "open file " << index_name << " error!" << std::endl;
      exit(-1);
    }
  }

  /**
   * save the partition result
   * @tparam T
   * @param filename
   * @param partition
   */
  void save_partition(const char *filename) {
    // 正确性测试
    // re_id2pid();
    std::ofstream writer(filename, std::ios::binary | std::ios::out);
    std::cout << "writing bin: " << filename << std::endl;
    writer.write((char *)&C, sizeof(_u64));
    writer.write((char *)&_partition_number, sizeof(_u64));
    writer.write((char *)&_nd, sizeof(_u64));
    std::cout << "_partition_num: " << _partition_number << " C: " << C << " _nd: " << _nd << std::endl;
  
    for (unsigned i = 0; i < _partition_number; i++) {
      auto& p = _partition[i];
      unsigned s = _sizes[i].load();
    
      writer.write((char *)&s, sizeof(unsigned));
      writer.write((char *)p.data(), sizeof(unsigned) * s);
    }
    std::vector<unsigned> id2pidv(_nd);
    for(size_t i = 0; i < _nd; i++) {
      id2pidv[i] = _labels[i];
    }
    writer.write((char *)id2pidv.data(), sizeof(unsigned) * _nd);
  }

  /**
   * load partition from disk
   * @param filename
   */
  void load_partition(const char *filename) {
    std::ifstream reader(filename, std::ios::binary);
    reader.read((char *)&C, sizeof(_u64));
    reader.read((char *)&_partition_number, sizeof(_u64));
    reader.read((char *)&_nd, sizeof(_u64));
    std::cout << "load partition _partition_num: " << _partition_number << ", C: " << C << std::endl;
    _partition.resize(_partition_number);
    _sizes = (std::atomic<unsigned> *)malloc(sizeof(std::atomic<unsigned>) * _partition_number);
    new (_sizes) std::atomic<unsigned>[ _partition_number ];
    for (unsigned i = 0; i < _partition_number; i++) {
      unsigned c;
      reader.read((char *)&c, sizeof(unsigned));
      _partition[i].resize(c);
      _sizes[i].store(c);
      reader.read((char *)_partition[i].data(), sizeof(unsigned) * c);
    }
    _labels.resize(_nd);
    reader.read((char *)_labels.data(), sizeof(unsigned) * _nd);
  }
  void re_id2pid() {
    id2pid.clear();
    for (unsigned i = 0; i < _partition_number; i++) {
      for (unsigned j = 0; j < _partition[i].size(); j++) {
        id2pid[_partition[i][j]] = i;
      }
    }
  }
  /**
   * count the id overlap according to the graph partitioning
   */
  void partition_statistic() {
    std::vector<unsigned> overlap(_nd, 0);
    std::vector<unsigned> blk_neighbor_overlap(_partition_number, 0);
    double overlap_ratio = 0;

#pragma omp parallel for schedule(dynamic, 100) reduction(+ : overlap_ratio)
    for (size_t i = 0; i < _partition_number; i++) {
      std::unordered_set<unsigned> neighbors;
      unsigned blk_neighbor_num = 0;
      for (size_t j = 0; j < _partition[i].size(); j++) {
        blk_neighbor_num += full_graph[_partition[i][j]].size();
        std::unordered_set<unsigned> ne;
        for (unsigned &x : full_graph[_partition[i][j]]) {
          neighbors.insert(x);
          ne.insert(x);
        }
        blk_neighbor_overlap[i] = blk_neighbor_num - neighbors.size();
        for (size_t z = 0; z < _partition[i].size(); z++) {
          if (_partition[i][j] == _partition[i][z]) continue;
          if (ne.find(_partition[i][z]) != ne.end()) {
            overlap[_partition[i][j]]++;
          }
        }
        overlap_ratio +=
            (_partition[i].size() == 1 ? 0 : (1.0 * overlap[_partition[i][j]] / (_partition[i].size() - 1)));
      }
    }
    unsigned max_overlaps = 0;
    unsigned min_overlaps = std::numeric_limits<unsigned>::max();
    double ave_overlap_ratio = 0;
    std::map<unsigned, unsigned> overlap_count;
    for (size_t i = 0; i < _nd; i++) {
      if (overlap_count.count(overlap[i])) {
        overlap_count[overlap[i]]++;
      } else {
        overlap_count[overlap[i]] = 1;
      }
      if (overlap[i] > max_overlaps) max_overlaps = overlap[i];
      if (overlap[i] < min_overlaps) min_overlaps = overlap[i];
    }
    ave_overlap_ratio = overlap_ratio / (double)_nd;
    for (auto &it : overlap_count) {
      std::cout << "each id, overlap number " << it.first << ", count: " << it.second << std::endl;
    }
    std::cout << "each id, max overlaps: " << max_overlaps << std::endl;
    std::cout << "each id, min overlaps: " << min_overlaps << std::endl;
    std::cout << "each id, average overlap ratio: " << ave_overlap_ratio << std::endl;
  }

  void partition_statistic_v2() {
    std::vector<unsigned> overlap(_nd, 0);
    std::vector<unsigned> blk_neighbor_overlap(_partition_number, 0);
    double overlap_ratio = 0;
    std::map<unsigned, unsigned> neighbors_set;
    /* for (size_t i = 0; i < _partition_number; i++) { */
    /*   for (size_t j = 0; j < _sizes[i].load(); j++) { */
    /*     if (neighbors_set.find(_partition[i][j]) != neighbors_set.end()) { */
    /*       std::cout << "error!" << std::endl; */
    /*       std::cout << "id: " << _partition[i][j] << " partition: " << i << std::endl; */
    /*       std::cout << "set has this id in blk " << neighbors_set[_partition[i][j]] << std::endl; */
    /*       exit(-1); */
    /*     } */
    /*     neighbors_set.insert({_partition[i][j], i}); */
    /*   } */
    /* } */

#pragma omp parallel for schedule(dynamic, 100) reduction(+ : overlap_ratio)
    for (size_t i = 0; i < _partition_number; i++) {
      std::unordered_set<unsigned> neighbors;
      unsigned blk_neighbor_num = 0;
      if (_sizes[i].load() > C) {
        std::cout << "partition " << i << " size: " << _sizes[i].load() << std::endl;
        std::cout << "error!" << std::endl;
        exit(-1);
      }
      for (size_t j = 0; j < _sizes[i].load(); j++) {
        auto id = _partition[i][j];
        blk_neighbor_num += full_graph[_partition[i][j]].size();
        std::unordered_set<unsigned> ne;
        for (unsigned &x : full_graph[_partition[i][j]]) {
          neighbors.insert(x);
          ne.insert(x);
        }
        blk_neighbor_overlap[i] = blk_neighbor_num - neighbors.size();
        for (size_t z = 0; z < _sizes[i].load(); z++) {
          if (_partition[i][j] == _partition[i][z]) continue;
          if (ne.find(_partition[i][z]) != ne.end()) {
            overlap[_partition[i][j]]++;
          }
        }
        overlap_ratio += (_partition[i].size() == 1 ? 0 : (1.0 * overlap[_partition[i][j]] / (_sizes[i].load() - 1)));
      }
    }
    unsigned max_overlaps = 0;
    unsigned min_overlaps = std::numeric_limits<unsigned>::max();
    double ave_overlap_ratio = 0;
    std::map<unsigned, unsigned> overlap_count;
    for (size_t i = 0; i < _nd; i++) {
      if (overlap_count.count(overlap[i])) {
        overlap_count[overlap[i]]++;
      } else {
        overlap_count[overlap[i]] = 1;
      }
      if (overlap[i] > max_overlaps) max_overlaps = overlap[i];
      if (overlap[i] < min_overlaps) min_overlaps = overlap[i];
    }
    ave_overlap_ratio = overlap_ratio / (double)_nd;
    for (auto &it : overlap_count) {
      std::cout << "each id, overlap number " << it.first << ", count: " << it.second << std::endl;
    }
    std::cout << "each id, max overlaps: " << max_overlaps << std::endl;
    std::cout << "each id, min overlaps: " << min_overlaps << std::endl;
    std::cout << "each id, average overlap ratio: " << ave_overlap_ratio << std::endl;
  }
  unsigned select_partition(unsigned i) {
#pragma omp atomic
    select_nums++;

    float maxn = 0.0;
    unsigned res = INF;
    std::unordered_map<unsigned, unsigned> pcount;
    unsigned tpid = 0;
    for (auto n : direct_graph[i]) {
      unsigned pid = id2pid[n];
      if (pid == INF) continue;
      pcount[pid] = pcount[pid] + 1;
      if (tpid < pid) {
        tpid = pid;
      }
    }
    for (auto n : reverse_graph[i]) {
      unsigned pid = id2pid[n];
      if (pid == INF) continue;
      pcount[pid] = pcount[pid] + 1;
      if (tpid < pid) {
        tpid = pid;
      }
    }
    for (auto c : pcount) {
      unsigned pid = c.first;
      float cnt = c.second;
      std::lock_guard<std::mutex> lock(*pmutex[pid]);
      double s = _partition[pid].size();
      cnt *= (1 - s / C);
      if (cnt > maxn && _partition[pid].size() < C) {
        res = pid;
        maxn = cnt;
      }
    }
    pcount.clear();
    if (res == INF) {
#pragma omp atomic
      select_free++;
      res = getUnfilled();
    }
    return res;
  }

  unsigned getUnfilled() {
#pragma omp atomic
    getUnfilled_nums++;
    unsigned res;
    do {
      free_q.pop(res);
    } while (_partition[res].size() == C);
    return res;
  }

  unsigned getUnfilled_v2() {
#pragma omp atomic
    getUnfilled_nums++;
    unsigned res;
    do {
      free_q.pop(res);
    } while (_sizes[res].load() == C);
    return res;
  }
  // graph partition
  void graph_partition(const char *filename, int k, int lock_nums = 0) {
    _labels.resize(_nd, INF);

    // uninform random unsinged generator from 0 to _nd -1
    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 g(seed);
    // distribute random gener can produce a random number in range [0, _nd -1]
    std::uniform_int_distribution<> dis(0, _partition_number - 1);

    _partition = std::vector<std::vector<unsigned>>(_partition_number, std::vector<unsigned>(C, INF));
    /* #pragma omp parallel for schedule(guided) */
    /*     for (unsigned i = 0; i < _nd; i++) { */
    /* auto random_label = dis(g); */
    /* _labels[i] = random_label; */
    /* } */
    // malloc _sizes
    _sizes = (std::atomic<unsigned> *)malloc(sizeof(std::atomic<unsigned>) * _partition_number);
    // place new
    new (_sizes) std::atomic<unsigned>[ _partition_number ];
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::uniform_int_distribution<> dis1(0, _nd - 1);
    std::uniform_int_distribution<> dis2(0, _partition_number - 1);
    _fixed.resize(_fixed_sz);

    for (size_t i = 0; i < _fixed_sz; i++) {
      auto id = i;
      auto pid = dis2(g);
      while (_sizes[pid].load() == C) {
        pid = dis2(g);
      }
      _fixed[id] = pid;
      _sizes[pid].fetch_add(1);
    }
    std::cout << "init over." << std::endl;

    for (int i = 0; i < k; i++) {
      select_free = 0;
      graph_partition_LDG();
      std::cout << "select free: " << (double)select_free / _partition_number << std::endl;
      /* partition_statistic_v2(); */
      /* auto ivf_file_name = std::string(filename) + std::string(".ivf") + std::to_string(i + 1); */
      std::cout << "total ivf time: " << ivf_time << std::endl;
      /* save_partition(ivf_file_name.c_str()); */
    }
    save_partition(filename);
    std::cout << "select pid nums" << select_nums << " get unfilled partition nums: " << getUnfilled_nums << std::endl;
    std::cout << "total ivf time: " << ivf_time << std::endl;
  }
  void graph_partition_LDG() {
    free_q.clear();

#pragma omp parallel for schedule(guided)
    for (unsigned i = 0; i < _partition_number; i++) {
      _sizes[i].store(0);
      free_q.push(i);
    }
#pragma omp parallel for schedule(guided)
    for (unsigned i = 0; i < _fixed_sz; i++) {
      unsigned id = i;
      unsigned pid = _fixed[id];
      _partition[pid][_sizes[pid].fetch_add(1)] = id;
      _labels[id] = pid;
    }
    cur = 0;
    std::cout << "start" << std::endl;
    /* std::vector<unsigned> stream(_nd); */
    /* std::iota(stream.begin(), stream.end(), 0); */
    /* auto rng = std::default_random_engine{}; */
    /* std::shuffle(std::begin(stream), std::end(stream), rng); */
    auto start = omp_get_wtime();
#pragma omp parallel for schedule(guided)
    for (unsigned i = _fixed_sz; i < _nd; i++) {
      size_t n = i;
      while (sync_v2(n) == false) {
        continue;
      }
      cout_step();
    }
    auto end = omp_get_wtime();
    std::cout << "ivf time: " << end - start << " round: " << round << std::endl;
    ivf_time += end - start;
    round++;
  }
  unsigned sync(unsigned i) {
    unsigned pid = select_partition(i);
    pmutex[pid]->lock();

    while (_partition[pid].size() == C) {
      pmutex[pid]->unlock();
      pid = select_partition(i);
      pmutex[pid]->lock();
    }
    _partition[pid].emplace_back(i);
    id2pid[i] = pid;
    unsigned s = _partition[pid].size();
    pmutex[pid]->unlock();

    if (s != C) {
      free_q.push(pid);
    }

    return pid;
  }

  unsigned select_partition_v2(unsigned u) {
    unsigned res = INF;
    std::unordered_map<unsigned, unsigned> pcount;
    for (auto n : direct_graph[u]) {
      unsigned label = _labels[n];
      if (label == INF) continue;
      pcount[label] = pcount[label] + 1;
    }
    for (auto n : reverse_graph[u]) {
      unsigned label = _labels[n];
      if (label == INF) continue;
      pcount[label] = pcount[label] + 1;
    }
    float max_gain = 0.0;
    for (auto c : pcount) {
      unsigned label = c.first;
      float cnt = c.second;
      unsigned sz = _sizes[label].load();
      double gain = (double)cnt * ((double)1 - (double)sz / C);
      if (gain > max_gain && sz < C) {
        res = label;
        max_gain = gain;
      }
    }
    return res;
  }
  bool sync_v2(unsigned u) {
    unsigned select_label = select_partition_v2(u);
    unsigned label = select_label == INF ? getUnfilled_v2() : select_label;
    auto old_sz = _sizes[label].load(std::memory_order_relaxed);
    while (old_sz != C && !_sizes[label].compare_exchange_weak(old_sz, old_sz + 1, std::memory_order_seq_cst)) {
      continue;
    }
    if (old_sz == C) {
      return false;
    }
    _labels[u] = label;
    _partition[label][old_sz] = u;
    if (old_sz + 1 != C and select_label == INF) {
      free_q.push(label);
    }
    return true;
  }

 private:
  size_t _dim;  // vector dimension
  _u64 _nd;     // vector number
  _u64 _max_node_len;
  unsigned _width;                                  // max out-degree
  unsigned _ep;                                     // seed vertex id
  std::vector<std::vector<unsigned>> direct_graph;  // neighbor list
  std::vector<std::vector<unsigned>> full_graph;
  unsigned select_free;
  _u64 C;                                                  // partition size threshold
  _u64 _partition_number = 0;                              // the number of partitions
  std::vector<std::vector<unsigned>> _partition{1000000};  // each partition set
  std::vector<std::unique_ptr<std::mutex>> pmutex;
  int cur = 0;
  std::vector<std::vector<unsigned>> reverse_graph;
  std::vector<std::vector<unsigned>> undirect_graph;
  std::unordered_map<unsigned, unsigned> id2pid;
  std::unordered_map<unsigned, unsigned> id2ratio;
  int round = 0;
  double ivf_time = 0.0;
  bool _visual = false;
  unsigned cursize = 10000;
  uint64_t select_nums = 0;
  uint64_t getUnfilled_nums = 0;
  _u64 E;
  std::uniform_real_distribution<> *_dis;
  std::mt19937 *_gen;
  std::random_device *_rd;
  concurrent_queue free_q;

  std::vector<puu> _freq_list;
  vpu _freq_nei_list;
  std::vector<bool> _lock_nodes;
  std::vector<bool> _lock_pids;

  std::atomic<unsigned> *_sizes;
  std::vector<unsigned> _labels;
  std::vector<unsigned> _fixed;
  unsigned _fixed_sz = 0;
};
}  // namespace GP
