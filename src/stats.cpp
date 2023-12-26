#include "partitioner.h"

int main (int argc, char *argv[])
{
  std::string index_file, data_type, gp_file;

  if (argc != 4) {
    std::cout << "Usage: " << argv[0] << " <index_file> <data_type> <gp_file>" << std::endl;
    return -1;
  }
  index_file = argv[1];
  data_type = argv[2];
  gp_file = argv[3];

  auto stats = GP::graph_partitioner(index_file.c_str(), data_type.c_str());
  stats.load_partition(gp_file.c_str());
  stats.partition_statistic_v2();
  return 0;
}
