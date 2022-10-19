# SSD_Based_Plan

# clone

```
git clone --recursive git@github.com:SonglinLife/SSD_BASED_PLAN.git 
```

## compile

``` sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target partitioner -j8
```

## useage

`gp.sh`
``` sh
bigann_10K="/data/wsl/SSD_Based_Plan/bigann/BIGANN_10K_R48_L128_B0.0003_disk.index"
bigann_1M="/data/wsl/SSD_Based_Plan/bigann/BIGANN1M_R48_L128_B0.03_disk.index"
bigann_10M="/data/wsl/SSD_Based_Plan/bigann/BIGANN10M_R48_L128_B0.03_M32_disk.index"
bigann_10M_R48="/data/wsl/SSD_Based_Plan/bigann/BIGANN_10M_R48_L128_B0.3_disk.index"
bigann_10M_R32="/data/wsl/SSD_Based_Plan/bigann/BIGANN_10M_R32_L96_B0.3_disk.index"
deep_10M="/data/wsl/diskann_re/build/tests/index/DEEP_10M_R48_L128_B0.225_disk.index_cp"
ssnpp_10M="/data/wsl/SSD_Based_Plan/bigann/SSNPP_ORG_R48_L128_B0.6_M32_disk.index"
sift_10M_R32="/data/wsl/SSD_Based_Plan/bigann/SIFT_10M_R32_L96_B0.3_disk.index"

index=$deep_10M
dim=96
data_type=float
nd=10000000
mkdir -p ./gpfile
gp_file=./gpfile/deep_10m_r48_rb1
lock_nums=1000000
./build/partitioner --data_type $data_type --index_file $index --gp_file $gp_file -N $nd -D $dim -L 16 --lock_nums $lock_nums
#>> ${gp_file}.log
```

## args


```
./partitioner --h

-h [ --help ]                print information
--data_type arg              data type <int8/uint8/float>
--index_file arg             diskann diskann index or mem index
--gp_file arg                output gp file
--freq_file arg              freq_file
-N [ --npts ] arg            data size
-D [ --dim ] arg             data vector dim
--lock_nums arg (=0)         lock node nums
-B [ --block_size ] arg (=1) block size for one partition
-L [ --ldg_times ] arg (=4)  exec ldg partition alg times
-A [ --alg ] arg (=0)        algorithm for graph Partition
--use_disk arg (=1)          Use 1 for use disk index (default), 0 for 
                             DiskANN mem index
--visual arg (=0)            see real time progress of
```
