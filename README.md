# SSD_Based_Plan

# clone

```
git clone --recursive git@github.com:SonglinLife/SSD_BASED_PLAN.git 
```

## compile

``` sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

## use

``` sh
./SSD_Based_Plan dim nd index data_type gp_file thread_num LDG_times 0 DEBUG blocksize 
```

`dim` vector dim, like 128

`nd` total vector number, like 1000000

`index` diskann index, like ./index/BIGANN_R48_disk.index

`data_type` vector type , float or uint8. int8 is equal to uint8, just use uint8

`gp_file` the partition file name

`thread_num` thread number

`LDG_times` the algorithm will execute , set 0 will output unpartitioned result

`0` set it 0

`DEBUG` see DEBUG message, set it 0 for no debug

`blocksize` set 1, the parition block size would be 4KB, set 2 would be 8KB, and so on

use like 

```
./SSD_Based_Plan 96 10000000 /data/wsl/diskann_re/build/tests/index/DEEP_10M_R48_L128_B0.225_disk.index_cp float ./gpfile/deep_10m_r48_rb1 64 16 0 1 1
```
