#!/bin/bash
#PBS -q q-s
###PBS -l select=1:ncpus=36::mpiprocs=36:mem=360gb
#PBS -l select=1:ncpus=36:mpiprocs=36:mem=360gb
#PBS -l walltime=48:00:00
#PBS -N dSFMTSearch
###PBS -o search_test.out.txt 不要か
#PBS -j oe
#PBS -m abe
#PBS -M sai10@hiroshima-u.ac.jp

## 自動並列化でこれを指定せよとの記述あり gcc ならなし
## export OMP_NUM_THREADS=6

### mpiexec ではいらない　dplace では必須
### export KMP_AFFINITY=disabled

cd $PBS_O_WORKDIR # 現ディレクトリに変更 ほぼ必須

### openMP/MPI ハイブリッドでは omplace を指定せよとのこと。（たぶん自動並列化も）
### gcc で openMP を使わないのでなし
mpiexec_mpt ./dSFMTdc --file dSFMTdc.19937 --fixed-sl1 19 --fixed-pos1 117 --seed 100 --count 100 19937
