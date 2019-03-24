#!/bin/bash
#PBS -q q-s
#PBS -l select=1:ncpus=36:mpiprocs=36:mem=360gb
#PBS -N dSFMTSearch
#PBS -o out_file
#PBS -j oe
#PBS -m abe
#PBS -M sai10@hiroshima-u.ac.jp

export KMP_AFFINITY=disabled

cd $PBS_O_WORKDIR # 現ディレクトリに変更 ほぼ必須

mipexec_mpt ./dSFMTdc
