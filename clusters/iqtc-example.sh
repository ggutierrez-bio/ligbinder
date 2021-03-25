#!/bin/bash
#$ -N ligbinder  
#$ -pe gpu 1
#$ -q fartorgpu.q
#$ -S /bin/bash
#$ -cwd
#$ -o ligbinder.out
#$ -e ligbinder.err
#$ -m e

# load environment variables and modules
. /etc/profile
module use ~/modules
module load ligbinder
export CUDA_VISIBLE_DEVICES=`cat $TMPDIR/.gpus`


persistent_folder=$PWD

cp -r * $TMPDIR/.
cd $TMPDIR
ligbinder -l ligbinder.log
cp -r * $persistent_folder/.
exit
