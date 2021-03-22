#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH --nodes=1

module load amber/20
module load ligbinder

ligbinder
