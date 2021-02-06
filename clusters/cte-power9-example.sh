#!/bin/bash
#!/bin/bash
#SBATCH --job-name=ligbinder
#SBATCH -D .
#SBATCH --output=ligbinder.o
#SBATCH --error=ligbinder.e
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=40

export OMP_NUM_THREADS=1

module load gcc/7.3.0 cuda/9.1 openmpi/3.0.0 boost/1.69.0 pnetcdf/1.11.2 amber/20 python/3.6.5

# remember to activate your virtual environment if the ligbinder package is not installed locally
# this is an example, modify the path for your own environment
# source ~/environments/ligbinder/bin/activate

ligbinder
