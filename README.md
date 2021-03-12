# LigBinder 0.1a

## installation:

dependencies automatic installation is disabled as BSC rejects outgoing traffic originated from within  

make sure your environment has a working installation of:

python packages:
* pyyaml
* pytraj
* parmed

binaries:
* pmemd.cuda
* sander (optional)

in order to be able to use ambertools20 python bindings, it is required to use python version 3.6.5  
make sure you loaded the necessary modules before installing this package.  
this must be done by running this:

```
module load gcc/7.3.0
module load boost/1.69.0
module load pnetcdf/1.11.2
module load amber/20
module load python/3.6.5
```

create a python virtual environment in a convenient folder by:

```
virtualenv <path> --prompt <prompt>
```

for example:

```
virtualenv .venv --prompt "(ligbinder) "
```
activate this environment by running
```
source .venv/bin/activate
```

with the virtual environment active, copy the repository to a folder of your choice, get in and run:  
```
pip install .
```

if everything went fine you should be able to run ligbinder by executing:

```
ligbinder
```

you can also install it into your user general environment instead by (not creating and activating the virtual environment and) running

```
pip install --user .
```

this way you won't need to source any environment in the submission script.


## quickstart guide

1) create a new folder for the ligbinder dynamics and get in:

2) create a folder named `data`. it will contain the relevant files

3) copy your equilibrated pose into `data/pose.rst7`. it must be in rst7 format (amber binary restart file)

4) copy the corresponding topology file into `data/top.prmtop`. it must be in amber parm file format

5) copy the target pose into `data/ref.crd`. (amber binary format again)

6) run `ligbinder` and go for pop-corn

```
mkdir lb_test && cd lb_test
mkdir data
cp ../my_dynamic/eq.rst7 data/pose.rst7
cp ../my_dynamic/my_topology.prmtop data/top.prmtop
cp ../my_dynamic/bound_state.crd data/ref.crd
ligbinder
```
## detailed guide and configuration

configuration can be specified by providing a yaml file.

the currently available parameters are:
```
tree:
    max_depth: 20                   # maximum number of successive iterations
    max_children: 5                 # maximum number of dynamics spawn from a single intermediate iteration
    tolerance: 0.2                  # rmsd tolerance to accept a node as a hit
    max_nodes: 500                  # maximum number of nodes in total
    min_relative_improvement: 0.1   # minimum rmsd improvement in relation to current rmsd. default is 10% (0.1)
    min_absolute_improvement: 0.1   # minimum rmsd improvement in absolute value
    max_steps_back: 5               # maximum difference between maximum depth reached and current depth

data_files:
    crd_file: "./data/pose.rst7"    # default location for project starting pose
    top_file: "./data/top.prmtop"   # default location for project topology file
    ref_file: "./data/ref.crd"      # default location for project reference file

system:
    protein_mask: "@CA"     # mask for system alignment
    ligand_mask: ":LIG"     # mask for ligand rmsd calculation

md:
    crd_file: "initial.rst7"    # default name for node starting file
    top_file: "top.prmtop"      # default name for node topology file
    trj_file: "traj.nc"         # default name for node trajectory file
    rst_file: "final.rst7"      # default name for node final structure file
    ref_file: "ref.crd"         # default name for node reference file
    log_file: "md.log"          # default name for md logging file
    inp_file: "md.in"           # default name for md input file
    use_hmr: true               # use hydrogen mass repartitioning. default is true
    tstep: 4                    # integration time of 4fs. reduce to 2 if use_hmr is set to false
    iteration_steps: 250000     # md steps. 250k * 4fs = 1 nano per node
    use_gpu: true               # specify md engine version to use. use pmemd.cuda by default.
                                # if use_gpu is set to false, sander is used instead.
                                # it is useful for debugging, but too slow for production.
```

## running in cte-power9

it will be probably already installed over there, but in case it is not, just follow the instructions to install it.  

there's an example submission script at `clusters/cte-power9-example.sh`

## running in kraken

same story as cte-power9.  

the example submissino script is `clusters/kraken-example.sh`
