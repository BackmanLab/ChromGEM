# ChromGEM
By Ruyi Gong and Rivaan Kakkaramadam

## Overview

The **Chromatin Geometry Empowered Model (ChromGEM)** is a three-layer computational framework that models how the ground-truth organization of chromatin is transformed by fluorescence labeling and single-molecule localization microscopy (SMLM) measurement physics. The framework begins with a heterogeneous distribution of nucleosomes, simulates the accessibility and binding of fluorescent probes, and then simulates the stochastic localization events detected by an SMLM experiment. The resulting data can be compared directly with the underlying nucleosome distribution to quantify how probe geometry and fluorophore photokinetics influence the apparent structure of chromatin.

The three layers are connected sequentially:

1. **Layer 1 — SR-EV Chromatin Generation:** generates ground-truth 3D nucleosome coordinates representative of electron microscopy-observed chromatin in a cell line of choice.
3. **Layer 2 — Molecular Dynamics Simulation of Label Probes:** places mobile EdU- or BrdU-associated probes around the fixed nucleosome configuration, simulates probe diffusion and binding, and removes unbound probes in a simulated wash.
4. **Layer 3 — Monte Carlo Simulation of Fluorophore Emissions:** converts the post-wash probe coordinates into stochastic SMLM localization events using experimentally informed blinking and localization-uncertainty distributions.

The downstream workflows operate on voxel-based 3D TIFF volumes rather than coordinate lists. The nucleosome, probe, or simulated-localization coordinates are therefore binned into a 3D grid and smoothed with a 3D Gaussian filter before domain identification and analysis. The workflows then identify packing-domain centers, extract individual three-dimensional domains, calculate their radial density and mass-scaling properties, and generate the plots used in Figure 3d–e.

```text
Layer 1: SR-EV nucleosome coordinates
              |
              v
Layer 2: accessible, bound-probe coordinates
              |
              v
Layer 3: simulated SMLM localization coordinates
              |
              v
Coordinate binning + 3D Gaussian blur (radius = 5 pixels)
              |
              v
Voxel-based 3D TIFF volumes
              |
              v
3D domain identification and domain analysis
              |
              v
Radial-density plots
```

## Repo Contents

* [BrdU_Execution_Scripts](https://github.com/BackmanLab/SMLM-SREV/tree/490e57f64faff08d02154158d3eefdb6bcfec7a4/BrdU_Execution_Scripts): The scripts used to simulate BrdU labeling (Layer 2 of ChromGEM) with an example script provided for running through a cluster.
* [EdU_Execution_Scripts](https://github.com/BackmanLab/SMLM-SREV/tree/490e57f64faff08d02154158d3eefdb6bcfec7a4/EdU_Execution_Scripts): The scripts used to simulate EdU labeling (Layer 2 of ChromGEM) with an example script provided for running through a cluster.

## Inputs/Outputs

### Input

The inputs into this code are model chromatin configurations outputted by the Self-Returning Excluded Volume (SR-EV) polymer model). The computational procedure for SR-EV model generation was obtained from published work [Carignano, Kroeger, et al., eLife 2024] (https://elifesciences.org/articles/97604). The code for SR-EV generation can be made available from the co-first authors of the prior publication upon request. Two user-selected parameters are inputted into SR-EV: alpha, which regulates the clustering behavior in the SR-EV polymer, and phi, which modulates the volume fraction occupied by the model chromatin in the simulation box (i.e., the number of nucleosomes in the simulation). Here, we inputted 1.15 for $\alpha$ and 0.12 for $\phi$ to generate six model HCT116 chromatin distributions. They are provided in the [Inputs](https://github.com/BackmanLab/SMLM-SREV/tree/3700d2414475549b374b2fe3bf918d11cdbcc6eb/Inputs) folder. 

## System Requirements

### Hardware Requirements

* It is recommended to allocate at least 30GB and 16 processors to each simulation. Each simulation will take up to 48 hours to run. ChromGEM can be run in serial (simulated SMLM in one chromatin configuration at a time) or in parallel (simulated SMLM in multiple chromatin configurations at once).

### Software Requirements

The following dependencies are required to run ChromGEM:
1.	LAMMPS Molecular Dynamics Software (Mar 3, 2020 version) built with the ASPHERE, MPI, MPIIO, OpenCL, and OpenMP packages. The instructions for installation with these packages are provided in the LAMMPS User Guide here: https://docs.lammps.org/Packages.html 
2.	OpenMPI (version 4.1.4). The release can be downloaded from the Open MPI site here: https://www.open-mpi.org/software/ompi/v4.1/
3.	Fortran 90. The instructions for obtaining Fortran can be found here: https://fortran-lang.org/learn/os_setup/install_gfortran/ 

## Workflow

The order of operations is the following:

1. Go to the directory containing the SREV configuration (ideally should have 1 SREV configuration per folder).

2. Copy over the fluorophore diffusivity simulation scripts to that directory.

3. Load the LAMMPS Molecular Dynamics and MPI packages.

4. Enter SREV configuration filenames into the fluorophore diffusivity scripts.

5. Initialize the random positions of the dye labels within the simulation box containing the SREV configuration.

5. Run a brief LAMMPS simulation to generate a version of the SREV configuration that lacks overlapping nucleosomes.

6. Run a brief LAMMPS simulation to output the number of SREV nucleosomes that will be used in the simulation. This is vital for properly setting the temperature of the mobile objects at the desired level.

7. Run the first iteration of the simulation. This simulation run will output a trajectory of all 10000 dye labels over the course of this simulation iteration.

For *n* number of user-specified iterations, loop through Steps #8-10 below:

> 8. Input the dye label trajectory file to perform the distance-based colocalization check + appropriate tagging of all dye labels and generate a new dye label input file containing those tags.
> 9. Run a brief simulation to obtain output files that give the number of immobilized dye labels and stuck SREV nucleosomes in the simulation. This is important for recomputing the parameter needed to properly set the temperature of the mobile objects at the desired level.
> 10. Run the next iteration of the simulation with the new dye label input file. This simulation run will output a trajectory of all 10000 dye labels over the course of this simulation iteration.

11. Simulate a "wash" step of the simulated labeled chromatin sample by filtering out all dye labels that have not bound to an SREV nucleosome by that time. Output the coordinates of these dye labels still present in the model chromatin post-wash. This output is your post-wash labeled model chromatin sample. 


## BEFORE YOU BEGIN:

Set up your directory of SREV configurations in the following manner:


- Have a production directory that is for containing all of your SMLM simulations in your set of input SREV configurations. This directory should ONLY contain SREV configurations (configs) that were made with the same input parameters (ex: all the SREV configs in this directory being made with input alpha = 1.15 & phi = 0.12 is good, some SREV configs in this directory being made with input alpha = 1.15 & phi = 0.12 and other SREV configs in this directory being made with input alpha = 1.15 & phi = 0.16 is bad).
  

- For the *n* number of SREV configs in this production directory, have *n* directories (one for each SREV config). Each SREV config should be assigned its own index that is between 1 to *n*. The SREV config file name should be named as "config-$_i_.dump", where $_i_ is the index of that SREV config. Each of the directories within the production directory should also be assigned its own index that is between 1 to *n*. Each of these directories within the production directory should be named as "Config_$_i_", where $_i_ is the assigned index of that directory. Each of these directories should only contain the SREV config that is of the same index as the directory.

<ins>**Specify the simulation box dimensions.**</ins> 


- This should match the SREV config file dimensions. All these scripts are in units of reduced units or "ru", where 1 ru = 10 nm.
Go to InitWalkers.f90, ImmobilConversion.f90, and FinalWash.f90, and ensure all xlength, ylength, and zlength variables are consistent with the box dimensions indicated in the SREV config file and that in the InitWalkers.f90 file all xlowerbound, ylowerbound, and zlowerbound variables are consistent with the lowerbound dimensions stipulated in lines 6-8 of the SREV config file. To identify the xlowerbound, ylowerbound, and zlowerbound information in the SREV config file, understand that line 6 specifies the lower and upper limits of the box in the x direction, line 7 specifies the lower and upper limits of the box in the y direction, and line 8 specifies the lower and upper limits of the box in the z direction. In lines 6-8, the order of the values is lower limit in that dimension and then upper limit in that dimension. Understand that these box dimensions are given in ru and in scientific notation. So, if line 7 says " -1.0000000000000000e+02 1.0000000000000000e+02", that means the box dimension in the y dimension is -100 ru to 100 ru. This means you should go to the lines that define the ylowerbound and ylength variables and set ylowerbound as -100 and ylength as 200 (as y_upper_limit - y_lower_limit = ylength).


<ins>**Update the bash file.**</ins>


Open the execute_sim_in_array.sh file in EdU_Execution_Scripts and in BrdU_Execution_Scripts. Make the following updates:
- In line #3: change the range specified in the definition of array to match the range of SREV configs you wish to simulate SMLM labeling in. So, if you would like to simulate SMLM labeling in 100 SREV configs, change the "_#SBATCH --array=1-10_" line to be "_#SBATCH --array=1-100_".

  
- In line #13: change the production directory for the location of the directory containing your SREV configs and utilize the SLURM_ARRAY_TASK_ID variable to tell the supercomputer how to locate each SREV config. So, if you have 100 SREV configs in a directory with the name "_/home/my_SMLM_simulations/SREV_configs_made_with_phi_0.12_alpha_1.15/_" and in this directory, each SREV configuration of index _i_ is in a directory with a name of "_SREV_Configi_" (ex: "SREV_Config1", "SREV_Config2", "SREV_Config3", etc.), change the "_dirname="/projects/p31375/SMLM-SREV-main/BrdU_Phi08/Config_${SLURM_ARRAY_TASK_ID}"_" line to be "_dirname="/home/my_SMLM_simulations/SREV_configs_made_with_phi_0.12_alpha_1.15/SREV_Config_${SLURM_ARRAY_TASK_ID}"_".

  
- In line 19: change the source directory line for the location of the directory containing the BrdU and EdU SMLM simulation execution scripts. These are the locations of the BrdU_Execution_Scripts folder and EdU_Execution_Scripts folder on your machine/workspace. So, if the BrdU_Execution_Scripts folder and EdU_Execution_Scripts folder are in a directory with the name _"/home/SMLM_Execution_Scripts/BrdU_Execution_Scripts/"_ and _"/home/SMLM_Execution_Scripts/EdU_Execution_Scripts/"_ respectively, in the execute_sim_in_array.sh file located in BrdU_Execution_Scripts, change the "_cp /projects/p31375/SMLM-SREV-main/BrdU_Execution_Scripts/* ._" line to "_cp /home/SMLM_Execution_Scripts/BrdU_Execution_Scripts/* ._" and in the execute_sim_in_array.sh file located in EdU_Execution_Scripts, change the "_cp /projects/p31375/SMLM-SREV-main/EdU_Execution_Scripts/* ._" line to "_cp /home/SMLM_Execution_Scripts/EdU_Execution_Scripts/* ._" . Make sure to preserve the "_* ._" at the end!

  
- In line 2 and 9: change the line "_#SBATCH --account=p31375_" to specify instead the project you would be running this through on Quest and change the "_#SBATCH --output=log.%j_" to your preferred name for output logs.

## HOW TO EXECUTE:

Go to the directory containing the execution scripts (so EdU_Execution_Scripts or BrdU_Execution_Scripts). After you complete all the steps in this README's "BEFORE YOU BEGIN:" section, type "sbatch execute_sim_in_array.sh". This will submit the SMLM simulation jobs over all your inputted SREV configs to Quest.

## Conceptual Notes

### Why is redefining the temperature parameters important?

Over the course of all the simulation iterations, it is important that you are collecting molecule localizations that occur at the same temperature. When specifying the input target temperature parameters, LAMMPS takes that to specify the temperature of the whole simulation. Since the vast majority of the simulation particles are SREV nucleosome ellipsoids that are kept stuck in place, this complicates the actual energy that gets distributed among the mobile particles to have them move at the specified temperature. Therefore, the input temperature parameters that are given to the LAMMPS simulator program should take this into account so that the mobile particles are actually diffusing at your desired target temperature (in this case the script is set to simulate the objects at a target temperature T = 2.50*). Since the number of particles the simulation should mark as immobilized gets updated with each iteration, the temperature to input is re-computed each iteration to account for that updated number.

The equation that computes the temperature to input to produce the desired temperature is the following: number_of_mobile_particles / (number_of_mobile_particles + number_of_immobile_particles) * desired_temp = temp_to_input

## References and Acknowledgements

The simulations scripted here were performed using the LAMMPS Molecular Dynamics Simulator<sup>1</sup>. For more on LAMMPS, consult https://www.lammps.org . Interactions between the particles were modeled using LAMMPS' *pair gayberne* framework for modeling Gay-Berne interaction potentials between ellipsoidal particles<sup>2</sup>.

1. LAMMPS - A flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales, Comp. Phys. Comm. 271, 108171 (2022)

2. Brown, W. M., Petersen, M. K., Plimpton, S. J., & Grest, G. S. (2009). Liquid crystal nanodroplets in solution. *The Journal of Chemical Physics, 130*(4), 044901. https://doi.org/10.1063/1.3058435

We acknowledge the contributions of Dr. Marcelo Carignano in determining the formula for obtaining the desired simulation temperature based on the number of mobile particles and immobilized particles.
