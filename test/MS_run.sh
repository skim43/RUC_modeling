#!/bin/bash
# make sure your current directory is  "***/gen_PGFem3D_input_from_Rocpack/test"
# make sure working_dir.txt file shows "***/gen_PGFem3D_input_from_Rocpack/test"
module load matlab/2021b

#matlab -nodisplay -nodesktop -nosplash -r "MS_generate_particles; exit"
matlab -nodesktop -nosplash -r "MS_generate_particles_rand; exit"
#matlab -nodisplay -nodesktop -nosplash -r "MS_generate_packs; exit"
#matlab -nodisplay -nodesktop -nosplash -r "MS_generate_packs_rand; exit"

matlab -nodisplay -nodesktop -nosplash -r "MS_generate_PGFem3D_inputs; exit"


