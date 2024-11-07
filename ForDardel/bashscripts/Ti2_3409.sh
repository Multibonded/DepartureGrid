#!/bin/bash -l                                                                                                                                                                                                                                                

#SBATCH -A pdc-bus-2022-4                                                                                                                                                                                                                                     
#SBATCH -p eggnog                                                                                                                                                                                                                                             
#SBATCH -n 1                                                                                                                                                                                                                                                  
#SBATCH -t 06:00:00                                                                                                                                                                                                                                         
#SBATCH -J 3409
#SBATCH --ntasks-per-node=128
        # Email notifications, BEGIN, END, FAIL, REQUEUE, ALL                                                                                                                                                                                                         
#SBATCH --mail-user=jack.mallinson@astro.su.se                                                                                                                                                                                                                
#SBATCH --mail-type=ALL
module load PDC
#module load anaconda3/2022.05
module load gcc/11.2.0
#module load all-spack-modules/0.16.3
module load parallel/20230422
python ./reset_all.py "3409"
for line in "3409"
do
	seq 0 73371 | parallel -j 128 --joblog "joblog_Ti1_${3409}_NLTE.log" python ./create_grid_files_bash.py {}  "$line"
done
