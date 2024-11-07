#!/bin/bash -l                                                                                                                                                                                                                                                

#SBATCH -A pdc-bus-2022-4                                                                                                                                                                                                                                     
#SBATCH -p eggnog                                                                                                                                                                                                                                             
#SBATCH -n 1                                                                                                                                                                                                                                                  
#SBATCH -t 05:00:00                                                                                                                                                                                                                                         
#SBATCH -J 4758
#SBATCH --ntasks-per-node=128
        # Email notifications, BEGIN, END, FAIL, REQUEUE, ALL                                                                                                                                                                                                         
module load PDC
#module load anaconda3/2022.05                                                                                                                                                                                                                                
module load gcc/11.2.0
#module load all-spack-modules/0.16.3                                                                                                                                                                                                                         
module load parallel/20230422
python ./reset_all.py 4758
cd /cfs/klemming/projects/snic/pdc-bus-2022-4/Jack/
source my-venv-dardel/bin/activate
cd Departures/DepartureGrid/CoG/
python ./reset_all.py "4758"
for line in "4758"
do
	seq 0 73371 | parallel -j 128 --joblog "joblog_Ti1_${4758}_NLTE.log" python ./create_grid_files_bash.py {}  "$line"
done
