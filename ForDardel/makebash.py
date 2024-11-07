

ion = "1"
lines= ["4758", "4759", "4778", "4781", "4797", "4801", "4820", "5022", "5129", "5689", "5716", "5720", "5739", "5866", "6716", "6599", "7852"] # ti1
#lines = ['4444', '4501', '4572'] #ti2
#lines = ["4443", "4468", "3409", "4501", "4571"] # 2
for line in lines:
    with open("bashscripts/Ti"+ion+"_"+line+".pysme", "w") as f:
        print("""cd /cfs/klemming/projects/snic/pdc-bus-2022-4/Jack/
source my-venv-dardel/bin/activate
cd Departures/DepartureGrid/CoG/
sbatch Ti"""+ion+"""_"""+line+""".sh""", file=f)
    with open("bashscripts/Ti"+ion+"_" + line + ".sh", "w") as f:
        print("""#!/bin/bash -l                                                                                                                                                                                                                                                

#SBATCH -A pdc-bus-2022-4                                                                                                                                                                                                                                     
#SBATCH -p eggnog                                                                                                                                                                                                                                             
#SBATCH -n 1                                                                                                                                                                                                                                                  
#SBATCH -t 05:00:00                                                                                                                                                                                                                                         
#SBATCH -J """+line,  file=f)
        print("""#SBATCH --ntasks-per-node=128
        # Email notifications, BEGIN, END, FAIL, REQUEUE, ALL                                                                                                                                                                                                         
module load PDC
#module load anaconda3/2022.05                                                                                                                                                                                                                                
module load gcc/11.2.0
#module load all-spack-modules/0.16.3                                                                                                                                                                                                                         
module load parallel/20230422""", file=f)
        print("python ./reset_all.py "+str(line), file=f)
        print("""cd /cfs/klemming/projects/snic/pdc-bus-2022-4/Jack/
source my-venv-dardel/bin/activate
cd Departures/DepartureGrid/CoG/""", file=f)

        print('python ./reset_all.py "'+line+'"', file=f)

        print('for line in "'+line+'"', file=f)
        print("do", file=f)


        print('\tseq 0 73371 | parallel -j 128 --joblog "joblog_Ti1_${'+line+'}_NLTE.log" python ./create_grid_files_bash.py {}  "$line"', file=f)
        print("done", file=f)
