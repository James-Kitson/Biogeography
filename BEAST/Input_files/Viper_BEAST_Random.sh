#!/bin/sh -e

#!/bin/bash
#SBATCH -J Random_BEAST
#SBATCH -N 1
#SBATCH -o %A-%a.log
#SBATCH -e %A-%a.err
#SBATCH -p compute

export PATH=$PATH:~/applications/beast2/bin

module add libbeagle/2.1.2 java/jdk1.8.0_102


for i in ${SLURM_ARRAY_TASK_ID}
do
	cp Random_Cratopus_BEAST.xml Random_Cratopus_BEAST_${SLURM_ARRAY_TASK_ID}.xml
	sed -i -e "s/Random_Cratopus_BEAST/Random_Cratopus_BEAST_${SLURM_ARRAY_TASK_ID}/g" Random_Cratopus_BEAST_${SLURM_ARRAY_TASK_ID}.xml
done

beast "Random_Cratopus_BEAST_${SLURM_ARRAY_TASK_ID}".xml

