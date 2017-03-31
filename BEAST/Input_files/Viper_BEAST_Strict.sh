#!/bin/sh -e

#!/bin/bash
#SBATCH -J Strict_BEAST
#SBATCH -N 1
#SBATCH -o %A-%a.log
#SBATCH -e %A-%a.err
#SBATCH -p compute

export PATH=$PATH:~/applications/beast2/bin

module add libbeagle/2.1.2 java/jdk1.8.0_102


for i in ${SLURM_ARRAY_TASK_ID}
do
	cp Strict_Cratopus_BEAST.xml Strict_Cratopus_BEAST_${SLURM_ARRAY_TASK_ID}.xml
	sed -i -e "s/Strict_Cratopus_BEAST/Strict_Cratopus_BEAST_${SLURM_ARRAY_TASK_ID}/g" Strict_Cratopus_BEAST_${SLURM_ARRAY_TASK_ID}.xml
done

beast "Strict_Cratopus_BEAST_${SLURM_ARRAY_TASK_ID}".xml

