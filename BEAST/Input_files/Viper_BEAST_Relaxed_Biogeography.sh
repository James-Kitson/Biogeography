#!/bin/sh -e

#!/bin/bash
#SBATCH -J Biogeography
#SBATCH -N 1
#SBATCH -o %A-%a.log
#SBATCH -e %A-%a.err
#SBATCH -p compute

export PATH=$PATH:~/applications/beast2/bin

module add libbeagle/2.1.2 java/jdk1.8.0_102


for i in ${SLURM_ARRAY_TASK_ID}
do
	cp Relaxed_Biogeography_Cratopus_BEAST.xml Relaxed_Biogeography_Cratopus_BEAST_${SLURM_ARRAY_TASK_ID}.xml
	sed -i -e "s/Relaxed_Biogeography_Cratopus_BEAST/Relaxed_Biogeography_Cratopus_BEAST_${SLURM_ARRAY_TASK_ID}/g" Relaxed_Biogeography_Cratopus_BEAST_${SLURM_ARRAY_TASK_ID}.xml
done

beast "Relaxed_Biogeography_Cratopus_BEAST_${SLURM_ARRAY_TASK_ID}".xml

