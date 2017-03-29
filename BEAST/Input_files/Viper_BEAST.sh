#!/bin/sh -e

#!/bin/bash
#SBATCH -J BEAST_array
#SBATCH -N 1
#SBATCH -o %A-%a.log
#SBATCH -e %A-%a.err
#SBATCH -p compute

export PATH=$PATH:~/applications/beast2/bin

module add libbeagle/2.1.2


for i in $(seq 1 10)
do
	cp Cratopus_BEAST.xml Cratopus_BEAST_$i.xml
	sed -i -e "s/Cratopus_BEAST/Cratopus_BEAST_${i}/g" Cratopus_BEAST_$i.xml
done

beast "Cratopus_BEAST_${SLURM_ARRAY_TASK_ID}".xml

