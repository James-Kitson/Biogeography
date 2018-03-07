#!/bin/sh -e

#!/bin/bash
#SBATCH -J Relaxed_BEAST
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -o %A-%a.log
#SBATCH -e %A-%a.err
#SBATCH -p compute

module add R/3.4.1 gcc/6.3.0


for i in ${SLURM_ARRAY_TASK_ID}
do
	cp BEAST_array_unstratified_DEC.R BEAST_array_unstratified_DEC_${SLURM_ARRAY_TASK_ID}.R
	sed -i -e "s/array_tag/${SLURM_ARRAY_TASK_ID}/g" BEAST_array_unstratified_DEC_${SLURM_ARRAY_TASK_ID}.R
done

R CMD BATCH "BEAST_array_unstratified_DECJX_${SLURM_ARRAY_TASK_ID}".R

