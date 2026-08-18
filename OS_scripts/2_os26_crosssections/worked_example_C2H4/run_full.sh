#!/bin/bash
#SBATCH --job-name=os26_C2H4_full
#SBATCH --time=96:00:00
#SBATCH --partition=astro2_long
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=16G
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL

module load astro

CWD=$(pwd)
RUNDIR=${CWD}/full_run
mkdir -p ${RUNDIR}
cp ${CWD}/os.input ${RUNDIR}/os.input
cd ${RUNDIR}

${CWD}/os26

echo "Full run done."
