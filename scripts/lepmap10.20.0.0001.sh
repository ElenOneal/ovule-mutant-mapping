#!/bin/bash
#
#$ -S /bin/bash
#SBATCH --get-user-env
#SBATCH --job-name=10,.20
#SBATCH --output=10.20.out
#SBATCH --error=10.20.err
#SBATCH --cpus-per-task=8
#SBATCH --chdir=/work/
#SBATCH --mem=16G

java -cp /hpc/home/eo22/lepmap3/bin ParentCall2 data=ovule.ped.transpose.txt vcfFile=ovule_f2s.min60.max3.thinned.vcf removeNonInformative=1 > data.ovule_f2s.min60.max3.call

java -cp /hpc/home/eo22/lepmap3/bin Filtering2 data=data.ovule_f2s.min60.max3.call dataTolerance=0.0001 removeNonInformative=1 outputHWE=0 missingLimit=0.30 > data.ovule_f2s.min60.max3.0.0001.call

java -cp /hpc/home/eo22/lepmap3/bin/ SeparateChromosomes2 numThreads=8 grandparentPhase=1 data=data.ovule_f2s.min60.max3.0.0001.call lodLimit=10 theta=.20 > map.10.20.0.0001.ovule_f2s.min60.max3.txt

java -cp /hpc/home/eo22/lepmap3/bin/ OrderMarkers2 numThreads=8 grandparentPhase=1 data=data.ovule_f2s.min60.max3.0.0001.call map=map.10.20.0.0001.ovule_f2s.min60.max3.txt outputPhasedData=1 identicalLimit=0.01 > order.10.20.0.0001.id.ovule_f2s.min60.max3.txt

