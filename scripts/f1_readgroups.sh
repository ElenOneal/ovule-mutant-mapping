#!/bin/bash
set -euo pipefail

# Add read groups to pseudo-F1 parent BAMs in the 03_align directory.
# Assumes this script is launched from the scripts/ directory.

cd 03_align
mkdir -p temp/PseudoF1

java -Xmx8g -jar /hpc/home/picard.jar AddOrReplaceReadGroups \
  RGLB=DHRO22_S23 RGPL=ILLUMINA RGPU=DHRO22 RGSM=PseudoF1 \
  I=DHRO22_S23.RG.bam O=DHRO22_S23.F1.RG.bam \
  VALIDATION_STRINGENCY=LENIENT USE_JDK_DEFLATER=true USE_JDK_INFLATER=true TMP_DIR=temp/PseudoF1

java -Xmx8g -jar /hpc/home/picard.jar AddOrReplaceReadGroups \
  RGLB=DHRO22_S28 RGPL=ILLUMINA RGPU=DHRO22 RGSM=PseudoF1 \
  I=DHRO22_S8.RG.bam O=DHRO22_S8.F1.RG.bam \
  VALIDATION_STRINGENCY=LENIENT USE_JDK_DEFLATER=true USE_JDK_INFLATER=true TMP_DIR=temp/PseudoF1

java -Xmx8g -jar /hpc/home/picard.jar AddOrReplaceReadGroups \
  RGLB=GMR2_S22 RGPL=ILLUMINA RGPU=GMR2 RGSM=PseudoF1 \
  I=GMR2_S22.RG.bam O=GMR2_S22.F1.RG.bam \
  VALIDATION_STRINGENCY=LENIENT USE_JDK_DEFLATER=true USE_JDK_INFLATER=true TMP_DIR=temp/PseudoF1

java -Xmx8g -jar /hpc/home/picard.jar AddOrReplaceReadGroups \
  RGLB=GMR2_S7 RGPL=ILLUMINA RGPU=GMR2 RGSM=PseudoF1 \
  I=GMR2_S7.RG.bam O=GMR2_S7.F1.RG.bam \
  VALIDATION_STRINGENCY=LENIENT USE_JDK_DEFLATER=true USE_JDK_INFLATER=true TMP_DIR=temp/PseudoF1

samtools index DHRO22_S23.F1.RG.bam
samtools index DHRO22_S8.F1.RG.bam
samtools index GMR2_S22.F1.RG.bam
samtools index GMR2_S7.F1.RG.bam
