#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 2:00:00
#SBATCH --mail-user=<email>
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --job-name="sepchrRef"
#SBATCH --output=~/bbcar_project/out/210109_separate_chr_geneanno_refGene.out

dn="/projects/b1122/saya/bbcar_project/geneanno_refs"
fin="hg38_refGene.txt"

chrname="chr1"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr2"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr3"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr4"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr5"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr6"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr7"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr8"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr9"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr10"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr11"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr12"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr13"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr14"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr15"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr16"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr17"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr18"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr19"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr20"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr21"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chr22"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"

chrname="chrX"
grep "${chrname}\s" ${dn}/${fin} > ${dn}/"hg38_refGene_"${chrname}".txt"
