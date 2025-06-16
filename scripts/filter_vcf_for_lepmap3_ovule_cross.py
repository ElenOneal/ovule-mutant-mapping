'''This program filters vcf file with parents and offscript for importation into LepMap3.'''

#This script was written under Python 3.8.5

##########################DEPENDENCIES#####################

import sys #version 3.7.6
import argparse
import getopt
#import pysam
import os
import doctest
import re
import collections
import tempfile
import shutil
import numpy as np
import warnings
from timeit import default_timer as timer
import time
import logging
import itertools


###########################################################

logger = logging.getLogger() #use to figure out how long it takes to run
logger.setLevel(logging.DEBUG)

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
        
def get_args():

	parser=argparse.ArgumentParser(description='Filter vcf file for LepMap3')
	parser.print_help()
	parser.add_argument('--invcf', type=str, required=True, help="input VCF v 4.1 file")
	parser.add_argument('--ParentList', type=str, required=True, help="input file with parent names, crossing direction (mom or dad), inbred or F1, and minimum and max depth for filtering")
	parser.add_argument('--SampleList', type=str, required=True, help="tab-delimited file with sample Names")
#	parser.add_argument('--IgnoreSamples', type=str, required=True, help="comma separated list of samples to exclude")
	parser.add_argument('--out', type=str, required=True, help = 'Prefix for output file')
	parser.add_argument('--minMapQ', type=int, required=True, help = 'Minimum phred-scaled Mapping Quality')
	parser.add_argument('--filterMissingParents', type=str2bool, nargs='?', const=True, default=False, help="Activate filtering out sites with missing parents.")
	parser.add_argument('--filterInbredParents', type=str2bool, nargs='?', const=True, default=False, help="Activate filtering out heterozygous inbred parents.")
	parser.add_argument('--filterforhets', type=str2bool, nargs='?', const=True, default=False, help="Activate filtering homozygous F1s.")
	parser.add_argument('--minparentdepth',type=str2bool, nargs='?', const=True, default=False, help = 'use depth info from parent file to filter by minimum depth')
	parser.add_argument('--maxparentdepth', type=str2bool, nargs='?', const=True, default=False,help = 'use depth info from parent file to filter by maximum depth')
	parser.add_argument('--FractionF2', type = int, default = 0, help = 'proportion offspring required to pass filter')
	parser.add_argument('--minF2depth', type=int, required=True, default = 0, help = 'minimum depth of F2s to be included')
	parser.add_argument('--missingList', type=str, required=True, help = 'File with fraction of missing data for each backcross offspring')
	parser.add_argument('--MaxMissingData', type = float, required=True, help = 'Maximum missing data allowed for backcross offspring; set to 0 to allow all samples')

	args=parser.parse_args()
	return args

missing = ["./.",".",".|."]
heterozygous = ["0/1","0|1","1/2","1|2","0/2","0|2"]
homozygous = ["1/1", "1|1", "0/0", "0|0", "2/2", "2|2"]
bad_hets = ['0/2','0|2']

               
def open_parent_list(ParFile):

# ParentName	Duplicate	Type	Mindepth	MaxDepth	Direction
# DHRO22	Y	Inbred	10	500	GP1
# GMR2	Y	Inbred	10	500	GP2
# PseudoF1	Y	F1	20	1000	Mom


	parent_names = []
	parent_duplicate = []
	parent_type = []
	min_depth = []
	max_depth = []
	gP1 = []
	gP2 = []
	
	for line in ParFile:

		if line.startswith('ParentName') : continue #this should be header
		else:
			cols = line.replace('\n', '').split('\t')

			name = cols[0]
			dup = cols[1] #crossing direction, mom or dad
			type = cols[2]
			min = int(cols[3])
			max = int(cols[4])
			gpinfo = cols[5]

			parent_names.append(name)
			parent_duplicate.append(dup)
			parent_type.append(type)
			min_depth.append(min)
			max_depth.append(max)
			if gpinfo == 'GP1':
				gP1.append(name)
			if gpinfo == 'GP2':
				gP2.append(name)
			
	parent_info = []
	parent_info.append(parent_names)
	parent_info.append(parent_duplicate)
	parent_info.append(parent_type)
	parent_info.append(min_depth)
	parent_info.append(max_depth)
	parent_info.append(gP1)
	parent_info.append(gP2)
	

	return parent_info


def get_sample_index(vcf,parents,f1,inbred,duplicate,gP1s,gP2s):

	parent_list = list(parents)
	#print(parent_list)
	f1_list = list(f1)
	inbred_list = list(inbred)
	dup_list = list(duplicate)
	gP1s = list(gP1s)
	gP2s = list(gP2s)
	parent_index = []
	f1_index = []
	inbred_index = []
	dup_index = []
	f2_index = []
	Gp1_index = []
	Gp2_index = []

	headerInfo = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] 
	
	with open(vcf, 'r') as vF:
		for line in vF:
			if line.startswith("##"):continue#this should be header
			elif line.startswith('#CHROM'): #this is header with info fields and sample info
				info = line.replace('\n', '').split('\t')
				#print(info)
				info = [i for i in info]
				for p in parent_list:
					if p in info:
						parent_index.append(info.index(p))
				for gP in gP1s:
					if gP in info:
						Gp1_index.append(info.index(gP))
				for gP in gP2s:
					if gP in info:
						Gp2_index.append(info.index(gP))
				for f1 in f1_list:
					if f1 in info:
						f1_index.append(info.index(f1))
				for inbred in inbred_list:
					if inbred in info:
						inbred_index.append(info.index(inbred))
				for dup in dup_list:
					if dup in info:
						dup_index.append(info.index(dup))
				for i in info:
					if i not in headerInfo:
						if i not in parents:
							f2_index.append(info.index(i))
			else: break
	
	return(parent_index,f1_index,inbred_index,f2_index,dup_index,Gp1_index,Gp2_index)


def main():

	args = get_args()
	print(args)

	outprefix = args.out

	fh = logging.FileHandler(outprefix+'.log') #create log file that will record the time to completion
	fh.setLevel(logging.DEBUG) # or any level you want
	logger.addHandler(fh)

	start_time = time.time()
	
	headerInfo = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] 
	
	minMQ = args.minMapQ
	
	fracF2 = float(args.FractionF2)
	
	minF2depth = args.minF2depth
	
	vcffile = args.invcf
	pF = args.ParentList
	samples = args.SampleList

# ParentName	Duplicate	Type	Mindepth	MaxDepth	Direction
# DHRO22	Y	Inbred	10	500	GP1
# GMR2	Y	Inbred	10	500	GP2
# PseudoF1	Y	F1	20	1000	Mom

	ParentFile = open(pF) #File with information about parents
	parentInfo = open_parent_list(ParentFile)
#	print(parentInfo)
	parents = parentInfo[0]
	duplicate = parentInfo[1]
	Type = parentInfo[2]
	mindepthlist = parentInfo[3]
	maxdepthlist = parentInfo[4]
	gp1 = parentInfo[5]
	gp2 = parentInfo[6]

	fracMissing = args.MaxMissingData

	f1_type = [i for i, element in enumerate(Type) if element == 'F1']
	inbred_type = [i for i, element in enumerate(Type) if element == 'Inbred']
	dup_type = [i for i, element in enumerate(duplicate) if element == 'Y']
	gp1_type = [i for i, element in enumerate(parents) if element in gp1]
	gp2_type = [i for i, element in enumerate(parents) if element in gp2]
	#print(f1_type)
	#print(inbred_type)
	#print(dup_type)
	#print(gp1_type)
	#print(gp2_type)

	f1_list = [parents[i] for i in f1_type]
	inbred_list = [parents[i] for i in inbred_type]
	dup_list = [parents[i] for i in dup_type]
	gp1_list = [parents[i] for i in gp1_type]
	gp2_list = [parents[i] for i in gp2_type]
	 
# 	print(f1_list)
# 	print(inbred_list)
	#print(dup_list)
# 	print(gp1_list)
# 	print(gp2_list)
	#print(duplicate,dup_type,dup_list)
	
	outFile = open(outprefix+'.filtered.vcf','w')

#Things to do:
#1a. remove missing parent sites
#1b. remove heterozygous inbred parents
#1c. remove homozygous F1 sites #and missing sites
#2. filter by minimum depth for parents
#3. filter by maximum depth for parents
#4. filter by minimum #F2s
#5. filter by mean F2 depth

# 	missing = ["./.",".",".|."]
# 	heterozygous = ["0/1","0|1","1/2","1|2","0/2","0|2"]
# 	homozygous = ["1/1", "1|1", "0/0", "0|0", "2/2", "2|2"]
# 	bad_hets = ['0/2','0|2']

	si = get_sample_index(vcffile,parents,f1_list,inbred_list,dup_list,gp1_list,gp2_list)
	#print(si)
	parent_index = si[0]
	f1_index = si[1]
	inbred_index = si[2]
	dup_index = si[4]
	f2_index = si[3]
	gp1_index = si[5]
	gp2_index = si[6]
	print(parents,parent_index)
	print(inbred_list,inbred_index)
	print(f1_list,f1_index)
	print(gp1_list,gp1_index)
	print(gp2_list,gp2_index)
#	print(duplicate,dup_type,dup_list,dup_index)
	#print(f2_index)
	print(dup_index)
	

	missing_samples =[]
	
	with open(args.missingList, 'r') as mS:
		for line in mS:
			if line.startswith('INDV'):
				continue
			else:
				cols = line.split('\t')
				name = cols[0]
				fmiss = float(cols[4])
				if fmiss >= fracMissing:
					missing_samples.append(name)
				else: continue

	mS.close()
	
	tab = '\t'
	
	neg_samples = headerInfo + parents + missing_samples 
	
	#get sample index, excluding samples with too much missing data

	sampleList = []
	
	with open(args.invcf, 'r') as vF:
		for line in vF:
			if line.startswith("##"):continue#this should be header
			elif line.startswith('#CHROM'): #this is header with info fields and sample info
				info = line.replace('\n', '').split('\t') #split into tab delimited fields
				info = [i for i in info] #create list of header
				sNames = list(set(info).symmetric_difference(neg_samples))
				sampleList = sNames
				#print(sNames)
			else: break
	 	
	vF.close() #close file

	sampleSize = len(sampleList)
	#print(parentName,parent_index) #should check this before proceeding
#	print(sampleList)
	
	
	dup_names = []
	#print(dup_names)
	
	no_snps = 0	
	no_sites = 0
	good_mapping_qual = 0
	missing_parents_sites = 0
	heterozygous_inbred_sites = 0
	homozygous_f1_sites = 0
	low_parent_depth = 0
	high_parent_depth = 0
	missing_f2_sites = 0
	low_f2_depth = 0
	good_sites = 0
	lowMinF2 = 0
	low_median_f2 = 0
	low_meanf2_depth = 0
	high_genotypes = 0
	low_genotypes = 0
	discordant_alleles = 0
	same_grandparents = 0
	fewer_genotypes = 0
	min_f2s = 0
		
	with open(args.invcf, 'r') as variantFile:
		for line in variantFile:
			if line.startswith("##"): #continue#this should be header and #CHROM line
				outFile.write(line)
			elif line.startswith('#CHROM'): #continue #this is header with info fields and sample info
				line = line.strip('\n')
				cols = line.split('\t')
				dup_info = [cols[i] for i in dup_index]
				dup_string1 = '_D1'
				dup_string2 = '_D2'
				duplicate_names1 = [s + dup_string1 for s in dup_info]
				duplicate_names1 = '\t'.join(duplicate_names1)
				duplicate_names2 = [s + dup_string2 for s in dup_info]
				duplicate_names2 = '\t'.join(duplicate_names2)
				#outFile.write(line+'\t'+duplicate_names1+'\t'+duplicate_names2+'\n')
				outFile.write(line+'\t'+duplicate_names1+'\n')
				print(line+'\t'+duplicate_names1+'\n')
				#print(line+'\t'+duplicate_names1+'\t'+duplicate_names2+'\n')
				#outFile.write(line+'\n')
			#else: break
			else: #continue
				line = line.strip('\n')
				no_sites += 1
				
				cols = line.replace('\n', '').split('\t')
				contig = cols[0]

				pos = cols[1]
				ref = cols[3]
				alt = str(cols[4])
				qual = ''.join([i.strip('MQ=') for i in cols[7].split(";") if re.match('MQ=',i)])
				#print(qual)
				filter = str(cols[6])
				format = str(cols[8])

				#indexEnd = dup_index + dup_index
				indexEnd = dup_index
				#print(indexEnd)
				#indexEnd = dup_index + inbred_index
				indexData = [cols[i] for i in indexEnd]
				#print(indexEnd,indexData)

				#outFile.write(line+'\tCompF1_1\tComp_CSS4_1\tCompF1_2\tComp_CSS4_2\t'+duplicate_names1+'\n')

				if len(ref) > 1: continue #exclude deletions
				if len(alt) > 3: continue #exclude large insertions
				if len(alt) == 3 and ',' not in alt: continue #exclude small insertions
				if len(alt) == 2: continue #exclude small insertions
				if qual == '.': continue
				
				#no_sites += 1

				if alt == '.': continue #exclude invariant sites
# 				if ',' in alt:
# 					print(contig,pos,ref,alt)

				no_snps += 1
				
				if float(qual) < minMQ: continue #these are sites with all missing data
								
				good_mapping_qual += 1
				
				parent_data = [cols[i] for i in parent_index]
				
				missing_parents = 0
				heterozygous_inbreds = 0
				homozygous_f1s = 0
			
				parent_genotypes = []
				
				zipped_parents = tuple(zip(parents,parent_data))
				
				formatList = [i for i in format.split(':')]
				dp = formatList.index("DP")

				if len(dup_index) > 0:
					duplicate_data = '\t'.join([cols[i] for i in dup_index])
				else:
					duplicate_data = ''
#				print(duplicate_data)

				inbred_data = [cols[i] for i in inbred_index]
				f1_data = [cols[i] for i in f1_index]
				gp1_data = [cols[i] for i in gp1_index]
				gp2_data = [cols[i] for i in gp2_index]
				f2_data = [cols[i] for i in f2_index]
				f2_depths = []
				f2_genotypes = []
				
				#print(f2_index)
				#print(f2_data)
# 				print(inbred_data)
# 				print(f1_data)

				f1_geno = [i.split(':')[0] for i in f1_data if i.split(':')[0] not in missing]
				gp1_geno = [i.split(':')[0] for i in gp1_data if i.split(':')[0] not in missing]
				gp2_geno = [i.split(':')[0] for i in gp2_data if i.split(':')[0] not in missing]
				
				bad_site = 0
				
				if args.filterInbredParents:
					for inbred in inbred_data: #get genotypes of inbred parents
						inbred_geno = inbred.split(":")[0]
						if inbred_geno in heterozygous:
							heterozygous_inbreds += 1
				
					if heterozygous_inbreds > 0: #if an inbred parent is heterozygous, exclude site
						#print(contig,pos,inbred_data,zipped_parents,'Heterozygous inbred line')
						heterozygous_inbred_sites += 1
						bad_site += 1
				
				if args.filterforhets: #filter out sites where F1 is homozygous
					f1_data = [cols[i] for i in f1_index]
					for f1 in f1_data:
						f1_geno = f1.split(":")[0]
						if f1_geno in homozygous:
							homozygous_f1s += 1
				
					if homozygous_f1s > 0:
						homozygous_f1_sites += 1
						bad_site += 1

				if args.filterMissingParents:
				
					for parent in parent_data:
						parent_geno = parent.split(":")[0]
						if parent_geno in missing:
							missing_parents += 1
				
					if missing_parents > 0: #remove sites where one parent is missing (noninformative)
						#print(contig,pos,parent_data,zipped_parents,'Missing Parents')
						missing_parents_sites += 1
						bad_site += 1

				if args.minparentdepth or args.maxparentdepth: #filter out minimum or maximum depths based on user preference
					parent_depths = []
					for d in parent_data:
						depth = d.split(":")[dp]
						try:
							depth = int(depth)
						except ValueError:
							depth = 0
						parent_depths.append(depth)
					good_site_min = [i for i, j in zip(parent_depths,mindepthlist) if i >= j and i!= 0]
					good_site_max = [i for i, j in zip(parent_depths,maxdepthlist) if i <= j]

	
				if args.minparentdepth and len(good_site_min) < len(parent_depths):
					#print('Minimum parent depth',contig,pos,parent_depths,mindepthlist,zipped_parents)
					low_parent_depth += 1
					bad_site += 1
					
					
				if args.maxparentdepth and len(good_site_max) < len(parent_depths):
					#print('Maximum parent depth',contig,pos,parent_depths,maxdepthlist,zipped_parents)
					high_parent_depth += 1
					bad_site += 1

				
				if gp1_geno == gp2_geno:
					#print(contig,pos,gp1_geno,gp2_geno)
					same_grandparents += 1
					bad_site += 1		
			
				for i in f2_data:
					f2_geno = i.split(':')[0]
					if f2_geno not in missing:
						f2_depth = i.split(':')[dp]
						try:
							f2_depth = int(f2_depth)
						except ValueError:
							f2_depth = 0
						if f2_depth >= minF2depth:
							f2_depths.append(f2_depth)
							f2_genotypes.append(f2_geno)
				
				unique_f2s = set(f2_genotypes)
				
				if len(unique_f2s) > 3:
					high_genotypes += 1
					bad_site += 1
					#continue
				
				no_f2s = len(f2_genotypes)
				
				frac_offspring = float(no_f2s/len(sampleList))*100

				if frac_offspring < fracF2:
					#print(no_f2s,len(sampleList),frac_offspring)
					missing_f2_sites += 1
					bad_site += 1
					#continue
				
				if bad_site > 0:
					continue
												
				else: 
					good_sites += 1
					indexData = '\t'.join(indexData)
					#print(indexData)
					#print(line+'\t'+indexData+'\n')
					outFile.write(line+'\t'+indexData+'\n')
				

	logger.info('All done.')
	minute_time = float(time.time()/60)
	minute_start = float(start_time/60)
	logger.info('Arguments =  %s' % (args))
	logger.info('--- %d minutes to complete task ---' % (minute_time - minute_start))
	logger.info('There are %d sites' % (no_sites))
	logger.info('There are %d variant SNP sites' % (no_snps))
	logger.info('There are %d sites with MQ >= %d' % (good_mapping_qual,minMQ))
	logger.info('There are %d sites with missing parents, %d with het inbreds, and %d with homozygous F1s' % (missing_parents_sites,heterozygous_inbred_sites,homozygous_f1_sites))
	logger.info('There are %d sites with low parental depth' % (low_parent_depth))
	logger.info('There are %d sites with high parental depth' % (high_parent_depth))
	logger.info('There are %d sites where grandparents have the same genotype' %(same_grandparents))
	logger.info('There are %d sites with fewer than %d percent F2s' % (missing_f2_sites,fracF2))
	logger.info('There are %d sites with more than 3 F2 genotypes' % (high_genotypes))
	logger.info('There are %d final sites' % (good_sites))

	
	variantFile.close()		

p = main()
print(p)
