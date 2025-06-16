'''This program prepares a mapping order file from LepMap3 for import into R'''

##########################DEPENDENCIES#####################

import sys #version 3.7.6
###########################################################
file = sys.argv[1] #file with mapping order prepared by LepMap3
outprefix = sys.argv[2] #prefix for output file
minorscafs = sys.argv[3] #separate minor scaffolds True or False
###########################################################

minorscafs = minorscafs

outFile1 = open(outprefix+'_mappingorder.txt','w')
#outFile2 = open(outprefix+'_scaffolds_mappingorder.txt','w')

with open(file, 'r') as mF:
	for line in mF:
		if line.startswith('#java'): continue
		if line.startswith('#marker_number'): continue
		if line.startswith('#***'):
			lG = line.strip('\n')
			sep = ' likelihood'
			lG = lG.split(sep,1)[0]
			lG = lG.replace('#*** LG = ','')
		else:
			cols=line.strip('\n').split('\t')
			contig = cols[0]
			pos = cols[1]
			distance = float(cols[2])
			if minorscafs:
	# 			outFile2 = open(outprefix+'_scaffolds_mappingorder.txt','w')
# 				if 'scaffold' in contig: #continue
# 					outFile2.write(lG+'\t'+contig+'\t'+str(pos)+'\t'+str(distance)+'\n')
# 				else:
				outFile1.write(lG+'\t'+contig+'\t'+str(pos)+'\t'+str(distance)+'\n')
			else:
				if 'scaffold' in contig:
					continue
				else:
					outFile1.write(lG+'\t'+contig+'\t'+str(pos)+'\t'+str(distance)+'\n')