with open("hgmd_pro_allmut_2013.4") as file:
	out = open('out.vcf', 'w')
	out.write("##fileformat=VCFv4.1\n")
	out.write('\t'.join(['#CHROM', 'POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','FIRST'])+'\n')
	for line in file:
		if not line: continue
		if line[0] == '#': continue
		info = line.split('\t')
		if not info[2]: continue
		if not info[3]: continue
		if info[1] == '41243454':continue
		#if info[0][3:] == '17':continue
		out.write(info[0][3:] + '\t' + info[1] + '\t' + '.' + '\t' + info[2] + '\t' + info[3] + '\t' + '50' + '\t' + 'PASS' + '\t' + '.' + '\t' + 'GT' + '\t' + './.' +'\t\n')
 
			 	
			
