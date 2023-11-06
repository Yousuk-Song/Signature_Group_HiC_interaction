

import sys
import os

group1 = ['HNT-025', 'HNT-029', 'HNT-047', 'HNT-084', 'HNT-085']
group2 = ['HNT-024', 'HNT-031', 'HNT-037', 'HNT-042', 'HNT-089', 'HNT-104', 'HNT-114', 'HNT-121', 'HNT-134', 'HNT-140', 'HNT-142']
group3 = ['HNT-044', 'HNT-060', 'HNT-064', 'HNT-086', 'HNT-098', 'HNT-112']
group4 = ['HNT-045', 'HNT-051', 'HNT-053', 'HNT-062', 'HNT-070', 'HNT-073', 'HNT-074', 'HNT-075', 'HNT-078','HNT-082', 'HNT-087', 'HNT-090', 'HNT-091', 'HNT-095', 'HNT-097', 'HNT-101', 'HNT-103', 'HNT-108']
group5 = ['HNT-054', 'HNT-055', 'HNT-056', 'HNT-059', 'HNT-065', 'HNT-066', 'HNT-072', 'HNT-080', 'HNT-081', 'HNT-088', 'HNT-092', 'HNT-093', 'HNT-094', 'HNT-096', 'HNT-100', 'HNT-107']

all_group = group1 + group2 + group3 + group4 + group5

group_up_gene_dict = {'g1':{},'g2':{}, 'g3':{}, 'g4':{}, 'g5':{}}
def run(cn, sample):
	### FPKM ###
	gene_pos = open('Gencode_Protein_Coding_POS_list.txt') 
	fpkm = '200805.NCC.ICGC.merged.RNA-seq.fpkm.log2.txt'
	open_fpkm = open(fpkm)
	fpkm_dict = {}
	n = 0
	for line in open_fpkm:
		cols = line.rstrip().split()
		if n == 0:
			header = [i.replace('-T', '').replace('.FIXED', '') for i in cols]
			n = 1
			continue
		fpkm_line = dict(zip(header, cols))
		fpkm_dict[fpkm_line['gene']] = fpkm_line[sample]
	### Gene promoter region ###
	n = 0
	tss_dict = {}
	gene_pos_dict = {}
	for line in gene_pos:
		cols = line.rstrip().split()
		if n == 0:
			header = cols
			n = 1
			continue
		gene_line = dict(zip(header, cols))
		gene = gene_line['gene_name']
		if gene not in fpkm_dict:
			continue
		fpkm = float(fpkm_dict[gene])
		strand = gene_line['strand']
		chrom = gene_line['chrom']
		start = int(gene_line['start'])
		end = int(gene_line['end'])
		if strand == '+':
			tss = (chrom, start)
		elif strand == '-':
			tss = (chrom, end)
		gene_pos_dict[gene] = tss
		if tss in tss_dict:
			tss_dict[tss].append((gene, fpkm))
		else:
			tss_dict[tss] = [(gene, fpkm)]
	### CNV region ###
	cn_dict = {'2':[], '3':[], 'over_4':[]}
	for tss in tss_dict:
		tss_chrom = tss[0]
		tss_pos = tss[1]
		for (gene, fpkm) in tss_dict[tss]:
			if fpkm < 2:
				continue
			overlap_with_cn = 0
			n = 0
			open_cn = open(cn)
			for line in open_cn:
				if line[0] == '#':
					continue
				cols = line.rstrip().split()
				seg_chrom = 'chr' + cols[0]
				seg_start = int(cols[1])-1
				seg_end = int(cols[1]) -1 + 100000
				CNt = cols[4]
				if CNt == 'NA':
					CNt = 0
				else:
					CNt = float(CNt)
				if seg_chrom == tss_chrom:
					if seg_start <= tss_pos <= seg_end:
						overlap_with_cn = 1
						if 1.5 < CNt < 2.5:
							cn_dict['2'] += tss_dict[tss]
						elif 2.5 <= CNt < 4:
							cn_dict['3'] += tss_dict[tss]
						elif CNt >= 4:
							cn_dict['over_4'] += tss_dict[tss]
				del line
			open_cn.close()
#			if not overlap_with_cn:
#				cn_dict['2'] += tss_dict[tss]
	for cn in cn_dict:
		if cn not in ['2']:
			continue
		for (gene, fpkm) in cn_dict[cn]:
			fpkm = float(fpkm)
			if fpkm < 2:
				continue
			(gene_chrom, gene_tss) = gene_pos_dict[gene]
			print(gene,chrom, gene_tss)
			if sample in group1:
				if gene not in group_up_gene_dict['g1']:
					group_up_gene_dict['g1'][gene] = [fpkm]
				else:
					group_up_gene_dict['g1'][gene] += [fpkm]
			elif sample in group2:
				if gene not in group_up_gene_dict['g2']:
					group_up_gene_dict['g2'][gene] = [fpkm]
				else:
					group_up_gene_dict['g2'][gene] += [fpkm]
			elif sample in group3:
				if gene not in group_up_gene_dict['g3']:
					group_up_gene_dict['g3'][gene] = [fpkm]
				else:
					group_up_gene_dict['g3'][gene] += [fpkm]
			elif sample in group4:
				if gene not in group_up_gene_dict['g4']:
					group_up_gene_dict['g4'][gene] = [fpkm]
				else:
					group_up_gene_dict['g4'][gene] += [fpkm]
			elif sample in group5:
				if gene not in group_up_gene_dict['g5']:
					group_up_gene_dict['g5'][gene] = [fpkm]
				else:
					group_up_gene_dict['g5'][gene] += [fpkm]
	gene_pos.close()
	open_fpkm.close()

def resource(name):
	return [f'{name}.tumor.mpileup.100kbcov.absCN.gen_fi', name]

for name in all_group:
	[cn, sample] = resource(name)
	if os.path.exists(cn):
		run(cn, sample)

fo = open(f'Up_Gene_without_CNV.recurrent.csv', 'w')
fo.write('group,gene,average_fpkm,sample_count,group_count\n')

all_gene_list = []
for group in group_up_gene_dict:
	for gene in group_up_gene_dict[group]:
		all_gene_list.append(gene)


g_len = {'g1':len(group1), 'g2':len(group2), 'g3':len(group3), 'g4':len(group4), 'g5':len(group5)}
for group in group_up_gene_dict:
	for gene in group_up_gene_dict[group]:
		count = len(group_up_gene_dict[group][gene])
		average_fpkm = sum(group_up_gene_dict[group][gene])/len(group_up_gene_dict[group][gene])
		fo.write(f"{group},{gene},{round(average_fpkm,2)},{count},{g_len[group]}\n")
		
fo.close()




