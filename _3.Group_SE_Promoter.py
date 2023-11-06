

import sys
import os

group1 = ['HNT-025', 'HNT-029', 'HNT-047', 'HNT-084', 'HNT-085']
group2 = ['HNT-024', 'HNT-031', 'HNT-037', 'HNT-042', 'HNT-089', 'HNT-104', 'HNT-114', 'HNT-121', 'HNT-134', 'HNT-140', 'HNT-142']
group3 = ['HNT-044', 'HNT-060', 'HNT-064', 'HNT-086', 'HNT-098', 'HNT-112']
group4 = ['HNT-045', 'HNT-051', 'HNT-053', 'HNT-062', 'HNT-070', 'HNT-073', 'HNT-074', 'HNT-075', 'HNT-078','HNT-082', 'HNT-087', 'HNT-090', 'HNT-091', 'HNT-095', 'HNT-097', 'HNT-101', 'HNT-103', 'HNT-108']
group5 = ['HNT-054', 'HNT-055', 'HNT-056', 'HNT-059', 'HNT-065', 'HNT-066', 'HNT-072', 'HNT-080', 'HNT-081', 'HNT-088', 'HNT-092', 'HNT-093', 'HNT-094', 'HNT-096', 'HNT-100', 'HNT-107']

group_dict = {'g1':['HNT-025', 'HNT-029', 'HNT-047', 'HNT-084', 'HNT-085'], 'g2':['HNT-024', 'HNT-031', 'HNT-037', 'HNT-042', 'HNT-089', 'HNT-104', 'HNT-114', 'HNT-121', 'HNT-134', 'HNT-140', 'HNT-142'], 'g3':['HNT-044', 'HNT-060', 'HNT-064', 'HNT-086', 'HNT-098', 'HNT-112'], 'g4':['HNT-045', 'HNT-051', 'HNT-053', 'HNT-062', 'HNT-070', 'HNT-073', 'HNT-074', 'HNT-075', 'HNT-078','HNT-082', 'HNT-087', 'HNT-090', 'HNT-091', 'HNT-095', 'HNT-097', 'HNT-101', 'HNT-103', 'HNT-108'], 'g5':['HNT-054', 'HNT-055', 'HNT-056', 'HNT-059', 'HNT-065', 'HNT-066', 'HNT-072', 'HNT-080', 'HNT-081', 'HNT-088', 'HNT-092', 'HNT-093', 'HNT-094', 'HNT-096', 'HNT-100', 'HNT-107']}

open_gene = open(f'{sys.path[0]}/Gencode_Protein_Coding_POS_list.txt')

gene_promoter_dict = {}
n = 0
for line in open_gene:
	if n == 0:
		n = 1
		continue
	cols = line.rstrip().split()
	gene = cols[0]
	strand = cols[1]
	chrom = cols[2]
	pos1 = int(cols[3])
	pos2 = int(cols[4])
	if strand == '+':
		start = max(0, pos1 - 1000)
		end = pos1
	elif strand == '-':
		start = pos1
		end = pos1 + 1000
	gene_promoter_dict[gene] = [chrom, start, end]
open_gene.close()

SE = open(f'{sys.path[0]}/SE_HNSC_hg19.bed')
SE_dict = {}
n = 0
for line in SE:
	cols = line.rstrip().split()
	if n == 0:
		header = cols
		n = 1
		continue
	se_info = dict(zip(header, cols))
	chrom = se_info['se_chr']
	start = se_info['se_start']
	end = se_info['se_end']
	se_id = se_info['se_id']
	if chrom not in SE_dict:
		SE_dict[chrom] = [[start, end, se_id]]
	else:
		SE_dict[chrom] += [[start, end, se_id]]
SE.close()

def find_overlap(tuple1, tuple2):
	start1, end1 = tuple1
	start2, end2 = tuple2
	if end1 < start2 or start1 > end2:
		return False
	else:
		overlap_start = max(start1, start2)
		overlap_end = min(end1, end2)
	return (overlap_start, overlap_end)

recurrent_group_gene_dict = {}
n = 0
for line in open('recurrent_group_gene_list.csv'):
	cols = line.rstrip().split(',')
	if n == 0:
		header = cols
		n = 1
		continue
	gene_info = dict(zip(header, cols))
	group = gene_info['group']
	gene = gene_info['gene']
	rate = gene_info['recurrent_rate(%)']
	sample_per_group_count = gene_info['sample_count;group_count']
	p = gene_info['p-value']
	gene_info = [gene] + gene_promoter_dict[gene] + [rate, sample_per_group_count, p]
	if group not in recurrent_group_gene_dict:
		recurrent_group_gene_dict[group] = [gene_info]
	else:
		recurrent_group_gene_dict[group] += [gene_info]

interaction_dict = {'g1':{}, 'g2':{}, 'g3':{}, 'g4':{}, 'g5':{}}

#loop_file_5kb = f"{sample.replace('HNT-', 'HNT_')}.Tumor.HiC.ICE.genome-loops.0.9.5000bp.bedpe"
#loop_file_10kb = f"{sample.replace('HNT-', 'HNT_')}.Tumor.HiC.ICE.genome-loops.0.9.10000bp.bedpe"
#loop_file_25kb = f"{sample.replace('HNT-', 'HNT_')}.Tumor.HiC.ICE.genome-loops.0.9.25000bp.bedpe"

loop_padding = 10000
def group_se_promoter(interaction_dict, res):
	for group in group_dict:
		for sample in group_dict[group]:
			recurrent_gene_list = recurrent_group_gene_dict[group]
			loop_file = f"{sample.replace('HNT-', 'HNT_')}.Tumor.HiC.ICE.genome-loops.0.9.{res}bp.bedpe"
			if os.path.exists(loop_file):
				open_loop_file = open(loop_file)
				for loop in open_loop_file:
					loop = loop.rstrip().split()
					[loop_chrom1, loop_start1, loop_end1] = [loop[0], int(loop[1]), int(loop[2])] 
					[loop_chrom2, loop_start2, loop_end2] = [loop[3], int(loop[4]), int(loop[5])]
					loop_start1 = max(0, loop_start1 - loop_padding)
					loop_start2 = max(0, loop_start2 - loop_padding)
					loop_end1 = loop_end1 + loop_padding
					loop_end2 = loop_end2 + loop_padding

					for gene_info in recurrent_gene_list:
						[gene, gene_chrom, gene_start, gene_end, rate, sample_per_group_count, p] = gene_info
						gene = (gene, rate, sample_per_group_count, p)
						loop1_promoter_list = []
						if loop_chrom1 == gene_chrom:
							if find_overlap((loop_start1, loop_end1), (gene_start, gene_end)):
								loop1_promoter_list.append(gene)
						loop2_promoter_list = []
						if loop_chrom2 == gene_chrom:
							if find_overlap((loop_start2, loop_end2), (gene_start, gene_end)):
								loop2_promoter_list.append(gene)
						if not loop1_promoter_list == loop2_promoter_list == []:
							loop1_se_list = []
							loop2_se_list = []
							for se_chrom in SE_dict:
								if loop_chrom1 == se_chrom:
									for [se_start, se_end, se_id] in SE_dict[se_chrom]:
										se_start = int(se_start)
										se_end = int(se_end)
										if find_overlap((loop_start1,loop_end1), (se_start,se_end)):
											loop1_se_list.append([se_chrom, se_start, se_end, se_id, loop_chrom1, loop_start1, loop_end1, loop_chrom2, loop_start2, loop_end2])
								if loop_chrom2 == se_chrom:
									for [se_start, se_end, se_id] in SE_dict[se_chrom]:
										se_start = int(se_start)
										se_end = int(se_end)
										if find_overlap((loop_start2,loop_end2), (se_start,se_end)):
											loop2_se_list.append([se_chrom, se_start, se_end, se_id, loop_chrom1, loop_start1, loop_end1, loop_chrom2, loop_start2, loop_end2])
							if not loop1_se_list == loop2_se_list == []:
								if [] not in [loop1_se_list, loop2_promoter_list]:
									if sample not in interaction_dict[group]:
										interaction_dict[group][sample] = {}
										for gene in loop2_promoter_list:
											interaction_dict[group][sample][gene] = loop1_se_list
									else:
										for gene in loop2_promoter_list:
											if gene not in interaction_dict[group][sample]:
												interaction_dict[group][sample][gene] = loop1_se_list
											else:
												interaction_dict[group][sample][gene] += loop1_se_list
								if [] not in [loop1_promoter_list, loop2_se_list]:
									if sample not in interaction_dict[group]:
										interaction_dict[group][sample] = {}
										for gene in loop1_promoter_list:
											interaction_dict[group][sample][gene] = loop2_se_list
									else:
										for gene in loop1_promoter_list:
											if gene not in interaction_dict[group][sample]:
												interaction_dict[group][sample][gene] = loop2_se_list
											else:
												interaction_dict[group][sample][gene] += loop2_se_list
group_se_promoter(interaction_dict, 5000)
group_se_promoter(interaction_dict, 10000)
group_se_promoter(interaction_dict, 25000)

fo = open(f'Promoter_SE_interaction_fpkm_up_gene_list.csv', 'w')
fo.write('group,gene,interacting_samples,recurrent_rate(%),sample_count;group_count,p-value,se_info\n')

for group in interaction_dict:
	se_dict = {}
	group_gene_dict = {}
	for sample in interaction_dict[group]:
		gene_list = list(interaction_dict[group][sample].keys())
		for gene_rate in gene_list:
			gene = gene_rate[0]
			se_info = sample + ';' + ';'.join([':'.join([str(i) for i in single_se[:4]]) for single_se in interaction_dict[group][sample][gene_rate]])
			if gene not in se_dict:
				se_dict[gene] = [se_info]
			else:
				if se_info not in se_dict[gene]:
					se_dict[gene].append(se_info)
			if gene_rate not in group_gene_dict:
				group_gene_dict[gene_rate] = [sample]
			else:
				group_gene_dict[gene_rate] += [sample]
	for gene_rate in sorted(group_gene_dict.keys()):
		(gene, rate, sample_per_group_count, p) = gene_rate
		fo.write(group + ',' +  gene + ',' + ';'.join(group_gene_dict[gene_rate]) + ',' + rate + ',' +  sample_per_group_count + ',' + p + ',' + '['.join(se_dict[gene]) + '\n')
fo.close()



				









