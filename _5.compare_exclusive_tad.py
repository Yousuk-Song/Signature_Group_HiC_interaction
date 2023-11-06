#!/home/juyoung/anaconda3/envs/py311/bin/python3.11

import sys
import os
import math
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

group1 = ['HNT-025', 'HNT-029', 'HNT-047', 'HNT-084', 'HNT-085']
group2 = ['HNT-024', 'HNT-031', 'HNT-037', 'HNT-042', 'HNT-089', 'HNT-104', 'HNT-114', 'HNT-121', 'HNT-134', 'HNT-140', 'HNT-142']
group3 = ['HNT-044', 'HNT-060', 'HNT-064', 'HNT-086', 'HNT-098', 'HNT-112']
group4 = ['HNT-045', 'HNT-051', 'HNT-053', 'HNT-062', 'HNT-070', 'HNT-073', 'HNT-074', 'HNT-075', 'HNT-078','HNT-082', 'HNT-087', 'HNT-090', 'HNT-091', 'HNT-095', 'HNT-097', 'HNT-101', 'HNT-103', 'HNT-108']
group5 = ['HNT-054', 'HNT-055', 'HNT-056', 'HNT-059', 'HNT-065', 'HNT-066', 'HNT-072', 'HNT-080', 'HNT-081', 'HNT-088', 'HNT-092', 'HNT-093', 'HNT-094', 'HNT-096', 'HNT-100', 'HNT-107']

group_dict = {'g1':group1,'g2':group2,'g3':group3,'g4':group4,'g5':group5}

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

open_gene = open('Gencode_Protein_Coding_POS_list.txt') 
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

def find_overlap(tuple1, tuple2):
	start1, end1 = tuple1
	start2, end2 = tuple2
	if end1 < start2 or start1 > end2:
		return False
	else:
		overlap_start = max(start1, start2)
		overlap_end = min(end1, end2)
		return (overlap_start, overlap_end)


loop_padding = 10000

def group_se_promoter(group_list, interaction_dict, gene_info, res):
	[gene, gene_chrom, gene_start, gene_end] = gene_info
	for group in group_list:
		for sample in group_dict[group]:
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
					loop1_promoter_list = []
					loop2_promoter_list = []
					if loop_chrom1 == gene_chrom:
						if find_overlap((loop_start1, loop_end1), (gene_start, gene_end)):
							loop1_promoter_list.append(gene)
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
											if loop1_se_list not in interaction_dict[group][sample][gene]:
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
											if loop2_se_list not in interaction_dict[group][sample][gene]:
												interaction_dict[group][sample][gene] += loop2_se_list

def get_tad_len(tad):
	value = 2 * ((tad[1] - tad[2])**2)
	return round(math.sqrt(value), 2)


box_plot_numb_dict = {}
box_plot_len_dict = {}
for line in open('Group_Gene_TAD_list.csv'):
	cols = line.rstrip().split(',')
	group = cols[0]
	gene = cols[1]
	gene_info = [gene] + gene_promoter_dict[gene]
	tad_numb_list = [[i.split(':')[0] for i in cols[2].split(';')].count(j) for j in set([i.split(':')[0] for i in cols[2].split(';')])]
	tad_len_list = [float(i.split(':')[1]) for i in cols[2].split(';')]
	tad_numb = len(tad_len_list)
	other_group_list = [i for i in group_dict if i != group]
	interaction_dict = {'g1':{}, 'g2':{}, 'g3':{}, 'g4':{}, 'g5':{}}
	del interaction_dict[group]

	box_plot_numb_dict[gene] = {group:tad_numb_list}
	box_plot_len_dict[gene] = {group:tad_len_list}
#	print(gene)
#	print(group, tad_numb_list)
#	print(group, tad_len_list)
	group_se_promoter(other_group_list, interaction_dict, gene_info, 5000)
	group_se_promoter(other_group_list, interaction_dict, gene_info, 10000)
	group_se_promoter(other_group_list, interaction_dict, gene_info, 25000)
	for other_group in interaction_dict:
		other_tad_numb_list = []
		other_tad_len_list = []
		for other_sample in interaction_dict[other_group]:
			for other_gene in interaction_dict[other_group][other_sample]:
				other_tad_list = []
				for loop in interaction_dict[other_group][other_sample][other_gene]:
					loop = loop[4:]
					tad_chrom = loop[0]
					tad_pos1 = float(loop[1])
					tad_pos2 = float(loop[4])
					tad = [tad_chrom, tad_pos1, tad_pos2]
					if tad not in other_tad_list:
						other_tad_list.append(tad)
				other_tad_len_list += [get_tad_len(tad) for tad in other_tad_list]
				other_tad_numb = len(other_tad_list)
				other_tad_numb_list.append(other_tad_numb)
		box_plot_numb_dict[gene][other_group] = other_tad_numb_list
		box_plot_len_dict[gene][other_group] = other_tad_len_list
#		print(other_group, other_tad_numb_list)
#		print(other_group, other_tad_len_list)
#	print()


def box_plot(data, gene, type):
	if not os.path.exists('Boxplot'):
		os.system(f'mkdir Boxplot')
	df = pd.DataFrame.from_dict(data, orient='index').T
	plt.figure(figsize=(4, 8))
	palette = {key: 'r' if key == list(data.keys())[0] else 'g' for key in data.keys()}
	sns.set(style="whitegrid")  # Set the style
	sns.boxplot(data=df, palette=palette, width=0.5)
	plt.title(f'{gene} TAD {type} distribution')
	plt.xlabel('Groups')
	plt.ylabel(f'TAD {type}')
	plt.savefig(f'Boxplot/{gene}_TAD_{type}.dist.boxplot.png')
	plt.close()

def hist(data, gene):
	if not os.path.exists('Histogram'):
		os.system(f'mkdir Histogram')
	df = pd.DataFrame.from_dict(data, orient='index', columns=['Values'])
	plt.figure(figsize=(4, 8))
	values = list(data.values())
	palette = {key: 'r' if key == list(data.keys())[0] else 'g' for key in data.keys()}
	sns.set(style="whitegrid")

	ax = sns.barplot(data=df, x=df.index, y='Values', palette=palette, width=0.5)
	labels = list(data.keys())
	plt.xlabel('Group')
	plt.ylabel('TAD Length')
	plt.xticks(rotation=45, ha='right')
	plt.title('Histogram of smallest TAD length')
	plt.tight_layout()
	plt.savefig(f'Histogram/{gene}_Smallest_TAD_Length.hist.png')
	plt.close()

for gene in box_plot_numb_dict:
	data = box_plot_numb_dict[gene]
	box_plot(data, gene, 'Number')


for gene in box_plot_len_dict:
	data = box_plot_len_dict[gene]
	smallest_len_dict = {}
	for group in data:
		smallest_tad_len = sorted(data[group])[0]
		smallest_len_dict[group] = smallest_tad_len
	hist(smallest_len_dict, gene)
	box_plot(data, gene, 'Length')





