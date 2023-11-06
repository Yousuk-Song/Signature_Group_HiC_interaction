#!/home/juyoung/anaconda3/envs/py311/bin/python3.11

import sys
import math
import os


def get_tad_len(tad):
	[tad_chrom, (tad_x, tad_x), (tad_x, tad_y,), (tad_y, tad_y)] = tad
	value = 2 * ((tad_x - tad_y)**2)
	return round(math.sqrt(value), 2)

	
gene_loop_file = open('Promoter_SE_info.csv')

group_loop_dict = {}
n = 0
for line in gene_loop_file:
	cols = line.rstrip().split(',')
	if n == 0:
		header = cols
		n = 1
		continue
	line_dict = dict(zip(header, cols))
	group = line_dict['group']
	gene = line_dict['gene']
	sample = line_dict['sample']
	loop_chrom1 = line_dict['loop_chrom1']
	loop_start1 = line_dict['loop_start1']
	loop_end1 = line_dict['loop_end1']
	loop_chrom2 = line_dict['loop_chrom2']
	loop_start2 = line_dict['loop_start2']
	loop_end2= line_dict['loop_end2']
	single_loop = [loop_chrom1, loop_start1, loop_end1, loop_chrom2, loop_start2, loop_end2]
	if group not in group_loop_dict:
		group_loop_dict[group] = {(sample, gene):[single_loop]}
	else:
		if (sample, gene) not in group_loop_dict[group]:
			group_loop_dict[group][(sample, gene)] = [single_loop]
		else:
			if single_loop not in group_loop_dict[group][(sample, gene)]:
				group_loop_dict[group][(sample, gene)].append(single_loop)

group_tad_dict = {}
for group in group_loop_dict:
	for sample_gene in group_loop_dict[group]:
		for single_loop in group_loop_dict[group][sample_gene]:
			[loop_chrom1, loop_start1, loop_end1, loop_chrom2, loop_start2, loop_end2] = single_loop
			loop_start1 = int(loop_start1)
			loop_end1 = int(loop_end1)
			loop_start2 = int(loop_start2)
			loop_end2 = int(loop_end2)
			tad = [loop_chrom1, (loop_start1, loop_start1), (loop_start1, loop_start2), (loop_start2, loop_start2)]
			if group not in group_tad_dict:
				group_tad_dict[group] = {sample_gene:[tad]}
			else:
				if sample_gene not in group_tad_dict[group]:
					group_tad_dict[group][sample_gene] = [tad]
				else:
					if tad not in group_tad_dict[group][sample_gene]:
						group_tad_dict[group][sample_gene].append(tad)

if not os.path.exists('Promoter_SE_TAD_Dir'):
	os.system('mkdir Promoter_SE_TAD_Dir')

recurrent_gene_file = open('Promoter_SE_interaction_fpkm_up_recurrent_gene.csv')

group_tad_file = open('Group_Gene_TAD_list.csv', 'w')
recur_tad_out_dict = {}
n = 0
for line in recurrent_gene_file:
	cols = line.rstrip().split(',')
	if n == 0:
		header = cols
		n = 1
		continue
	line_dict = dict(zip(header, cols))
	group = line_dict['Group']
	gene = line_dict['Gene']
	recur_numb = int(line_dict['Recurrent_numb'])
	recur_sample = line_dict['Recurrent_samples'].split(';')
	if group not in ['g4', 'g5']:
		continue
	if recur_numb < 7:
		continue
	tad_list = []
	for sample in recur_sample:
		sample_gene = (sample, gene)
		for tad in group_tad_dict[group][sample_gene]:
			tad_len = get_tad_len(tad)
			tad_list.append(sample + ':' + str(tad_len))
			key = (group, gene, sample)
			if key not in recur_tad_out_dict:
				recur_tad_out_dict[key] = [tad]
			else:
				recur_tad_out_dict[key].append(tad)
	group_tad_file.write(','.join([group, gene, ';'.join(tad_list)]) + '\n')
group_tad_file.close()
for key in recur_tad_out_dict:
	(group, gene, sample) = key
	if not os.path.exists(f'Promoter_SE_TAD_Dir/{group}'):
		os.system(f'mkdir Promoter_SE_TAD_Dir/{group}')
	if not os.path.exists(f'Promoter_SE_TAD_Dir/{group}/{gene}'):
		os.system(f'mkdir Promoter_SE_TAD_Dir/{group}/{gene}')
	fo = open(f'Promoter_SE_TAD_Dir/{group}/{gene}/{sample}_{gene}_tad.bedpe', 'w')
	for tad in recur_tad_out_dict[key]:
		tad_chrom = tad[0]
		tad_pos1 = str(tad[1][0])
		tad_pos2 = str(tad[3][0])
		fo.write('\t'.join([tad_chrom, tad_pos1, tad_pos2, tad_chrom, tad_pos1, tad_pos2]) + '\n')
	fo.close()







