#!/home/juyoung/anaconda3/envs/py311/bin/python3.11

import sys
from scipy.stats import chi2_contingency
from scipy.stats import kruskal
import matplotlib.pyplot as plt
import numpy as np


group1 = ['HNT-025', 'HNT-029', 'HNT-047', 'HNT-084', 'HNT-085']
group2 = ['HNT-024', 'HNT-031', 'HNT-037', 'HNT-042', 'HNT-089', 'HNT-104', 'HNT-114', 'HNT-121', 'HNT-134', 'HNT-140', 'HNT-142']
group3 = ['HNT-044', 'HNT-060', 'HNT-064', 'HNT-086', 'HNT-098', 'HNT-112']
group4 = ['HNT-045', 'HNT-051', 'HNT-053', 'HNT-062', 'HNT-070', 'HNT-073', 'HNT-074', 'HNT-075', 'HNT-078','HNT-082', 'HNT-087', 'HNT-090', 'HNT-091', 'HNT-095', 'HNT-097', 'HNT-101', 'HNT-103', 'HNT-108']
group5 = ['HNT-054', 'HNT-055', 'HNT-056', 'HNT-059', 'HNT-065', 'HNT-066', 'HNT-072', 'HNT-080', 'HNT-081', 'HNT-088', 'HNT-092', 'HNT-093', 'HNT-094', 'HNT-096', 'HNT-100', 'HNT-107']

group_len = {'g1':len(group1), 'g2':len(group2), 'g3':len(group3), 'g4':len(group4), 'g5':len(group5)}

D = {}
n = 0
for line in open('Up_Gene_without_CNV.recurrent.csv'):
	cols = line.rstrip().split(',')
	if n == 0:
		header = cols
		n = 1
		continue
	gene_dict = dict(zip(header, cols))
	group = gene_dict['group']
	gene = gene_dict['gene']
	sample_count = int(gene_dict['sample_count'])
	group_count = int(gene_dict['group_count'])
	if gene not in D:
		D[gene] = {'g1':[0,5], 'g2':[0,11], 'g3':[0,6], 'g4':[0,18], 'g5':[0,16]}
		D[gene][group] = [sample_count, group_count - sample_count]
	else:
		D[gene][group] = [sample_count, group_count - sample_count]

group_list = ['g1','g2','g3','g4','g5']

def plot_rate(L):
	L = sorted(L)
	x = np.arange(len(L))
	plt.plot(L, linestyle='-', color='blue', label='Recurrent Rate(%)')
	plt.axhline(y=35, color='red', linestyle='--', label='Threshold=35%')
	plt.xlabel('Gene')
	plt.ylabel('Recurrent Rate(%)')
	plt.legend()
	plt.grid(True)
	plt.savefig('rate_plot.png')
	
rate_list = []
group_recurrent_dict = {'g1':[],'g2':[],'g3':[],'g4':[],'g5':[]}
for gene in D:
	observed_data = [[], []]
	for sublist in D[gene].values():
		for idx, value in enumerate(sublist):
			observed_data[idx].append(value)
	total_counts = list(group_len.values())
	chi2, p, dof, expected = chi2_contingency(observed_data)

	if p < 0.05:
		total_count_list = [5, 11, 6, 18, 16]
		condition_counts_list = observed_data[0]
		non_condition_counts_list = observed_data[1]

		count_ratio_list = [condition / (condition + non_condition) for condition, non_condition in zip(condition_counts_list, non_condition_counts_list)]
		most_frequent_list = [group_list[i] for i in range(len(condition_counts_list)) if condition_counts_list[i] == max(condition_counts_list)]	

		recurrent_rate = round(max(count_ratio_list) * 100, 2)
		rate_list.append(recurrent_rate)
		if recurrent_rate < 35:
			continue
#		print(gene)
#		print(condition_counts_list)
#		print(non_condition_counts_list)
#		print(round(chi2, 2), round(p, 2), round(dof))
#		print()
		highest_recur_group = group_list[count_ratio_list.index(max(count_ratio_list))]
		if highest_recur_group in most_frequent_list:
			recurrent_group = highest_recur_group
			recurrent_count = max(condition_counts_list)
		else:
			group_ratio_dict = dict(zip(group_list, [round(i*100, 2) for i in count_ratio_list]))
			continue
			count_dict = dict(zip(group_list, condition_counts_list))
			ratio_dict = dict(zip(group_list, count_ratio_list))

			case1 = highest_recur_group
			case1_score = (ratio_dict[case1]+0.5)*count_dict[case1]
			(case2, case2_score) = sorted([(case, (ratio_dict[case]+0.5)*count_dict[case]) for case in most_frequent_list], key = lambda x:x[-1])[0]
			if case2_score > case1_score:
				recurrent_group = case2
				recurrent_count = count_dict[case1]
				recurrent_rate = group_ratio_dict[recurrent_group]
			elif case2_score < case1_score:
				recurrent_group = case1
				recurrent_count = count_dict[case2]
				recurrent_rate = group_ratio_dict[recurrent_group]
		group_recurrent_dict[recurrent_group].append([gene, recurrent_rate, recurrent_count,group_len[recurrent_group], p])

fo = open(f'recurrent_group_gene_list.csv', 'w')
fo.write('group,gene,recurrent_rate(%),sample_count;group_count,p-value\n')
for group in group_recurrent_dict:
	for i in group_recurrent_dict[group]:
		[gene, rate, sample_count, group_count, p] = i
		fo.write(f'{group},{gene},{rate},{sample_count};{group_count},{p}\n')
fo.close()

plot_rate(rate_list)



