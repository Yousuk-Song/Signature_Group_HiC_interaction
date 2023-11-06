

import sys
import os
from tadlib.visualize.heatmaps import *

mcool = sys.argv[1]
loop1 = sys.argv[2]
loop2 = loop1.replace('10000bp', '25000bp')


#os.system(f'domaincaller --uri {mcool}::/resolutions/{res} -O {mcool.split("/")[-1].replace(".mcool", ".tad.bed")} --DI-output {mcool.split("/")[-1].replace(".mcool", ".DIs.bedGraph")}')

'''
# prepare metadata
res=25000
metadata = mcool.split('/')[-1].replace('.mcool', '.metadata.txt')
fo = open(metadata, 'w')
fo.write(f'res:{res}\n')
fo.write(f' rep1:{mcool}::/resolutions/{res}')
fo.close()

hitad_out = metadata.replace(".metadata.txt", f".tad.{res}.txt")
hitad_log = metadata.replace(".metadata.txt", ".hitad.log")
#os.system(f'/home/juyoung/HNSC/1.Scripts/TADLib/build/scripts-3.7/hitad -O {hitad_out} -d {metadata} --logFile {hitad_log}')
'''

res=5000
chrom = 'chr6'
plot_window = 800000
tss = 134639250
start = tss - plot_window
end = tss + plot_window
print(f'{chrom}:{start}-{end}')

vis = Triangle(f'{mcool}::/resolutions/{res}', chrom, start, end)
vis.matrix_plot()
vis.plot_loops(loop1)
vis.plot_loops(loop2)
vis.outfig(mcool.split('/')[-1].replace('.mcool', f'.vis.{res}.png'))





