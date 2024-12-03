import sys
import os
import numpy as np
def main():
	num = sys.argv[1]

	inp = './raw.UMAP.for.MT.filter/whole.stat.'+num+'.for.raw.txt'
	f1 = open(inp,'r')
	MT_dict = {}
	for line1 in f1:
		info1 = line1.strip().split('\t')

		if info1[0] == 'orig.ident':
			continue

		cluster_num = int(info1[-1])

		try:
			MT_dict[cluster_num].append(float(info1[4]))
		except KeyError as reason:
			MT_dict[cluster_num] = []
			MT_dict[cluster_num].append(float(info1[4]))
	f1.close()

	oup = './raw.UMAP.for.MT.filter/whole.res.'+num+'.median.MT.txt'
	out = open(oup,'w')

	for cluster_num in sorted(MT_dict.keys()):
		Median_num = np.median(MT_dict[cluster_num])
		info_list = [cluster_num,Median_num]
		out.write('{0}\n'.format('\t'.join(map(str,info_list))))
	out.close()

if __name__ == '__main__':
	main()
