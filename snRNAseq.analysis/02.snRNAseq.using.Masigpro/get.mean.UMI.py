import sys
import os
import numpy as np
def main():
	sample = sys.argv[1]
	celltype = sys.argv[2]

	oup = './result.data/'+sample+'.'+celltype+'.mean.UMI.txt'
	out = open(oup,'w')

	inp = './result.data/'+sample+'.'+celltype+'.whole.RNA.gene.data.csv'
	f1 = open(inp,'r')
	k = 0
	for line1 in f1:
		info1 = line1.strip().split(',')

		if k == 0:
			k += 1
			continue

		gid = info1[0]
		meanUMI = np.mean([float(i) for i in info1[1:]])
		info_list = [gid,meanUMI]
		out.write('{0}\n'.format('\t'.join(map(str,info_list))))
	f1.close()
	out.close()



if __name__ == '__main__':
	main()
