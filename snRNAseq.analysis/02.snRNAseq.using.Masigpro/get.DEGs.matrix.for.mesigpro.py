import sys
import os
def main():
	celltype = sys.argv[1]
	sample_list = ['G5-D3-1', 'G5-D3-2', 'G5-D9-1', 'G5-D9-2', 'G5-D13-1', 'G5-D13-2']

	inp = 'Mm39.110.PCG.no.MT.genes.txt'
	f1 = open(inp,'r')
	rep_dict = {}
	for line1 in f1:
		info1 = line1.strip().split('\t')
		rep_dict[info1[3]] = info1[-2]
	f1.close()

	inp = './DEGs.overlap.annotated/'+celltype+'.DEGs.overlap.between.times.txt'
	f1 = open(inp,'r')
	DEG_list = []
	for line1 in f1:
		info1 = line1.strip().split('\t')

		DEG_list.append(info1[0])
	f1.close()

	UMI_dict = {}
	for sample in sample_list:
		UMI_dict[sample] = {}
		inp = './result.data/'+sample+'.'+celltype+'.mean.UMI.txt'
		f1 = open(inp,'r')
		for line1 in f1:
			info1 = line1.strip().split('\t')

			UMI_dict[sample][info1[0]] = info1[1]
		f1.close()

	oup = celltype+'.DEGs.mean.UMI.matrix.txt'
	out = open(oup,'w')
	info_list = ['gid']+sample_list
	out.write('{0}\n'.format('\t'.join(map(str,info_list))))

	for gid in DEG_list:
		info_list = [rep_dict[gid]]
		for sample in sample_list:
			info_list.append(UMI_dict[sample][gid])
		out.write('{0}\n'.format('\t'.join(map(str,info_list))))
	out.close()



if __name__ == '__main__':
	main()
