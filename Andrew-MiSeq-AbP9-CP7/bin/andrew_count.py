#!/usr/bin/python
import sys

def main():
	header = ['TOTAL', 'TOTAL_ALIGNED_READS', 'ALIGN_TARGET_ONLY', 'ALIGN_primer_as_ONLY', 'ALIGN_primer_ss_ONLY', 'ALIGN_TO_ALL_THREE', 'ALIGN_TO_TWO_PRIMERS', 'ALIGN_TO_VWF_THEN_PRIMER', 'OVERLAP_LESS_THAN_0', 'OVERLAP_GREATER_THAN_2', 'OVERLAP_IS_0', 'OVERLAP_IS_1', 'OVERLAP_IS_2', 'SUM']
	#print '\t%s' % '\t'.join(header[:-1])
	stat, total = sys.argv[1:3]
	samples = {}
	for s in sys.argv[3:]:
		samples[s] = {}
		for ln in open('%s/%s_count.txt' % (stat, s)):
			fs = ln[:-1].split('\t')
			if len(fs) != 2: raise Exception('Err1')
			if fs[0] not in header: raise Exception('Err2:' + fs[0])
			if fs[0] in samples[s]: raise Exception('Err3')
			samples[s][fs[0]] = int(fs[1])
		samples[s]['TOTAL'] = int(open('%s/%s.txt' % (total, s)).read())
	lst = []
	for v in samples.values(): lst += v.keys()
	lst = set(lst)
	header = [col for col in header if col in lst and col != 'SUM']
	print '\t%s' % '\t'.join(header)
	for s in sys.argv[3:]:
		fs = [s]
		for col in header:
			fs.append(samples[s].get(col, 0))
		print '\t'.join([str(s) for s in fs])


main()
#line1 = None
#for txt in sys.argv[1:]:
#	sample = txt.split('/')[1].replace('_count.txt', '')
#	fs = open(txt).read().split()
#	for i in xrange(0, len(fs), 2):
#		if fs[i] not in header: raise Exception('Err1: %s' % fs[i])
#	if line1 == None:
#		line1 = [fs[i] for i in xrange(0, len(fs), 2)]
#		print '\t%s' % '\t'.join(line1[:-1])
#	if int(fs[1]) * 2 != int(fs[-1]):
#		raise Exception('Err2')
#	print '%s\t%s' % (sample, '\t'.join([fs[i] for i in xrange(1, len(fs)-2, 2)]))
