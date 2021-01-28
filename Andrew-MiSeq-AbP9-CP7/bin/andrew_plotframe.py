#!/usr/bin/python
import sys
#echo -e "set origin 0,0\nset size 1,0.125\nplot '$(<D)/9170.txt' title '9170' with lines"; 
#echo -e "set origin 0,0.125\nset size 1,0.125\nplot '$(<D)/9171.$$E.1.txt' title '9171.$$E.1' with lines, '$(<D)/9171.$$E.2.txt' title '9171.$$E.2' with lines, '$(<D)/9171.$$E.3.txt' title '9171.$$E.3' with lines"; \
def plot(f_list, ss_as, column):
	fout = open('%s_%s.plot' % (sys.argv[1], ss_as), 'w')
	print >> fout, "set terminal png size 1000,6000"
	print >> fout, "set output '%s_%s.png'" % (sys.argv[1], ss_as)
	print >> fout, "set multiplot"
	n, height = 0, 1.0 / len(f_list)
	for txt in f_list:
		id = txt.replace('.txt', '').split('/')[1]
		print >> fout, "set origin 0,%f" % (height*n)
		print >> fout, "set size 1,%f" % height
		print >> fout, "set title '%s'" % id
		print >> fout, "plot",
		for i in xrange(1, 4):
			if i > 1: print >> fout, ',',
			print >> fout, "'%s' using 1:%d title 'frame_%s' with lines" % (txt, i*2 + column - 2, i),
		n += 1
		print >> fout

if __name__ == '__main__':
	f_list = [txt[:-1] for txt in sys.stdin]
	plot(f_list, 'sense', 2)
	plot(f_list, 'antisense', 3)
