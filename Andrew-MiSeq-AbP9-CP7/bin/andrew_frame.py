#!/usr/bin/python
import sys
shift = int(sys.argv[1])
report = {}
for ln in sys.stdin:
	if ln[0] == '#': continue
	loc, cnt_ss, cnt_as = [int(s) for s in ln.split()]
	loc_frame = ((loc - shift) / 3 + 1, (loc - shift) % 3 + 1) 
	report[loc_frame] = [cnt_ss, cnt_as]

locs = [loc_frame[0] for loc_frame in report.keys()]

print '#\tFRAME1_SS\tFRAME1_AS\tFRAME2_SS\tFRAME2_AS\tFRAME3_SS\tFRAME3_AS'
for loc in xrange(min(locs), max(locs) + 1):
	ln = [loc]
	for frame in xrange(1, 4):
		ln += report.get((loc, frame), [0, 0])
	print '\t'.join([str(n) for n in ln])
	
