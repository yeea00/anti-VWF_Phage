#!/usr/bin/python
import sys
#_locs = range(222, 8730)
#_locs = range(1, 8833)
def print_out(fread, id, lst, report, frame):
	if id == None: return
	key = 'a:TOTAL_ALIGNED_READS'
	report[key] = report.get(key, 0) + 1
	filtered = []
	prev = None
	for fs in sorted(lst, key=lambda fs: (fs[0], fs[4]-fs[3])):
		if fs[8] < fs[9]:
			fs += [fs[8], fs[9], '+']
		else:
			fs += [fs[9], fs[8], '-']
		if fs[0] != prev:
			filtered.append(fs)
			prev = fs[0]
	if len(filtered) == 1: 
		key = 'b:ALIGN_TARGET_ONLY'
		report[key] = report.get(key, 0) + 1
		return
	filtered = sorted(filtered, key=lambda fs: fs[12])
	if len(filtered) == 3: 
		key = 'c:ALIGN_TO_ALL_THREE'
		report[key] = report.get(key, 0) + 1
		return
	if len(filtered) > 3: raise Exception('Err2')
	if filtered[0][0].startswith('primer_') and filtered[1][0].startswith('primer_'):
		key = 'd:ALIGN_TO_TWO_PRIMERS'
		report[key] = report.get(key, 0) + 1
		return
	if not filtered[0][0].startswith('primer_'):
		key = 'e:ALIGN_TO_VWF_THEN_PRIMER'
		report[key] = report.get(key, 0) + 1
		return
	overlap = filtered[0][13] + 1 - filtered[1][12] 
	if overlap < 0:
		key = 'f:OVERLAP_LESS_THAN_0'
		report[key] = report.get(key, 0) + 1
		return
	if overlap > 2:
		key = 'g:OVERLAP_GREATER_THAN_2'
		report[key] = report.get(key, 0) + 1
		return
	key = 'h:OVERLAP_IS_%d' % overlap
	report[key] = report.get(key, 0) + 1
#	for fs in filtered:
#		print fs
#	print overlap
	if filtered[1][14] == '+':
		loc = filtered[1][6] + overlap
	elif filtered[1][14] == '-':
		loc = filtered[1][7] - overlap
	else:
		raise Exception('Err6')
	t = frame[filtered[0][0]]
	t[loc] = t.get(loc, 0) + 1
	print >> fread, filtered[0][1]
def main():
	report, frame = {}, {'primer_ss':{}, 'primer_as':{}}
	prev_read, lst = None, None
	fread = open('%s_read.txt' % sys.argv[2], 'w')
	for ln in sys.stdin:
		fs = ln.split()
		if len(fs) != 12: raise Exception('Err1')
		if float(fs[2]) < float(sys.argv[1]): continue
		for i in xrange(3, 10): fs[i] = int(fs[i])
		if fs[1] != prev_read:
			print_out(fread, prev_read, lst, report, frame)
			prev_read, lst = fs[1], []
		lst.append(fs)
	print_out(fread, prev_read, lst, report, frame)
	if len(report) > 12: raise Exception('Err3')
	fout = open('%s_count.txt' % sys.argv[2], 'w')
	for key in sorted(report.keys()):
		print >> fout, '%s\t%d' % (key[2:], report[key])
	print >> fout, 'SUM\t%d' % sum(report.values())
	loc_range = range(min(min(frame['primer_ss'].keys()), min(frame['primer_as'].keys())), max(max(frame['primer_ss'].keys()), max(frame['primer_as'].keys())) + 1)
#	if min(min(frame['primer_ss'].keys()), min(frame['primer_as'].keys())) < min(_locs):
#		print >> sys.stderr, min(frame['primer_ss'].keys()), min(frame['primer_as'].keys()), min(_locs)
#		raise Exception('Err4')
#	if max(max(frame['primer_ss'].keys()), max(frame['primer_as'].keys())) > max(_locs):
#		print >> sys.stderr, max(frame['primer_ss'].keys()), max(frame['primer_as'].keys()), max(_locs)
#		raise Exception('Err5')
	fout = open('%s_location.txt' % sys.argv[2], 'w')
	print >> fout, '#LOC\tSS\tAS'
	for loc in loc_range:
		print >> fout, '%d\t%d\t%d' % (loc, frame['primer_ss'].get(loc, 0), frame['primer_as'].get(loc, 0))
if __name__ == '__main__': main()
