#!/usr/bin/python3
#19500101 initial version daimon:/etc/Model/sh
import sys, os, argparse, subprocess, hashlib, tempfile
def get_args():
	parser = argparse.ArgumentParser()
#	parser.add_argument("square", help="display a square of a given number", type=int, nargs='+', action='append', required=True)
	parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true", required=False, default=False)
	parser.add_argument("-q", "--query-sequence", help="query sequence", required=True)
	parser.add_argument("-p", "--pident-of-blast", help="higt that has less than this pident is ignored", type=float, default=97)
	parser.add_argument("fastq", help="fastq files", nargs='+')
	args = parser.parse_args()
	return args
def err(code, msg):
	print('ERR-%03d:' % code, msg, file=sys.stderr)
	sys.exit(code)
def checkexe(exe, md5):
	proc = subprocess.run([exe, '-help'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	if hashlib.md5(proc.stdout).hexdigest() != md5:
		yn = input(f'Your {exe} version is not \'2.11.0+\', which is this program is tested with. Do you want to continue anyway? (Y/n)')
		if yn.lower().strip() not in ['', 'y']: sys.exit(1)

def makeblastdb(verbose, filenames, query):
	cnt, db = 0, tempfile.NamedTemporaryFile().name
	proc = subprocess.Popen(['makeblastdb', '-in', '-', '-dbtype', 'nucl', '-title', db, '-out', db], stdin=subprocess.PIPE, stdout=subprocess.DEVNULL)
	for fn in filenames:
		if verbose: print('#makeblastdb', fn, file=sys.stderr)
		if fn.endswith('.gz'):
			import gzip
			fq = gzip.open(fn, 'rt')
		else:
			fq = open(fn)
		while True:
			lns = [fq.readline() for i in range(4)]
			if lns[0] == '': break
			lns = [ln.rstrip() for ln in lns]
			if lns[0][0] != '@': err(101, 'wrong fastq line ' + lns[0])
			if len(lns[1]) != len(lns[3]): err(103, 'wrong fastq line ' + lns[0])
			if lns[2] != '+': err(102, 'wrong fastq line ' + lns[0])
#			id_seq = '>%s\n%s\n' % (lns[0][1:].replace(' ', '-'), lns[1])
			cnt += 1
			id_seq = '>s%d\n%s\n' % (cnt, lns[1])
			proc.stdin.write(str.encode(id_seq))
	proc.stdin.close()
	if proc.wait() != 0: err(106, f'failed to run makeblastdb')
	return db, cnt

def blastn(verbose, db, cnt, query, pident):
	proc = subprocess.Popen(['blastn', '-outfmt', '6', '-word_size', '7', '-ungapped', '-max_target_seqs', str(cnt), '-db', db, '-query', query], stdin=subprocess.DEVNULL, stdout=subprocess.PIPE)
	blastrtn = []
	for ln in proc.stdout:
		fs = ln.decode().split()
		if float(fs[2]) < pident: continue
		fs[0], fs[1] = fs[1], fs[0]
		blastrtn.append(fs)
		if verbose and len(blastrtn) % 10000 == 0:
			print('#blastn', len(blastrtn), file=sys.stderr)
	if proc.wait() != 0: err(107, f'failed to run blastn')
	if verbose: print('#blastn', len(blastrtn), '\n#sort', file=sys.stderr)
	blastrtn.sort()
	return blastrtn

def countprint(sid, lst, frame):
	if sid == None: return
	filtered = []
	prev = None
	for fs in sorted(lst, key=lambda fs: (fs[1], fs[4]-fs[3])):
		if fs[8] < fs[9]:
			fs += [fs[8], fs[9], '+']
		else:
			fs += [fs[9], fs[8], '-']
		if fs[1] != prev:
			filtered.append(fs)
			prev = fs[1]
	if len(filtered) == 1: 
		key = 'b:ALIGN_TARGET_ONLY'
		return
	filtered = sorted(filtered, key=lambda fs: fs[12])
	if len(filtered) == 3: 
		key = 'c:ALIGN_TO_ALL_THREE'
		return
	if len(filtered) > 3: raise Exception(f'Err2 {filtered}')
	if filtered[0][1].startswith('primer_') and filtered[1][1].startswith('primer_'):
		key = 'd:ALIGN_TO_TWO_PRIMERS'
		return
	if not filtered[0][1].startswith('primer_'):
		key = 'e:ALIGN_TO_VWF_THEN_PRIMER'
		return
	overlap = filtered[0][13] + 1 - filtered[1][12] 
	if overlap < 0:
		key = 'f:OVERLAP_LESS_THAN_0'
		return
	if overlap > 2:
		key = 'g:OVERLAP_GREATER_THAN_2'
		return
	key = 'h:OVERLAP_IS_%d' % overlap
#	for fs in filtered:
#		print fs
#	print overlap
	if filtered[1][14] == '+':
		loc = filtered[1][6] + overlap
	elif filtered[1][14] == '-':
		loc = filtered[1][7] - overlap
	else:
		raise Exception('Err6')
	t = frame[filtered[0][1]]
	t[loc] = t.get(loc, 0) + 1
def counter(blast):
	frame = {'primer_ss':{}, 'primer_as':{}}
	prev_read, lst = None, None
	for fs in blast:
		if len(fs) != 12: raise Exception('Err1')
		for i in range(3, 10): fs[i] = int(fs[i])
		if fs[0] != prev_read:
			countprint(prev_read, lst, frame)
			prev_read, lst = fs[0], []
		lst.append(fs)
	countprint(prev_read, lst, frame)
	loc_range = range(min(min(frame['primer_ss'].keys()), min(frame['primer_as'].keys())), max(max(frame['primer_ss'].keys()), max(frame['primer_as'].keys())) + 1)
	print('#LOC\tSS\tAS')
	for loc in loc_range:
		print(loc, frame['primer_ss'].get(loc, 0), frame['primer_as'].get(loc, 0), sep='\t')
		
def main():
	args = get_args()
	checkexe('blastn', '18199e01fff784af84fd4d9c102fcabb')
	checkexe('makeblastdb', 'e958f13a0efb06c9f140b3065abd9ab6')
	db, cnt = makeblastdb(args.verbose, args.fastq, args.query_sequence)
	blastrtn = blastn(args.verbose, db, cnt, args.query_sequence, args.pident_of_blast)
#	for ln in blastrtn: print('\t'.join(ln))
#	sys.exit(9)
#	print('daimh1', len(blastrtn))
	counter(blastrtn)
	for suffix in [ 'db', 'hr', 'in', 'ot', 'sq', 'tf', 'to']:
		os.remove(f'{db}.n{suffix}')
if __name__ == '__main__': main()
