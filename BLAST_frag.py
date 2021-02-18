#!/usr/bin/python3
#20200205 --output, --kilo-reads-per-blast
#20200204 initial version
import sys, os, argparse, subprocess, hashlib, tempfile
def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true", required=False, default=False)
	parser.add_argument("-q", "--query-sequence", help="query sequence", required=True)
	parser.add_argument("-p", "--pident-of-blast", help="only count blast result that has a pident higher than N", type=float, default=97)
	parser.add_argument("-w", "--wordsize-of-blast", help="wordsize in blast", type=int, default=11)
	parser.add_argument("-e", "--evalue-of-blast", help="evalue in blast", type=float, default=10)
	parser.add_argument("-o", "--output", help="output file, - means stdout", required=True)
	parser.add_argument("-k", "--kilo-reads-per-blast", help="make a blastdb for N thousands of reads", type=int, required=True)
	parser.add_argument("fastq", help="sequencing reads files", nargs='+')
	args = parser.parse_args()
	return args

def err(code, msg):
	print('ERR-%03d:' % code, msg, file=sys.stderr)
	sys.exit(code)

def checkexe(exe, md5):
	try:
		proc = subprocess.run([exe, '-help'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	except:
		err(104, f'Please make sure command "{exe}" is in your PATH!')
	if hashlib.md5(proc.stdout).hexdigest() != md5:
		yn = input(f'Your {exe} version is not \'2.11.0+\', which is this program is tested with. Do you want to continue anyway? (Y/n)')
		if yn.lower().strip() not in ['', 'y']:
			sys.exit(1)

def blast_frag(verbose, filenames, query_sequence, reads_per_blast, pident_of_blast, wordsize_of_blast, evalue_of_blast, output):
	loop, cnt, db = 0, 0, tempfile.NamedTemporaryFile().name
	proc = None
	frame = {'primer_ss':{}, 'primer_as':{}}
	for fn in filenames:
		if verbose:
			print(fn, file=sys.stderr)
		if fn.endswith('.gz'):
			import gzip
			fq = gzip.open(fn, 'rt')
		else:
			fq = open(fn)
		while True:
			lns = [fq.readline() for i in range(4)]
			if lns[0] == '':
				break
			lns = [ln.rstrip() for ln in lns]
			if lns[0][0] != '@':
				err(101, f'fastq {lns[0]} does not start with "@"')
			if len(lns[1]) != len(lns[3]):
				err(102, f'fastq {lns[0]} has different length of sequence and quality')
			if lns[2] != '+':
				err(103, f'fastq {lns[0]} misses "+"')
			id_seq = '>%s\n%s\n' % (lns[0][1:].replace(' ', '-'), lns[1])
			if proc is None:
				if verbose:
					loop += 1
					print('#makeblastdb: loop', loop, file=sys.stderr)
				proc = subprocess.Popen(['makeblastdb', '-in', '-', '-dbtype', 'nucl', '-title', db, '-out', db], stdin=subprocess.PIPE, stdout=subprocess.DEVNULL)
			proc.stdin.write(str.encode(id_seq))
			cnt += 1
			if (cnt % reads_per_blast == 0):
				blastn_and_count(verbose, db, reads_per_blast, query_sequence, pident_of_blast, wordsize_of_blast, evalue_of_blast, proc, frame)
				proc = None
	blastn_and_count(verbose, db, reads_per_blast, query_sequence, pident_of_blast, wordsize_of_blast, evalue_of_blast, proc, frame)
	if output == '-':
		fout = sys.stdout
	else:
		fout = open(output, 'w')
	print('#LOC\tSS\tAS', file=fout)
	loc_range = range(min(min(frame['primer_ss'].keys()), min(frame['primer_as'].keys())), max(max(frame['primer_ss'].keys()), max(frame['primer_as'].keys())) + 1)
	for loc in loc_range:
		print(loc, frame['primer_ss'].get(loc, 0), frame['primer_as'].get(loc, 0), sep='\t', file=fout)
	fout.close()

def blastn_and_count(verbose, db, cnt, query, pident, wordsize, evalue, proc, frame):
	if proc is None: return
	proc.stdin.close()
	rtn = proc.wait()
	if rtn != 0:
		err(106, f'makeblastdb returned error code {rtn}')
	proc = subprocess.Popen(['blastn', '-outfmt', '6', '-word_size', str(wordsize), '-ungapped', '-max_target_seqs', str(cnt), '-evalue', str(evalue), '-db', db, '-query', query], stdin=subprocess.DEVNULL, stdout=subprocess.PIPE)
	blastrtn = []
	print('#blastn', file=sys.stderr)
	for ln in proc.stdout:
		fs = ln.decode().split()
		fs[2] = float(fs[2])
		if fs[2] < pident:
			continue
		for i in range(3, 10):
			fs[i] = int(fs[i])
		for i in range(10, 12):
			fs[i] = float(fs[i])
		blastrtn.append(fs)
#		if fs[1] == '1:M03079:25:000000000-AJJW2:1:1110:5239:9582': print('###daimh2', fs)
		if verbose and len(blastrtn) % 100000 == 0:
			print(len(blastrtn), file=sys.stderr)
	rtn = proc.wait()
	if rtn != 0:
		err(107, f'blastn returned error code {rtn}')
	if verbose:
		print(len(blastrtn), file=sys.stderr)
		print('#sort', file=sys.stderr)
	blastrtn.sort(key=lambda fs : (fs[1], fs[0]))
	counter(blastrtn, frame)
	for suffix in [ 'db', 'hr', 'in', 'ot', 'sq', 'tf', 'to']:
		try:
			os.remove(f'{db}.n{suffix}')
		except FileNotFoundError:
			pass

def countprint(sid, lst, frame):
	if sid == None:
		return
#	if sid == '1:M03079:25:000000000-AJJW2:1:1110:5239:9582': print('###daimh3', sid, lst)
	one_per_seq = []
	prev = None
	for fs in sorted(lst, key=lambda fs : (fs[0], (fs[4] - fs[3]))):
		if fs[8] < fs[9]:
			fs += [fs[8], fs[9], '+']
		else:
			fs += [fs[9], fs[8], '-']
		if fs[0] != prev:
			one_per_seq.append(fs)
			prev = fs[0]
	if len(one_per_seq) == 1:
		key = 'b:ALIGN_TARGET_ONLY'
		return
	elif len(one_per_seq) == 3:
		key = 'c:ALIGN_TO_ALL_THREE'
		return
	elif len(one_per_seq) != 2:
		err(108, f'contact developer with this error message please! {one_per_seq}')
	one_per_seq.sort(key=lambda fs : fs[12])
	if one_per_seq[0][0].startswith('primer_') and one_per_seq[1][0].startswith('primer_'):
		key = 'd:ALIGN_TO_TWO_PRIMERS'
		return
	if not one_per_seq[0][0].startswith('primer_'):
		key = 'e:ALIGN_TO_VWF_THEN_PRIMER'
		return
	overlap = one_per_seq[0][13] + 1 - one_per_seq[1][12] 
	if overlap < 0:
		key = 'f:OVERLAP_LESS_THAN_0'
		return
	if overlap > 2:
		key = 'g:OVERLAP_GREATER_THAN_2'
		return
	key = 'h:OVERLAP_IS_%d' % overlap
	if one_per_seq[1][14] == '+':
		loc = one_per_seq[1][6] + overlap
	elif one_per_seq[1][14] == '-':
		loc = one_per_seq[1][7] - overlap
	else:
		err(109, f'contact developer with this error message please! {one_per_seq}')
#	if loc == 1801: print('####daimh1', one_per_seq)
	loc_cnt = frame[one_per_seq[0][0]]
	loc_cnt[loc] = loc_cnt.get(loc, 0) + 1

def counter(blast, frame):
	print('#count', file=sys.stderr)
	prev_read, lst = None, None
	for fs in blast:
		if fs[1] != prev_read:
			countprint(prev_read, lst, frame)
			prev_read, lst = fs[1], []
		lst.append(fs)
	countprint(prev_read, lst, frame)
	if len(frame['primer_ss']) < 2 or len(frame['primer_as']) < 2:
		err(105, 'there are not enough reads in the fastq file(s)')

def checkquery(filename):
	ids = []
	for ln in open(filename):
		if ln[0] == '>':
			ids.append(ln[1:].strip())
	if len(ids) != 3 or len(set(ids)) != 3 or  "primer_ss" not in ids or "primer_as" not in ids:
		err(110, f'query sequence file "{filename}" should have three sequences "primer_ss", "primer_as" and one VWF sequence')
		
def main():
	args = get_args()
	checkquery(args.query_sequence)
	checkexe('blastn', '18199e01fff784af84fd4d9c102fcabb')
	checkexe('makeblastdb', 'e958f13a0efb06c9f140b3065abd9ab6')
	blast_frag(args.verbose, args.fastq, args.query_sequence, 1000 * args.kilo_reads_per_blast, args.pident_of_blast, args.wordsize_of_blast, args.evalue_of_blast, args.output)

if __name__ == '__main__':
	main()
