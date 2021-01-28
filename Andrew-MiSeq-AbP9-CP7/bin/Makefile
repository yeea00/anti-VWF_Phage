SAMPLES=$(shell ls INPUT/fq.gz/|awk -F "_R[12]" '{print $$1}' |sort -u)

all : OUTPUT/plotframe OUTPUT/plotvwf OUTPUT/count.xls
	echo "andrew analysis is done" | mail daimh@umich.edu
	touch $@
list : 
	echo $(SAMPLES)| sed -e "s/ /\n/g" | nl
OUTPUT/plotframe : $(addprefix step1.5.frame/,$(SAMPLES))
	pqcheck -r SGE -n 30 $^
	mkdir -p $(@D)
	ls $(<D)/*.txt | bin/andrew_plotframe.py $@
	gnuplot < $@_sense.plot
	gnuplot < $@_antisense.plot
	rm $@_sense.plot $@_antisense.plot
	touch $@
OUTPUT/plotvwf : $(addprefix step1.4.stat/,$(SAMPLES))
	pqcheck -r SGE -n 30 $^
	mkdir -p $(@D)
	ls step1.4.stat/*_location.txt|bin/andrew_plotvwf.py $@
	gnuplot < $@_sense.plot
	gnuplot < $@_antisense.plot
	rm $@_sense.plot $@_antisense.plot
	touch $@
OUTPUT/count.xls : $(addprefix step1.4.stat/,$(SAMPLES)) $(addprefix step1.6.total/,$(SAMPLES))
	pqcheck -r SGE -n 30 $^
	mkdir -p $(@D)
	bin/andrew_count.py step1.4.stat step1.6.total $(SAMPLES) > $@.tmp
	mv $@.tmp $@

step1.6 : $(addprefix step1.6.total/,$(SAMPLES))
step1.5 : $(addprefix step1.5.frame/,$(SAMPLES))
step1.4 : $(addprefix step1.4.stat/,$(SAMPLES))
step1.3 : $(foreach S,$(SAMPLES),$(foreach R,1 2,step1.3.blastn/$S_R$R))
step1.2 : $(foreach S,$(SAMPLES),$(foreach R,1 2,step1.2.blastdb/$S_R$R))
step1.1 : $(foreach S,$(SAMPLES),$(foreach R,1 2,step1.1.split/$S_R$R))
define TMPL_SAMPLES
step1.6.total/$1 : step1.1.split/$1_R1
	pqsub -r SGE -t $$@ -c "cat $$<.???? | wc -l > $$@.txt" $$^
step1.5.frame/$1 : step1.4.stat/$1
	pqsub -r SGE -t $$@ -c "bin/andrew_frame.py 251 < $$<_location.txt > $$@.txt" $$^
step1.4.stat/$1 : step1.3.blastn/$1_R1 step1.3.blastn/$1_R2
	pqsub -r SGE -t $$@ -c "(cat $$(<D)/$1_R1.????.txt|sed -e 's/\t/\t1:/'; cat $$(<D)/$1_R2.????.txt|sed -e 's/\t/\t2:/') | bin/andrew_stat.py 97 $$@" $$^
endef
$(foreach S,$(SAMPLES),$(eval $(call TMPL_SAMPLES,$S)))

NUM=4000000
define TMPL_SAMPLE_ENDS
step1.3.blastn/$1_R$2 : step1.2.blastdb/$1_R$2 INPUT/fa/andrew.fa
	for I in `cat step1.1.split/$1_R$2.lst`; do pqsub -r SGE -t $$@.$$$$I -c "blastn -word_size 7 -ungapped -max_target_seqs $(NUM) -query INPUT/fa/andrew.fa -db $$<.$$$$I -outfmt 6 | Csort -S 7g -t '	' -k 2,2 -k 1,1 > $$@.$$$$I.txt" $$<.$$$$I || exit 1; done
	cat $$@.???? > $$@
step1.2.blastdb/$1_R$2 : step1.1.split/$1_R$2
	pqcheck -r SGE $$^
	for I in `cat step1.1.split/$1_R$2.lst`; do pqsub -r SGE -t $$@.$$$$I -c "makeblastdb -in $$<.$$$$I -dbtype nucl -out $$@.$$$$I" || exit 1; done
	cat $$@.???? > $$@
step1.1.split/$1_R$2 : INPUT/fq.gz/$1_R$2_001.fastq.gz
	pqsub -r SGE -t $$@ -c 'gunzip -c $$^ | sed -n "s/^@/>/; s/ .*//; 1~4p; 2~4p" | split -a 4 -l $(NUM) --verbose - $$@. |rev |cut -c 2-5 | rev > $$@.lst'
endef
$(foreach S,$(SAMPLES),$(foreach R,1 2,$(eval $(call TMPL_SAMPLE_ENDS,$S,$R))))

test :
	echo $(SAMPLES)
