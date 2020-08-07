#################################################################

# Uses new alignments and non-mappable regions 
# Finds reads that map to both haplotypes and reads that map to just one
# caplist file should have list of each file prefix to go through

#Time for demo:19m27.568s

#################################################################
# Variables to assign

caplist=list/caplist.min.txt
fastq_location=../fastq/

#################################################################

mkdir -p rawdata/bam@hap1/
mkdir -p rawdata/bam@hap2/

STARTTIME=$(date +%s)
DONE=0
REMAIN=$(cat $caplist | wc -l)
while read NAME ID; do
	if false; then
		if [ $ID != "GM18505" ]; then
			continue
		fi
	fi
	if [ ! -d hapbed/$ID ]; then
		continue;
	fi
	echo $ID
	if true; then
	# Alignment and UMI collapse 
	gunzip -c $fastq_location"${NAME}".fastq.gz > tmp/raw.fastq
	bowtie -p8 -v0 -m1 -S ebwt/$ID.TSSc.hap1 tmp/raw.fastq > tmp/map1.sam 2>tmp/log2.txt
	bowtie -p8 -v0 -m1 -S ebwt/$ID.TSSc.hap2 tmp/raw.fastq > tmp/map2.sam 2>>tmp/log2.txt
	rm tmp/raw.fastq
	head -3 tmp/map1.sam > tmp/map1.sorted.sam
	awk 'BEGIN{OFS="\t"}NF>3&&$3!="*"{print substr($1,length($1)-7,8)"\t"$0}' tmp/map1.sam |\
        sort -u -k4,4 -k5,5n -k1,1 --parallel=4 -S 50% | cut -f2- >> tmp/map1.sorted.sam 
	rm tmp/map1.sam
	samtools view -Sb tmp/map1.sorted.sam > tmp/map1.bam
	rm tmp/map1.sorted.sam
	head -3 tmp/map2.sam > tmp/map2.sorted.sam
	awk 'BEGIN{OFS="\t"}NF>3&&$3!="*"{print substr($1,length($1)-7,8)"\t"$0}' tmp/map2.sam |\
		sort -u -k4,4 -k5,5n -k1,1 --parallel=4 -S 50% | cut -f2- >> tmp/map2.sorted.sam 
	rm tmp/map2.sam
	samtools view -Sb tmp/map2.sorted.sam > tmp/map2.bam
	rm tmp/map2.sorted.sam
	cp tmp/map1.bam rawdata/bam@hap1/$ID.TSSc.hap1.bam
	cp tmp/map2.bam rawdata/bam@hap2/$ID.TSSc.hap2.bam
	fi
	cp rawdata/bam@hap1/$ID.TSSc.hap1.bam tmp/map1.bam
	cp rawdata/bam@hap2/$ID.TSSc.hap2.bam tmp/map2.bam
	# Count reads at TSSc regions
	bedtools intersect -v -f 1.0 -a tmp/map1.bam -b hapbed/$ID/nonmap.$ID.hap1.bed3 -sorted > tmp/map1.map.bam 
	bedtools coverage -a hapbed/$ID/TSSc.$ID.hap1.pl.bed6 -b tmp/map1.map.bam -s > tmp/pl1.txt
	bedtools coverage -a hapbed/$ID/TSSc.$ID.hap1.mn.bed6 -b tmp/map1.map.bam -s > tmp/mn1.txt
	paste <(cut -f4,7 tmp/pl1.txt) <(cut -f7 tmp/mn1.txt) |\
		awk '{split($1,a,":"); print a[1]"\t"a[2]"\t"$2+$3"\t"$2"\t"$3}' > readcount/$ID.hap1.txt
	bedtools intersect -v -f 1.0 -a tmp/map2.bam -b hapbed/$ID/nonmap.$ID.hap2.bed3 -sorted > tmp/map2.map.bam 
	bedtools coverage -a hapbed/$ID/TSSc.$ID.hap2.pl.bed6 -b tmp/map2.map.bam -s > tmp/pl2.txt
	bedtools coverage -a hapbed/$ID/TSSc.$ID.hap2.mn.bed6 -b tmp/map2.map.bam -s > tmp/mn2.txt
	paste <(cut -f4,7 tmp/pl2.txt) <(cut -f7 tmp/mn2.txt) |\
		awk '{split($1,a,":"); print a[1]"\t"a[2]"\t"$2+$3"\t"$2"\t"$3}' > readcount/$ID.hap2.txt
	if true; then
	# extract ase reads
	bedtools intersect -a tmp/map1.bam -b hapbed/$ID/snp.het.$ID.hap1.bed3 -sorted > tmp/map1.ase.bam
	bedtools intersect -v -f 1.0 -a tmp/map1.ase.bam -b hapbed/$ID/snp.TSSc.$ID.nonmap.indiv.hap1.bed3 > tmp/map1.ase2.bam
	bedtools intersect -v -f 1.0 -a tmp/map1.ase2.bam -b hapbed/$ID/ref.unmap.$ID.hap1.bed3 > asebam/$ID.ase.hap1.bam
	
	bedtools intersect -a tmp/map2.bam -b hapbed/$ID/snp.het.$ID.hap2.bed3 -sorted > tmp/map2.ase.bam
	bedtools intersect -v -f 1.0 -a tmp/map2.ase.bam -b hapbed/$ID/snp.TSSc.$ID.nonmap.indiv.hap2.bed3 > tmp/map2.ase2.bam
	bedtools intersect -v -f 1.0 -a tmp/map2.ase2.bam -b hapbed/$ID/ref.unmap.$ID.hap2.bed3 > asebam/$ID.ase.hap2.bam
	fi
	CURTIME=$(date +%s)
	ELAPSED=$(($CURTIME-$STARTTIME))
	DONE=$(($DONE+1))
	REMAIN=$(($REMAIN-1))
	ESTIMATED=$(($ELAPSED*$REMAIN/$DONE))
	HR=$((ESTIMATED/3600))
	MIN=$(($ESTIMATED/60-$HR*60))
	SEC=$(($ESTIMATED%60))
	echo Estimated $HR:$MIN:$SEC left
done < $caplist
