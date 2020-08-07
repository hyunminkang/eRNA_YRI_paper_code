#################################################################

# Uses new alignments and non-mappable regions 
# Finds reads that map to both haplotypes and reads that map to just one
# rep_caplist file should have list of each file prefix to go through
# Same as 3 but for replicates

# Time for demo: 10m47.659s

#################################################################
# Variables to assign

rep_caplist=list/caplist.min.rep.txt # as described previously specifically with replicates
fastq_location=../fastq/

#################################################################

STARTTIME=$(date +%s)
DONE=0
REMAIN=$(cat $rep_caplist | wc -l)
while read NAME SID; do
	if false; then
		if [ $SID != "GM18520r" ]; then
			continue
		fi
	fi
	ID=${SID:0:7}
	echo $SID
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
	cp tmp/map1.bam rawdata/bam@hap1/$SID.TSSc.hap1.bam
	cp tmp/map2.bam rawdata/bam@hap2/$SID.TSSc.hap2.bam
	# Count reads at TSSc regions
	bedtools intersect -v -f 1.0 -a tmp/map1.bam -b hapbed/$ID/nonmap.$ID.hap1.bed3 -sorted > tmp/map1.map.bam 
	bedtools coverage -a hapbed/$ID/TSSc.$ID.hap1.pl.bed6 -b tmp/map1.map.bam -s > tmp/pl1.txt
	bedtools coverage -a hapbed/$ID/TSSc.$ID.hap1.mn.bed6 -b tmp/map1.map.bam -s > tmp/mn1.txt
	paste <(cut -f4,7 tmp/pl1.txt) <(cut -f7 tmp/mn1.txt) |\
		awk '{split($1,a,":"); print a[1]"\t"a[2]"\t"$2+$3"\t"$2"\t"$3}' > readcount/$SID.hap1.txt
	bedtools intersect -v -f 1.0 -a tmp/map2.bam -b hapbed/$ID/nonmap.$ID.hap2.bed3 -sorted > tmp/map2.map.bam 
	bedtools coverage -a hapbed/$ID/TSSc.$ID.hap2.pl.bed6 -b tmp/map2.map.bam -s > tmp/pl2.txt
	bedtools coverage -a hapbed/$ID/TSSc.$ID.hap2.mn.bed6 -b tmp/map2.map.bam -s > tmp/mn2.txt
	paste <(cut -f4,7 tmp/pl2.txt) <(cut -f7 tmp/mn2.txt) |\
		awk '{split($1,a,":"); print a[1]"\t"a[2]"\t"$2+$3"\t"$2"\t"$3}' > readcount/$SID.hap2.txt
	# extract ase reads
	bedtools intersect -a tmp/map1.bam -b hapbed/$ID/snp.het.$ID.hap1.bed3 -sorted > tmp/map1.ase.bam
	bedtools intersect -v -f 1.0 -a tmp/map1.ase.bam -b hapbed/$ID/snp.TSSc.$ID.nonmap.indiv.hap1.bed3 > tmp/map1.ase2.bam
	bedtools intersect -v -f 1.0 -a tmp/map1.ase2.bam -b hapbed/$ID/ref.unmap.$ID.hap1.bed3 > asebam/$SID.ase.hap1.bam
	
	bedtools intersect -a tmp/map2.bam -b hapbed/$ID/snp.het.$ID.hap2.bed3 -sorted > tmp/map2.ase.bam
	bedtools intersect -v -f 1.0 -a tmp/map2.ase.bam -b hapbed/$ID/snp.TSSc.$ID.nonmap.indiv.hap2.bed3 > tmp/map2.ase2.bam
	bedtools intersect -v -f 1.0 -a tmp/map2.ase2.bam -b hapbed/$ID/ref.unmap.$ID.hap2.bed3 > asebam/$SID.ase.hap2.bam
	fi
	CURTIME=$(date +%s)
	ELAPSED=$(($CURTIME-$STARTTIME))
	DONE=$(($DONE+1))
	REMAIN=$(($REMAIN-1))
	ESTIMATED=$(($ELAPSED*$REMAIN/$DONE))
	HR=$((ESTIMATED/3600))
	MIN=$(($ESTIMATED/60-$HR*60))
	SEC=$(($ESTIMATED%60))
	if [ $REMAIN -gt 0 ]; then echo Estimated $HR:$MIN:$SEC left; fi
done < $rep_caplist
