###################################################################

# Takes previously made tiling sequence for each haplotype and aligns each to haplotype genomes and h19 genome
# Makes list of non-mappable regions and creates list of non-mappable SNPs
# Time for demo: 2m25.178s

###################################################################
# Variables to set

caplist=list/caplist.min.txt # list of files with file name in first column and corresponding sample name in second
bowtie_genome=hg19/hg19

##########################################################



STARTTIME=$(date +%s)
DONE=0
REMAIN=$(cat $caplist | wc -l)

# Map mini-haploid reference fastq files and extract non-mappable SNPs
# Loop through all fastq files
if true; then
while read FILE ID; do
	if false; then
		if [ $ID != "GM18499" ]; then
			continue
		fi
	fi

	echo Processing $ID snp reads fastq file
	if [ ! -f fastq/snp.TSSc.reads.$ID.hap1.fastq ]; then
		echo No fastq file
		continue
	fi
	# Make bowtie reference
	bowtie-build -q fasta/$ID.TSSc.hap1.fa ebwt/$ID.TSSc.hap1
	bowtie-build -q fasta/$ID.TSSc.hap2.fa ebwt/$ID.TSSc.hap2

	# Make list of het snps in the individual
	awk '{split($1,a,"|"); print a[2];getline;getline;getline}' fastq/snp.TSSc.reads.$ID.hap1.fastq | uniq > tmp/het.snp.txt
	awk 'BEGIN{OFS="\t"}NR==FNR{split($4,a,"|"); pos[a[1]]=$3;next}
		{print "hap1\t"pos[$1]-1"\t"pos[$1]}' \
		hapbed/$ID/snp.TSSc.$ID.hap1.bed6 tmp/het.snp.txt | sort -k2,2n -k3,3n > hapbed/$ID/snp.het.$ID.hap1.bed3
	awk 'BEGIN{OFS="\t"}NR==FNR{split($4,a,"|"); pos[a[1]]=$3;next}
		{print "hap2\t"pos[$1]-1"\t"pos[$1]}' \
		hapbed/$ID/snp.TSSc.$ID.hap2.bed6 tmp/het.snp.txt | sort -k2,2n -k3,3n > hapbed/$ID/snp.het.$ID.hap2.bed3

	# Map hap1 reads to hap1 reference 
	bowtie -p8 -v0 -m1 -S ebwt/$ID.TSSc.hap1 fastq/snp.TSSc.reads.$ID.hap1.fastq > tmp/snp.reads.hap1.sam 2> tmp/log.txt
	# Split unique and multimap reads
	awk '$3=="*"{print "@"$1 > "tmp/hap1.multi.txt";n++;next}
		$3=="hap1"{print "@"$1"\n"$10"\n+\n"$11 > "tmp/hap1.uniq.fastq"}
		END{print "\nhap1 reads multimapped to hap1: "n"\n";}' \
		tmp/snp.reads.hap1.sam >> tmp/log.txt
	# Map unique hap1 reads to hap2 reference
	bowtie -p8 -v0 -a -S ebwt/$ID.TSSc.hap2 tmp/hap1.uniq.fastq > tmp/snp.reads.hap12.sam 2>> tmp/log.txt
	# Remove any mapped reads
	awk '$3=="hap2"{print "@"$1 >> "tmp/hap1.multi.txt";n++;next}
		$3=="*"{print "@"$1"\n"$10"\n+\n"$11 > "tmp/hap12.uniq.fastq"}
		END{print "\nhap1 reads mapped to hap2: "n"\n";}' \
		tmp/snp.reads.hap12.sam >> tmp/log.txt
	# Map to hg19 genome (allowing 2 mismatches)
	bowtie -p8 -v2 -m1 -S --max tmp/hap1hg19.multi.fastq $bowtie_genome tmp/hap12.uniq.fastq \
		> tmp/snp.reads.hap1hg19.sam 2>> tmp/log.txt
	# Remove multi-mapped reads
	awk '{print $1; getline; getline; getline;n++}END{print "\nhap1 reads multimapped to hg19: "n"\n" >> "tmp/log.txt"}' \
		tmp/hap1hg19.multi.fastq >> tmp/hap1.multi.txt
	# Remoce single mapped reads that are not mapped to the SNP site
	awk 'substr($1,1,1)!="@"&&$3!="*"{split($1,a,/[:|]/); if($3!="chr"a[2]||$4<a[3]-60||$4>a[3]+30) {print "@"$1;++n}}
		END{print "\nhap1 reads mismapped to hg19: "n"\n" >> "tmp/log.txt"}' \
		tmp/snp.reads.hap1hg19.sam >> tmp/hap1.multi.txt
	awk '{split($1,a,/[:|]/); print "chr"a[2]"\t"a[3]"\t"a[2]":"a[3]"\t"a[6]}' tmp/hap1.multi.txt | sort -k1,1 -k2,2n -k4,4n |\
		uniq > tmp/hap1.nonmap.snp.txt 

	# Map hap2 reads to hap2 reference 
	bowtie -p8 -v0 -m1 -S ebwt/$ID.TSSc.hap2 fastq/snp.TSSc.reads.$ID.hap2.fastq > tmp/snp.reads.hap2.sam 2>> tmp/log.txt
	# Split unique and multimap reads
	awk '$3=="*"{print "@"$1 > "tmp/hap2.multi.txt";n++;next}
		$3=="hap2"{print "@"$1"\n"$10"\n+\n"$11 > "tmp/hap2.uniq.fastq"}
		END{print "\nhap2 reads multimapped to hap2: "n"\n";}' \
		tmp/snp.reads.hap2.sam >> tmp/log.txt
	# Map unique hap2 reads to hap1 reference
	bowtie -p8 -v0 -a -S ebwt/$ID.TSSc.hap1 tmp/hap2.uniq.fastq > tmp/snp.reads.hap21.sam 2>> tmp/log.txt
	# Remove any mapped reads
	awk '$3=="hap1"{print "@"$1 >> "tmp/hap2.multi.txt";n++;next}
		$3=="*"{print "@"$1"\n"$10"\n+\n"$11 > "tmp/hap21.uniq.fastq"}
		END{print "\nhap2 reads mapped to hap1: "n"\n";}' \
		tmp/snp.reads.hap21.sam >> tmp/log.txt
	# Map to hg19 genome (allowing 2 mismatches)
	bowtie -p8 -v2 -m1 -S --max tmp/hap2hg19.multi.fastq $bowtie_genome tmp/hap21.uniq.fastq \
		> tmp/snp.reads.hap2hg19.sam 2>> tmp/log.txt
	# Remove multi-mapped reads
	awk '{print $1; getline; getline; getline;n++}END{print "\nhap2 reads multimapped to hg19: "n"\n" >> "tmp/log.txt"}' \
		tmp/hap2hg19.multi.fastq >> tmp/hap2.multi.txt
	# Remove single mapped reads that are not mapped to the SNP site
	awk 'substr($1,1,1)!="@"&&$3!="*"{split($1,a,/[:|]/); if($3!="chr"a[2]||$4<a[3]-60||$4>a[3]+30) {print "@"$1;++n}}
		END{print "\nhap2 reads mismapped to hg19: "n"\n" >> "tmp/log.txt"}' \
		tmp/snp.reads.hap2hg19.sam >> tmp/hap2.multi.txt
	awk '{split($1,a,/[:|]/); print "chr"a[2]"\t"a[3]"\t"a[2]":"a[3]"\t"a[6]}' tmp/hap2.multi.txt | sort -k1,1 -k2,2n -k4,4n |\
		uniq > tmp/hap2.nonmap.snp.txt 
	
	# Combine non-mappable snps in both haplotypes
	cat tmp/hap1.nonmap.snp.txt tmp/hap2.nonmap.snp.txt | sort -k1,1 -k2,2n | uniq > tmp/nonmap.snp.txt
	cp tmp/nonmap.snp.txt nonmap/snps.$ID.nonmap.txt

	# Make bed files of non-mappable snps for each haplotype in individuals
	awk 'BEGIN{OFS="\t"}NR==FNR{split($4,a,"|"); pos[a[1]]=$3;next}
		{print "hap1\t"pos[$3]+$4"\t"pos[$3]+$4+29}' \
		hapbed/$ID/snp.TSSc.$ID.hap1.bed6 tmp/nonmap.snp.txt | sort -k2,2n -k3,3n > tmp/hap1.nonmap.snp.bed3
	bedtools merge -i tmp/hap1.nonmap.snp.bed3 > hapbed/$ID/snp.TSSc.$ID.nonmap.indiv.hap1.bed3
	awk 'BEGIN{OFS="\t"}NR==FNR{split($4,a,"|"); pos[a[1]]=$3;next}
		{print "hap2\t"pos[$3]+$4"\t"pos[$3]+$4+29}' \
		hapbed/$ID/snp.TSSc.$ID.hap2.bed6 tmp/nonmap.snp.txt | sort -k2,2n -k3,3n > tmp/hap2.nonmap.snp.bed3
	bedtools merge -i tmp/hap2.nonmap.snp.bed3 > hapbed/$ID/snp.TSSc.$ID.nonmap.indiv.hap2.bed3
	mv tmp/log.txt alignstat/log.$ID.txt
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
fi

if true; then
awk '' $caplist > tmp/nonmap.snp.all.txt

cat nonmap/snps.GM*.txt > tmp/nonmap.snp.all.txt
sort -k1,1 -k2,2n tmp/nonmap.snp.all.txt | uniq > nonmap/snps.all.nonmap.txt

while read FILE ID; do
	if [ ! -f nonmap/snps.$ID.nonmap.txt ]; then
		continue
	fi
	# Make bed files of non-mappable snps for all haplotypes
	awk 'BEGIN{OFS="\t"}NR==FNR{split($4,a,"|"); pos[a[1]]=$3;next}
		{print "hap1\t"pos[$3]+$4"\t"pos[$3]+$4+29}' \
		hapbed/$ID/snp.TSSc.$ID.hap1.bed6 nonmap/snps.all.nonmap.txt | sort -k2,2n -k3,3n \
		> hapbed/$ID/snp.TSSc.$ID.nonmap.all.hap1.bed3
	cat hapbed/$ID/snp.TSSc.$ID.nonmap.all.hap1.bed3 hapbed/$ID/ref.unmap.$ID.hap1.bed3 | \
		sort -k2,2n -k3,3n > tmp/unmap.bed3
	bedtools merge -i tmp/unmap.bed3 > hapbed/$ID/nonmap.$ID.hap1.bed3
	awk 'BEGIN{OFS="\t"}NR==FNR{split($4,a,"|"); pos[a[1]]=$3;next}
		{print "hap2\t"pos[$3]+$4"\t"pos[$3]+$4+29}' \
		hapbed/$ID/snp.TSSc.$ID.hap2.bed6 nonmap/snps.all.nonmap.txt | sort -k2,2n -k3,3n \
		> hapbed/$ID/snp.TSSc.$ID.nonmap.all.hap2.bed3
	cat hapbed/$ID/snp.TSSc.$ID.nonmap.all.hap2.bed3 hapbed/$ID/ref.unmap.$ID.hap2.bed3 | \
		sort -k2,2n -k3,3n > tmp/unmap.bed3
	bedtools merge -i tmp/unmap.bed3 > hapbed/$ID/nonmap.$ID.hap2.bed3
done < $caplist
fi
