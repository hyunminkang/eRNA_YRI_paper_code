# Generate mini-ref fasta of TSSc

###################################################################

# Code takes TSS region file (contains TRE center, readcounts for + and - strand TSS, TSS distance from center, and dominant TSS strand)
# Makes tiled 30mer of region and aligns reads to genome to find multimapping regions in the reference genome
# Takes each individual and creates individualized haplotypes and individualized TRE positions in haplotype coordinates
# Extracts tiling sequences for mappability test from each haplotype
# Time for demo: 19m3.259s

###################################################################

# File locations to fill in

region_list=data/procap.TSSc.bed6 # regions of interest as a bed file
genotype_file=genotype/genotype.imputed.vcf.gz # genotypes of all samples to be used
genome_fasta=hg19/hg19.fa
ID_list=list/idlist.min.txt # list of sample names
bowtie_genome=hg19/hg19

###################################################################

mkdir -p tmp
mkdir -p alignstat
mkdir -p asebam
mkdir -p aserc
mkdir -p bed
mkdir -p ebwt
mkdir -p fasta
mkdir -p fastq
mkdir -p hapbed
mkdir -p nonmap
mkdir -p rawdata
mkdir -p readcount
mkdir -p testbam #
mkdir -p tiQTL #

if true; then
# Merge overlapping regions
awk '{s=$2-250; if(s<1) s=1; print $1"\t"s"\t"$2+250}' $region_list > tmp/ref.TSSc.bed3
bedtools merge -i tmp/ref.TSSc.bed3 > tmp/ref.TSSc.merge.bed3
# Extract new TSSc positions
bedtools intersect -a $region_list -b tmp/ref.TSSc.merge.bed3 -wb |\
	awk '{print $1"\t"$2"\t"$7"\t"$8"\t"$2-$8}' > tmp/procap.TSSc.ref.pos.txt
# Make mini-ref padded fasta file
bedtools getfasta -fi $genome_fasta -bed tmp/ref.TSSc.merge.bed3 -tab -fo tmp/ref.txt
awk '{split($1,a,/[:\-]/); print a[1]"\t"a[2]"\t"a[3]"\t"0+s; s+=a[3]-a[2]+100}' tmp/ref.txt > tmp/ref.pos.txt
cut -f1-3 tmp/ref.pos.txt > tmp/ref.pos.bed3
awk '{printf $2; for(i=1;i<=100;++i) printf "N"; printf "\n"}' tmp/ref.txt > tmp/ref.fa
# TSSc positions on TSSc reference
awk 'NR==FNR{line[$1"\t"$2]=FNR;next}
	line[$3"\t"$4]{
		id=$3"\t"$4; print $1"\t"$2"\t"line[id]"\t"$5"\t"$5-250"\t"$5+250;}'\
	tmp/ref.pos.txt tmp/procap.TSSc.ref.pos.txt > tmp/procap.TSSc.rel.pos.txt
# Make tiled 30mer reads for reference mappability test
awk '{	l=length($1);
		for(i=1;i<=l-30;++i) {
			s=substr($1,i,30);
			if(substr(s,30,1)=="N") break;
			print "@"abpos+i"|"NR"|"i"\n"s"\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"; }
		abpos+=l; }' tmp/ref.fa > fastq/TSSc.reads.ref.fastq		
fi
if true; then
# Align tiled 30mer reads to hg19 genome and collect multimapped regions
bowtie -p6 -v2 -m1 -S --max tmp/ref.multi.fastq  $bowtie_genome fastq/TSSc.reads.ref.fastq \
	> tmp/ref.sam 2> tmp/log.txt
awk '{split($1,a,"|");print a[2]"\t"a[3]"\t"a[3]+29;getline;getline;getline}' tmp/ref.multi.fastq | \
	sort -k1,1n -k2,2n -k3,3n > tmp/ref.unmap.rel.pos.bed3
bedtools merge -i tmp/ref.unmap.rel.pos.bed3 > tmp/ref.unmap.rel.pos.merge.bed3
fi

if true; then
# Make snp position file
gzip -cd genotype/genotype.imputed.vcf.gz | head -n 1 > tmp/snps.TSSc.vcf
tabix $genotype_file -R tmp/ref.TSSc.merge.bed3 |\
	grep -vF '>' >> tmp/snps.TSSc.vcf
fi

if true; then
# Loop through all individuals
while read ID; do
	# get genotype bed file of the cellID
	awk -v id=$ID 'NR==1{id=substr(id,3,5); for(i=10;i<=NF;++i) if(substr($i,3,5)==id) col=i; if(col==0) exit 99;next}
		{print $1"\t"$2-1"\t"$2"\t"$4"|"$5"|"$col"|"$3"\t0\t+"}' tmp/snps.TSSc.vcf > tmp/SNP.bed
	if [ $? -eq 99 ]; then
		echo "No sample match in vcf"
		continue
	fi
	if [ ! -d hapbed/$ID ]; then
		mkdir -p hapbed/$ID
	fi

	bedtools intersect -a tmp/ref.pos.bed3 -b tmp/SNP.bed -wa -wb > tmp/OSNP.txt
	# Make relative snp position in merged fasta file
	# RSNP.txt : line pos ref alt hap1 hap2 snpName|genotype
	awk 'NR==FNR{pos[$1"\t"$2]=FNR;next}
		{	split($7,a,"|"); h1=a[a[3]+1];h2=a[a[4]+1];
			print pos[$1"\t"$2]"\t"$5-$2+1"\t"a[1]"\t"a[2]"\t"h1"\t"h2"\t"a[5]"|"a[3]"|"a[4]}' \
		tmp/ref.pos.bed3 tmp/OSNP.txt > tmp/RSNP.txt
		
	# replace padded fasta with phased haplotypes and output fasta file
	awk 'NR==FNR{hap1[FNR]=$1;offset1[FNR]=0;hap2[FNR]=$1;offset2[FNR]=0;n++;next}
		{	reflen=length($3);
			hap1[$1]=substr(hap1[$1],1,$2-1+offset1[$1]) $5 substr(hap1[$1],$2+reflen+offset1[$1]);
			if(length($5)!=reflen) hap1offset[$1]=hap1offset[$1]$2+offset1[$1]"|"offset1[$1]+length($5)-reflen"|";
			offset1[$1]=offset1[$1]+length($5)-reflen;
			hap2[$1]=substr(hap2[$1],1,$2-1+offset2[$1]) $6 substr(hap2[$1],$2+reflen+offset2[$1]);
			if(length($6)!=reflen) hap2offset[$1]=hap2offset[$1]$2+offset2[$1]"|"offset2[$1]+length($6)-reflen"|";
			offset2[$1]=offset2[$1]+length($6)-reflen; }
		END{
			print ">hap1" > "tmp/hap1.fa";
			for(i=1;i<=n;++i){
				print hap1[i] > "tmp/hap1.fa";
				print i"\t"hap1offset[i] > "tmp/hap1offset.txt";
			}
			print ">hap2" > "tmp/hap2.fa";
			for(i=1;i<=n;++i){
				print hap2[i] > "tmp/hap2.fa";
				print i"\t"hap2offset[i] > "tmp/hap2offset.txt";
			}
		}' \
		tmp/ref.fa tmp/RSNP.txt
	if true; then
	mv tmp/hap1.fa fasta/$ID.TSSc.hap1.fa
	mv tmp/hap2.fa fasta/$ID.TSSc.hap2.fa
	# New reference positions on modified haploid references
	paste tmp/ref.pos.bed3 <(awk 'NR>1{print 0+s;s+=length($1)}' fasta/$ID.TSSc.hap1.fa) > tmp/ref1.pos.txt
	paste tmp/ref.pos.bed3 <(awk 'NR>1{print 0+s;s+=length($1)}' fasta/$ID.TSSc.hap2.fa) > tmp/ref2.pos.txt

	# New TSSc coordinates on haploid references
	awk 'BEGIN{OFS="\t"}
		NR==FNR{offset[$1]=$2;next}
		{	n=(split(offset[$3],o,"|")-1)/2;
			for(i=1;i<=n;++i) {
				if($4>o[i*2-1]) $4+=o[i*2];
				if($5>o[i*2-1]) $5+=o[i*2];
				if($6>o[i*2-1]) $6+=o[i*2]; }
			print }' \
		tmp/hap1offset.txt tmp/procap.TSSc.rel.pos.txt > tmp/TSSc.rel.hap1.txt
	awk 'BEGIN{OFS="\t"}
		NR==FNR{offset[$1]=$2;next}
		{	n=(split(offset[$3],o,"|")-1)/2;
			for(i=1;i<=n;++i) {
				if($4>o[i*2-1]) $4+=o[i*2];
				if($5>o[i*2-1]) $5+=o[i*2];
				if($6>o[i*2-1]) $6+=o[i*2]; }
			print }' \
		tmp/hap2offset.txt tmp/procap.TSSc.rel.pos.txt > tmp/TSSc.rel.hap2.txt
	awk 'BEGIN{OFS="\t"}
		NR==FNR{ refpos[FNR]=$4; next }
		{	pos=refpos[$3];
			print "hap1\t"pos+$5"\t"pos+$4"\t"$1":"$2"\t0\t-" > "tmp/hap1.mn.bed6";
			print "hap1\t"pos+$4"\t"pos+$6"\t"$1":"$2"\t0\t+" > "tmp/hap1.pl.bed6"; }' \
		tmp/ref1.pos.txt tmp/TSSc.rel.hap1.txt
	awk 'BEGIN{OFS="\t"}
		NR==FNR{ refpos[FNR]=$4; next }
		{	pos=refpos[$3];
			print "hap2\t"pos+$5"\t"pos+$4"\t"$1":"$2"\t0\t-" > "tmp/hap2.mn.bed6";
			print "hap2\t"pos+$4"\t"pos+$6"\t"$1":"$2"\t0\t+" > "tmp/hap2.pl.bed6"; }' \
		tmp/ref2.pos.txt tmp/TSSc.rel.hap2.txt
	sort -k1,1 -k2,2n -k3,3n tmp/hap1.pl.bed6 > hapbed/$ID/TSSc.$ID.hap1.pl.bed6
	sort -k1,1 -k2,2n -k3,3n tmp/hap1.mn.bed6 > hapbed/$ID/TSSc.$ID.hap1.mn.bed6
	sort -k1,1 -k2,2n -k3,3n tmp/hap2.pl.bed6 > hapbed/$ID/TSSc.$ID.hap2.pl.bed6
	sort -k1,1 -k2,2n -k3,3n tmp/hap2.mn.bed6 > hapbed/$ID/TSSc.$ID.hap2.mn.bed6
	# New SNP coordinates	
	awk '{	reflen=length($3);
			haplen1=length($5);
			haplen2=length($6);
			print $1"\t"$2+os1[$1]"\t"$2+os2[$1]"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7;
			os1[$1]+=haplen1-reflen;
			os2[$1]+=haplen2-reflen;
		}' tmp/RSNP.txt > tmp/NSNP.txt
	awk 'BEGIN{OFS="\t"}
		NR==FNR{ refpos[FNR]=$4; next }
		{	pos=refpos[$1]
			print "hap1\t"pos+$2-1"\t"pos+$2"\t"$8"\t0\t+" > "tmp/hap1.snp.bed6"; }' \
		tmp/ref1.pos.txt tmp/NSNP.txt
	awk 'BEGIN{OFS="\t"}
		NR==FNR{ refpos[FNR]=$4; next }
		{	pos=refpos[$1]
			print "hap2\t"pos+$3-1"\t"pos+$3"\t"$8"\t0\t+" > "tmp/hap2.snp.bed6"; }' \
		tmp/ref2.pos.txt tmp/NSNP.txt
	sort -k1,1 -k2,2n -k3,3n tmp/hap1.snp.bed6 > hapbed/$ID/snp.TSSc.$ID.hap1.bed6	
	sort -k1,1 -k2,2n -k3,3n tmp/hap2.snp.bed6 > hapbed/$ID/snp.TSSc.$ID.hap2.bed6
	fi
	# New reference nonmappable region coordinates
	awk 'BEGIN{OFS="\t"}
		NR==FNR{offset[$1]=$2;next}
		{	n=(split(offset[$1],o,"|")-1)/2;
			for(i=1;i<=n;++i) {
				if($2>o[i*2-1]) $2+=o[i*2];
				if($3>o[i*2-1]) $3+=o[i*2]; }
			if($3>$2) print }' \
		tmp/hap1offset.txt tmp/ref.unmap.rel.pos.merge.bed3 > tmp/ref.rel.unmap.hap1.bed3
	awk 'BEGIN{OFS="\t"}
		NR==FNR{ refpos[FNR]=$4; next }
		{	pos=refpos[$1];
			print "hap1\t"pos+$2"\t"pos+$3}' \
		tmp/ref1.pos.txt tmp/ref.rel.unmap.hap1.bed3 | sort -k2,2n -k3,3n > tmp/ref.unmap.hap1.bed3
	bedtools merge -i tmp/ref.unmap.hap1.bed3 > hapbed/$ID/ref.unmap.$ID.hap1.bed3
	awk 'BEGIN{OFS="\t"}
		NR==FNR{offset[$1]=$2;next}
		{	n=(split(offset[$1],o,"|")-1)/2;
			for(i=1;i<=n;++i) {
				if($2>o[i*2-1]) $2+=o[i*2];
				if($3>o[i*2-1]) $3+=o[i*2]; }
			if($3>$2) print }' \
		tmp/hap2offset.txt tmp/ref.unmap.rel.pos.merge.bed3 > tmp/ref.rel.unmap.hap2.bed3
	awk 'BEGIN{OFS="\t"}
		NR==FNR{ refpos[FNR]=$4; next }
		{	pos=refpos[$1];
			print "hap2\t"pos+$2"\t"pos+$3}' \
		tmp/ref2.pos.txt tmp/ref.rel.unmap.hap2.bed3 | sort -k2,2n -k3,3n > tmp/ref.unmap.hap2.bed3
	bedtools merge -i tmp/ref.unmap.hap2.bed3 > hapbed/$ID/ref.unmap.$ID.hap2.bed3

	if true; then
	# Make tiling sequences at the SNPs for mappability test
	awk 'NR==FNR{
			if(substr($1,1,1)==">") { hap++; l=1; }
			else { seq[hap"\t"l]=$1; ++l; } }
		$6!=$7{
			start=$2-29;
			if(start<1) start=1;
			for(i=start;i<=$2+length($6)-1;++i) {
				s=substr(seq["1\t"$1],i,30);
				if(length(s)==30&&substr(s,length(s),1)!="N") print "@hap1|"$8"|"i-$2"\n"s"\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"; }
		}' fasta/$ID.TSSc.hap1.fa tmp/NSNP.txt > tmp/reads.hap1.tmp
	awk 'NR==FNR{
			if(substr($1,1,1)==">") { hap++; l=1; }
			else { seq[hap"\t"l]=$1; ++l; } }
		$6!=$7{
			start=$3-29;
			if(start<1) start=1;
			for(i=start;i<=$3+length($7)-1;++i) {
				s=substr(seq["1\t"$1],i,30);
				if(length(s)==30&&substr(s,length(s),1)!="N") print "@hap2|"$8"|"i-$3"\n"s"\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"; }
		}' fasta/$ID.TSSc.hap2.fa tmp/NSNP.txt > tmp/reads.hap2.tmp
	mv tmp/reads.hap1.tmp fastq/snp.TSSc.reads.$ID.hap1.fastq	
	mv tmp/reads.hap2.tmp fastq/snp.TSSc.reads.$ID.hap2.fastq	
	fi
done < $ID_list

fi
