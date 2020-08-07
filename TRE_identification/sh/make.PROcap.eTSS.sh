####################################################################################################

# This script generates TREs from bedgraph files and generates background regions with read counts,
# then calls an R script to calculate p-values for the TREs and indicate the P-value and FDR of an
# 0.5 RPM cut off

# Demo time: 38m24.851s

####################################################################################################
# Variables

bedgraph_pl=data/procap.pl.sorted.bg.gz
bedgraph_mn=data/procap.mn.sorted.bg.gz
chrom_sizes=data/hg19.chrom.sizes
total_reads=1378818883
reference_annotation=data/refFlat-hg19.bed

####################################################################################################

mkdir -p bed/
mkdir -p bedgraph/
mkdir -p readsumfile/
mkdir -p tmp/
mkdir -p results/

echo "Generating 100 bp bins"

# generate 100 bp bins of the procap bedgraph files
awk 'NR==1{prevChr=$1;prevStart=int($2/100);s=$4;next}
($1!=prevChr||int($2/100)!=prevStart){print prevChr"\t"prevStart*100"\t"prevStart*100+100"\t"s;prevChr=$1;prevStart=int($2/100);s=$4;next}
{s+=$4}' $bedgraph_pl > bedgraph/procap.100bp.pl.bedgraph

awk 'NR==1{prevChr=$1;prevStart=int($2/100);s=$4;next}
($1!=prevChr||int($2/100)!=prevStart){print prevChr"\t"prevStart*100"\t"prevStart*100+100"\t"s;prevChr=$1;prevStart=int($2/100);s=$4;next}
{s+=$4}' $bedgraph_mn > bedgraph/procap.100bp.mn.bedgraph

echo "Identifying top expressed bins"

# Select top 5% of the genomic regions with highest plus and minus read counts

LC_ALL="C" sort -k4,4nr bedgraph/procap.100bp.pl.bedgraph | head -1500000 | LC_ALL="C" sort -k1,1 -k2,2n | cut -f1-3 > bed/procap.top5.pl.bed
LC_ALL="C" sort -k4,4n bedgraph/procap.100bp.mn.bedgraph | head -1500000 | LC_ALL="C" sort -k1,1 -k2,2n | cut -f1-3 > bed/procap.top5.mn.bed

echo "Extending bins"

# Filter PRO-cap reads from top5 bins +- 100 bp
awk '{if($2<100) $2=100; print $1"\t"$2-100"\t"$3+100}' bed/procap.top5.pl.bed > tmp/procap.top5ext.pl.tmp
awk '{if($2<100) $2=100; print $1"\t"$2-100"\t"$3+100}' bed/procap.top5.mn.bed > tmp/procap.top5ext.mn.tmp

bedtools intersect -sorted -wb -a $bedgraph_pl -b tmp/procap.top5ext.pl.tmp > tmp/procap.top5ext.pl.tmp
bedtools intersect -sorted -wb -a $bedgraph_mn -b tmp/procap.top5ext.mn.tmp > tmp/procap.top5ext.mn.tmp

echo "Finding maximum count position in bins"

# Find the maximum of the 300 bp regions by sorting by read counts (filter by 5 reads)
LC_ALL="C" sort -k5,5 -k6,6n -k4,4nr tmp/procap.top5ext.pl.tmp | LC_ALL="C" sort -u -k5,5 -k6,6n | awk '{if($4>=5) print $1"\t"$2"\t"$3"\t"$4}' | LC_ALL="C" sort -u -k1,1 -k2,2n > tmp/procap.top5peak.pl.tmp
LC_ALL="C" sort -k5,5 -k6,6n -k4,4n tmp/procap.top5ext.mn.tmp | LC_ALL="C" sort -u -k5,5 -k6,6n | awk '{if($4<=-5) print $1"\t"$2"\t"$3"\t"$4}' | LC_ALL="C" sort -u -k1,1 -k2,2n > tmp/procap.top5peak.mn.tmp

echo "Finding upstream peak"

# Make upstream region from -300 to 0 from peaks
awk '$2>300{print $1"\t"$2-300"\t"$2"\t"$4}' tmp/procap.top5peak.pl.tmp > tmp/procap.top5peak.pl.tmp
awk '{print $1"\t"$3"\t"$3+300"\t"$4}' tmp/procap.top5peak.mn.tmp > tmp/procap.top5peak.mn.tmp

bedtools intersect -sorted -wb -a $bedgraph_mn -b tmp/procap.top5peak.pl.tmp > tmp/procap.top5as.pl.tmp
bedtools intersect -sorted -wb -a $bedgraph_pl -b tmp/procap.top5peak.mn.tmp > tmp/procap.top5as.mn.tmp


# Find the peak positions from opposite strand
LC_ALL="C" sort -k5,5 -k6,6n -k4,4n tmp/procap.top5as.pl.tmp | sort -u -k5,5 -k6,6n | awk '{if($4<=-5) print $1"\t"$2"\t"$3"\t"$4","$8"\t"$7"\t+"}' > tmp/procap.top5aspeak.pl.tmp
LC_ALL="C" sort -k5,5 -k6,6n -k4,4nr tmp/procap.top5as.mn.tmp | sort -u -k5,5 -k6,6n | awk '{if($4>=5) print $1"\t"$2"\t"$3"\t"$4","$8"\t"$6"\t-"}' > tmp/procap.top5aspeak.mn.tmp
cat tmp/procap.top5aspeak.pl.tmp tmp/procap.top5aspeak.mn.tmp | sort -u -k1,1 -k2,2n -k5,5n > tmp/procap.top5aspeak.tmp
awk '{if($2<$5) {a=$2;b=$5} else {a=$5;b=$2} pos=int((b+a)/2); d=int((b-a)/2); print $1"\t"pos"\t"pos+1"\t"$4"\t"d"\t"$6"\t"int(a/2)"\t"int(b/2)"\t"int((a+1)/2)"\t"int((b+1)/2)}' tmp/procap.top5aspeak.tmp | sort -u -k1,1 -k7,7n -k8,8n | sort -u -k1,1 -k9,9n -k10,10n | awk -F'[\t,]' '{if($4<0) {a=$5;b=$4} else {a=$4;b=$5} if(a+b>=0) $7="+"; else $7="-"; print $1"\t"$2"\t"$3"\t"a","b"\t"$6"\t"$7}' > bed/procap.cenPeak.bed6

echo "Getting TRE read counts"


# make tss positions for read sumA
awk '{st=$2; up=$2-200; dn=$2+200; if(up<1) up=1; pl=$1"\t"st"\t"dn"\t"$4"\t"st"\t"$6; mn=$1"\t"up"\t"st"\t"$4"\t"st"\t"$6; if($6=="+") {print pl"\tss">"tmp/tss.pl.tmp";print mn"\tas">"tmp/tss.mn.tmp"}else{print pl"\tas">"tmp/tss.pl.tmp";print mn"\tss">"tmp/tss.mn.tmp"}}' bed/procap.cenPeak.bed6
LC_ALL="C" sort -k1,1 -k2,2n tmp/tss.pl.tmp > tmp/tss.pl.sort.tmp
LC_ALL="C" sort -k1,1 -k2,2n tmp/tss.mn.tmp > tmp/tss.mn.sort.tmp

# get read counts
bedtools map -a tmp/tss.pl.sort.tmp -b $bedgraph_pl -c 4 -o sum -sorted > tmp/procap.pl.tmp
bedtools map -a tmp/tss.mn.sort.tmp -b $bedgraph_mn -c 4 -o sum -sorted > tmp/procap.mn.tmp
cat readsumfile/procap.pl.tmp readsumfile/procap.mn.tmp | awk '{if($8<0) $8=-$8; if($7=="ss") print $1"\t"$5"\t"$4"\t"$8>"tmp/procap.ss.tmp";else print $1"\t"$5"\t"$4"\t"$8>"tmp/procap.as.tmp"}'
LC_ALL="C" sort -k1,1 -k2,2n tmp/procap.ss.tmp > readsumfile/procap.ss.txt
LC_ALL="C" sort -k1,1 -k2,2n tmp/procap.as.tmp > readsumfile/procap.as.txt


awk 'NR==FNR{if($4>=1417/2) pos[$1$2]=1;next}pos[$1$2]{print}' readsumfile/procap.ss.txt bed/procap.cenPeak.bed6 > bed/procap.ati.bed6

awk '{print $1"\t"$2"\t"$4}' readsumfile/procap.ss.txt > bed/procap.ati.ss.txt
awk '{print $1"\t"$2"\t"$4}' readsumfile/procap.as.txt > bed/procap.ati.as.txt

awk '{for(i=3;i<=NF;++i) s[i-2]=$i; getline < "bed/procap.ati.as.txt"; printf $1"\t"$2; for(i=3;i<=NF;++i) printf "\t"s[i-2]+$i; printf "\n"}' bed/procap.ati.ss.txt > bed/procap.ati.txt

awk -v rc="$total_reads" '{print $1"\t"$2"\t"$3/(rc/1000000)}' bed/procap.ati.txt > bed/procap.rpm.txt

echo "Finding local maxima"

# select local maximum within +-150 bp
bin=150
awk -v bin="${bin}" 'BEGIN{OFS="\t"}NR==1{prevMaxChr=$1;prevMaxPos=$2;prevMaxCount=$3;next}{if($1!=prevMaxChr) {print prevMaxChr"\t"prevMaxPos;prevMaxChr=$1;prevMaxPos=$2;prevMaxCount=$3;} else if($2>prevMaxPos+bin) {print prevMaxChr"\t"prevMaxPos;prevMaxChr=$1;prevMaxPos=$2;prevMaxCount=$3;} else if($3>prevMaxCount) {prevMaxPos=$2;prevMaxCount=$3}}' bed/procap.ati.txt > tmp/procap.ati.lmf.tmp
sort -k1,1 -k2,2nr bed/procap.ati.txt | awk -v bin="${bin}" 'BEGIN{OFS="\t"}NR==1{prevMaxChr=$1;prevMaxPos=$2;prevMaxCount=$3;next}{if($1!=prevMaxChr) {print prevMaxChr"\t"prevMaxPos;prevMaxChr=$1;prevMaxPos=$2;prevMaxCount=$3;} else if($2<prevMaxPos-bin) {print prevMaxChr"\t"prevMaxPos;prevMaxChr=$1;prevMaxPos=$2;prevMaxCount=$3;} else if($3>prevMaxCount) {prevMaxPos=$2;prevMaxCount=$3}}' | sort -k1,1 -k2,2n > tmp/procap.ati.lmr.tmp
awk 'NR==FNR{a[$0]=1;next}a[$0]{print}'  tmp/procap.ati.lmf.tmp tmp/procap.ati.lmr.tmp | sort -k1,1 -k2,2n > tmp/procap.ati.lm.tmp

echo "Filtering by PRM"

# filter RPM data by local maxima
awk 'BEGIN{OFS="\t"}NR==FNR{pos[$1$2]=1;next}pos[$1$2]{gsub(/[ \t]+/,"\t",$0);print $0}' tmp/procap.ati.lm.tmp bed/procap.rpm.txt > bed/procap.rpm.lm.txt



awk 'BEGIN{OFS="\t"}NR==FNR{pos[$1$2]=1;next}pos[$1$2]{gsub(/[ \t]+/,"\t",$0);print $0}' tmp/procap.ati.lm.tmp bed/procap.ati.bed6 > results/procap.atilm.bed6


echo "Making null distribution"
# making the null distribution

awk '{if ($6="+"&&$2-1000>0) {print $1"\t"$2-1000"\t"$2+1000} else if ($6="-"&&$3-1000>0) {print $1"\t"$3-1000"\t"$3+1000}}' $reference_annotation > bed/ref_TSS.bed


bedtools random -l 1 -n 2000000 -g $chrom_sizes | awk '{if ($2-200>0) print;}' > tmp/background.tmp
bedtools intersect -v -a tmp/background.tmp -b bed/ref_TSS.bed | shuf -n 1000000 | sort -k1,1 -k2,2n > bed/background.bed

awk '{st=$2; up=$2-200; dn=$2+200; if(up<=0) up=1; pl=$1"\t"st"\t"dn"\t"$4"\t"st"\t"$6; mn=$1"\t"up"\t"st"\t"$4"\t"st"\t"$6; if($6=="+") {print pl"\tss">"tmp/background.pos.pl.tmp";print mn"\tas">"tmp/background.pos.mn.tmp"}else{print pl"\tas">"tmp/background.pos.pl.tmp";print mn"\tss">"tmp/background.pos.mn.tmp"}}' bed/background.bed


bedtools map -o sum -c 4 -a tmp/background.pos.pl.tmp -b $bedgraph_pl > bedgraph/background_counts.pl.bed
bedtools map -o sum -c 4 -a tmp/background.pos.mn.tmp -b $bedgraph_mn > bedgraph/background_counts.mn.bed

cat bedgraph/background_counts.pl.bed bedgraph/background_counts.mn.bed | awk '{if($8<0) $8=-$8; if($7=="ss") print $1"\t"$5"\t"$4"\t"$8>"tmp/background.ss.tmp";else print $1"\t"$5"\t"$4"\t"$8>"tmp/background.as.tmp"}'
LC_ALL="C" sort -k1,1 -k2,2n tmp/background.ss.tmp > bedgraph/background.ss.txt
LC_ALL="C" sort -k1,1 -k2,2n tmp/background.as.tmp > bedgraph/background.as.txt

awk '{print $1"\t"$2"\t"$4}' bedgraph/background.ss.txt > bed/background.ati.ss.txt
awk '{print $1"\t"$2"\t"$4}' bedgraph/background.as.txt > bed/background.ati.as.txt

awk '{for(i=3;i<=NF;++i) s[i-2]=$i; getline < "bed/background.ati.as.txt"; printf $1"\t"$2; for(i=3;i<=NF;++i) printf "\t"s[i-2]+$i; printf "\n"}' bed/background.ati.ss.txt > bed/background.ati.txt

awk -v rc="$total_reads" '{print $1"\t"$2"\t"$3/(rc/1000000)}' bed/background.ati.txt > bed/background_rpm.txt

echo "Filtering by p-value"

Rscript --vanilla rscript/pval_filter.R

rm tmp/*.tmp
