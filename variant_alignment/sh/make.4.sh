#################################################################

# Finds haplotype specific read counts using allele specific reads

# Time for demo: 0m12.239s

#################################################################
# Variables to assign

ID_list=list/idlist.min.txt
region_list=data/procap.TSSc.bed6

#################################################################

STARTTIME=$(date +%s)
DONE=0
REMAIN=$(cat $ID_list | wc -l)
while read ID; do
	if true; then
	bedtools coverage -a hapbed/$ID/snp.TSSc.$ID.hap1.bed6 -b asebam/$ID.ase.hap1.bam -sorted -s > tmp/RC1P.txt
	bedtools coverage -a hapbed/$ID/snp.TSSc.$ID.hap1.bed6 -b asebam/$ID.ase.hap1.bam -sorted -S > tmp/RC1M.txt
	bedtools coverage -a hapbed/$ID/snp.TSSc.$ID.hap2.bed6 -b asebam/$ID.ase.hap2.bam -sorted -s > tmp/RC2P.txt
	bedtools coverage -a hapbed/$ID/snp.TSSc.$ID.hap2.bed6 -b asebam/$ID.ase.hap2.bam -sorted -S > tmp/RC2M.txt
		
	paste <(cut -f4,7 tmp/RC1P.txt) <(cut -f7 tmp/RC1M.txt) <(cut -f7 tmp/RC2P.txt) <(cut -f7 tmp/RC2M.txt) > tmp/RC.txt
	awk 'BEGIN{OFS="\t"}{split($1,a,"|"); if(a[2]==a[3]) print $1"\t0\t0\t0\t0"; else print}' tmp/RC.txt > aserc/$ID.ase.txt
	awk '{split($1,a,/[:|]/); print "chr"a[1]"\t"a[2]-1"\t"a[2]"\t"$1"\t"$2":"$3":"$4":"$5"\t+"}' aserc/$ID.ase.txt > tmp/snp.bed
	bedtools window -a $region_list -b tmp/snp.bed -w 250 > tmp/snprc.TSSc.txt
	awk '{split($10,a,"|"); split($11,b,":"); print $1"\t"$2"\t"a[1]"\t"$9-$2"\t"a[2]"|"a[3]"\t"b[1]"\t"b[2]"\t"b[3]"\t"b[4]}' \
		tmp/snprc.TSSc.txt > tmp/snprc.TSSc.rel.txt
	bedtools coverage -a hapbed/$ID/TSSc.$ID.hap1.pl.bed6 -b asebam/$ID.ase.hap1.bam -sorted -s > tmp/RC1P.txt
	bedtools coverage -a hapbed/$ID/TSSc.$ID.hap1.mn.bed6 -b asebam/$ID.ase.hap1.bam -sorted -s > tmp/RC1M.txt
	bedtools coverage -a hapbed/$ID/TSSc.$ID.hap2.pl.bed6 -b asebam/$ID.ase.hap2.bam -sorted -s > tmp/RC2P.txt
	bedtools coverage -a hapbed/$ID/TSSc.$ID.hap2.mn.bed6 -b asebam/$ID.ase.hap2.bam -sorted -s > tmp/RC2M.txt
	
	paste <(cut -f4,7 tmp/RC1P.txt) <(cut -f7 tmp/RC1M.txt) <(cut -f7 tmp/RC2P.txt) <(cut -f7 tmp/RC2M.txt) > tmp/RC.txt
	awk '{split($1,a,":"); print a[1]"\t"a[2]"\t"$2"\t"$3"\t"$4"\t"$5}' tmp/RC.txt > readcount/$ID.ast.txt 
	fi
	awk 'FNR==1{++fid}fid==1{hap1[FNR]=$3;next}fid==2{hap2[FNR]=$3;next}{print $1"\t"$2"\t"int((hap1[FNR]+hap2[FNR]+$3+$4+$5+$6)/2)}' \
		readcount/$ID.hap1.txt readcount/$ID.hap2.txt readcount/$ID.ast.txt > readcount/$ID.txt
	awk 'FNR==1{++fid}fid==1{hap1[FNR]=$4;next}fid==2{hap2[FNR]=$4;next}{print $1"\t"$2"\t"int((hap1[FNR]+hap2[FNR]+$3+$5)/2)}' \
		readcount/$ID.hap1.txt readcount/$ID.hap2.txt readcount/$ID.ast.txt > readcount/$ID.pl.txt
	awk 'FNR==1{++fid}fid==1{hap1[FNR]=$5;next}fid==2{hap2[FNR]=$5;next}{print $1"\t"$2"\t"int((hap1[FNR]+hap2[FNR]+$4+$6)/2)}' \
		readcount/$ID.hap1.txt readcount/$ID.hap2.txt readcount/$ID.ast.txt > readcount/$ID.mn.txt
	awk 'NR==FNR{rc[$1]=$2"\t"$3"\t"$4"\t"$5;next}{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"rc[$1":"$2]}' \
		tmp/RC.txt tmp/snprc.TSSc.rel.txt > aserc/$ID.asTSSr.txt
	paste <(cut -f1-3 readcount/$ID.pl.txt) <(cut -f3 readcount/$ID.mn.txt) | \
		awk '$3+$4==0{print $1"\t"$2"\t0";next}{print $1"\t"$2"\t"$3/($3+$4)}' > tmp/RR.txt

	if [ $DONE -eq 0 ]; then
		cut -f1,2 readcount/$ID.txt > readcount/procap.readcount.txt
		cut -f1,2 readcount/$ID.txt > readcount/procap.readratio.txt
	fi
	paste readcount/procap.readcount.txt <(cut -f3 readcount/$ID.txt) > tmp/RC.txt
	mv tmp/RC.txt readcount/procap.readcount.txt

	paste readcount/procap.readratio.txt <(cut -f3 tmp/RR.txt) > tmp/RRC.txt
	mv tmp/RRC.txt readcount/procap.readratio.txt
	CURTIME=$(date +%s)
	ELAPSED=$(($CURTIME-$STARTTIME))
	DONE=$(($DONE+1))
	REMAIN=$(($REMAIN-1))
	ESTIMATED=$(($ELAPSED*$REMAIN/$DONE))
	HR=$((ESTIMATED/3600))
	MIN=$(($ESTIMATED/60-$HR*60))
	SEC=$(($ESTIMATED%60))
	echo Estimated $HR:$MIN:$SEC left
done < $ID_list
