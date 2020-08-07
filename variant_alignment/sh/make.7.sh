#################################################################

# Finds haplotype specific read counts using allele specific reads
# Same as 4 but processes replicates

# Time for demo: 0m10.185s

#################################################################
# Variables to assign

rep_caplist=list/caplist.min.rep.txt
region_list=data/procap.TSSc.bed6

#################################################################

STARTTIME=$(date +%s)
DONE=0
REMAIN=$(cat $rep_caplist | wc -l)
while read FILE SID; do
	ID=${SID:0:7}
	echo $SID
	if true; then
	bedtools coverage -a hapbed/$ID/snp.TSSc.$ID.hap1.bed6 -b asebam/$SID.ase.hap1.bam -sorted -s > tmp/RC1P.txt
	bedtools coverage -a hapbed/$ID/snp.TSSc.$ID.hap1.bed6 -b asebam/$SID.ase.hap1.bam -sorted -S > tmp/RC1M.txt
	bedtools coverage -a hapbed/$ID/snp.TSSc.$ID.hap2.bed6 -b asebam/$SID.ase.hap2.bam -sorted -s > tmp/RC2P.txt
	bedtools coverage -a hapbed/$ID/snp.TSSc.$ID.hap2.bed6 -b asebam/$SID.ase.hap2.bam -sorted -S > tmp/RC2M.txt
		
	paste <(cut -f4,7 tmp/RC1P.txt) <(cut -f7 tmp/RC1M.txt) <(cut -f7 tmp/RC2P.txt) <(cut -f7 tmp/RC2M.txt) > tmp/RC.txt
	awk 'BEGIN{OFS="\t"}{split($1,a,"|"); if(a[2]==a[3]) print $1"\t0\t0\t0\t0"; else print}' tmp/RC.txt > aserc/$SID.ase.txt
	awk '{split($1,a,/[:|]/); print "chr"a[1]"\t"a[2]-1"\t"a[2]"\t"$1"\t"$2":"$3":"$4":"$5"\t+"}' aserc/$SID.ase.txt > tmp/snp.bed
	bedtools window -a $region_list -b tmp/snp.bed -w 250 > tmp/snprc.TSSc.txt
	awk '{split($10,a,"|"); split($11,b,":"); print $1"\t"$2"\t"a[1]"\t"$9-$2"\t"a[2]"|"a[3]"\t"b[1]"\t"b[2]"\t"b[3]"\t"b[4]}' \
		tmp/snprc.TSSc.txt > tmp/snprc.TSSc.rel.txt
	bedtools coverage -a hapbed/$ID/TSSc.$ID.hap1.pl.bed6 -b asebam/$SID.ase.hap1.bam -sorted -s > tmp/RC1P.txt
	bedtools coverage -a hapbed/$ID/TSSc.$ID.hap1.mn.bed6 -b asebam/$SID.ase.hap1.bam -sorted -s > tmp/RC1M.txt
	bedtools coverage -a hapbed/$ID/TSSc.$ID.hap2.pl.bed6 -b asebam/$SID.ase.hap2.bam -sorted -s > tmp/RC2P.txt
	bedtools coverage -a hapbed/$ID/TSSc.$ID.hap2.mn.bed6 -b asebam/$SID.ase.hap2.bam -sorted -s > tmp/RC2M.txt
	
	paste <(cut -f4,7 tmp/RC1P.txt) <(cut -f7 tmp/RC1M.txt) <(cut -f7 tmp/RC2P.txt) <(cut -f7 tmp/RC2M.txt) > tmp/RC.txt
	awk '{split($1,a,":"); print a[1]"\t"a[2]"\t"$2"\t"$3"\t"$4"\t"$5}' tmp/RC.txt > readcount/$SID.ast.txt 
	awk 'FNR==1{++fid}fid==1{hap1[FNR]=$3;next}fid==2{hap2[FNR]=$3;next}{print $1"\t"$2"\t"int((hap1[FNR]+hap2[FNR]+$3+$4+$5+$6)/2)}' \
		readcount/$SID.hap1.txt readcount/$SID.hap2.txt readcount/$SID.ast.txt > readcount/"${SID}".txt
	awk 'FNR==1{++fid}fid==1{hap1[FNR]=$4;next}fid==2{hap2[FNR]=$4;next}{print $1"\t"$2"\t"int((hap1[FNR]+hap2[FNR]+$3+$5)/2)}' \
		readcount/$SID.hap1.txt readcount/$SID.hap2.txt readcount/$SID.ast.txt > readcount/$SID.pl.txt
	awk 'FNR==1{++fid}fid==1{hap1[FNR]=$5;next}fid==2{hap2[FNR]=$5;next}{print $1"\t"$2"\t"int((hap1[FNR]+hap2[FNR]+$4+$6)/2)}' \
		readcount/$SID.hap1.txt readcount/$SID.hap2.txt readcount/$SID.ast.txt > readcount/$SID.mn.txt
	awk 'NR==FNR{rc[$1]=$2"\t"$3"\t"$4"\t"$5;next}{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"rc[$1":"$2]}' \
		tmp/RC.txt tmp/snprc.TSSc.rel.txt > aserc/$SID.asTSSr.txt
	fi
	if [ $DONE -eq 0 ]; then
		cut -f1,2 readcount/$SID.txt > readcount/procap.rep.readcount.txt
		cut -f1,2 readcount/$SID.pl.txt > readcount/procap.rep.readcount.pl.txt
		cut -f1,2 readcount/$SID.mn.txt > readcount/procap.rep.readcount.mn.txt
	fi
	paste readcount/procap.rep.readcount.txt <(cut -f3 readcount/"${SID}".txt) > tmp/RC.txt
	mv tmp/RC.txt readcount/procap.rep.readcount.txt
	paste readcount/procap.rep.readcount.pl.txt <(cut -f3 readcount/"${SID}".pl.txt) > tmp/RC.txt
	mv tmp/RC.txt readcount/procap.rep.readcount.pl.txt
	paste readcount/procap.rep.readcount.mn.txt <(cut -f3 readcount/"${SID}".mn.txt) > tmp/RC.txt
	mv tmp/RC.txt readcount/procap.rep.readcount.mn.txt
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
