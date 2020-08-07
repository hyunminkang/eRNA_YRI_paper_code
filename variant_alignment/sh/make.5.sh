#################################################################

# Merges all read counts together

# Time for demo: 0m2.299s

#################################################################
# Variables to assign

ID_list=list/idlist.min.txt

#################################################################

if true; then
awk '{split($1,a,"|"); print a[1]"\t0\t0\t0\t0"}' aserc/GM18520.ase.txt > tmp/merge.txt

while read ID; do
	awk 'NR==FNR{snp[FNR]=$1; rp[FNR]=$2; ap[FNR]=$3; rm[FNR]=$4; am[FNR]=$5; next} {split($1,a,"|")}
		{	g=a[2]"|"a[3];
			if(g=="1|0") {rp[FNR]+=$2; ap[FNR]+=$4; rm[FNR]+=$3; am[FNR]+=$5}
			else if(g=="0|1") {rp[FNR]+=$4; ap[FNR]+=$2; rm[FNR]+=$5; am[FNR]+=$3}
			print snp[FNR]"\t"rp[FNR]"\t"ap[FNR]"\t"rm[FNR]"\t"am[FNR];
		}' tmp/merge.txt aserc/$ID.ase.txt > tmp/m2.txt
	mv tmp/m2.txt tmp/merge.txt
done < $ID_list
mv tmp/merge.txt aserc/all.ase.txt
fi
awk '{print $1"\t"$2"\t"$3"\t"$4"\t0\t0\t0\t0"}' aserc/GM18520.asTSSr.txt > tmp/merge.txt

while read ID; do
	awk 'NR==FNR{rp[FNR]=$5; ap[FNR]=$6; rm[FNR]=$7; am[FNR]=$8; next}
		{	g=$5;
			if(g=="1|0") {rp[FNR]+=$6; ap[FNR]+=$8; rm[FNR]+=$7; am[FNR]+=$9}
			else if(g=="0|1") {rp[FNR]+=$8; ap[FNR]+=$6; rm[FNR]+=$9; am[FNR]+=$7}
			print $1"\t"$2"\t"$3"\t"$4"\t"rp[FNR]"\t"ap[FNR]"\t"rm[FNR]"\t"am[FNR];
		}' tmp/merge.txt aserc/$ID.asTSSr.txt > tmp/m2.txt
	mv tmp/m2.txt tmp/merge.txt
done < $ID_list
mv tmp/merge.txt aserc/all.aseTSSr.txt
