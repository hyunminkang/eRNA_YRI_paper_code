#################################################################

# Combines all readcounts from all replicates together

# Time for demo: 0m1.414s

#################################################################
# Variables to assign

rep_caplist=list/caplist.min.rep.txt
caplist_all=list/caplist.min.all.txt
ID_list=list/idlist.min.txt

#################################################################

mkdir -p bed/

# Readcount for all regions
paste readcount/procap.readcount.txt <(cut -f3- readcount/procap.rep.readcount.txt) > readcount/procap.all.readcount.txt
cat $ID_list <(cut -f2 $rep_caplist) > $caplist_all

# average RPM for all nTSS
awk '{s=0;for(i=3;i<=NF;++i) s+=$i;print $1"\t"$2"\t"s/1315.892015}' readcount/procap.all.readcount.txt > tmp/procap.sum.rpm.txt
awk '$3>0.5{print $1"\t"$2"\t"$2+1}' tmp/procap.sum.rpm.txt > bed/procap.ambr.bed3

# Readcounts for all regions with ambr
awk 'NR==FNR{a[$1":"$2]=1;next}a[$1":"$2]{print}' bed/procap.ambr.bed3 readcount/procap.all.readcount.txt > readcount/procap.all.ambr.readcount.txt

# Makes the normalization for readcounts
sh/normalize.sh readcount/procap.all.ambr.readcount.txt list/rpm.all.txt

