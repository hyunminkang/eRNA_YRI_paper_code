# Input file example
# First line should be column identifiers
#	chr		pos		reads
#	OR
#	chr		name	reads

# make Read per million normalization
# if total read count file is supplied as the second argument, use the total read count
# otherwise use the read count sum
if [ $# -gt 1 ]; then
	awk 'NR==FNR{rc[NR]=$1;next}
		FNR==1{print;next}
		{printf $1"\t"$2; for(i=3;i<=NF;++i) printf "\t%.6f",$i/rc[i-2]*1000000; printf "\n"}' \
		$2 $1 > rpm.tmp
else
	awk 'NR==FNR{ if(FNR>1) 
					for(i=3;i<=NF;++i) s[i]+=$i;
				next }
		FNR==1{print;next}
		{printf $1"\t"$2; for(i=3;i<=NF;++i) printf "\t%.6f",$i/s[i]*1000000; printf "\n"}' \
		$1 $1 > rpm.tmp
fi

# make rank normalization
Rscript --vanilla \
	-e 'a=read.table("rpm.tmp",header=T,stringsAsFactors=F)' \
	-e 'b=apply(a[,-(1:2)],2,rank);b=b/dim(b)[1]' \
	-e 'write.table(format(data.frame(a[,1:2],b),digits=6,scientific=F),file="quantile.tmp",quote=F,sep="\t",col.names=T,row.names=F);'

# make DEseq normalization
Rscript --vanilla \
	-e 'a=read.table("rpm.tmp",header=T,stringsAsFactors=F)' \
	-e 'b=apply(a[,-(1:2)],1,function(x) exp(mean(log(x))))' \
	-e 'c=a[b>0,-(1:2)]/b[b>0];d=t(t(a[,-(1:2)])/apply(c,2,median))' \
	-e 'write.table(format(data.frame(a[,1:2],d),digits=4,scientific=F),file="deseq.tmp",quote=F,sep="\t",col.names=T,row.names=F);'

# make quantile normalization to DEseq approximation
Rscript --vanilla \
	-e 'a=read.table("deseq.tmp",header=T,stringsAsFactors=F); b=a[,-(1:2)]' \
	-e 'c=apply(apply(b,2,rank),2,function(x) apply(apply(b,2,sort),1,mean)[x])' \
	-e 'write.table(data.frame(a[,1:2],c),file="norm.tmp",quote=F,sep="\t",col.names=T,row.names=F);'

mv norm.tmp ${1%.*}.norm.txt
rm *.tmp
