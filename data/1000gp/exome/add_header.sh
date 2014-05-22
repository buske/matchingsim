#for i in {1..22} X 
for i in 1 2
do
	zcat $i.vcf.gz | cat vcf_header.txt - | gzip >> n$i.vcf.gz
	echo Done $i!
done
