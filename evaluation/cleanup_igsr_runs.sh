cd igsr_data/populations
rm *.txt
cd ../vcf
for i in {1..22} 
do
 ./delete_chromosome_vcf_by_index.sh $i
done
cd ../../igsr_runs
for i in {1..22} 
do
 rm -r chr$i
done
