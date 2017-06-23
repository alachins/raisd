AA=1
NUMBER=$1
v1=$(grep SuccessRate d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=7 '{print $7s}')
echo $AA $v1 > success_table.txt

AA=2
NUMBER=$2
v1=$(grep SuccessRate d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=7 '{print $7s}')
echo $AA $v1 >> success_table.txt

AA=3
NUMBER=$3
v1=$(grep SuccessRate d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=7 '{print $7s}')
echo $AA $v1 >> success_table.txt

AA=4
NUMBER=$4
v1=$(grep SuccessRate d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=7 '{print $7s}')
echo $AA $v1 >> success_table.txt

AA=5
NUMBER=$5
v1=$(grep SuccessRate d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=7 '{print $7s}')
echo $AA $v1 >> success_table.txt

AA=6
NUMBER=$6
v1=$(grep SuccessRate d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=7 '{print $7s}')
echo $AA $v1 >> success_table.txt

AA=7
NUMBER=$7
v1=$(grep SuccessRate d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=7 '{print $7s}')
echo $AA $v1 >> success_table.txt

AA=8
NUMBER=$8
v1=$(grep SuccessRate d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=7 '{print $7s}')
echo $AA $v1 >> success_table.txt

AA=9
NUMBER=$9
v1=$(grep SuccessRate d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=7 '{print $7s}')
echo $AA $v1 >> success_table.txt

AA=10
NUMBER="${10}"
v1=$(grep SuccessRate d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=7 '{print $7s}')
echo $AA $v1 >> success_table.txt

Rscript create_suc_rate_barplot.R $@
rm success_table.txt

mv Rplots.pdf average_success_barplot.pdf


