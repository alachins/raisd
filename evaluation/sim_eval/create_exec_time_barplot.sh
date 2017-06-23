AA=1
NUMBER=$1
v1=$(grep ExecutionNeutral d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=15 '{print $15s}')
echo $AA $v1 > exectime_table.txt
v2=$(grep ExecutionNeutral d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=19 '{print $19s}')
echo $AA $v2 >> exectime_table.txt

AA=2
NUMBER=$2
v1=$(grep ExecutionNeutral d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=15 '{print $15s}')
echo $AA $v1 >> exectime_table.txt
v1=$(grep ExecutionNeutral d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=19 '{print $19s}')
echo $AA $v1 >> exectime_table.txt

AA=3
NUMBER=$3
v1=$(grep ExecutionNeutral d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=15 '{print $15s}')
echo $AA $v1 >> exectime_table.txt
v1=$(grep ExecutionNeutral d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=19 '{print $19s}')
echo $AA $v1 >> exectime_table.txt

AA=4
NUMBER=$4
v1=$(grep ExecutionNeutral d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=15 '{print $15s}')
echo $AA $v1 >> exectime_table.txt
v1=$(grep ExecutionNeutral d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=19 '{print $19s}')
echo $AA $v1 >> exectime_table.txt

AA=5
NUMBER=$5
v1=$(grep ExecutionNeutral d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=15 '{print $15s}')
echo $AA $v1 >> exectime_table.txt
v1=$(grep ExecutionNeutral d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=19 '{print $19s}')
echo $AA $v1 >> exectime_table.txt

AA=6
NUMBER=$6
v1=$(grep ExecutionNeutral d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=15 '{print $15s}')
echo $AA $v1 >> exectime_table.txt
v1=$(grep ExecutionNeutral d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=19 '{print $19s}')
echo $AA $v1 >> exectime_table.txt

AA=7
NUMBER=$7
v1=$(grep ExecutionNeutral d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=15 '{print $15s}')
echo $AA $v1 >> exectime_table.txt
v1=$(grep ExecutionNeutral d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=19 '{print $19s}')
echo $AA $v1 >> exectime_table.txt

AA=8
NUMBER=$8
v1=$(grep ExecutionNeutral d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=15 '{print $15s}')
echo $AA $v1 >> exectime_table.txt
v1=$(grep ExecutionNeutral d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=19 '{print $19s}')
echo $AA $v1 >> exectime_table.txt

AA=9
NUMBER=$9
v1=$(grep ExecutionNeutral d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=15 '{print $15s}')
echo $AA $v1 >> exectime_table.txt
v1=$(grep ExecutionNeutral d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=19 '{print $19s}')
echo $AA $v1 >> exectime_table.txt

AA=10
NUMBER="${10}"
v1=$(grep ExecutionNeutral d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=15 '{print $15s}')
echo $AA $v1 >> exectime_table.txt
v1=$(grep ExecutionNeutral d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=19 '{print $19s}')
echo $AA $v1 >> exectime_table.txt

gcc -o calc_acc_2 calc_acc_2.c
size=$(tail -n 1 exectime_table.txt | awk -v N=1 '{print $1s}')
./calc_acc_2 $size

mv exectime_table2.txt exectime_table.txt
Rscript create_exec_time_barplot.R $@
rm exectime_table.txt

mv Rplots.pdf average_exectime_barplot.pdf


