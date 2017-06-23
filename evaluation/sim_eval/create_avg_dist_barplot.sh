AA=1
NUMBER=$1
v1=$(grep DistanceError d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=3 '{print $3s}')
echo $AA $v1 > distance_table.txt

AA=2
NUMBER=$2
v1=$(grep DistanceError d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=3 '{print $3s}')
echo $AA $v1 >> distance_table.txt

AA=3
NUMBER=$3
v1=$(grep DistanceError d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=3 '{print $3s}')
echo $AA $v1 >> distance_table.txt

AA=4
NUMBER=$4
v1=$(grep DistanceError d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=3 '{print $3s}')
echo $AA $v1 >> distance_table.txt

AA=5
NUMBER=$5
v1=$(grep DistanceError d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=3 '{print $3s}')
echo $AA $v1 >> distance_table.txt

AA=6
NUMBER=$6
v1=$(grep DistanceError d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=3 '{print $3s}')
echo $AA $v1 >> distance_table.txt

AA=7
NUMBER=$7
v1=$(grep DistanceError d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=3 '{print $3s}')
echo $AA $v1 >> distance_table.txt

AA=8
NUMBER=$8
v1=$(grep DistanceError d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=3 '{print $3s}')
echo $AA $v1 >> distance_table.txt

AA=9
NUMBER=$9
v1=$(grep DistanceError d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=3 '{print $3s}')
echo $AA $v1 >> distance_table.txt

AA=10
NUMBER="${10}"
v1=$(grep DistanceError d$NUMBER/d$NUMBER\_evaluation_report.txt | awk -v N=3 '{print $3s}')
echo $AA $v1 >> distance_table.txt

Rscript create_avg_dist_barplot.R $@
rm distance_table.txt

mv Rplots.pdf average_distance_barplot.pdf


