NUMBER=$1 # dataset index
LENGTH=$2
TARGET=$3
THRESHOLD=$4

TOOL=OmegaPlus
tool=omegaplus

echo "$TOOL Accuracy for Dataset $NUMBER"
v1=$TARGET 
echo "Target selection location: $v1"

v2=$(grep Location ../sim_runs/d$NUMBER/$TOOL\_Info.d$NUMBER\_selection\_$tool | awk -v N=2 '{print $2s}' | perl -e 'print sort{$a<=>$b}<>')
#IFS=$'\n' sorted=($(sort -n <<<"${v2[*]}"))
IFS=$'\n' sorted=($(cat <<<"${v2[*]}"))
printf "%s\n" "${#sorted[@]}" > sorted_list.txt
printf "%s\n" "$v1" >> sorted_list.txt
printf "%s\n" "$LENGTH" >> sorted_list.txt
printf "%s\n" "$THRESHOLD" >> sorted_list.txt
printf "%s\n" "${sorted[@]}" >> sorted_list.txt

v2=$(grep Max ../sim_runs/d$NUMBER/$TOOL\_Info.d$NUMBER\_neutral\_$tool | awk -v N=3 '{print $3s}' | perl -e 'print sort{$a<=>$b}<>')
#IFS=$'\n' sorted=($(sort -n <<<"${v2[*]}"))
IFS=$'\n' sorted=($(cat <<<"${v2[*]}"))
printf "%s\n" "${#sorted[@]}" >> sorted_list2.txt
printf "50\n" >> sorted_list2.txt
printf "%s\n" "${sorted[@]}" >> sorted_list2.txt

v2=$(grep Max ../sim_runs/d$NUMBER/$TOOL\_Info.d$NUMBER\_selection\_$tool | awk -v N=3 '{print $3s}' | perl -e 'print sort{$a<=>$b}<>')
#IFS=$'\n' sorted=($(sort -n <<<"${v2[*]}"))
IFS=$'\n' sorted=($(cat <<<"${v2[*]}"))
printf "%s\n" "${#sorted[@]}" >> sorted_list3.txt
printf "%s\n" "${sorted[@]}" >> sorted_list3.txt

gcc -o calc_acc calc_acc.c
./calc_acc d$NUMBER/d$NUMBER\_evaluation_report.txt $TOOL

v5=$(grep time ../sim_runs/d$NUMBER/$TOOL\_Info.d$NUMBER\_neutral\_$tool | awk -v N=4 '{print $4s}')
echo "Execution time (neutral): $v5"
V4=$v5
v5=$(grep time ../sim_runs/d$NUMBER/$TOOL\_Info.d$NUMBER\_selection\_$tool | awk -v N=4 '{print $4s}')
echo "Execution time (selection): $v5"
V5=$v5

echo "ExecutionNeutral: $V4 s | ExecutionSelection: $V5 s" >> d$NUMBER/d$NUMBER\_evaluation_report.txt

rm sorted_list.txt
rm sorted_list2.txt
rm sorted_list3.txt





