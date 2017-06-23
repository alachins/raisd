NUMBER=$1 # dataset index
LENGTH=$2
TARGET=$3
THRESHOLD=$4

TOOL=SweepFinder2
tool=sweepfinder2

echo "$TOOL Accuracy for Dataset $NUMBER"
v1=$TARGET 
echo "Target selection location: $v1"

v2=$(grep maxsweep ../sim_runs/d$NUMBER/d$NUMBER\_selection\_$tool\_all/$TOOL\_Info.d$NUMBER\_selection\_$tool\_*.txt | awk -v N=3 '{print $3s}' | cut -d'=' -f 2 | perl -e 'print sort{$a<=>$b}<>') # this is selection
#IFS=$'\n' sorted=($(sort -n <<<"${v2[*]}"))
IFS=$'\n' sorted=($(cat <<<"${v2[*]}"))
printf "%s\n" "${#sorted[@]}" >> sorted_list.txt
printf "%s\n" "$v1" >> sorted_list.txt
printf "%s\n" "$LENGTH" >> sorted_list.txt
printf "%s\n" "$THRESHOLD" >> sorted_list.txt
printf "%s\n" "${sorted[@]}" >> sorted_list.txt

v2=$(grep maxsweep ../sim_runs/d$NUMBER/d$NUMBER\_neutral\_$tool\_all/$TOOL\_Info.d$NUMBER\_neutral\_$tool\_*.txt | awk -v N=2 '{print $2s}' | cut -d'=' -f 2 | perl -e 'print sort{$a<=>$b}<>')
#IFS=$'\n' sorted=($(sort -n <<<"${v2[*]}"))
IFS=$'\n' sorted=($(cat <<<"${v2[*]}"))
printf "%s\n" "${#sorted[@]}" >> sorted_list2.txt
printf "50\n" >> sorted_list2.txt
printf "%s\n" "${sorted[@]}" >> sorted_list2.txt

v2=$(grep maxsweep ../sim_runs/d$NUMBER/d$NUMBER\_selection\_$tool\_all/$TOOL\_Info.d$NUMBER\_selection\_$tool\_*.txt | awk -v N=2 '{print $2s}' | cut -d'=' -f 2 | perl -e 'print sort{$a<=>$b}<>') # this is selection
#IFS=$'\n' sorted=($(sort -n <<<"${v2[*]}"))
IFS=$'\n' sorted=($(cat <<<"${v2[*]}"))
printf "%s\n" "${#sorted[@]}" >> sorted_list3.txt
printf "%s\n" "${sorted[@]}" >> sorted_list3.txt

gcc -o calc_acc calc_acc.c
./calc_acc d$NUMBER/d$NUMBER\_evaluation_report.txt $TOOL

v5=$(grep real ../sim_runs/d$NUMBER/d$NUMBER\_neutral\_$tool\_all/$TOOL\_Info.d$NUMBER\_neutral\_$tool\_*.txt | awk -v N=2 '{print $2s}' | cut -d'.' -f 1 | cut -d'm' -f 2)
total_seconds=0
for i in ${v5[@]}; do
  let total_seconds+=$i
done
v5=$(grep real ../sim_runs/d$NUMBER/d$NUMBER\_neutral\_$tool\_all/$TOOL\_Info.d$NUMBER\_neutral\_$tool\_*.txt | awk -v N=2 '{print $2s}' | cut -d'm' -f 1)
total_minutes=0
for i in ${v5[@]}; do
  let total_minutes+=$i
done
v5=$(grep real ../sim_runs/d$NUMBER/d$NUMBER\_neutral\_$tool\_all/$TOOL\_Info.d$NUMBER\_neutral\_$tool\_*.txt | awk -v N=2 '{print $2s}' | cut -d'.' -f 2 | cut -d "s" -f 1)
total_seconds_fraction=0
for i in ${v5[@]}; do
  total_seconds_fraction=$(bc <<<"scale=5; $total_seconds_fraction + $i")
done
total_seconds_fraction=$(bc <<<"scale=5; $total_seconds_fraction / 1000.0")
total_time=$total_minutes
total_time=$(bc <<<"scale=5; $total_minutes * 60.0")
total_time=$(bc <<<"scale=5; $total_time + $total_seconds")
total_time=$(bc <<<"scale=5; $total_time + $total_seconds_fraction")
echo "Execution time (neutral): $total_time"
V4=$total_time

v5=$(grep real ../sim_runs/d$NUMBER/d$NUMBER\_selection\_$tool\_all/$TOOL\_Info.d$NUMBER\_selection\_$tool\_*.txt | awk -v N=2 '{print $2s}' | cut -d'.' -f 1 | cut -d'm' -f 2) #selection
total_seconds=0
for i in ${v5[@]}; do
  let total_seconds+=$i
done
v5=$(grep real ../sim_runs/d$NUMBER/d$NUMBER\_selection\_$tool\_all/$TOOL\_Info.d$NUMBER\_selection\_$tool\_*.txt | awk -v N=2 '{print $2s}' | cut -d'm' -f 1) #selection
total_minutes=0
for i in ${v5[@]}; do
  let total_minutes+=$i
done
v5=$(grep real ../sim_runs/d$NUMBER/d$NUMBER\_selection\_$tool\_all/$TOOL\_Info.d$NUMBER\_selection\_$tool\_*.txt | awk -v N=2 '{print $2s}' | cut -d'.' -f 2 | cut -d "s" -f 1) #selection
total_seconds_fraction=0
for i in ${v5[@]}; do
  total_seconds_fraction=$(bc <<<"scale=5; $total_seconds_fraction + $i")
done
total_seconds_fraction=$(bc <<<"scale=5; $total_seconds_fraction / 1000.0")
total_time=$total_minutes
total_time=$(bc <<<"scale=5; $total_minutes * 60.0")
total_time=$(bc <<<"scale=5; $total_time + $total_seconds")
total_time=$(bc <<<"scale=5; $total_time + $total_seconds_fraction")
echo "Execution time (selection): $total_time"
V5=$total_time

echo "ExecutionNeutral: $V4 s | ExecutionSelection: $V5 s" >> d$NUMBER/d$NUMBER\_evaluation_report.txt

rm sorted_list.txt
rm sorted_list2.txt
rm sorted_list3.txt



