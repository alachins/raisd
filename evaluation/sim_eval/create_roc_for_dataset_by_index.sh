NUMBER=$1 # dataset index

TOOL=RAiSD
tool=raisd
#perl -e 'print sort { $a<=>$b } <>' < input-file
v1=$(grep MuStat ../sim_runs/d$NUMBER/$TOOL\_Info.d$NUMBER\_neutral\_$tool | awk -v N=28 '{print $28s}' | perl -e 'print sort{$a<=>$b}<>') 
IFS=$'\n' sorted=($(cat <<<"${v1[*]}"))
#sorted=$v1
printf "%s\n" "${#sorted[@]}" > sorted_list.txt
printf "%s\n" "${sorted[@]}" >> sorted_list.txt

v1=$(grep MuStat ../sim_runs/d$NUMBER/$TOOL\_Info.d$NUMBER\_selection\_$tool | awk -v N=28 '{print $28s}' | perl -e 'print sort{$a<=>$b}<>')
#IFS=$'\n' sorted=($(sort -g <<<"${v1[*]}"))
IFS=$'\n' sorted=($(cat <<<"${v1[*]}"))
printf "%s\n" "${#sorted[@]}" >> sorted_list.txt
printf "%s\n" "${sorted[@]}" >> sorted_list.txt

TOOL=OmegaPlus
tool=omegaplus

v1=$(grep "Max Omega" ../sim_runs/d$NUMBER/$TOOL\_Info.d$NUMBER\_neutral\_$tool | awk -v N=3 '{print $3s}' | perl -e 'print sort{$a<=>$b}<>' )
#IFS=$'\n' sorted=($(sort -n <<<"${v1[*]}"))
IFS=$'\n' sorted=($(cat <<<"${v1[*]}"))
printf "%s\n" "${#sorted[@]}" >> sorted_list.txt
printf "%s\n" "${sorted[@]}" >> sorted_list.txt

v1=$(grep "Max Omega" ../sim_runs/d$NUMBER/$TOOL\_Info.d$NUMBER\_selection\_$tool | awk -v N=3 '{print $3s}' | perl -e 'print sort{$a<=>$b}<>')
#IFS=$'\n' sorted=($(sort -n <<<"${v1[*]}"))
IFS=$'\n' sorted=($(cat <<<"${v1[*]}"))
printf "%s\n" "${#sorted[@]}" >> sorted_list.txt
printf "%s\n" "${sorted[@]}" >> sorted_list.txt

TOOL=SweeD
tool=sweed

v1=$(grep Likelihood ../sim_runs/d$NUMBER/$TOOL\_Info.d$NUMBER\_neutral\_$tool | awk -v N=2 '{print $2s}' | perl -e 'print sort{$a<=>$b}<>')
#IFS=$'\n' sorted=($(sort -n <<<"${v1[*]}"))
IFS=$'\n' sorted=($(cat <<<"${v1[*]}"))
printf "%s\n" "${#sorted[@]}" >> sorted_list.txt
printf "%s\n" "${sorted[@]}" >> sorted_list.txt

v1=$(grep Likelihood ../sim_runs/d$NUMBER/$TOOL\_Info.d$NUMBER\_selection\_$tool | awk -v N=2 '{print $2s}' | perl -e 'print sort{$a<=>$b}<>')
#IFS=$'\n' sorted=($(sort -n <<<"${v1[*]}"))
IFS=$'\n' sorted=($(cat <<<"${v1[*]}"))
printf "%s\n" "${#sorted[@]}" >> sorted_list.txt
printf "%s\n" "${sorted[@]}" >> sorted_list.txt

TOOL=SweepFinder2
tool=sweepfinder2

v2=$(grep maxsweep ../sim_runs/d$NUMBER/d$NUMBER\_neutral\_$tool\_all/$TOOL\_Info.d$NUMBER\_neutral\_$tool\_*.txt | awk -v N=2 '{print $2s}' | cut -d'=' -f 2 | perl -e 'print sort{$a<=>$b}<>') # this is selection
#IFS=$'\n' sorted=($(sort -n <<<"${v2[*]}"))
IFS=$'\n' sorted=($(cat <<<"${v2[*]}"))
printf "%s\n" "${#sorted[@]}" >> sorted_list.txt
printf "%s\n" "${sorted[@]}" >> sorted_list.txt

v2=$(grep maxsweep ../sim_runs/d$NUMBER/d$NUMBER\_selection\_$tool\_all/$TOOL\_Info.d$NUMBER\_selection\_$tool\_*.txt | awk -v N=2 '{print $2s}' | cut -d'=' -f 2 | perl -e 'print sort{$a<=>$b}<>') # this is selection
#IFS=$'\n' sorted=($(sort -n <<<"${v2[*]}"))
IFS=$'\n' sorted=($(cat <<<"${v2[*]}"))
printf "%s\n" "${#sorted[@]}" >> sorted_list.txt
printf "%s\n" "${sorted[@]}" >> sorted_list.txt

gcc -o calc_roc_data calc_roc_data.c -lm
./calc_roc_data
rm calc_roc_data

gcc -o trim_string trim_string.c -lm

rm -r tmp
mkdir tmp
cd tmp
wget 139.91.162.50/raisd_data/d$NUMBER.tar.gz
tar -xvzf d$NUMBER.tar.gz
headerline=$(head -n 1 d$NUMBER/msneutral$NUMBER.out | cut -d' ' -f4-)
rm -r d$NUMBER 
echo $headerline > command_line.txt
cp ../trim_string .
./trim_string
headerline=$(head -n 1 command_line2.txt)
cd ..
rm -r tmp

#exit

Rscript create_roc_for_dataset_by_index.R $NUMBER $headerline

mv Rplots.pdf d$NUMBER\_roc.pdf 

rm sorted_list.txt
rm tpr_table.txt
rm command_line2.txt



