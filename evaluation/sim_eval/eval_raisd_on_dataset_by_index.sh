NUMBER=$1 # dataset index
LENGTH=$2

TOOL=RAiSD
tool=raisd

echo "RAiSD Accuracy for Dataset $NUMBER"
v1=$(grep Target ../sim_runs/d$NUMBER/RAiSD\_Info.d$NUMBER\_selection\_$tool | tail -n -1 | head -n 1 | awk -v N=4 '{print $4s}')
v1=${v1::-1}
echo "Target selection location: $v1"

v2=$(grep MuStat ../sim_runs/d$NUMBER/RAiSD\_Info.d$NUMBER\_selection\_$tool | tail -n -2 | head -n 1 | awk -v N=2 '{print $2s}')
v3=$(bc <<<"scale=5; $v2 / $LENGTH")
v4=$(bc <<<"scale=5; $v3 * 100.0")
echo "Average reported distance to Target: $v2 ($v4%)"
V1=$v4

v2=$(grep Distance ../sim_runs/d$NUMBER/RAiSD\_Info.d$NUMBER\_selection\_$tool | tail -n -1 | head -n 1 | awk -v N=4 '{print $4s}')
v2=${v2::-1}
v3=$(bc <<<"scale=5; $v2 / $LENGTH")
v4=$(bc <<<"scale=5; $v3 * 100.0")
echo "Error distance threshold: $v2 ($v4%)"

v2=$(grep MuStat ../sim_runs/d$NUMBER/RAiSD\_Info.d$NUMBER\_selection\_$tool | tail -n -1 | head -n 1 | awk -v N=2 '{print $2s}')
v4=$(bc <<<"scale=5; $v2 * 100.0")
echo "Success rate (reported locations at distance to selection target < threshold): $v2 ($v4%)"
V2=$v4

v2=$(grep TPR ../sim_runs/d$NUMBER/RAiSD\_Info.d$NUMBER\_selection\_$tool | tail -n -1 | head -n 1 | awk -v N=2 '{print $2s}')
v2=$(bc <<<"scale=5; $v2 * 100.0")
echo "True Positive Rate for 5% False Positive Rate: $v2%"
V3=$v2

v5=$(grep time ../sim_runs/d$NUMBER/$TOOL\_Info.d$NUMBER\_neutral\_$tool | awk -v N=4 '{print $4s}')
echo "Execution time (neutral): $v5"
V4=$v5
v5=$(grep time ../sim_runs/d$NUMBER/$TOOL\_Info.d$NUMBER\_selection\_$tool | awk -v N=4 '{print $4s}')
echo "Execution time (selection): $v5"
V5=$v5

echo "" >> d$NUMBER/d$NUMBER\_evaluation_report.txt
echo "[RAiSD]         DistanceError: $V1 % | SuccessRate: $V2 % | TPR: $V3 % | ExecutionNeutral: $V4 s | ExecutionSelection: $V5 s" >> d$NUMBER/d$NUMBER\_evaluation_report.txt
 





