NUMBER=$1 # dataset index
REGIONLENGTH=$2 # region length

DATASET=d$NUMBER
TOOLFOLDER=raisd
TOOL=RAiSD

../../tools/$TOOLFOLDER/$TOOL -n $DATASET\_neutral\_$TOOLFOLDER -I ../../sim_data/$DATASET/msneutral$NUMBER.out -f -L $REGIONLENGTH -k 0.05
../../tools/$TOOLFOLDER/$TOOL -n $DATASET\_selection\_$TOOLFOLDER -I ../../sim_data/$DATASET/msselection$NUMBER.out -f -L $REGIONLENGTH -l $(grep "FPR Threshold" $TOOL\_Info.$DATASET\_neutral\_$TOOLFOLDER | awk -v N=3 '{print $3}') -T $(($REGIONLENGTH/2)) -d $(($REGIONLENGTH/100))

