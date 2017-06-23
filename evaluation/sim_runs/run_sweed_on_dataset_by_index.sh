NUMBER=$1 # dataset index
REGIONLENGTH=$2 # region length
GRID=$3 # grid

DATASET=d$NUMBER
TOOLFOLDER=sweed
TOOL=SweeD
../../tools/$TOOLFOLDER/$TOOL -name $DATASET\_neutral\_$TOOLFOLDER -input ../../sim_data/$DATASET/msneutral$NUMBER.out -length $2 -grid $3
../../tools/$TOOLFOLDER/$TOOL -name $DATASET\_selection\_$TOOLFOLDER -input ../../sim_data/$DATASET/msselection$NUMBER.out -length $2 -grid $3

