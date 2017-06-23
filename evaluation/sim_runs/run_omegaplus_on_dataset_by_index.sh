NUMBER=$1 # dataset index
MINWIN=$2 # minwin
MAXWIN=$3 # maxwin
REGIONLENGTH=$4 # region length
GRID=$5 # grid

DATASET=d$NUMBER
TOOLFOLDER=omegaplus
TOOL=OmegaPlus
../../tools/$TOOLFOLDER/$TOOL -name $DATASET\_neutral\_$TOOLFOLDER -input ../../sim_data/$DATASET/msneutral$NUMBER.out -minwin $2 -maxwin $3 -length $4 -grid $5
../../tools/$TOOLFOLDER/$TOOL -name $DATASET\_selection\_$TOOLFOLDER -input ../../sim_data/$DATASET/msselection$NUMBER.out -minwin $2 -maxwin $3 -length $4 -grid $5

