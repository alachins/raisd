NUMBER=$1 # dataset index
REGIONLENGTH=$2 # region length
GRID=$3 # grid

DATASET=d$NUMBER
TOOLFOLDER=sweepfinder2
TOOL=SweepFinder2


../../tools/sweed/SweeD -name $DATASET\_neutral\_$TOOLFOLDER -input ../../sim_data/$DATASET/msneutral$NUMBER.out -length $REGIONLENGTH -osf $DATASET\_neutral\_$TOOLFOLDER.txt -reports > monitor_$DATASET\_neutral.txt
rm monitor_$DATASET\_neutral.txt
rm SweeD_Info.$DATASET\_neutral\_$TOOLFOLDER
rm SweeD_Warnings.$DATASET\_neutral\_$TOOLFOLDER
rm $DATASET\_neutral\_$TOOLFOLDER.txt
for i in {1..1000}
do
    (time ../../tools/$TOOLFOLDER/$TOOL -s $GRID $DATASET\_neutral\_$TOOLFOLDER.txt.$i $TOOL\_Report.$DATASET\_neutral\_$i) &> $TOOL\_Info.$DATASET\_neutral\_$TOOLFOLDER\_$i.txt
    rm $DATASET\_neutral\_$TOOLFOLDER.txt.$i
    rm SweeD_Report.$DATASET\_neutral\_$TOOLFOLDER.$i
done

mkdir $DATASET\_neutral\_$TOOLFOLDER\_all
mv $TOOL\_Report.* $DATASET\_neutral\_$TOOLFOLDER\_all
mv $TOOL\_Info.* $DATASET\_neutral\_$TOOLFOLDER\_all


../../tools/sweed/SweeD -name $DATASET\_selection\_$TOOLFOLDER -input ../../sim_data/$DATASET/msselection$NUMBER.out -length $REGIONLENGTH -osf $DATASET\_selection\_$TOOLFOLDER.txt -reports > monitor_$DATASET\_selection.txt
rm monitor_$DATASET\_selection.txt
rm SweeD_Info.$DATASET\_selection\_$TOOLFOLDER
rm SweeD_Warnings.$DATASET\_selection\_$TOOLFOLDER
rm $DATASET\_selection\_$TOOLFOLDER.txt
for i in {1..1000}
do
    (time ../../tools/$TOOLFOLDER/$TOOL -s $GRID $DATASET\_selection\_$TOOLFOLDER.txt.$i $TOOL\_Report.$DATASET\_selection\_$i) &> $TOOL\_Info.$DATASET\_selection\_$TOOLFOLDER\_$i.txt
    rm $DATASET\_selection\_$TOOLFOLDER.txt.$i
    rm SweeD_Report.$DATASET\_selection\_$TOOLFOLDER.$i
done

mkdir $DATASET\_selection\_$TOOLFOLDER\_all
mv $TOOL\_Report.* $DATASET\_selection\_$TOOLFOLDER\_all
mv $TOOL\_Info.* $DATASET\_selection\_$TOOLFOLDER\_all


#../../tools/$TOOLFOLDER/$TOOL -name $DATASET\_selection\_$TOOLFOLDER -input ../../sim_data/bottleneck/$DATASET/msselection$NUMBER.out -length $2 -grid $3

