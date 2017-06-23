NUMBER=$1 # dataset index
GRID=$2 # grid

DATASET=chr$NUMBER
TOOLFOLDER=sweepfinder2
TOOL=SweepFinder2

PATH_TO_TOOL=../../tools/$TOOLFOLDER/$TOOL
CHR=chr$NUMBER
PATH_TO_DATASET=../../igsr_data/vcf/chr$NUMBER/ALL.chr$NUMBER.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf

POPULATION=acb
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=asw
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=beb
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=cdx
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=ceu
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=chb
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=chd
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=chs
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=clm
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=esn
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=fin
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=gbr
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=gih
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=gwd
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=gwf
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=gwj
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=gww
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=ibs
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=itu
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=jpt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=khv
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=lwk
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=msl
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=mxl
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt


POPULATION=pel
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=pjl
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=pur
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=stu
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=tsi
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt

POPULATION=yri
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
../../tools/sweed/SweeD -name $CHR\_$POPULATION\_$TOOLFOLDER -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic -osf $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
(time $PATH_TO_TOOL -s $GRID $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt $TOOL\_Report.$DATASET) &> $TOOL\_Info.$DATASET
rm $CHR\_$POPULATION\_$TOOLFOLDER.sf.txt


