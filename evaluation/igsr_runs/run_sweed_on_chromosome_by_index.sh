NUMBER=$1 # dataset index
GRID=$2 # grid

DATASET=chr$NUMBER
TOOLFOLDER=sweed
TOOL=SweeD-P
THREADS=2

PATH_TO_TOOL=../../tools/$TOOLFOLDER/$TOOL
CHR=chr$NUMBER
PATH_TO_DATASET=../../igsr_data/vcf/chr$NUMBER/ALL.chr$NUMBER.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf


POPULATION=acb
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=asw
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=beb
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=cdx
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=ceu
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=chb
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=chd
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=chs
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=clm
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=esn
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=fin
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=gbr
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=gih
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=gwd
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=gwf
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=gwj
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=gww
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=ibs
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=itu
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=jpt
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=khv
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=lwk
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=msl
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=mxl
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=pel
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=pjl
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=pur
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=stu
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=tsi
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

POPULATION=yri
rm SweeD_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm SweeD_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt -strictPolymorphic

