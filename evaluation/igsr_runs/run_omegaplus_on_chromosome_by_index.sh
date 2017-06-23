NUMBER=$1 # dataset index
MINWIN=$2 # minwin
MAXWIN=$3 # maxwin
GRID=$4 # grid

DATASET=chr$NUMBER
TOOLFOLDER=omegaplus
TOOL=OmegaPlus-C
THREADS=2

PATH_TO_TOOL=../../tools/$TOOLFOLDER/$TOOL
CHR=chr$NUMBER
PATH_TO_DATASET=../../igsr_data/vcf/chr$NUMBER/ALL.chr$NUMBER.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf

POPULATION=acb
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=asw
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=beb
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=cdx
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=ceu
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=chb
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=chd
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=chs
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=clm
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=esn
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=fin
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=gbr
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=gih
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=gwd
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=gwf
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=gwj
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=gww
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=ibs
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=itu
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=jpt
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=khv
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=lwk
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=msl
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=mxl
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=pel
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=pjl
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=pur
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=stu
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=tsi
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt

POPULATION=yri
rm OmegaPlus_Info.$CHR\_$POPULATION\_$TOOLFOLDER
rm OmegaPlus_Report.$CHR\_$POPULATION\_$TOOLFOLDER
$PATH_TO_TOOL -name $CHR\_$POPULATION\_$TOOLFOLDER -threads $THREADS -grid $GRID -minwin $MINWIN -maxwin $MAXWIN -seed 12345 -input $PATH_TO_DATASET -sampleList ../../igsr_data/populations/$POPULATION.txt
