INDEX=$1
OMEGAPLUS_MINWIN=$2
OMEGAPLUS_MAXWIN=$3
GRID=$4

rm -r chr$INDEX
mkdir chr$INDEX
cd chr$INDEX
../run_raisd_on_chromosome_by_index.sh $INDEX
../run_omegaplus_on_chromosome_by_index.sh $INDEX $OMEGAPLUS_MINWIN $OMEGAPLUS_MAXWIN $GRID
../run_sweed_on_chromosome_by_index.sh $INDEX $GRID
../run_sweepfinder2_on_chromosome_by_index.sh $INDEX $GRID
