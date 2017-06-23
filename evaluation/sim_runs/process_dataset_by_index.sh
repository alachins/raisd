INDEX=$1
LENGTH=$2
OMEGAPLUS_MINWIN=$3
OMEGAPLUS_MAXWIN=$4
GRID=$5

rm -r d$INDEX
mkdir d$INDEX
cd d$INDEX
../run_raisd_on_dataset_by_index.sh $INDEX $LENGTH
../run_omegaplus_on_dataset_by_index.sh $INDEX $OMEGAPLUS_MINWIN $OMEGAPLUS_MAXWIN $LENGTH $GRID
../run_sweed_on_dataset_by_index.sh $INDEX $LENGTH $GRID
../run_sweepfinder2_on_dataset_by_index.sh $INDEX $LENGTH $GRID
