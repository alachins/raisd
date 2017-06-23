INDEX=$1
LENGTH=$2
OMEGAPLUS_MINWIN=$3
OMEGAPLUS_MAXWIN=$4
GRID=$5

cd ../sim_data/bottleneck
./delete_dataset_by_index.sh $INDEX
./fetch_dataset_by_index.sh $INDEX
cd ../../sim_runs
./process_dataset_by_index.sh $INDEX $LENGTH $OMEGAPLUS_MINWIN $OMEGAPLUS_MAXWIN $GRID
cd ../sim_data/bottleneck
./delete_dataset_by_index.sh $INDEX

