#Execution parameters for real data
OMEGAPLUS_MINWIN=15000
OMEGAPLUS_MAXWIN=100000
GRID=10

cd igsr_data/populations
./generate_population_samples.sh

cd ../../igsr_runs
./execute_all_on_chromosome_by_index.sh 22 $OMEGAPLUS_MINWIN $OMEGAPLUS_MAXWIN $GRID
