#Execution parameters for simulated data
LENGTH=1000000
OMEGAPLUS_MINWIN=10000
OMEGAPLUS_MAXWIN=50000
GRID=10

cd sim_runs

./execute_all_on_dataset_by_index.sh 1 $LENGTH $OMEGAPLUS_MINWIN $OMEGAPLUS_MAXWIN $GRID > run_1_monitor.txt
./execute_all_on_dataset_by_index.sh 2 $LENGTH $OMEGAPLUS_MINWIN $OMEGAPLUS_MAXWIN $GRID > run_2_monitor.txt
./execute_all_on_dataset_by_index.sh 3 $LENGTH $OMEGAPLUS_MINWIN $OMEGAPLUS_MAXWIN $GRID > run_3_monitor.txt
./execute_all_on_dataset_by_index.sh 4 $LENGTH $OMEGAPLUS_MINWIN $OMEGAPLUS_MAXWIN $GRID > run_4_monitor.txt
./execute_all_on_dataset_by_index.sh 5 $LENGTH $OMEGAPLUS_MINWIN $OMEGAPLUS_MAXWIN $GRID > run_5_monitor.txt

./execute_all_on_dataset_by_index.sh 6 $LENGTH $OMEGAPLUS_MINWIN $OMEGAPLUS_MAXWIN $GRID > run_6_monitor.txt
./execute_all_on_dataset_by_index.sh 7 $LENGTH $OMEGAPLUS_MINWIN $OMEGAPLUS_MAXWIN $GRID > run_7_monitor.txt
./execute_all_on_dataset_by_index.sh 8 $LENGTH $OMEGAPLUS_MINWIN $OMEGAPLUS_MAXWIN $GRID > run_8_monitor.txt
./execute_all_on_dataset_by_index.sh 9 $LENGTH $OMEGAPLUS_MINWIN $OMEGAPLUS_MAXWIN $GRID > run_9_monitor.txt
./execute_all_on_dataset_by_index.sh 10 $LENGTH $OMEGAPLUS_MINWIN $OMEGAPLUS_MAXWIN $GRID > run_10_monitor.txt


