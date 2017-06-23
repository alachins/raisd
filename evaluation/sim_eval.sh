LENGTH=1000000
TARGET=500000
THRESHOLD=10000

cd sim_eval

./evaluate_dataset_by_index.sh 2 $LENGTH $TARGET $THRESHOLD
./evaluate_dataset_by_index.sh 5 $LENGTH $TARGET $THRESHOLD
./evaluate_dataset_by_index.sh 11 $LENGTH $TARGET $THRESHOLD
./evaluate_dataset_by_index.sh 12 $LENGTH $TARGET $THRESHOLD
./evaluate_dataset_by_index.sh 20 $LENGTH $TARGET $THRESHOLD
./evaluate_dataset_by_index.sh 21 $LENGTH $TARGET $THRESHOLD
./evaluate_dataset_by_index.sh 24 $LENGTH $TARGET $THRESHOLD
./evaluate_dataset_by_index.sh 45 $LENGTH $TARGET $THRESHOLD
./evaluate_dataset_by_index.sh 59 $LENGTH $TARGET $THRESHOLD
./evaluate_dataset_by_index.sh 60 $LENGTH $TARGET $THRESHOLD

./evaluate_dataset_by_index.sh 61 $LENGTH $TARGET $THRESHOLD
./evaluate_dataset_by_index.sh 62 $LENGTH $TARGET $THRESHOLD
./evaluate_dataset_by_index.sh 63 $LENGTH $TARGET $THRESHOLD
./evaluate_dataset_by_index.sh 64 $LENGTH $TARGET $THRESHOLD
./evaluate_dataset_by_index.sh 65 $LENGTH $TARGET $THRESHOLD
./evaluate_dataset_by_index.sh 66 $LENGTH $TARGET $THRESHOLD
./evaluate_dataset_by_index.sh 67 $LENGTH $TARGET $THRESHOLD
./evaluate_dataset_by_index.sh 68 $LENGTH $TARGET $THRESHOLD
./evaluate_dataset_by_index.sh 69 $LENGTH $TARGET $THRESHOLD
./evaluate_dataset_by_index.sh 70 $LENGTH $TARGET $THRESHOLD

./evaluate_dataset_by_index.sh 71 $LENGTH $TARGET $THRESHOLD
./evaluate_dataset_by_index.sh 72 $LENGTH $TARGET $THRESHOLD

./create_avg_dist_barplot.sh 2 5 11 12 20 21 24 45 59 60
./create_suc_rate_barplot.sh 2 5 11 12 20 21 24 45 59 60
./create_exec_time_barplot.sh 2 5 11 12 20 21 24 45 59 60
./create_roc_for_dataset_by_index.sh 2
./create_roc_for_dataset_by_index.sh 5
./create_roc_for_dataset_by_index.sh 11
./create_roc_for_dataset_by_index.sh 12
./create_roc_for_dataset_by_index.sh 20
./create_roc_for_dataset_by_index.sh 21
./create_roc_for_dataset_by_index.sh 24
./create_roc_for_dataset_by_index.sh 45
./create_roc_for_dataset_by_index.sh 59
./create_roc_for_dataset_by_index.sh 60

mv average_distance_barplot.pdf average_distance_barplot.pdf.1
mv average_exectime_barplot.pdf average_exectime_barplot.pdf.1
mv average_success_barplot.pdf average_success_barplot.pdf.1


./create_avg_dist_barplot.sh 61 62 63 64 65 66 67 68 69 70
./create_suc_rate_barplot.sh 61 62 63 64 65 66 67 68 69 70
./create_exec_time_barplot.sh 61 62 63 64 65 66 67 68 69 70
./create_roc_for_dataset_by_index.sh 61
./create_roc_for_dataset_by_index.sh 62
./create_roc_for_dataset_by_index.sh 63
./create_roc_for_dataset_by_index.sh 64
./create_roc_for_dataset_by_index.sh 65
./create_roc_for_dataset_by_index.sh 66
./create_roc_for_dataset_by_index.sh 67
./create_roc_for_dataset_by_index.sh 68
./create_roc_for_dataset_by_index.sh 69
./create_roc_for_dataset_by_index.sh 70

mv average_distance_barplot.pdf average_distance_barplot.pdf.2
mv average_exectime_barplot.pdf average_exectime_barplot.pdf.2
mv average_success_barplot.pdf average_success_barplot.pdf.2

./create_avg_dist_barplot.sh 71 72 71 72 71 72 71 72 71 72
./create_suc_rate_barplot.sh 71 72 71 72 71 72 71 72 71 72
./create_exec_time_barplot.sh 71 72 71 72 71 72 71 72 71 72
./create_roc_for_dataset_by_index.sh 71
./create_roc_for_dataset_by_index.sh 72

mv average_distance_barplot.pdf average_distance_barplot.pdf.3
mv average_exectime_barplot.pdf average_exectime_barplot.pdf.3
mv average_success_barplot.pdf average_success_barplot.pdf.3




#evince Rplots.pdf


