INDEX=$1
LENGTH=$2
TARGET=$3
THRESHOLD=$4

rm -r d$INDEX
mkdir d$INDEX

echo "Dataset $INDEX Comparison Summary [Length: $LENGTH] [Target: $TARGET] [Threshold: $THRESHOLD] " > d$INDEX/d$INDEX\_evaluation_report.txt

./eval_raisd_on_dataset_by_index.sh $INDEX $LENGTH
./eval_omegaplus_on_dataset_by_index.sh $INDEX $LENGTH $TARGET $THRESHOLD
./eval_sweed_on_dataset_by_index.sh $INDEX $LENGTH $TARGET $THRESHOLD
./eval_sweepfinder2_on_dataset_by_index.sh $INDEX $LENGTH $TARGET $THRESHOLD
