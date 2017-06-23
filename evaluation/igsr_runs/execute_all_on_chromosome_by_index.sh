INDEX=$1
OMEGAPLUS_MINWIN=$2
OMEGAPLUS_MAXWIN=$3
GRID=$4

cd ../igsr_data/vcf
./delete_chromosome_vcf_by_index.sh $INDEX
./fetch_chromosome_vcf_by_index.sh $INDEX
./unzip_chromosome_vcf_by_index.sh $INDEX
cd ../../igsr_runs
./process_chromosome_by_index.sh $INDEX $OMEGAPLUS_MINWIN $OMEGAPLUS_MAXWIN $GRID

