#!/bin/bash -l
#SBATCH -A b2011026
#SBATCH -p node
#SBATCH -n 32
#SBATCH -t 7-00:00:00
#SBATCH -J coca_train_combiner

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#define default values
PREDICTION_RESULT_FILE_DEFAULT="/home/jessada/development/scilifelab/assignments/20121119_CombiVEP_publication/resources/prediction_result"

usage=$(
cat <<EOF
usage:
$0 [OPTION]
option:
-p FILE    specify path of prediction result file
EOF
)

die () {
    echo >&2 "[exception] $@"
    echo >&2 "$usage"
    exit 1
}

#parse param
while getopts "p:" OPTION; do
  case "$OPTION" in
    p)
      prediction_result_file="$OPTARG"
      ;;
    *)
      die "unrecognized option"
      ;;
  esac
done

#setting default values:
: ${prediction_result_file=$PREDICTION_RESULT_FILE_DEFAULT}

#display parameter
cat <<EOF
configuration:
prediction result file : $prediction_result_file
EOF

################################## start coding section ####################################
tmp_dir=$script_dir/tmp
task_list_file=$tmp_dir/task_list
tmp_condel_result=$tmp_dir/tmp_condel_result
tmp_file=$tmp_dir/tmp_file
parsed_condel_result_file=$tmp_dir/parsed_condel_result
parsed_existing_prediction_result=$tmp_dir/parsed_existing_prediction_result

#clear old data
if [ -e $tmp_condel_result ]
then
    rm $tmp_condel_result
fi
if [ -e $tmp_file ]
then
    rm $tmp_file
fi
if [ -e $parsed_condel_result_file ]
then
    rm $parsed_condel_result_file
fi
if [ -e $parsed_existing_prediction_result ]
then
    rm $parsed_existing_prediction_result
fi

pad_length=2
pad=$(printf '%0.1s' "0"{1..2})

while read line
do
    curl -X get $line | grep -v "^#" >> $tmp_condel_result
done < $task_list_file

nawk -F $'\t' '{printf("%s:%s:%s\n", $2, $3, $15)}' < $tmp_condel_result >> $tmp_file

grep "^[0-9]" $tmp_file | awk -F ':' '{printf("%d\t%d\t%s\t%.3f\t%02d|%012d|%s|%.3f\t%02d|%012d|%s\n", $1, $2, $3, $4, $1, $2, $3, $4, $1, $2, $3)}' | sort -k5 -r | uniq -f5 | sort -k5 >> $parsed_condel_result_file
grep -v "^[0-9]" $tmp_file | awk -F ':' '{printf("%s\t%d\t%s\t%.3f\t%s|%012d|%s|%.3f\t%s|%012d|%s\n", $1, $2, $3, $4, $1, $2, $3, $4, $1, $2, $3)}' | sort -k5 -r | uniq -f5 | sort -k5 >> $parsed_condel_result_file

grep -v "^#" $prediction_result_file | grep "^[0-9]" | awk -F'\t' '{ printf "%s\t%02d|%012d|%s\n", $0, $1, $2, $4}' | sort -k 13 >> $parsed_existing_prediction_result
grep -v "^#" $prediction_result_file | grep -v "^[0-9]" | awk -F'\t' '{ printf "%s\t%s|%012d|%s\n", $0, $1, $2, $4}' | sort -k 13 >> $parsed_existing_prediction_result

echo -e "#CHROM\tPOS\tREF\tALT\tACTUAL_DELETERIOUS_EFFECT\tPREDICTED_DELETERIOUS_PROBABILITY\tPHYLOP_SCORE\tSIFT_SCORE\tPP2_SCORE\tLRT_SCORT\tMT_SCORE\tGERP_SCORE\tCONDEL_SCORE"
join -t $'\t' -1 6 -2 13 -o 2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,1.4 $parsed_condel_result_file $parsed_existing_prediction_result

