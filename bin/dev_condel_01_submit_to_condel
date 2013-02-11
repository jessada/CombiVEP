#!/bin/bash -l
#SBATCH -A b2011026
#SBATCH -p node
#SBATCH -n 32
#SBATCH -t 7-00:00:00
#SBATCH -J coca_train_combiner

#************************* DISCLAIMER **************************
# This script is not a part of CombiVEP application
# It is made for the ease of author to reproduce the result
# Please feel free to use.
# Please don't expect to see it run smoothly without knowing what
# this script actually do.


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
tmp_parsed_file=$tmp_dir/tmp_parsed_file
task_list_file=$tmp_dir/task_list

if [ ! -e $tmp_dir ]
then
    mkdir $tmp_dir
fi

#parse the predictiion result into condel format
grep -v "^#" $prediction_result_file | awk -F'\t' '{ printf "%s\t%s\t%s\t%s\n", $1, $2, $2, $4 }' > $tmp_parsed_file

function split_and_submit_to_condel {
    tmp_condel_file_name=$1"_"`printf "%05d" "$2"`"_"`printf "%05d" "$3"`
    echo $tmp_condel_file_name
    sed -n $2","$3"p" $tmp_parsed_file > $tmp_condel_file_name
    curl -X PUT -T $tmp_condel_file_name http://bg.upf.edu/condel/taskService >> $4
}  

#clear old data
if [ -e $task_list_file ]
then
    rm $task_list_file
fi

#split condel input and then submit it to condel
file_size="$( wc -l $tmp_parsed_file | awk '{print $1}')"
begin_position=1
end_position=300
round=1
echo $file_size
while :
do
    echo "round" $round
    if [ $end_position -ge $file_size ]
    then
        split_and_submit_to_condel $tmp_parsed_file $begin_position $file_size $task_list_file
        break
    fi
    split_and_submit_to_condel $tmp_parsed_file $begin_position $end_position $task_list_file

    round=$(( round + 1 ))
    begin_position=$(( begin_position + 300 ))
    end_position=$(( end_position + 300 ))
done



