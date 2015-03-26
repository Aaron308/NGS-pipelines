###
#code to make it work on osx and linux
if
[[ $OSTYPE == darwin* ]]
then
readlink=$(which greadlink)
scriptdir="$(dirname $($readlink -f $0))"
else
scriptdir="$(dirname $(readlink -f $0))"
fi
#

usage="USAGE:
03-runner.sh <number of threads> <min length> <max length>"

######### Setup ################
threads=$1
min=$2
max=$3
# kefile format: (tab seperated)
#Ordinal Sample <factor1_name> [<factor2_name>]
if [ "$#" -lt "3" ]
then
echo $usage
exit -1
else
echo "initiating $1 parallel cutadapt adapter removal jobs, min/max length after filtering $min, $max"
fi
########## Run #################


#user defined variables that could be changed:
workingdir=./
script=$scriptdir/03-cutadapt.sh
outdir=reads_noadapt_cutadapt
###

function findSamples () {
find reads/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
}

mkdir $outdir
timestamp=$(date +%Y%m%d-%H%M%S)

logdir="./logs/${outdir}.${timestamp}"
mkdir $logdir

cat $script > "$logdir/script.log"
cat $0 > "$logdir/runner.log"
cat $script

findSamples | parallel bash $script {} $min $max \>logs/${outdir}.${timestamp}/{}.log 2\>\&1

#To run, got to directory containing reads directory and call:
#bash ~/path_to/02-runner.sh

