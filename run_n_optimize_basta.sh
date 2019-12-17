#!/usr/bin/env bash
# This code run a grid search for the best parameters of a basta sequence call
# using a mock community list of species
#set -e
## Commands
OLDIFS=${IFS}
basta_db=$1
blast=$2
out_prefix=$3
config=$4
true_file=$5
basta=$6
## Set the grid
evalues=( 1E-80 1E-40 1E-20 1E-10 1E-5 1E-2 )
p_ids=( 70 80 90 95 100 )
m_hits={1..10}
p_hits=( 60 70 80 90 99 )
n_hits={0..10}
string=`IFS=,;eval echo "{${evalues[*]}}_{${p_ids[*]}}_${m_hits[*]}_{${p_hits[*]}}_${n_hits}"`
IFS=" ";perms=(${string})
IFS=${OLDIFS}
# set variuable for hash table
declare -A results
## Set functions
prog() {
    # Progress bar, courtesy of ilkkachu (stackoverflow)
    local w=80 p=$1;  shift
    printf -v dots "%*s" "$(( $p*$w/100 ))" ""; dots=${dots// /.};
    printf "\r\e[K|%-*s| %3d %% %s" "$w" "$dots" "$p" "$*";
}

run_basta(){
# Run basta with parameters given by variable 1
# 1) "_"-delimited string of basta parameters
# 2) blast file
# 3) mapping (gb or the one you constructed)
# 4) config file
# 5) File with list of True lables
# 6) workerid
#echo "running run_basta $@"
local IFS="_"; params=(${1})
IFS=${OLDIFS}
local IFS=${OLDIFS}
local evalue="${params[0]}"
local p_id="${params[1]}"
local m_hit="${params[2]}"
local p_hit="${params[3]}"
local n_hit="${params[4]}"
local blast=${2}
local mapping=${3}
local config_file=${4}
local true_lables=$(sort -u ${5})
local workerid=${6}
#echo "running basta"
python2 ${basta} sequence "${blast}" basta${workerid}.out ${mapping} \
-e ${evalue} -i ${p_id} -m ${m_hit} -n ${n_hit} -p ${p_hit} \
-c ${config_file} 2>/dev/null
#echo "Extracting metrics"
pred=$(cut -f 2 basta.out| cut -d ';' -f 7- | sed -e '/^$/d' -e 's/;//' -e 's/_/ /g' | sort -u)
local false_negatives=`comm -23 <(echo "${true_lables[@]}") <(echo "${pred[@]}") | wc -l`
local false_positives=`comm -13 <(echo "${true_lables[@]}") <(echo "${pred[@]}") | wc -l`
local true_positives=`comm -12 <(echo "${true_lables[@]}") <(echo "${pred[@]}"}) | wc -l`
#echo "Computing F1"
intF1=$(( (true_positives * 100) / (true_positives + false_negatives + false_positives) ))
#echo "F1 is ${intF1}"
#echo "Removing basta file"
echo -e "${param}\t$intF1" >> results.basta
rm basta${workerid}.out
}

export -f run_basta
echo '' > results.basta # if already exist overwrite
total="${#perms[@]}"
counter=0
for param in "${perms[@]}"
  do
    let counter++
    prog $(( (counter * 100) / total ))
    run_basta "${param}" "${blast}" "${basta_db}" "${config}" ${true_file}
#    echo "out F1 ${intF1}"
#    results["${param}"]="${intF1}"
#    unset intF1
  done


#for k in "${!results[@]}"
#  do
#    echo "${k} ${results[$k]}" >> results.basta
#  done
IFS="_"; read -ra best <<< `cat results.basta| sort -rn -k2| head -n 1 | cut -f 1 -d ' '`
IFS=${OLDIFS}
echo "basta sequence "${blast}" ${out_prefix}.out ${basta_db} \\" > best.param
echo "-e ${best[0]} -i ${best[1]} -m ${best[2]} -n ${best[3]} -p ${best[4]}"  >> best.param