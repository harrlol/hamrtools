#!/bin/bash
set -u

# GO Analysis: mod based

if [ "$#" -lt 2 ]; then
echo "Missing arguments!"
echo "USAGE: " 
echo "EXAMPLE:"
exit 1
fi

# assign project directory
dir=$1

# assign panther API folder
py=$2

#current dir
curdir=$(dirname $0)

echo "generating genelist from mod table..."
# produce gene lists for all GMUCT (for now) groups
Rscript $curdir/produceGenelist.R \
    $dir/mod_long.csv \
    $dir/go/genelists

echo "sending each gene list to panther for overrepresentation analysis..."
# Send each gene list into panther API and generate a overrepresentation result file in another folter
for f in $dir/go/genelists/*.txt
do
n=$(basename $f)
echo "$n"
python $py/pthr_go_annots.py \
    --service enrich \
    --params_file $py/params/enrich.json \
    --seq_id_file $f \
    > $dir/go/pantherout/$n
done

echo "producing heatmap..."
# Run the R script that scavenges through a directory for result files and produce heatmap from it
Rscript $curdir/panther2heatmap.R \
    $dir/go/pantherout \
    $dir



