#!/bin/bash
mkdir -p ../Emma_im/

#!/bin/bash

path="/work/FAC/FBM/CIG/dgatfiel/smg6_kinnex/DG_kinnexmRNA_DEMUXED_Full_Length_Non_Concensus_HQ_transcript"

EmmaPath="/work/FAC/FBM/CIG/dgatfiel/smg6_kinnex/Emma_im/"

file1="collapse_isoforms.flnc_count.txt"
file2="pigeon.classification.txt"

# Copy files
cp "$path/$file1" "$EmmaPath/"
cp "$path/$file2" "$EmmaPath/"

# Count lines
lines1=$(wc -l < "$path/$file1")
lines2=$(wc -l < "$path/$file2")

echo "In $file1, there are $lines1 lines."
echo "In $file2, there are $lines2 lines."

sed 's/,/\t/g' /work/FAC/FBM/CIG/dgatfiel/smg6_kinnex/DG_kinnexmRNA_DEMUXED_Full_Length_Non_Concensus_HQ_transcript/collapse_isoforms.flnc_count.txt > /work/FAC/FBM/CIG/dgatfiel/smg6_kinnex/Emma_im/collapse_isoforms.flnc_count_tab.txt

cut -f1 "$EmmaPath/$file1" > "$EmmaPath/flnc_ids.txt"
cut -f1 "$EmmaPath/$file2" > "$EmmaPath/class_ids.txt"
echo "PB IDs from reads and classifications are stored in 'flnc_ids.txt' and 'class_ids.txt'.”

comm -12 <(sort "$EmmaPath/flnc_ids.txt") <(sort "$EmmaPath/class_ids.txt") > "$EmmaPath/common_ids.txt"

echo "two files are sorted and filtered by common PB IDs stored in 'common_ids.txt'."

awk 'FNR==NR {ids[$1]; next} $1 in ids' "$EmmaPath/common_ids.txt" "$EmmaPath/$file1" > "$EmmaPath/filtered_flnc.txt" 
awk 'FNR==NR {ids[$1]; next} $1 in ids' "$EmmaPath/common_ids.txt" "$EmmaPath/$file2” > "$EmmaPath/filtered_class.txt" 

echo "in classification and flnc, only lines that contains the common ids are extracted and stored in the filtered_class.txt and filtered_flnc.txt"

sort "$EmmaPath/filtered_flnc.txt"  > "$EmmaPath/sorted_flnc.txt" 
sort "$EmmaPath/filtered_class.txt"  > "$EmmaPath/sorted_class.txt" 
join -t $'\t' "$EmmaPath/sorted_flnc.txt"  "$EmmaPath/sorted_class.txt"  > "$EmmaPath/joined.txt" 

echo "The sorted joined counts and classification are stored in joined.txt" 



