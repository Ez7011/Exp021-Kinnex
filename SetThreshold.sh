# access the data on the server
ssh -X wzhang3@curnagl.dcsr.unil.ch
cd /work/FAC/FBM/CIG/dgatfiel/smg6_kinnex/Emma_try1


# we will use collapse_isoform.fl_count file. 
# but first we need to join the classification file and this isoform.fl file
# the two files has a common first colume: pacBio id (PB1.1, PB 2.1 PB 3.1 ...) 
# I counte the lines
wc -l collapse_isoforms.flnc_count.txt 
# 800968 collapse_isoforms.flnc_count.txt
wc -l ../DG_kinnexmRNA_DEMUXED_Full_Length_Non_Concensus_HQ_transcript/pigeon.classification.txt
# 800523 ../DG_kinnexmRNA_DEMUXED_Full_Length_Non_Concensus_HQ_transcript/pigeon.classification.txt
# we see they don't match
# so I want to see which lines didn't match
# to do so, i conver the collapse_isoforms.flnc_count.txt to a tab seprated file first
sed 's/,/\t/g' collapse_isoforms.flnc_count.txt > collapse_isoforms.flnc_count_tab.txt
# then i compared the two file:
awk 'NR==FNR {a[$1]; next} !($1 in a)' ../DG_kinnexmRNA_DEMUXED_Full_Length_Non_Concensus_HQ_transcript/pigeon.classification.txt collapse_isoforms.flnc_count_tab.txt
# here is what it showed:
# id	BioSample_1	BioSample_2	BioSample_3	BioSample_4	BioSample_5	BioSample_6
# PB.13029.1	0	0	0	2	0	0
# PB.13031.1	0	1	1	0	0	1
# PB.13031.2	1	0	0	1	0	0
# PB.13031.3	0	2	0	0	0	0 ...
# i asked it to give me the count of lines that it printed out:
awk 'NR==FNR {a[$1]; next} !($1 in a)' ../DG_kinnexmRNA_DEMUXED_Full_Length_Non_Concensus_HQ_transcript/pigeon.classification.txt collapse_isoforms.flnc_count_tab.txt | wc -l
# this give: 446, this is correct

# from here I only want to include ID that exist in both files
# I first extracted the first columns (id) from the two files
cut -f1 collapse_isoforms.flnc_count_tab.txt > flnc_ids.txt
cut -f1 ../DG_kinnexmRNA_DEMUXED_Full_Length_Non_Concensus_HQ_transcript/pigeon.classification.txt > class_ids.txt
# then i sort them and got the common ids

comm -12 <(sort flnc_ids.txt) <(sort class_ids.txt) > common_ids.txt
wc -l common_ids.txt

#800522 common_ids.txt

#now i only want the collapse_isoforms.flnc_count_tab.txt to contains the transcripts in the common_ids.txt
#since i already converted everything to tab, i don't need to mind the comma
awk 'FNR==NR {ids[$1]; next} $1 in ids' common_ids.txt collapse_isoforms.flnc_count_tab.txt > filtered_flnc.txt

head filtered_flnc.txt 
# PB.1.1	2	1	0	0	0	0
# PB.2.1	2	0	0	0	0	0
# PB.3.1	4	0	1	0	0	0
# PB.3.2	0	0	1	2	0	0

# same for the class
awk ' FNR==NR {ids[$1]; next} $1 in ids' common_ids.txt ../DG_kinnexmRNA_DEMUXED_Full_Length_Non_Concensus_HQ_transcript/pigeon.classification.txt > filtered_class.txt

# filtering by threshold:
# i want to look at a summary statistics and distributions for these counts
# do this in R

module load r-light
R

df <- read.table("filtered_flnc.txt", header=FALSE, sep="\t")
data <- df[, 2:7]
summary(data)
for(i in 2:7) {
  col <- df[, i]
  hist(col, breaks=10, main=paste("Column", i), xlab=paste("Column", i))
}

write.table(summary(data), file="KinnexReads_stats.txt", sep="\t")

deciles <- apply(data, 2, function(x) round(quantile(x, probs=seq(0.1, 0.9, by=0.1)), 4))

# Append deciles to the existing file
write.table(deciles, file="KinnexReads_stats.txt", sep="\t", quote=FALSE, col.names=NA, append=TRUE)

# this table gives us most importantly the decitile:
#      V2      V3      V4      V5      V6      V7
# 0%      0       0       0       0       0       0
# 10%     0       0       0       0       0       0
# 20%     0       0       0       0       0       0
# 30%     0       0       0       0       0       0
# 40%     0       0       0       0       0       0
# 50%     1       1       1       0       1       1
# 60%     1       1       1       1       1       1
# 70%     2       2       2       1       2       2
# 80%     3       2       3       2       3       3
# 90%     6       6       6       5       6       7
# 100%    942746  865276  1080237 900677  758491  1056027

# now i joing the class and the reads and just to be safe i will sort the two files 
sort filtered_class.txt > sorted_class.txt
sort filtered_flnc.txt  >  sorted_flnc.txt
join -t $'\t' sorted_flnc.txt sorted_class.txt > joined.txt

# i want to keep all the lines that have no zero counts and at least one samples that has more than 5 reads
#   HERE IS WHERE YOU CHANGE THRESHOLDS  
awk -F'\t' '($2>5 || $3>5 || $4>5 || $5>5 || $6>5 || $7>5) && ($2*$3*$4*$5*$6*$7 !=0 ) ' joined.txt > joined_filtered.txt

# this will print out all the classes for these trasncripts:
awk -F'\t' '{print $12}' joined_filtered.txt | sort | uniq




# The further analysis is done in R see next script
