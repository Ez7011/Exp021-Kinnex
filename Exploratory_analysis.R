#module load r-light

library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(DESeq2)

#load file

data <- read.csv("joined_filtered.txt", sep = "\t", stringsAsFactors = FALSE, header = FALSE)

# put the header back

header <- paste(readLines("collapse_isoforms.flnc_count_tab.txt", n = 1),readLines("../DG_kinnexmRNA_DEMUXED_Full_Length_Non_Concensus_HQ_transcript/pigeon.classification.txt", n = 1), sep = '\t')
header <- strsplit(header, "\t")[[1]]
header <- header[-8]
colnames(data) <- header

length(unique(data[,"associated_gene"]))
# here is you filter by at least one read is 5, you get 12982, if you do 10 you have 11079
length(unique(data[,"associated_transcript"]))
# if by 5, 18676, if by 10 - 15187

df_raw <- data [, c("id", "BioSample_1", "BioSample_2" , "BioSample_3" , "BioSample_4" , "BioSample_5" , "BioSample_6" , "structural_category", "associated_gene", "associated_transcript" )]

counts <- df_raw[, 2:7]
rownames(counts) <- df_raw$id
colData <- data.frame(
  row.names = colnames(counts),
  condition = c("WT", "WT", "WT", "TX", "TX", "TX")
)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = colData,
  design = ~ condition
)

dds <- DESeq(dds)

df <- as.data.frame(counts(dds, normalized = TRUE))



# exploratory analysis, we will look at the distribution of the structural category across samples

category_df <- df[,c("id", "BioSample_1", "BioSample_2" , "BioSample_3" , "BioSample_4" , "BioSample_5" , "BioSample_6" , "structural_category")]
cdf_long <- category_df %>%  pivot_longer(cols = starts_with("BioSample_"), names_to = "Sample", values_to = "Count")

colors_scheme <- brewer.pal(12, "Paired")
stacked_bar_all <- ggplot(cdf_long, aes(x = Sample, y = Count, fill = structural_category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors_scheme) +
  labs(x = "Sample", y = "Count", fill = "Category") +
  theme_minimal()
ggsave("stacked_bar_class.png", plot = stacked_bar_all, width = 8, height = 6, dpi = 300)


# note: if threshold is set to 5, you get a category genic_intron, if set to 10, this category disappears



#I want to know the top expressed gene

gene_df <- aggregate(df[, 2:7],
                       by = list(associated_gene = df$associated_gene),
                       FUN = sum)
gene_df$wt_avg <- rowMeans(gene_df[,2:4])
gene_df$tx_avg <- rowMeans(gene_df[,5:7])

gene_df <- gene_df[order(gene_df$wt_avg, decreasing = TRUE), ]

#head(gene_df,20)
#      associated_gene BioSample_1 BioSample_2 BioSample_3 BioSample_4
#11859             Trf     1002529      926578     1145575      949962
#6177          mt-Cytb      306483      278039      297219      183503
#10162            Scd1      228892      277419      236730      138644
#10300       Serpina3k      260941      242694      218668      343951
#3541              Fgb      246618      239515      218470      153211
#6183           mt-Nd4      175299      147502      168746      179540
#3553              Fgg      169070      170902      145996      110733
#6188          mt-Rnr2      168068      139277      147315       73873
#783              Apoe      162344      144170      125526      117392
#5282             Kng1      151377      142908      135267      114927
#6271            Mup22      161584      134333      126721      125690
#2462          Cyp3a11      125340      133571      139173      110334
#3781               Gc      119639      115965      108581       83425
#6272             Mup3       93431      128450      111227       29145
#2452           Cyp2e1      114144      116697       94339       94488
#3924          Gm10925       98270      102220      112178       81602
#1254             Bhmt      120711       92081       90971       89920
#10246         Selenop      104670       96180       94527       57868
#10295       Serpina1c      108456       97070       73678       89967
#6178           mt-Nd1       99362       90738       84959      127426

gene_df <- gene_df[order(gene_df$tx_avg, decreasing = TRUE), ]

# > head(gene_df,20)
# associated_gene BioSample_1 BioSample_2 BioSample_3 BioSample_4
# 11859             Trf     1002529      926578     1145575      949962
# 10300       Serpina3k      260941      242694      218668      343951
# 6271            Mup22      161584      134333      126721      125690
# 3541              Fgb      246618      239515      218470      153211
# 6177          mt-Cytb      306483      278039      297219      183503
# 10162            Scd1      228892      277419      236730      138644
# 6183           mt-Nd4      175299      147502      168746      179540
# 2462          Cyp3a11      125340      133571      139173      110334
# 3553              Fgg      169070      170902      145996      110733
# 6267            Mup19        9932        7031        1658       52126
# 2453           Cyp2f2       67629       67299       59835      168476
# 5282             Kng1      151377      142908      135267      114927
# 783              Apoe      162344      144170      125526      117392
# 439              Adh1       79285       76301       76325      143456
# 1254             Bhmt      120711       92081       90971       89920
# 6178           mt-Nd1       99362       90738       84959      127426
# 3781               Gc      119639      115965      108581       83425
# 10295       Serpina1c      108456       97070       73678       89967
# 4059          Gm28437       68596       66427       87971       85877
# 2452           Cyp2e1      114144      116697       94339       94488

# now I want to know the top expressed trasncript
# to do so, I need to first merge all the counts that have the same associated_transcript
# so first, i need to create some new names for the transcript, because of the following
# when i created a index for transcript to gene:

transcript_index_gene <- unique(df[, c("associated_transcript", "associated_gene")])
head(transcript_index_gene)

# I see this "novel" showing up, if we do direct collapse, it will merge all the novel gene
#   associated_transcript associated_gene
# 1  ENSMUST00000027727.15         Adipor1
# 2   ENSMUST00000112237.2         Adipor1
# 10                 novel         Adipor1
# 42  ENSMUST00000117836.8          Ifnar2
# 43                 novel          Ifnar2
# 47 ENSMUST00000023693.14          Ifnar2

# so i have to create a new column that id all trasncript

df$new_tp_id <- paste(df$associated_gene,df$associated_transcript, sep = "_")

transcript_df <- aggregate(df[, 2:7],
                       by = list(new_tp_id = df$new_tp_id),
                       FUN = sum)
transcript_df$wt_avg <- rowMeans(transcript_df[,2:4])
transcript_df$tx_avg <- rowMeans(transcript_df[,5:7])
transcript_df <- transcript_df[order(transcript_df$wt_avg, decreasing = TRUE), ]


# head(transcript_df,20)
#                           new_tp_id BioSample_1 BioSample_2 BioSample_3
#23280      Trf_ENSMUST00000035158.16      994186      918948     1136497
#12725   mt-Cytb_ENSMUST00000082421.1      306294      277816      296971
#19819 Serpina3k_ENSMUST00000043058.5      258233      240042      216199
#7463        Fgb_ENSMUST00000048246.5      246383      239264      218247
#19492      Scd1_ENSMUST00000041331.4      216920      262384      224329
#12733    mt-Nd4_ENSMUST00000082414.1      175248      147464      168676
#12740   mt-Rnr2_ENSMUST00000082390.1      168068      139277      147315
#1636       Apoe_ENSMUST00000174064.9      157765      139750      121476
#7487       Fgg_ENSMUST00000048486.13      147669      145808      122724
#5178    Cyp3a11_ENSMUST00000035918.8      124199      131972      137577
#12931     Mup22_ENSMUST00000211875.2      137442      114399      107026
#7972        Gc_ENSMUST00000049209.13      118526      114990      107676
#10851     Kng1_ENSMUST00000039492.14      116905      110946      104348
#5155     Cyp2e1_ENSMUST00000026552.9      113546      116028       93823
#8255    Gm10925_ENSMUST00000189941.2       98269      102217      112172
#2654       Bhmt_ENSMUST00000099309.6      120135       91595       90569
#19809 Serpina1c_ENSMUST00000074051.6      107600       96377       73118
#19699  Selenop_ENSMUST00000159216.10       97877       89045       88593
#12727    mt-Nd1_ENSMUST00000082392.1       99359       90731       84955
#19811 Serpina1d_ENSMUST00000078869.6      103683       89696       69642


transcript_df <- transcript_df[order(transcript_df$tx_avg, decreasing = TRUE), ]



# head(transcript_df,20)
# new_tp_id BioSample_1 BioSample_2 BioSample_3
# 23280      Trf_ENSMUST00000035158.16      994186      918948     1136497
# 19819 Serpina3k_ENSMUST00000043058.5      258233      240042      216199
# 7463        Fgb_ENSMUST00000048246.5      246383      239264      218247
# 12725   mt-Cytb_ENSMUST00000082421.1      306294      277816      296971
# 12931     Mup22_ENSMUST00000211875.2      137442      114399      107026
# 19492      Scd1_ENSMUST00000041331.4      216920      262384      224329
# 12733    mt-Nd4_ENSMUST00000082414.1      175248      147464      168676
# 5178    Cyp3a11_ENSMUST00000035918.8      124199      131972      137577
# 7487       Fgg_ENSMUST00000048486.13      147669      145808      122724
# 5159    Cyp2f2_ENSMUST00000003100.10       66565       66215       58849
# 12921     Mup19_ENSMUST00000080606.9        9031        6361        1418
# 1636       Apoe_ENSMUST00000174064.9      157765      139750      121476
# 875       Adh1_ENSMUST00000004232.10       78017       75209       75234
# 2654       Bhmt_ENSMUST00000099309.6      120135       91595       90569
# 12727    mt-Nd1_ENSMUST00000082392.1       99359       90731       84955
# 10851     Kng1_ENSMUST00000039492.14      116905      110946      104348
# 7972        Gc_ENSMUST00000049209.13      118526      114990      107676
# 19809 Serpina1c_ENSMUST00000074051.6      107600       96377       73118
# 8419    Gm28437_ENSMUST00000190277.2       68596       66427       87971
# 12740   mt-Rnr2_ENSMUST00000082390.1      168068      139277      147315


## now I want to move on the differential expression


df$p_value <- apply(df, 1, function(row) {
  group1 <- as.numeric(row[2:4])
  group2 <- as.numeric(row[5:7])
  
  if(length(unique(group1)) == 1 || length(unique(group2)) == 1) {
    return(NA)
  } else {
    return(t.test(group1, group2)$p.value)
  }
})

df$log2FC <- log2(rowMeans(df[, c("BioSample_4", "BioSample_5", "BioSample_6")])/rowMeans(df[, c("BioSample_1", "BioSample_2", "BioSample_3")]))


# adjusted pval
# first i need to remove all the rows that has pval = NA  or =1 becuase that skews the p adj

df_clean <- df[!is.na(df$p_val), ]
df_clean <- df <- df[df$p_value != 1, ]
df_clean$padj <- p.adjust(df_clean$p_value, method = "fdr")

# volcano plot

volcano_df <- data.frame(
  gene = df_clean$id,
  log2FC = df_clean$log2FC,
  p_value = df_clean$p_value,
  padj = df_clean$padj)

n_up <- volcano_df %>%
  filter(log2FC > 1, p_value < 0.05) %>%
  nrow()

n_down <- volcano_df %>%
  filter(log2FC < -1, p_value < 0.05) %>%
  nrow()


volcano_df$significant <- (abs(volcano_df$log2FC) > 1) & (volcano_df$p_value < 0.05)
volcano <- ggplot(volcano_df, aes(x = log2FC, y = -log10(p_value))) +
  geom_point(aes(color = significant)) +
  scale_y_continuous(limits = c(0, 8))+
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="blue") +
  geom_vline(xintercept = c(-1,1), linetype="dashed", color="blue") +
annotate(
    "text",
    x = 7, 
    y = 5,
    label = paste0("n = ", n_up),
    hjust = 1,
    size = 4
  ) +
  annotate(
    "text",
    x = -7,
    y = 5, 
    label = paste0("n = ", n_down),
    hjust = 0,
    size = 4
  ) +
  theme_minimal()
ggsave("volcano.png", plot = volcano, width = 8, height = 6, dpi = 300)


## P Adjusted
# use df_clean

volcano_df <- data.frame(
  gene = df_clean$id,
  log2FC = df_clean$log2FC,
  p_value = df_clean$p_value,
  padj = df_clean$padj)


n_up <- volcano_df %>%
  filter(log2FC > 1, padj < 0.05) %>%
  nrow()

n_down <- volcano_df %>%
  filter(log2FC < -1, padj < 0.05) %>%
  nrow()


volcano_df$significant <- (abs(volcano_df$log2FC) > 1) & (volcano_df$padj < 0.05)
volcano <- ggplot(volcano_df, aes(x = log2FC, y = -log10(padj))) +
  geom_point(aes(color = significant)) +
  scale_y_continuous(limits = c(0, 8))+
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="blue") +
  geom_vline(xintercept = c(-1,1), linetype="dashed", color="blue") +
annotate(
    "text",
    x = 7, 
    y = 5,
    label = paste0("n = ", n_up),
    hjust = 1,
    size = 4
  ) +
  annotate(
    "text",
    x = -7,
    y = 5, 
    label = paste0("n = ", n_down),
    hjust = 0,
    size = 4
  ) +
  theme_minimal()


ggsave("volcano_adj.png", plot = volcano, width = 8, height = 6, dpi = 300)

# Expanding on the significantly changed data (i.e red dots on volcano plot right)
df_sig <- df[!is.na(df$padj) & df$padj < 0.05 & (df$log2FC > 1 |df$log2FC < -1), ]

dim(df_sig)
> dim(df_sig)
[1] 4776   13

# This is an extra step (maybe unecessary) to keep the variance low - meaning i want wt and tx group to uniformly express this transcript
df_sig_crude <- df_sig
df_sig$wt_avg <- rowMeans(df_sig[,2:4])
df_sig$tx_avg <- rowMeans(df_sig[,5:7])
df_sig$wt_sd <- apply(df_sig[, 2:4], 1, sd)
df_sig$tx_sd <- apply(df_sig[, 5:7], 1, sd)
df_sig$wt_cv <- apply(df_sig[, 2:4], 1, sd)/df_sig$wt_avg
df_sig$tx_cv <- apply(df_sig[, 5:7], 1, sd)/df_sig$tx_avg

#df_sig <- df_sig[df_sig$tx_cv<0.3 & df_sig$wt_cv<0.3,]
#dim(df_sig)
#[1] 1379   19



#now a stacked bar plot


df_sig_long <- df_sig %>%
  pivot_longer(
    cols = 2:7,
    names_to = "Sample",
    values_to = "Count"
  )

stacked_sig <- ggplot(df_sig_long, aes(x = Sample, y = Count, fill = structural_category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors_scheme) +
  labs(x = "Sample", y = "Count", fill = "Category") +
  theme_minimal()

ggsave("stacked_sig.png", plot = stacked_sig, width = 8, height = 6, dpi = 300)

# without the variance control:


df_sig_long <- df_sig_crude %>%
  pivot_longer(
    cols = 2:7,
    names_to = "Sample",
    values_to = "Count"
  )

stacked_sig <- ggplot(df_sig_long, aes(x = Sample, y = Count, fill = structural_category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors_scheme) +
  labs(x = "Sample", y = "Count", fill = "Category") +
  theme_minimal()

ggsave("stacked_sig_crude.png", plot = stacked_sig, width = 8, height = 6, dpi = 300)

# we want to look at specifically the up reg'ed ones 
colors_scheme <- brewer.pal(12, "Paired")
df_sig_up <- df[!is.na(df$padj) & df$padj < 0.05 & (df$log2FC > 1), ]

df_sig_long <- df_sig_up %>%
  pivot_longer(
    cols = 2:7,
    names_to = "Sample",
    values_to = "Count"
  )

stacked_sig <- ggplot(df_sig_long, aes(x = Sample, y = Count, fill = structural_category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors_scheme) +
  labs(x = "Sample", y = "Count", fill = "Category") +
  theme_minimal()

ggsave("stacked_sig_up.png", plot = stacked_sig, width = 8, height = 6, dpi = 300)


# here is a normalized to 100% 
stacked_sig <- ggplot(df_sig_long, aes(x = Sample, y = Count, fill = structural_category)) +
  geom_bar(stat = "identity",position = "fill") +
  scale_fill_manual(values = colors_scheme) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Sample", y = "Percentage", fill = "Category") +
  theme_minimal()

ggsave("stacked_sig_up_norm.png", plot = stacked_sig, width = 8, height = 6, dpi = 300)


# and down

df_sig_down <- df[!is.na(df$padj) & df$padj < 0.05 & (df$log2FC < -1), ]

df_sig_long <- df_sig_down %>%
  pivot_longer(
    cols = 2:7,
    names_to = "Sample",
    values_to = "Count"
  )

stacked_sig <- ggplot(df_sig_long, aes(x = Sample, y = Count, fill = structural_category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors_scheme) +
  labs(x = "Sample", y = "Count", fill = "Category") +
  theme_minimal()

ggsave("stacked_sig_down.png", plot = stacked_sig, width = 8, height = 6, dpi = 300)


stacked_sig <- ggplot(df_sig_long, aes(x = Sample, y = Count, fill = structural_category)) +
  geom_bar(stat = "identity",position = "fill") +
  scale_fill_manual(values = colors_scheme) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Sample", y = "Percentage", fill = "Category") +
  theme_minimal()

ggsave("stacked_sig_down_norm.png", plot = stacked_sig, width = 8, height = 6, dpi = 300)

