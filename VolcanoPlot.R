#scp -r wzhang3@curnagl.dcsr.unil.ch:/work/FAC/FBM/CIG/dgatfiel//smg6_kinnex/Emma_im/Phkg2.csv /Users/wzhang3/Desktop/Smg6/Exp\ 021\ Kinnex\ Data 
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(DESeq2)

setwd("/Users/wzhang3/Desktop/Smg6/Exp 021 Kinnex Data")
#load file

data <- read.csv("joined_filtered.txt", sep = "\t", stringsAsFactors = FALSE, header = FALSE)


# put the header back

header <- paste(readLines("collapse_isoforms.flnc_count_tab.txt", n = 1),readLines("pigeon.classification.txt", n = 1), sep = '\t')
header <- strsplit(header, "\t")[[1]]
header <- header[-8]
colnames(data) <- header

counts <- data [, c("BioSample_1", "BioSample_2" , "BioSample_3" , "BioSample_4" , "BioSample_5" , "BioSample_6" )]


condition = c("c","c","c","t","t","t")    
condition = data.frame(condition)

dds <- DESeqDataSetFromMatrix(countData=counts,colData=condition,tidy=FALSE,design= ~condition)
dds <- DESeq(dds)

res = results(dds, contrast=c("condition","t","c"))

df <- cbind (as.data.frame(res),data)

write.csv(df,"df.csv")

#df_filtered <- df %>%
#  filter(
#    (rowMeans(across(8:10)) > 5) |
#      (rowMeans(across(11:13)) > 5)
# 3  )

volcano_df <- data.frame(
  id = df$id,
  category = df$structural_category,
  log2FC = df$log2FoldChange,
  p_value = df$pvalue,
  padj = df$padj)



n_up <- volcano_df %>%
  filter(log2FC > 1, padj< 0.05) %>%
  nrow()

n_down <- volcano_df %>%
  filter(log2FC < -1, padj < 0.05) %>%
  nrow()


volcano_df$significant <- (abs(volcano_df$log2FC) > 1) & (volcano_df$padj < 0.05)
volcano <- ggplot(volcano_df, aes(x = log2FC, y = -log10(padj))) +
  geom_point(aes(color = significant)) +
  
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="blue") +
  geom_vline(xintercept = c(-1,1), linetype="dashed", color="blue") +
  annotate(
    "text",
    x = 7, 
    y = 42,
    label = paste0("n = ", n_up),
    hjust = 1,
    size = 4
  ) +
  annotate(
    "text",
    x = -7,
    y = 42, 
    label = paste0("n = ", n_down),
    hjust = 0,
    size = 4
  ) +
  theme_minimal()
ggsave("volcano_normalized_filtered.pdf", plot = volcano, width = 8, height = 6, dpi = 300)

# volcano plot coloured by categories

struct_colors <- c(
  'full-splice_match' = '#1f97b4',
  'incomplete-splice_match' = "orange",
  'novel_not_in_catalog' = '#e377c2',
  'novel_in_catalog' = '#bcbd22',
  'intergenic' = '#9060bd',
  'genic' = '#8c564b',
  'moreJunctions' = '#d62728',
  'fusion' = '#7f7f7f',
  'antisense' = '#1aa76c'
)

volcano_cat <- ggplot(volcano_df, aes(x = log2FC, y = -log10(padj))) +
  geom_point(aes(color = category)) +
  scale_color_manual(values = struct_colors) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="blue") +
  geom_vline(xintercept = c(-1,1), linetype="dashed", color="blue") +
  annotate(
    "text",
    x = 7, 
    y = 42,
    label = paste0("n = ", n_up),
    hjust = 1,
    size = 4
  ) +
  annotate(
    "text",
    x = -7,
    y = 42, 
    label = paste0("n = ", n_down),
    hjust = 0,
    size = 4
  ) +
  theme_minimal()
ggsave("volcano_normalized_cat.png", plot = volcano_cat, width = 8, height = 6, dpi = 300)


# volcano plot category ditribution
up <- df_filtered %>%
  filter(log2FoldChange > 1, padj< 0.05) 

down <- df_filtered %>%
  filter(log2FoldChange < -1, padj < 0.05)


total_dis <- df_filtered %>% dplyr::count(structural_category) %>% mutate(group = "Total")
up_dis <- up %>% dplyr::count(structural_category) %>% mutate(group = "Upregulated")
down_dis <- down %>% dplyr::count(structural_category) %>% mutate(group = "Downregulated")

distribution <- ggplot(down_dis,
       aes(x = group, y = n, fill = structural_category)) +
  geom_bar(stat = "identity", position = "fill")+
  scale_fill_manual(values = struct_colors)+
  coord_flip() +
  labs(x = "", y = "Proportion") +
  theme_minimal() +
  theme(aspect.ratio = 0.2,legend.position = "none")
distribution
ggsave("down_distribution.pdf", plot = distribution, width = 8, height = 6, dpi = 300)

bar_distribution <- bind_rows (total_dis,up_dis,down_dis)

distribution <- ggplot(bar_distribution,
       aes(x = group, y = n, fill = structural_category)) +
  geom_bar(stat = "identity", position = "fill")+
  scale_fill_manual(values = struct_colors)+
  coord_flip() +
  labs(x = "", y = "Proportion") +
  theme_minimal() +
  theme(aspect.ratio = 0.2,legend.title = element_text(size = 12), legend.text = element_text(size = 12)) 
ggsave("normalized_distribution.pdf", plot = distribution, width = 8, height = 6, dpi = 300)
