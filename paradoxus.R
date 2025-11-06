#conda info --envs
#conda create -n sequali_env python=3.10 -c conda-forge
#conda activate sequali_env
#conda install sequali=2.9.6 -c conda c bioforge
#sequali --outdir ~/sequali_reports '/Users/ipaz00/Downloads/BIOL4315_R/project/BIOL-4315-project/saccharomycesParadoxus.fastq'
#filtlong -h
#filtlong --min_length 140 '/Users/ipaz00/Downloads/BIOL4315_R/project/BIOL-4315-project/saccharomycesParadoxus.fastq' > paradoxus_filtered.fastq
#sequali --outdir ~/sequali_reports /Users/ipaz00/Downloads/BIOL4315_R/project/BIOL-4315-project/paradoxus_filtered.fastq
#flye --nano-hq /Users/ipaz00/Downloads/BIOL4315_R/project/BIOL-4315-project/paradoxus_filtered.fastq --out-dir flye_outputP/ --threads 8 --genome-size 12m
############

#Total length:   8868681
#Fragments:      729
#Fragments N50:  17034
#Largest frg:    71455
#Scaffolds:      0
#Mean coverage:  7

#############
#export PATH="$PATH:/Applications/Docker.app/Contents/Resources/bin/"
#type docker
#docker login
#docker run --rm \-v "/Users/ipaz00/Downloads/BIOL4315_R/project/BIOL-4315-project":/data \-w /data \ezlabgva/busco:v6.0.0_cv1 \busco -i assemblyP.fasta -o busco_outputP -l saccharomycetaceae_odb12 -m genome -c 8

#docker pull quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2

#docker run --rm \
#-v "/Users/ipaz00/Downloads/BIOL4315_R/project/BIOL-4315-project":/data \
#-w /data \
#quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2 \
#quast.py assemblyP.fasta \
#-r /data/paradoxusreferenceassembly.fna \
#-o quastP_output


#Funannotate
#make sure you're in the right directory!!!
#step 1 masking (need to do before annotation)
#docker run --rm \-v "/Users/ipaz00/Downloads/BIOL4315_R/project/BIOL-4315-project":/data \-w /data \--platform linux/amd64 \nextgenusfs/funannotate \funannotate mask -i assemblyP.fasta -o tr_masked_assemblyP.fasta
#step 2 prediction 
#docker run --rm \-v "/Users/ipaz00/Downloads/BIOL4315_R/project/BIOL-4315-project":/data \-w /data \--platform linux/amd64 \nextgenusfs/funannotate \funannotate predict -i tr_masked_assemblyP.fasta -o /data/outputP -s "Saccharomyces paradoxus" 
#Mapping 560,728 proteins
#Found 255,355 preliminary alignments
#[Oct 31 09:48 PM]: Exonerate finished in 2:10:03: found 5,154 alignments
#  Feature       Specificity   Sensitivity
#nucleotides   97.1%         86.9%      
# exons         53.8%         53.2%      
#  genes         72.0%         51.8% 

#docker run --rm \-v "/Users/ipaz00/Downloads/BIOL4315_R/project/BIOL-4315-project":/data \-w /data \--platform linux/amd64 \nextgenusfs/funannotate \funannotate annotate -i /data/outputP -o /data/annotation_outputP

#try giving it rna data -> reference genome 
#align the reference rna seqs against my assembly -> use for funannotate

library(readxl)
library(ggplot2)
library(dplyr)
library(scales)  # for percent formatting

# File paths
file_path1 <- "/Users/ipaz00/Downloads/paradoxus eggnog/out.emapper.annotations.xlsx"
file_path2 <- "/Users/ipaz00/Downloads/bayanus eggnog/out.emapper.annotations.xlsx"

# Function to process a file and return relative abundance per category
process_cog <- function(file_path, genome_name) {
  data <- read_excel(file_path, skip = 2)
  cog_data <- data$COG_category
  cog_data <- cog_data[!is.na(cog_data)]
  cog_data <- trimws(cog_data)
  cog_data <- cog_data[nchar(cog_data) == 1]
  
  category <- case_when(
    cog_data %in% c("C", "G") ~ "Fermentation",
    cog_data %in% c("V", "M", "W", "P", "Q") ~ "Resistance",
    TRUE ~ "Other"
  )
  
  category_counts <- as.data.frame(table(category))
  colnames(category_counts) <- c("Function", "Count")
  category_counts <- category_counts %>%
    mutate(RelAbundance = Count / sum(Count),
           Genome = genome_name)  # add genome column
  
  category_counts$Function <- factor(category_counts$Function, levels = c("Fermentation", "Resistance", "Other"))
  
  return(category_counts)
}

# Process both genomes
data1 <- process_cog(file_path1, "Paradoxus")
data2 <- process_cog(file_path2, "Bayanus")

# Combine into one data frame
combined_data <- bind_rows(data1, data2)

# Repeat rows so bubbles appear multiple times (fake category)
combined_data <- combined_data %>%
  slice(rep(1:n(), each = 6)) %>%  # repeat each row 6 times
  mutate(FakeY = 1)  # all on the same y-axis value

# Color-blind friendly palette
cb_palette <- c("Fermentation" = "#0072B2",   # blue
                "Resistance" = "#009E73",     # green
                "Other" = "#D55E00")          # reddish-orange

# Bubble plot with all bubbles on the same horizontal line
ggplot(combined_data, aes(x = Function, y = FakeY, size = RelAbundance, color = Function)) +
  geom_point(alpha = 0.7) +
  geom_text(aes(label = scales::percent(RelAbundance, accuracy = 0.1)), 
            vjust = -1.5, size = 5) +
  scale_size(range = c(5, 15)) +
  scale_color_manual(values = cb_palette) +
  labs(title = "COG Category Relative Abundance by Genome",
       x = NULL,
       y = NULL) +  # remove y-axis label
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),
    axis.text.y = element_blank(),  # remove y-axis text
    axis.ticks.y = element_blank(), # remove y-axis ticks
    axis.ticks.x = element_blank(), # optional: remove x-axis ticks
    plot.title = element_text(size = 20, face = "bold")
  ) +
  facet_wrap(~ Genome)
