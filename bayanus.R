#conda info --envs
#conda create -n sequali_env python=3.10 -c conda-forge
#conda activate sequali_env
#conda install sequali=2.9.6 -c conda c bioforge
#sequali --outdir ~/sequali_reports '/Users/ipaz00/Downloads/BIOL4315_R/project/BIOL-4315-project/Saccharomyces bayanus.fastq'
#filtlong -h
#filtlong --min_length 120 '/Users/ipaz00/Downloads/BIOL4315_R/project/BIOL-4315-project/Saccharomyces bayanus.fastq' > bayanus_filtered.fastq
#sequali --outdir ~/sequali_reports /Users/ipaz00/Downloads/BIOL4315_R/project/BIOL-4315-project/bayanus_filtered.fastq
#flye --nano-hq /Users/ipaz00/Downloads/BIOL4315_R/project/BIOL-4315-project/bayanus_filtered.fastq --out-dir flye_output/ --threads 8 --genome-size 12m
############

#Total length:   10163590
#Fragments:      1293
#Fragments N50:  11571
#Largest frg:    43665
#Scaffolds:      0
#Mean coverage:  9

#############
#export PATH="$PATH:/Applications/Docker.app/Contents/Resources/bin/"
#type docker
#docker login
#docker run --rm \-v "/Users/ipaz00/Downloads/BIOL4315_R/project/BIOL-4315-project":/data \-w /data \ezlabgva/busco:v6.0.0_cv1 \busco -i assembly.fasta -o busco_output -l saccharomycetaceae_odb12 -m genome -c 8

#docker pull quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2

#docker run --rm \
#-v "/Users/ipaz00/Downloads/BIOL4315_R/project/BIOL-4315-project":/data \
#-w /data \
#quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2 \
#quast.py assembly.fasta \
#-r /data/bayanusreferenceassembly.fna \
#-o quast_output


#Funannotate
#make sure you're in the right directory!!!
#step 1 masking (need to do before annotation)
#docker run --rm \-v "/Users/ipaz00/Downloads/BIOL4315_R/project/BIOL-4315-project":/data \-w /data \--platform linux/amd64 \nextgenusfs/funannotate \funannotate mask -i assembly.fasta -o tr_masked_assembly.fasta
#step 2 prediction 
#docker run --rm \-v "/Users/ipaz00/Downloads/BIOL4315_R/project/BIOL-4315-project":/data \-w /data \--platform linux/amd64 \nextgenusfs/funannotate \funannotate predict -i tr_masked_assembly.fasta -o /data/output -s "Saccharomyces bayanus" 
#Mapping 560,728 proteins
#Found 294,793 preliminary alignments
#[Oct 31 09:48 PM]: Exonerate finished in 2:10:03: found 5,154 alignments
#  Feature       Specificity   Sensitivity
#nucleotides   97.1%         86.9%      
 # exons         53.8%         53.2%      
#  genes         72.0%         51.8% 

#docker run --rm \-v "/Users/ipaz00/Downloads/BIOL4315_R/project/BIOL-4315-project":/data \-w /data \--platform linux/amd64 \nextgenusfs/funannotate \funannotate annotate -i /data/output -o /data/annotation_output

#try giving it rna data -> reference genome 
#align the reference rna seqs against my assembly -> use for funannotate

library(stringr)
gffFile <- readLines('/Users/ipaz00/Downloads/BIOL4315_R/project/BIOL-4315-project/output/predict_results/Saccharomyces_bayanus.gff3')
product_lines <- str_subset(gffFile, 'product=')
products <- str_match(product_lines, "product=([^;]+);")[,2]
hypo_count <- sum(str_detect(products, regex("^hypothetical", ignore_case = TRUE)))
filtered_products <- products[!str_detect(products, regex("^hypothetical", ignore_case = TRUE))]
summary_line <- paste(hypo_count, "hypothetical proteins")
final_list <- c(summary_line, filtered_products)
