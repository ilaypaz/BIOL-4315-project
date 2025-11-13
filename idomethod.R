library(GO.db)
library(AnnotationDbi)
library(dplyr)
library(stringr)
library(tidyr)
library(readxl)

para <- read_xlsx("/Users/ilaypaz/Downloads/paradoxus/out.emapper.annotations.xlsx", skip = 2)
baya <- read_xlsx("/Users/ilaypaz/Downloads/bayanus/out.emapper.annotations.xlsx", skip = 2)

process_go <- function(df) {
  df <- as.data.frame(df)
  names(df) <- str_trim(names(df))
  go_col <- names(df)[str_detect(names(df), regex("^GO", ignore_case = TRUE))][1]
  gene_col <- names(df)[str_detect(names(df), regex("^query$", ignore_case = TRUE))][1]
  names(df)[names(df) == go_col] <- "GO_IDs"
  names(df)[names(df) == gene_col] <- "Gene"
  
  df <- df[, c("Gene", "GO_IDs"), drop=FALSE]
  df$GO_IDs <- as.character(df$GO_IDs)
  df <- tidyr::separate_rows(df, GO_IDs, sep = ",|;|\\|")
  df <- df[!is.na(df$GO_IDs) & df$GO_IDs != "", ]
  
  keys <- unique(df$GO_IDs)
  keys <- as.character(keys)
  
  # This is the only select() call now
  go_info <- AnnotationDbi::select(x = GO.db,
                                   keys = keys,
                                   columns = c("GOID", "TERM", "ONTOLOGY"),
                                   keytype = "GOID")
  
  df <- merge(df, go_info, by.x = "GO_IDs", by.y = "GOID", all.x = TRUE)
  df$ONTOLOGY <- dplyr::case_when(
    df$ONTOLOGY == "MF" ~ "Molecular Function",
    df$ONTOLOGY == "CC" ~ "Cellular Component",
    df$ONTOLOGY == "BP" ~ "Biological Process",
    TRUE ~ NA_character_
  )
  
  return(df)
}

para_long <- process_go(para)
baya_long <- process_go(baya)
library(ggplot2)

# Use para_long for plotting
para_long <- para_long %>%
  mutate(ONTOLOGY = case_when(
    ONTOLOGY == "MF" ~ "Molecular Function",
    ONTOLOGY == "CC" ~ "Cellular Component",
    ONTOLOGY == "BP" ~ "Biological Process",
    TRUE ~ ONTOLOGY
  ))

# Summarize counts of genes per ontology and term
summary_counts <- para_long %>%
  group_by(ONTOLOGY, TERM) %>%
  summarise(Gene_Count = n(), .groups = "drop")

# Bubble plot of GO term counts
ggplot(summary_counts, aes(x = ONTOLOGY, y = TERM, size = Gene_Count, color = ONTOLOGY)) +
  geom_point(alpha = 0.7) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8)) +
  labs(title = "GO Annotation Summary", y = "GO Term", x = "Ontology", size = "Gene Count")

#remove NA's and plot by category 
#find associated words and then find thoise
