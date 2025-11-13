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
# Combine both species for joint analysis
combined_long <- bind_rows(
  para_long %>% mutate(Species = "Paradoxus"),
  baya_long %>% mutate(Species = "Bayanus")
)

# Keep only relevant ontologies
combined_long <- combined_long %>%
  filter(ONTOLOGY %in% c("Biological Process", "Molecular Function"))

# Define keyword groups
fermentation_keywords <- c("ferment", "glycolysis", "ethanol", "alcohol", "anaerobic", "carbon metabolism")
resistance_keywords   <- c("resistance", "drug", "toxin", "stress", "detox", "defense")

# Categorize GO terms based on TERM text
combined_long <- combined_long %>%
  mutate(
    Pathway = case_when(
      str_detect(str_to_lower(TERM), str_c(fermentation_keywords, collapse = "|")) ~ "Fermentation",
      str_detect(str_to_lower(TERM), str_c(resistance_keywords, collapse = "|")) ~ "Resistance",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Pathway))

summary_counts <- combined_long %>%
  group_by(Species, ONTOLOGY, Pathway, TERM) %>%
  summarise(Gene_Count = n(), .groups = "drop") %>%
  group_by(Species, Pathway) %>%
  mutate(Percent = Gene_Count / sum(Gene_Count))

ggplot(summary_counts, aes(x = Pathway, y = TERM, size = Percent, color = ONTOLOGY)) +
  geom_point(alpha = 0.7, position = position_jitter(width = 0.2, height = 0.3)) +
  facet_wrap(~ Species, scales = "free_y") +
  scale_size_continuous(range = c(3, 12), labels = scales::percent) +
  theme_bw(base_size = 12) +
  theme(axis.text.y = element_text(size = 8),
        legend.position = "bottom") +
  labs(
    title = "Fermentation and Resistance Pathways in GO Terms",
    x = "Pathway Category",
    y = "GO Term",
    color = "Ontology"
  ) +
  guides(size = "none")

#try ggsave
#try normalizing it by percent 
#get rid of gene count legend on the bottom 
#why are the bubbles overlapping

