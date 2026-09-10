library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(readr)

# ============================================================
DATA_DIR <- "C:/Users/user/OneDrive - CDC/Project/ww/data"
OUTPUT_DIR <- file.path(DATA_DIR, "version3")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUTPUT_DIR)
LINK_TABLE_PATH <- file.path(DATA_DIR, "SampleID_SRA.txt")
ACCESSION_NAME_PATH <- file.path(DATA_DIR, "accession_to_name.txt")
FASTP_QC_PATH <- file.path(DATA_DIR, "fastp_qc_summary_wide.tsv")
HUMAN_REMOVAL_PATH <- file.path(DATA_DIR, "human.reads.aligned.txt")
RESULTS_DIR <- file.path(DATA_DIR, "65")
SARS_ACC <- "NC_045512.2"
# Sample-level eligibility criterion (inclusive).
BREADTH_CUTOFF <- 0.5
# Call-level high-frequency criterion (strictly greater than 50%).
HIGH_FREQ_CUTOFF <- 0.5
# Q3/Q4.2 row-level detection criterion (inclusive). A virus group is detected
# in a library if any matching accession/segment reaches this threshold.
VIRUS_DETECTION_MIN_READ_PAIRS <- 50L
# ============================================================
# 1. Sample metadata and common input tables
# ============================================================
link_table <- read.table(
  LINK_TABLE_PATH,
  header = TRUE,
  sep = "\t",
  strip.white = TRUE
) %>%
  mutate(
    SRA = trimws(SRA),
    Enriched = if_else(
      str_detect(SampleID, "_INF_unenriched$"),
      "unenriched",
      "enriched"
    ),
    Plant = str_extract(SampleID, "^[A-Z]+"),
    Date = as.Date(
      str_extract(
        SampleID,
        "(?<=^[A-Z]{2,4}_)[0-9_]+(?=(_INF_unenriched)?$)"
      ),
      "%m_%d_%Y"
    ),
    Days_from_first_sample = as.numeric(Date - min(Date, na.rm = TRUE))
  )
sample_order <- link_table %>%
  arrange(Plant, Date, desc(Enriched)) %>%
  pull(SampleID)
# ============================================================
# 2. Q1: read quality and post-trim duplication
# ============================================================
qc_data <- read_tsv(FASTP_QC_PATH, show_col_types = FALSE) %>%
  left_join(
    link_table %>% select(SampleID, SRA),
    by = c(SRR = "SRA")
  ) %>%
  mutate(
    SampleID = factor(SampleID, levels = sample_order)
  )

p_q30 <- ggplot(qc_data) +
  geom_segment(
    aes(
      x = q30_raw,
      xend = q30_trim,
      y = SampleID,
      yend = SampleID
    ),
    linetype = "dashed",
    color = "gray50"
  ) +
  geom_point(
    aes(x = q30_raw, y = SampleID, color = "Raw"),
    size = 2
  ) +
  geom_point(
    aes(x = q30_trim, y = SampleID, color = "Trimmed"),
    size = 2
  ) +
  scale_color_manual(
    values = c(Raw = "darkorange", Trimmed = "forestgreen"),
    name = NULL
  ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Q30 rate (%)",
    y = NULL,
    title = "Q1: Read quality before vs after trimming"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 6),
    legend.position = "top"
  )

if (interactive()) print(p_q30)

ggsave(
  "q1.1_q30_before_after_trim.png",
  p_q30,
  width = 6,
  height = 12,
  dpi = 300
)

p_dup <- ggplot(
  qc_data,
  aes(x = dup_trim, y = SampleID)
) +
  geom_col(fill = "steelblue") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    x = "Duplication rate (after trimming)",
    y = NULL,
    title = "Q1: Post-trim duplication rate per sample"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6))

if (interactive()) print(p_dup)

ggsave(
  "q1.2_duplicate_rate_after_trim.png",
  p_dup,
  width = 6,
  height = 12,
  dpi = 300
)

# REVIEW - STEP 1: QUALITY CONTROL (DUPLICATION AND Q30)
# Q: Why do high duplication and low Q30 together threaten wastewater-sequencing interpretation?
# A. They only delay results.
# B. Duplicates inflate abundance, while low quality creates false sequence differences.
# C. They affect only GC content.
# D. High duplication always means human contamination.
# Answer: B.
# Insight: Duplication distorts abundance; low Q30 distorts identity and variant calls.
# ============================================================
# 3. Q2: human-read removal
# ============================================================
human_data <- read.table(
  HUMAN_REMOVAL_PATH,
  header = TRUE,
  sep = "\t"
) %>%
  mutate(
    human_aligned_rate = as.numeric(str_remove(human_aligned_rate, "%"))
  ) %>%
  left_join(
    link_table %>% select(SampleID, SRA),
    by = "SRA"
  ) %>%
  mutate(
    SampleID = factor(SampleID, levels = sample_order)
  )

p_human <- ggplot(
  human_data,
  aes(x = human_aligned_rate, y = SampleID)
) +
  geom_col(fill = "firebrick") +
  scale_x_continuous(limits = c(0, 100)) +
  labs(
    x = "Human-aligned rate (%)",
    y = NULL,
    title = "Q2: Fraction of reads mapping to human genome (removed)"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6))

if (interactive()) print(p_human)

ggsave(
  "q2_human_removal_rate.png",
  p_human,
  width = 6,
  height = 12,
  dpi = 300
)

# REVIEW - STEP 2: HOST FILTERING (BACKGROUND AND SENSITIVITY)
# Q: Why does removing human reads improve sensitivity for low-abundance microbes?
# A. Human reads are always viral.
# B. Host reads consume finite depth and bury rare microbial signals.
# C. Human reads stop the aligner.
# D. Host removal reveals taxa by changing GC content.
# Answer: B.
# Insight: Host background dilutes the read share available to rare targets.
# ============================================================
# 4. Virus-name lookup used by Q3/Q4
# ============================================================
name_lookup <- read.table(
  ACCESSION_NAME_PATH,
  sep = "\t",
  col.names = c("genome", "virus_name"),
  quote = "",
  fill = TRUE
)
clean_name <- function(name) {
  name <- str_remove(name, ",?\\s*complete genome.*$")
  name <- str_remove(name, ",?\\s*complete cds.*$")
  name <- str_remove(name, ",?\\s*complete sequence.*$")
  trimws(name)
}
# Family, genus, and species assignments for the ten RefSeq accessions
# displayed in the current Q3/Q4 legend.
# Verified against NCBI Nucleotide/Taxonomy records on 2026-09-10.
virus_taxonomy <- tibble::tribble(
  ~genome, ~legend_name, ~family, ~genus, ~species,
  "NC_028478.1", "Tomato brown rugose fruit virus", "Virgaviridae", "Tobamovirus", "Tobamovirus fructirugosum",
  "NC_003630.1", "Pepper mild mottle virus", "Virgaviridae", "Tobamovirus", "Tobamovirus capsici",
  "NC_001801.1", "Cucumber green mottle mosaic virus", "Virgaviridae", "Tobamovirus", "Tobamovirus viridimaculae",
  "NC_002692.1", "Tomato mosaic virus", "Virgaviridae", "Tobamovirus", "Tobamovirus tomatotessellati",
  "NC_022230.1", "Tomato mottle mosaic virus", "Virgaviridae", "Tobamovirus", "Tobamovirus maculatessellati",
  "NC_001556.1", "Tobacco mild green mosaic virus", "Virgaviridae", "Tobamovirus", "Tobamovirus mititessellati",
  "NC_030229.1", "Tropical soda apple mosaic virus", "Virgaviridae", "Tobamovirus", "Tobamovirus tropici",
  "NC_001504.1", "Melon necrotic spot virus", "Tombusviridae", "Gammacarmovirus", "Gammacarmovirus melonis",
  "NC_074405.1", "ssRNA phage SRR5466337_3", "Steitzviridae", "Kinglevirus", "Kinglevirus lutadaptatum",
  "NC_045512.2", "Severe acute respiratory syndrome coronavirus 2", "Coronaviridae",
    "Betacoronavirus", "Betacoronavirus pandemicum"
)

genome_info_files <- list.files(
  RESULTS_DIR,
  pattern = "_genome_info\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

snv_files <- list.files(
  RESULTS_DIR,
  pattern = "_SNVs\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)
if (length(genome_info_files) == 0) {
  stop("No genome_info.tsv files found. Check RESULTS_DIR.")
}
if (length(snv_files) == 0) {
  stop("No SNVs.tsv files found. Check RESULTS_DIR.")
}
combined <- lapply(genome_info_files, function(file_path) {
  genome_info <- read.table(
    file_path,
    header = TRUE,
    sep = "\t"
  )
  genome_info$SRA <- str_extract(
    basename(file_path),
    "SRR[0-9]+"
  )
  genome_info
}) %>%
  bind_rows() %>%
  left_join(link_table, by = "SRA") %>%
  left_join(name_lookup, by = "genome")
# ============================================================
# 5. Q3/Q4: paired enriched/unenriched viral composition
# ============================================================
pairing_check <- link_table %>%
  group_by(Plant, Date) %>%
  summarise(
    has_enriched = "enriched" %in% Enriched,
    has_unenriched = "unenriched" %in% Enriched,
    .groups = "drop"
  ) %>%
  mutate(
    status = case_when(
      has_enriched & has_unenriched ~ "paired",
      has_enriched ~ "enriched_only",
      has_unenriched ~ "unenriched_only"
    )
  )

pairing_summary <- pairing_check %>%
  count(Plant, status) %>%
  pivot_wider(
    names_from = status,
    values_from = n,
    values_fill = 0
  )

write.table(
  pairing_summary,
  "pairing_summary_by_plant.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

paired_dates <- pairing_check %>%
  filter(status == "paired") %>%
  select(Plant, Date)

if (nrow(paired_dates) == 0) {
  stop("No paired enriched/unenriched dates found.")
}

combined_paired <- combined %>%
  semi_join(paired_dates, by = c("Plant", "Date"))

top10 <- combined_paired %>%
  group_by(genome, virus_name) %>%
  summarise(
    total_coverage = sum(coverage, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_coverage)) %>%
  slice(1:10) %>%
  mutate(virus_name = clean_name(virus_name))

sars_row <- combined_paired %>%
  filter(genome == SARS_ACC) %>%
  group_by(genome, virus_name) %>%
  summarise(
    total_coverage = sum(coverage, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(virus_name = clean_name(virus_name))

if (!SARS_ACC %in% top10$genome && nrow(sars_row) > 0) {
  top10 <- bind_rows(top10 %>% slice(1:9), sars_row)
}

top10 <- top10 %>%
  left_join(virus_taxonomy, by = "genome") %>%
  mutate(
    legend_name = coalesce(legend_name, virus_name),
    taxonomy_label = if_else(
      !is.na(family),
      paste0(
        legend_name,
        "\n",
        "Family: ",
        family,
        " | ",
        "Genus: ",
        genus,
        " | ",
        "Species: ",
        species
      ),
      paste0(legend_name, "\nTaxonomy unavailable")
    )
  )

write.table(
  top10 %>%
    select(genome, legend_name, family, genus, species, total_coverage),
  "q3_q4_top10_taxonomy.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plot_data <- combined_paired %>%
  left_join(
    top10 %>% select(genome, taxonomy_label),
    by = "genome"
  ) %>%
  mutate(virus_group = coalesce(taxonomy_label, "Other")) %>%
  group_by(
    SRA,
    SampleID,
    Plant,
    Date,
    Enriched,
    Days_from_first_sample,
    virus_group
  ) %>%
  summarise(
    coverage = sum(coverage, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(SRA) %>%
  mutate(rel_abund = coverage / sum(coverage)) %>%
  ungroup() %>%
  mutate(
    virus_group = factor(
      virus_group,
      levels = c(top10$taxonomy_label, "Other")
    )
  )

sars_legend_label <- top10 %>%
  filter(genome == SARS_ACC) %>%
  pull(taxonomy_label)

other_viruses <- setdiff(top10$taxonomy_label, sars_legend_label)
n_other <- length(other_viruses)
base_colors <- RColorBrewer::brewer.pal(max(n_other, 3), "Set3")[1:n_other]

virus_colors <- setNames(
  c(base_colors, "red", "gray80"),
  c(other_viruses, sars_legend_label, "Other")
)

p_composition <- ggplot(
  plot_data,
  aes(
    x = rel_abund,
    y = reorder(
      format(Date, "%m.%d.%Y"),
      Days_from_first_sample
    ),
    fill = virus_group
  )
) +
  geom_bar(stat = "identity") +
  facet_grid(
    Plant ~ Enriched,
    scales = "free_y"
  ) +
  scale_x_continuous(
    labels = scales::percent_format(),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = virus_colors) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Virus taxonomy",
    title = paste("Q3/Q4: Top 10 taxa (SARS-CoV-2 highlighted)", "across plants, dates, enrichment status")
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    ),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 6.5, lineheight = 0.9),
    legend.key.height = grid::unit(0.75, "cm"),
    legend.key.width = grid::unit(
      0.5,
      "cm"
    ),
    legend.spacing.y = grid::unit(0.02, "cm")
  )

if (interactive()) print(p_composition)

ggsave(
  "q3_q4_top10_composition_enriched_vs_unenriched.png",
  p_composition,
  width = 15,
  height = 10.5,
  dpi = 300
)

# REVIEW - STEP 3: TAXONOMIC CLASSIFICATION
# Q: Why can respiratory viruses remain outside the top taxa after enrichment?
# A. Enrichment failed.
# B. Dominant plant viruses/phages still outrank respiratory viruses that began very rare.
# C. Respiratory viruses never survive wastewater.
# D. Classifiers reject animal viruses.
# Answer: B.
# Insight: Enrichment can improve detection without making a target taxonomically dominant.
# ============================================================
# 5b. Q3/Q4.2: selected-virus detection in paired libraries
# ============================================================
# Human norovirus is restricted to human GI/GII/GIV RefSeq records in the
# local accession lookup. Animal-associated noroviruses (including canine,
# marmot, primate, GIII, GV, and GIV.2 carnivore records) are excluded.
human_norovirus_accessions <- c(
  "NC_001959.2",
  "NC_029646.1",
  "NC_039475.1",
  "NC_039476.1",
  "NC_039477.1",
  "NC_039897.1",
  "NC_040876.1",
  "NC_044045.1",
  "NC_044046.1",
  "NC_044853.1",
  "NC_044854.1",
  "NC_044855.1",
  "NC_044856.1",
  "NC_044932.1"
)
# The four routinely recognized seasonal human coronaviruses. SARS-CoV-2 is
# deliberately excluded here and displayed as its own virus group.
seasonal_hcov_accessions <- c(
  "NC_002645.1",
  "NC_005831.2",
  "NC_006213.1",
  "NC_006577.2"
)

interest_virus_taxonomy <- tibble::tribble(
  ~virus_group, ~display_name, ~family, ~genus, ~species,
  "Pepper mild mottle virus", "Pepper mild mottle virus (reference control)", "Virgaviridae",
    "Tobamovirus", "Tobamovirus capsici",
  "Human norovirus", "Human norovirus", "Caliciviridae", "Norovirus", "Norovirus norwalkense",
  "Seasonal human coronavirus", "Seasonal human coronavirus", "Coronaviridae",
    "Alphacoronavirus / Betacoronavirus", "multiple species",
  "Human adenovirus", "Human adenovirus", "Adenoviridae", "Mastadenovirus", "multiple species",
  "Influenza A virus", "Influenza A virus", "Orthomyxoviridae", "Alphainfluenzavirus", "Alphainfluenzavirus influenzae",
  "SARS-CoV-2", "SARS-CoV-2", "Coronaviridae", "Betacoronavirus", "Betacoronavirus pandemicum"
) %>%
  mutate(
    taxonomy_label = paste0(
      display_name,
      "\n",
      "Family: ",
      family,
      " | Genus: ",
      genus,
      " | Species: ",
      species
    )
  )
# Assign each genome_info row to one requested group. A group is detected if at
# least one of its matching accession/segment rows reaches the >=50-pair rule;
# weak rows are not added together to manufacture a positive call.
interest_hits <- combined_paired %>%
  mutate(
    virus_group = case_when(
      genome == "NC_003630.1" ~ "Pepper mild mottle virus",
      genome == SARS_ACC ~ "SARS-CoV-2",
      genome %in% human_norovirus_accessions ~ "Human norovirus",
      genome %in% seasonal_hcov_accessions ~ "Seasonal human coronavirus",
      str_detect(
        coalesce(virus_name, ""),
        regex("^Human adenovirus\\b", ignore_case = TRUE)
      ) ~ "Human adenovirus",
      str_detect(
        coalesce(virus_name, ""),
        regex("^Influenza A virus\\b", ignore_case = TRUE)
      ) ~ "Influenza A virus",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(virus_group))

interest_calls_observed <- interest_hits %>%
  group_by(SRA, virus_group) %>%
  summarise(
    detected = any(
      filtered_read_pair_count >= VIRUS_DETECTION_MIN_READ_PAIRS,
      na.rm = TRUE
    ),
    total_filtered_read_pair_count = sum(
      filtered_read_pair_count,
      na.rm = TRUE
    ),
    max_reference_filtered_read_pair_count = max(
      c(0, filtered_read_pair_count),
      na.rm = TRUE
    ),
    matching_accession_count = n_distinct(genome),
    matching_accessions = paste(
      sort(unique(genome)),
      collapse = "; "
    ),
    qualifying_accession_count = n_distinct(
      genome[filtered_read_pair_count >= VIRUS_DETECTION_MIN_READ_PAIRS]
    ),
    qualifying_accessions = paste(
      sort(
        unique(
          genome[
            filtered_read_pair_count >= VIRUS_DETECTION_MIN_READ_PAIRS
          ]
        )
      ),
      collapse = "; "
    ),
    .groups = "drop"
  )
# Start from the full set of paired libraries so a missing genome_info record is
# represented as zero read pairs (Not detected), rather than disappearing.
paired_sample_metadata <- link_table %>%
  semi_join(
    paired_dates,
    by = c("Plant", "Date")
  ) %>%
  select(
    SRA,
    SampleID,
    Plant,
    Date,
    Enriched
  ) %>%
  distinct()

interest_heatmap_data <- tidyr::crossing(
  SRA = paired_sample_metadata$SRA,
  virus_group = interest_virus_taxonomy$virus_group
) %>%
  left_join(
    paired_sample_metadata,
    by = "SRA"
  ) %>%
  left_join(
    interest_calls_observed,
    by = c("SRA", "virus_group")
  ) %>%
  mutate(
    detected = coalesce(detected, FALSE),
    total_filtered_read_pair_count = coalesce(
      total_filtered_read_pair_count,
      0
    ),
    max_reference_filtered_read_pair_count = coalesce(
      max_reference_filtered_read_pair_count,
      0
    ),
    matching_accession_count = coalesce(
      matching_accession_count,
      0L
    ),
    matching_accessions = coalesce(
      matching_accessions,
      ""
    ),
    qualifying_accession_count = coalesce(
      qualifying_accession_count,
      0L
    ),
    qualifying_accessions = coalesce(
      qualifying_accessions,
      ""
    ),
    detection_status = factor(
      if_else(detected, "Detected", "Not detected"),
      levels = c("Not detected", "Detected")
    )
  ) %>%
  left_join(
    interest_virus_taxonomy,
    by = "virus_group"
  ) %>%
  mutate(
    Plant = factor(
      Plant,
      levels = c("HTP", "PL")
    ),
    Enriched = factor(
      Enriched,
      levels = c("enriched", "unenriched")
    ),
    date_label = factor(
      format(Date, "%m.%d.%Y"),
      levels = format(sort(unique(Date)), "%m.%d.%Y")
    ),
    taxonomy_label = factor(
      taxonomy_label,
      levels = interest_virus_taxonomy$taxonomy_label
    ),
    virus_axis_label = factor(
      display_name,
      levels = interest_virus_taxonomy$display_name
    ),
    detection_fill = factor(
      if_else(detected, as.character(taxonomy_label), "Not detected"),
      levels = c(interest_virus_taxonomy$taxonomy_label, "Not detected")
    ),
    detection_threshold = VIRUS_DETECTION_MIN_READ_PAIRS
  )

pmmov_check <- interest_heatmap_data %>%
  filter(virus_group == "Pepper mild mottle virus")

n_pmmov_total <- nrow(pmmov_check)
n_pmmov_detected <- sum(pmmov_check$detected, na.rm = TRUE)

if (n_pmmov_detected != n_pmmov_total) {
  warning(
    "PMMoV did not meet the Q3/Q4.2 threshold in every paired library: ",
    n_pmmov_detected,
    "/",
    n_pmmov_total
  )
}

interest_detection_summary <- interest_heatmap_data %>%
  mutate(
    taxonomy_label = as.character(taxonomy_label)
  ) %>%
  group_by(
    virus_group,
    display_name,
    family,
    genus,
    species,
    detection_threshold
  ) %>%
  summarise(
    detected_libraries = sum(detected, na.rm = TRUE),
    total_paired_libraries = n(),
    detection_percent = 100 * detected_libraries / total_paired_libraries,
    .groups = "drop"
  ) %>%
  mutate(
    virus_order = match(
      virus_group,
      interest_virus_taxonomy$virus_group
    )
  ) %>%
  arrange(virus_order) %>%
  select(-virus_order)

write.table(
  interest_heatmap_data %>%
    mutate(
      taxonomy_label = as.character(taxonomy_label),
      date_label = as.character(date_label),
      detection_status = as.character(detection_status)
    ) %>%
    arrange(
      Plant,
      Date,
      Enriched,
      match(virus_group, interest_virus_taxonomy$virus_group)
    ) %>%
    select(
      SRA,
      SampleID,
      Plant,
      Date,
      Enriched,
      virus_group,
      display_name,
      family,
      genus,
      species,
      total_filtered_read_pair_count,
      max_reference_filtered_read_pair_count,
      matching_accession_count,
      matching_accessions,
      qualifying_accession_count,
      qualifying_accessions,
      detection_threshold,
      detected,
      detection_status
    ),
  "q3_q4.2_viruses_of_interest_detection_calls.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  interest_detection_summary,
  "q3_q4.2_viruses_of_interest_detection_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
cat("\n=== Q3/Q4.2 VERSION 3: SELECTED-VIRUS DETECTION ===\n")
cat(
  "Detection rule: at least one genome_info row with filtered_read_pair_count >=",
  VIRUS_DETECTION_MIN_READ_PAIRS,
  "per virus group and library\n"
)
print(interest_detection_summary)

cat(
  "PMMoV reference control detected:",
  n_pmmov_detected,
  "of",
  n_pmmov_total,
  "paired libraries\n"
)

interest_virus_colors <- setNames(
  c(
    "#009E73",
    "#CC79A7",
    "#E69F00",
    "#56B4E9",
    "#0072B2",
    "#D62728"
  ),
  interest_virus_taxonomy$taxonomy_label
)

interest_heatmap_colors <- c(
  interest_virus_colors,
  `Not detected` = "#E6E6E6"
)

# Invisible points train the legend so virus colors remain visible even when a
# virus has zero qualifying detections at the selected threshold.
interest_legend_training <- interest_virus_taxonomy %>%
  transmute(
    virus_axis_label = factor(
      display_name,
      levels = interest_virus_taxonomy$display_name
    ),
    date_label = factor(
      format(
        min(paired_sample_metadata$Date),
        "%m.%d.%Y"
      ),
      levels = format(
        sort(unique(paired_sample_metadata$Date)),
        "%m.%d.%Y"
      )
    ),
    detection_fill = factor(
      taxonomy_label,
      levels = c(
        interest_virus_taxonomy$taxonomy_label,
        "Not detected"
      )
    ),
    Plant = factor(
      "HTP",
      levels = c("HTP", "PL")
    ),
    Enriched = factor(
      "enriched",
      levels = c("enriched", "unenriched")
    )
  )

p_interest_heatmap <- ggplot(
  interest_heatmap_data,
  aes(
    x = virus_axis_label,
    y = date_label,
    fill = detection_fill
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.35
  ) +
  geom_point(
    data = interest_legend_training,
    aes(
      x = virus_axis_label,
      y = date_label,
      fill = detection_fill
    ),
    inherit.aes = FALSE,
    shape = 22,
    size = 4,
    alpha = 0,
    stroke = 0,
    show.legend = TRUE
  ) +
  facet_grid(
    Plant ~ Enriched,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(
      Enriched = c(
        enriched = "Enriched",
        unenriched = "Unenriched"
      ),
      Plant = c(
        HTP = "HTP",
        PL = "PL"
      )
    )
  ) +
  scale_fill_manual(
    values = interest_heatmap_colors,
    breaks = names(interest_heatmap_colors),
    drop = FALSE,
    name = "Detected virus (Family | Genus | Species)"
  ) +
  scale_x_discrete(
    drop = FALSE,
    labels = function(x) str_wrap(x, width = 18),
    expand = expansion(add = 0.15)
  ) +
  scale_y_discrete(expand = expansion(add = 0.15)) +
  labs(
    x = "Virus",
    y = "Sampling date",
    title = "Q3/Q4.2: Detection of selected viruses in paired-end reads",
    subtitle = paste0(
      "Detected when >=1 matching genome_info row has filtered_read_pair_count >= ",
      VIRUS_DETECTION_MIN_READ_PAIRS,
      " per library; PMMoV reference control = ",
      n_pmmov_detected,
      "/",
      n_pmmov_total
    ),
    caption = paste(
      "Only plant/date combinations with both enriched and unenriched libraries are shown.",
      "Colored tiles are detected viruses; gray tiles are not detected.",
      "Seasonal human coronavirus excludes SARS-CoV-2."
    )
  ) +
  guides(
    fill = guide_legend(
      ncol = 2,
      byrow = TRUE,
      override.aes = list(
        fill = unname(interest_heatmap_colors),
        color = NA,
        alpha = 1,
        shape = 22,
        size = 4,
        stroke = 0
      )
    )
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      size = 8
    ),
    axis.text.y = element_text(size = 7.5),
    axis.title = element_text(size = 10),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(
      size = 6.5,
      lineheight = 0.9
    ),
    legend.key.height = grid::unit(0.55, "cm"),
    legend.key.width = grid::unit(0.55, "cm"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 9),
    plot.caption = element_text(size = 8, hjust = 0),
    panel.spacing = grid::unit(0.8, "lines")
  )

if (interactive()) print(p_interest_heatmap)

ggsave(
  "q3_q4.2_viruses_of_interest_detection_heatmap.png",
  p_interest_heatmap,
  width = 13,
  height = 15,
  dpi = 300
)

# REVIEW - STEP 4: TARGETED ABUNDANCE OVER TIME
# Q: A respiratory virus is present but absent from the top-taxa plot. What is the strongest next analysis?
# A. Enlarge the font.
# B. Conclude it is absent.
# C. Examine target-mapped reads, breadth/coverage, and detection over time.
# D. Change the random seed.
# Answer: C.
# Insight: Rank plots show dominance, not whether a specific low-abundance target is present.
# ============================================================
# 6. Q5a: SARS-CoV-2 detection, with high-breadth highlighting
# ============================================================
sars_genome_info <- combined %>%
  filter(genome == SARS_ACC) %>%
  select(
    SRA,
    SampleID,
    Plant,
    Date,
    Enriched,
    coverage,
    breadth,
    filtered_read_pair_count,
    SNV_count,
    SNS_count
  ) %>%
  mutate(
    detected = coverage > 0,
    high_breadth = !is.na(breadth) & breadth >= BREADTH_CUTOFF,
    SampleID = factor(SampleID, levels = sample_order)
  )

n_high_breadth <- sum(sars_genome_info$high_breadth, na.rm = TRUE)

cat("\n=== Q5 VERSION 2: SARS-CoV-2 detection ===\n")
cat("SARS genome-summary rows:", nrow(sars_genome_info), "\n")
cat(
  "Rows with coverage > 0:",
  sum(sars_genome_info$detected, na.rm = TRUE),
  "\n"
)
cat(
  "Rows with >=",
  BREADTH_CUTOFF * 100,
  "% breadth:",
  n_high_breadth,
  "\n"
)

write.table(
  sars_genome_info,
  "sars_cov2_sample_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sars_positive_samples <- sars_genome_info %>%
  filter(detected) %>%
  select(SampleID, SRA, Plant, Date, Enriched, coverage, breadth)

write.table(
  sars_positive_samples,
  "sars_cov2_positive_sampleIDs.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Retain continuous breadth fill and label every bar with its breadth percentage.
# Bold percentage labels identify samples eligible for the stringent SNV analysis.
p_sars_detect <- ggplot(
  sars_genome_info,
  aes(x = SampleID, y = coverage, fill = breadth)
) +
  geom_col() +
  geom_text(
    data = filter(sars_genome_info, !high_breadth),
    aes(label = scales::percent(breadth, accuracy = 0.1)),
    hjust = -0.15,
    size = 2.2,
    color = "gray35"
  ) +
  geom_text(
    data = filter(sars_genome_info, high_breadth),
    aes(label = scales::percent(breadth, accuracy = 0.1)),
    hjust = -0.15,
    size = 2.5,
    color = "black",
    fontface = "bold"
  ) +
  coord_flip(clip = "off") +
  scale_fill_viridis_c(
    name = "Genome breadth",
    limits = c(0, 1),
    labels = scales::percent_format()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    x = NULL,
    y = "SARS-CoV-2 coverage",
    title = "Q5: SARS-CoV-2 detection by sample",
    subtitle = paste0(
      "Breadth labels: bold values meet the >=",
      BREADTH_CUTOFF * 100,
      "% SNV-analysis threshold (n = ",
      n_high_breadth,
      ")"
    )
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 6),
    plot.margin = margin(5.5, 35, 5.5, 5.5)
  )

if (interactive()) print(p_sars_detect)

ggsave(
  "q5.1_sars_cov2_detection_by_sample.png",
  p_sars_detect,
  width = 8,
  height = 12,
  dpi = 300
)
# ============================================================
# 3. Strict subset and ORF annotation
# ============================================================
samples_stringent <- sars_genome_info %>%
  filter(high_breadth) %>%
  pull(SRA)

if (length(samples_stringent) == 0) {
  stop(
    "No samples meet the >=",
    BREADTH_CUTOFF * 100,
    "% breadth threshold; no Q5 SNV figures can be created."
  )
}

snv_combined <- lapply(snv_files, function(file_path) {
  snv_data <- read.table(
    file_path,
    header = TRUE,
    sep = "\t"
  )
  snv_data$SRA <- str_extract(
    basename(file_path),
    "SRR[0-9]+"
  )
  snv_data %>%
    filter(
      scaffold == SARS_ACC,
      class != "AmbiguousReference"
    )
}) %>%
  bind_rows() %>%
  left_join(link_table, by = "SRA") %>%
  filter(SRA %in% samples_stringent)

if (nrow(snv_combined) == 0) {
  stop("No SARS-CoV-2 SNVs remain after stringent sample filtering.")
}

write.table(
  snv_combined,
  "sars_cov2_snvs_stringent.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sars_orf_map <- tibble::tribble(
  ~ORF, ~start, ~end, ~function_note,
  "5'UTR", 1, 265, "Untranslated region",
  "ORF1ab", 266, 21555,
  "Replicase polyprotein: RNA synthesis, proofreading, immune evasion",
  "S", 21563, 25384, "Spike protein: receptor binding, cell entry",
  "ORF3a", 25393, 26220, "Viroporin; host inflammatory response",
  "E", 26245, 26472, "Envelope protein: viral assembly, budding",
  "M", 26523, 27191, "Membrane protein: viral assembly",
  "ORF6", 27202, 27387, "Interferon antagonist",
  "ORF7a", 27394, 27759, "Immune modulation",
  "ORF7b", 27756, 27887, "Immune modulation",
  "ORF8", 27894, 28259, "Immune evasion, MHC-I downregulation",
  "N", 28274, 29533, "Nucleocapsid protein: genome packaging, replication",
  "ORF10", 29558, 29674, "Uncertain function",
  "3'UTR", 29675, 29903, "Untranslated region"
)

assign_orf_vec <- function(positions) {
  index <- sapply(positions, function(position) {
    hits <- which(
      position >= sars_orf_map$start &
        position <= sars_orf_map$end
    )

    if (length(hits) == 0)
      NA_integer_
    else
      hits[1]
  })

  list(
    ORF = sars_orf_map$ORF[index],
    ORF_function = sars_orf_map$function_note[index]
  )
}

orf_result <- assign_orf_vec(snv_combined$position)

snv_annotated <- snv_combined %>%
  mutate(
    ORF = orf_result$ORF,
    ORF_function = orf_result$ORF_function,
    change = paste0(ref_base, position, var_base)
  ) %>%
  select(
    SRA,
    SampleID,
    Plant,
    Date,
    Enriched,
    position,
    ref_base,
    var_base,
    change,
    var_freq,
    position_coverage,
    class,
    ORF,
    ORF_function
  )

write.table(
  snv_annotated,
  "sars_cov2_snvs_orf_annotated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
# ============================================================
# 4. Explicit high-frequency SNV objects
# ============================================================
# These are call-level events from the high-breadth sample set above.
high_frequency_snvs <- snv_annotated %>%
  filter(
    !is.na(var_freq),
    var_freq > HIGH_FREQ_CUTOFF
  ) %>%
  distinct(
    SRA,
    SampleID,
    Plant,
    Date,
    Enriched,
    position,
    ref_base,
    var_base,
    change,
    var_freq,
    position_coverage,
    class,
    ORF,
    ORF_function,
    .keep_all = TRUE
  )

write.table(
  high_frequency_snvs,
  "q5_high_frequency_snvs.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
# ============================================================
# 5. Q5b: Panel A data for high-breadth sample versus SNV count
# ============================================================
# The left join is deliberate: it keeps high-breadth samples with zero
# high-frequency SNVs, making the 7-samples-versus-8-SNV-events relationship clear.
high_breadth_snv_counts <- sars_genome_info %>%
  filter(high_breadth) %>%
  transmute(
    SRA,
    SampleID = as.character(SampleID),
    Plant,
    Date,
    Enriched,
    coverage,
    breadth
  ) %>%
  left_join(
    high_frequency_snvs %>%
      distinct(SRA, position, change) %>%
      count(SRA, name = "n_high_frequency_snvs"),
    by = "SRA"
  ) %>%
  mutate(
    n_high_frequency_snvs = coalesce(
      n_high_frequency_snvs,
      0L
    ),
    SampleID = factor(SampleID, levels = sample_order)
  ) %>%
  arrange(Date, Plant, SampleID)

cat("\n=== High-breadth sample versus high-frequency-SNV count ===\n")
print(
  high_breadth_snv_counts %>%
    select(
      SampleID,
      Date,
      coverage,
      breadth,
      n_high_frequency_snvs
    )
)

# Panel A is embedded in the genome-map PNG rather than exported separately.
panel_a_table_data <- high_breadth_snv_counts %>%
  transmute(
    `Sample ID` = as.character(SampleID),
    `Sampling date` = format(Date, "%Y-%m-%d"),
    `Genome breadth` = scales::percent(
      breadth,
      accuracy = 0.1
    ),
    `High-frequency SNVs` = n_high_frequency_snvs
  ) %>%
  bind_rows(
    tibble::tibble(
      `Sample ID` = paste0(
        "Total (",
        nrow(high_breadth_snv_counts),
        " samples)"
      ),
      `Sampling date` = "-",
      `Genome breadth` = "-",
      `High-frequency SNVs` = sum(
        high_breadth_snv_counts$n_high_frequency_snvs
      )
    )
  )

panel_a_column_labels <- c(
  "Sample ID",
  "Sampling date",
  "Genome breadth",
  paste0(
    "High-frequency SNVs\n(>",
    HIGH_FREQ_CUTOFF * 100,
    "%)"
  )
)

panel_a_plot_data <- panel_a_table_data %>%
  mutate(
    table_row = rev(seq_len(n())),
    is_total = row_number() == n()
  ) %>%
  mutate(
    across(
      -c(table_row, is_total),
      as.character
    )
  ) %>%
  pivot_longer(
    cols = -c(table_row, is_total),
    names_to = "table_column",
    values_to = "display_value"
  ) %>%
  mutate(
    table_column = factor(
      table_column,
      levels = names(panel_a_table_data),
      labels = panel_a_column_labels
    )
  )

p_panel_a <- ggplot(
  panel_a_plot_data,
  aes(x = table_column, y = table_row)
) +
  geom_tile(
    aes(fill = is_total),
    color = "#D9D9D9",
    linewidth = 0.35
  ) +
  geom_text(
    aes(label = display_value),
    size = 3.5
  ) +
  scale_fill_manual(
    values = c(
      `FALSE` = "white",
      `TRUE` = "#E6EEF5"
    ),
    guide = "none"
  ) +
  scale_x_discrete(
    position = "top",
    expand = expansion(add = 0)
  ) +
  scale_y_continuous(
    breaks = NULL,
    expand = expansion(add = 0)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = paste0(
      "A. High-breadth samples (>=",
      BREADTH_CUTOFF * 100,
      "%) and high-frequency SNV counts"
    ),
    subtitle = paste0(
      nrow(high_breadth_snv_counts),
      " high-breadth samples; ",
      sum(high_breadth_snv_counts$n_high_frequency_snvs),
      " high-frequency SNV calls in total. Counts are distinct position/change calls per sample; ",
      "zero-count samples are retained."
    )
  ) +
  theme_void() +
  theme(
    axis.text.x = element_text(
      color = "black",
      face = "bold",
      size = 10,
      lineheight = 0.9,
      margin = margin(b = 5)
    ),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    plot.margin = margin(8, 10, 8, 10)
  )
# ============================================================
# 6. Q5c: high-frequency SNV genome map
# ============================================================
# The points and labels represent high-frequency SNV position/change events.
# Contributing SampleIDs are retained in the TSV, not used as point labels.
high_frequency_snv_positions <- high_frequency_snvs %>%
  group_by(position, ORF, ORF_function, change) %>%
  summarise(
    n_samples = n_distinct(SRA),
    samples = paste(
      sort(unique(SampleID)),
      collapse = "; "
    ),
    .groups = "drop"
  ) %>%
  mutate(
    snv_label = paste0(
      change,
      "\n[",
      coalesce(ORF, "Unannotated"),
      "]"
    )
  )

cat("\n=== High-frequency SNV positions ===\n")
print(
  high_frequency_snv_positions %>%
    select(position, change, ORF, n_samples, samples)
)

write.table(
  high_frequency_snv_positions,
  "q5_high_frequency_snv_position_sampleIDs.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

gene_track <- sars_orf_map %>%
  filter(!ORF %in% c("5'UTR", "3'UTR"))

p_high_frequency_genome_map <- ggplot() +
  geom_rect(
    data = gene_track,
    aes(
      xmin = start,
      xmax = end,
      ymin = -0.3,
      ymax = 0.3,
      fill = ORF
    ),
    alpha = 0.5,
    color = "black"
  ) +
  geom_point(
    data = high_frequency_snv_positions,
    aes(x = position, y = 0),
    shape = 21,
    size = 3.5,
    fill = "#D55E00",
    color = "black",
    stroke = 0.45
  ) +
  geom_text_repel(
    data = gene_track,
    aes(
      x = (start + end) / 2,
      y = 0.3,
      label = ORF
    ),
    ylim = c(0.6, 1.8),
    size = 3,
    fontface = "bold",
    segment.color = "gray40",
    segment.size = 0.4,
    min.segment.length = 0,
    direction = "x",
    force = 3,
    box.padding = 0.3
  ) +
  geom_text_repel(
    data = high_frequency_snv_positions,
    aes(
      x = position,
      y = 0,
      label = snv_label
    ),
    ylim = c(-2.8, -0.4),
    size = 2.5,
    color = "darkred",
    segment.color = "darkred",
    segment.size = 0.3,
    min.segment.length = 0,
    direction = "both",
    force = 2,
    box.padding = 0.3,
    point.padding = 0.15,
    max.overlaps = Inf,
    seed = 123,
    lineheight = 0.8
  ) +
  labs(
    x = "Genomic position (NC_045512.2)",
    y = NULL,
    title = paste0(
      "B. Genome map of high-frequency SARS-CoV-2 SNVs (variant frequency >",
      HIGH_FREQ_CUTOFF * 100,
      "%)"
    ),
    subtitle = paste0(
      "Each point is an SNV position/change from samples with >=",
      BREADTH_CUTOFF * 100,
      "% genome breadth"
    )
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  guides(fill = "none") +
  coord_cartesian(
    ylim = c(-3.2, 2),
    clip = "off"
  )

p_q5_genome_map_composite <- patchwork::wrap_plots(
  p_panel_a,
  p_high_frequency_genome_map,
  ncol = 1,
  heights = c(1.35, 2.65)
)

if (interactive()) print(p_q5_genome_map_composite)

ggsave(
  "q5.2_sars_cov_2_high_frequency_snv_genome_map.png",
  p_q5_genome_map_composite,
  width = 13,
  height = 11.5,
  dpi = 300
)
# ============================================================
# 7. Q5d: top variable SNVs over time, with ORF in legend
# ============================================================
# This panel intentionally uses all observed SNVs in high-breadth samples.
# It is a variability plot, not a high-frequency-only plot.
snv_over_time <- snv_annotated %>%
  mutate(ORF = coalesce(ORF, "Unannotated")) %>%
  group_by(Plant, Date, position, change, ORF) %>%
  summarise(
    overall_snv = mean(var_freq, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    position_label = paste0(position, " (", change, ")"),
    legend_label = paste0(
      position_label,
      " [",
      ORF,
      "]"
    )
  )

TOP_N_VARIABLE <- 10

top_variable_positions <- snv_over_time %>%
  group_by(position, change, ORF, legend_label) %>%
  summarise(
    sd_freq = sd(overall_snv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(sd_freq)) %>%
  slice(1:TOP_N_VARIABLE) %>%
  pull(legend_label)

snv_over_time_top <- snv_over_time %>%
  filter(legend_label %in% top_variable_positions) %>%
  mutate(
    legend_label = factor(
      legend_label,
      levels = top_variable_positions
    )
  )

p_snv_top_variable <- ggplot(
  snv_over_time_top,
  aes(
    x = Date,
    y = overall_snv,
    color = legend_label,
    group = legend_label
  )
) +
  geom_line(alpha = 0.8) +
  geom_point(size = 2) +
  facet_wrap(
    ~Plant,
    scales = "free_x"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  labs(
    y = "Variant frequency",
    x = "Sampling date",
    color = "Position (change) [ORF]",
    title = paste0("Q5: Top ", TOP_N_VARIABLE, " most variable SNV positions over time"),
    subtitle = paste0(
      "All observed SNVs in samples with >=",
      BREADTH_CUTOFF * 100,
      "% genome breadth; this panel is not high-frequency-only"
    )
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.text = element_text(size = 7)
  )

if (interactive()) print(p_snv_top_variable)

ggsave(
  "q5.4_sars_cov2_snp_top_variable_over_time.png",
  p_snv_top_variable,
  width = 10,
  height = 6,
  dpi = 300
)
# ============================================================
# 8. Q5e: high-frequency SNV burden by ORF over time
# ============================================================
snv_by_orf_time <- high_frequency_snvs %>%
  distinct(Plant, Date, position, ORF) %>%
  count(Plant, Date, ORF, name = "n_snvs")

orf_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_distinct(sars_orf_map$ORF)),
  sars_orf_map$ORF
)

p_snv_by_orf <- ggplot(
  snv_by_orf_time,
  aes(x = Date, y = n_snvs, fill = ORF)
) +
  geom_col() +
  facet_wrap(~Plant, scales = "free_x") +
  scale_fill_manual(values = orf_colors) +
  labs(
    x = "Sampling date",
    y = "Number of high-frequency SNVs",
    fill = "ORF",
    title = paste0(
      "Q5: High-frequency SNV burden by ORF over time (>",
      HIGH_FREQ_CUTOFF * 100,
      "% variant frequency)"
    ),
    subtitle = paste0(
      "Calls come from samples with >=",
      BREADTH_CUTOFF * 100,
      "% genome breadth"
    )
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

if (interactive()) print(p_snv_by_orf)

ggsave(
  "q5.3_sars_cov_high_frequency_snv_count_by_orf_over_time.png",
  p_snv_by_orf,
  width = 10,
  height = 6,
  dpi = 300
)

# REVIEW - STEP 5: SARS-CoV-2 SNVs OVER TIME
# Q: Why are some SNVs consistently high-frequency while others vary week to week?
# A. Stable SNVs are errors.
# B. Stable SNVs may define lineages; variable SNVs can reflect depth,
#    lineage shifts, and sampling noise.
# C. Any change proves a new variant.
# D. Variable SNVs mean the reference is wrong.
# Answer: B.
# Insight: Temporal stability helps separate population structure from
# low-depth noise in a mixed variant pool.
