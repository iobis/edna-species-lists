library(jsonlite)
library(glue)
library(purrr)
library(dplyr)
library(stringr)
library(robis)
library(worrms)
library(furrr)
library(rredlist)
library(tidyr)
library(jsonlite)

# TODO: remove contaminants and bacteria from reads/ASV statistics

include_dna <- TRUE
threatened_categories <- c("CR", "EN", "EW", "EX", "VU")
markers <- c("16s", "coi", "mifish", "mimammal", "teleo")
fish_classes <- c("Actinopteri", "Cladistii", "Coelacanthi", "Elasmobranchii", "Holocephali", "Myxini", "Petromyzonti", "Teleostei")
turtle_orders <- c("Testudines")
mammal_classes <- c("Mammalia")
bird_classes <- c("Aves")
mollusc_phyla <- c("Mollusca")
amphibia_classes <- c("Amphibia")
algae_phyla <- c("Chlorophyta", "Haptophyta", "Rhodophyta", "Ochrophyta", "Bacillariophyta")
echinoderm_phyla <- c("Echinodermata")
sponge_phyla <- c("Porifera")
cnidaria_phyla <- c("Cnidaria", "Ctenophora")
unicellular_phyla <- c("Cercozoa", "Amoebozoa", "Myzozoa")
fungi_phyla <- c("Ascomycota", "Oomycota")
worms_phyla <- c("Nemertea", "Gnathostomulida", "Annelida")
bryozoa_phyla <- c("Bryozoa")
phoronida_phyla <- c("Phoronida")
copepod_classes <- c("Copepoda")
crustacean_classes <- c("Malacostraca", "Thecostraca", "Branchiopoda")
arrowworm_phyla <- c("Chaetognatha")
ascidian_classes <- c("Ascidiacea")

# read OBIS species lists from https://github.com/iobis/mwhs-obis-species

sites <- fromJSON("https://raw.githubusercontent.com/iobis/mwhs-obis-species/master/lists/sites.json")

obis_species <- map(sites, function(site_name) {
  site_list <- fromJSON(glue("https://raw.githubusercontent.com/iobis/mwhs-obis-species/master/lists/{site_name}.json"))$species
  if (length(site_list) > 0) {
    site_list$site <- site_name
  }
  return(site_list)
}) %>%
  bind_rows()

# read eDNA species lists

dna_files <- list.files("output", "*DNADerivedData*", full.names = TRUE)
occurrence_files <- list.files("output", "*Occurrence*", full.names = TRUE)

dna <- map(dna_files, read.table, sep = "\t", quote = "", header = TRUE) %>%
  bind_rows() %>%
  mutate_if(is.character, na_if, "")

site_identifiers <- str_match(occurrence_files, "/([^\\/]*)\\_Occurrence.tsv")[,2]

occurrence <- map(occurrence_files, read.table, sep = "\t", quote = "", header = TRUE) %>%
  setNames(site_identifiers) %>%
  bind_rows(.id = "site") %>%
  mutate_if(is.character, na_if, "") %>%
  mutate(
    aphiaid = as.numeric(str_replace(scientificNameID, "urn:lsid:marinespecies.org:taxname:", ""))
  ) %>%
  left_join(dna, by = "occurrenceID")

# process annotations

occurrence <- occurrence %>% mutate(remove = FALSE, remove_reads = FALSE)

annotations_files <- list.files("edna-results/annotations", "*.json", full.names = TRUE)

for (annotations_file in annotations_files) {
  message(annotations_file)

  # contaminants
  
  if (annotations_file == "edna-results/annotations/contaminants.json") {
    # remove contaminants from species lists but also reads/asv statistics
    annotations <- fromJSON(annotations_file)
    occurrence <- occurrence %>%
      mutate(remove = ifelse(genus %in% na.omit(annotations$genus), TRUE, remove)) %>%
      mutate(remove_reads = ifelse(genus %in% na.omit(annotations$genus), TRUE, remove_reads)) %>%
      mutate(remove = ifelse(kingdom %in% na.omit(annotations$kingdom), TRUE, remove)) %>%
      mutate(remove_reads = ifelse(kingdom %in% na.omit(annotations$kingdom), TRUE, remove_reads)) %>%
      mutate(remove = ifelse(phylum %in% na.omit(annotations$phylum), TRUE, remove)) %>%
      mutate(remove_reads = ifelse(phylum %in% na.omit(annotations$phylum), TRUE, remove_reads))
    next
  }
  
  # site specific annotations
  
  site_name <- str_match(annotations_file, ".*/([^\\/]*)\\.json")[,2]
  annotations <- fromJSON(annotations_file) %>%
    mutate(remove = as.logical(remove))
  
  # remove
  
  aphiaid_remove <- annotations %>% filter(remove) %>% pull(AphiaID)
  message(glue("Removing {length(aphiaid_remove)} species"))
  occurrence <- occurrence %>%
    mutate(remove = ifelse(aphiaid %in% aphiaid_remove & site == site_name, TRUE, remove))
  
  # replace
  
  aphiaid_replace <- annotations %>%
    filter(!is.na(new_AphiaID)) %>%
    mutate(new_AphiaID = as.numeric(str_trim(new_AphiaID)))
  message(glue("Replacing {nrow(aphiaid_replace)} species"))
  
  aphiaid_replace_batches <- split(aphiaid_replace$new_AphiaID, as.integer((seq_along(aphiaid_replace$new_AphiaID) - 1) / 50))
  
  replacement_taxa <- map(aphiaid_replace_batches, wm_record) %>%
    bind_rows() %>%
    select(valid_aphiaid = AphiaID, kingdom, phylum, class, order, family, genus, scientificName = scientificname, taxonRank = rank, scientificNameID = lsid) %>%
    mutate(taxonRank = tolower(taxonRank)) %>%
    distinct()
  
  aphiaid_replace <- aphiaid_replace %>%
    select(old_AphiaID = AphiaID, aphiaid = new_AphiaID) %>%
    left_join(replacement_taxa, by = c("aphiaid" = "valid_aphiaid")) %>%
    mutate(site = site_name)
  
  occurrence <- occurrence %>%
    mutate(old_AphiaID = aphiaid) %>%
    rows_update(aphiaid_replace, by = c("site", "old_AphiaID"), unmatched = "ignore") %>%
    select(-old_AphiaID)

}

# resolve eDNA species to accepted names

aphiaids <- unique(occurrence$aphiaid)
aphiaid_batches <- split(aphiaids, as.integer((seq_along(aphiaids) - 1) / 50))
plan(multisession, workers = 3)
aphiaid_mapping <- future_map(aphiaid_batches, wm_record) %>%
  bind_rows() %>%
  select(aphiaid = AphiaID, valid_aphiaid = valid_AphiaID) %>%
  distinct() %>%
  filter(aphiaid != valid_aphiaid)

valid_aphiaids <- unique(aphiaid_mapping$valid_aphiaid)
valid_aphiaid_batches <- split(valid_aphiaids, as.integer((seq_along(valid_aphiaids) - 1) / 50))
valid_taxa <- map(valid_aphiaid_batches, wm_record) %>%
  bind_rows() %>%
  select(valid_aphiaid = AphiaID, scientificName = scientificname, scientificNameID = lsid, taxonRank = rank, kingdom, phylum, class, order, family, genus) %>%
  mutate(taxonRank = tolower(taxonRank))

occurrence <- occurrence %>%
  mutate(verbatimScientificName = scientificName) %>%
  left_join(aphiaid_mapping, by = "aphiaid") %>%
  rows_update(valid_taxa, by = "valid_aphiaid") %>%
  mutate(
    species = ifelse(taxonRank == "species", scientificName, NA),
    aphiaid = as.numeric(str_replace(scientificNameID, "urn:lsid:marinespecies.org:taxname:", ""))
  ) %>%
  select(-valid_aphiaid)

dna_species <- occurrence %>%
  filter(taxonRank == "species") %>%
  select(-species) %>%
  rename(species = scientificName) %>%
  mutate(AphiaID = as.numeric(str_extract(scientificNameID, "[0-9]+"))) %>%
  group_by(phylum, class, order, family, genus, species, AphiaID, site, remove) %>%
  summarize(
    target_gene = paste0(sort(unique(target_gene)), collapse = ";"),
    reads = sum(organismQuantity),
    asvs = n_distinct(DNA_sequence)
  ) %>%
  mutate(source_dna = TRUE)

# make list

if (include_dna) {
  combined_species <- bind_rows(obis_species, dna_species)
} else {
  combined_species <- obis_species %>% mutate(target_gene = NA, source_dna = NA, reads = NA, asvs = NA, remove = FALSE)
}

species <- combined_species %>%
  group_by(AphiaID, phylum, class, order, family, genus, species, site) %>%
  summarize(
    remove = any(remove),
    records = first(na.omit(records)),
    reads = first(na.omit(reads)),
    asvs = first(na.omit(asvs)),
    max_year = first(na.omit(max_year)),
    target_gene = first(na.omit(target_gene)),
    source_obis = first(na.omit(obis)),
    source_gbif = first(na.omit(gbif)),
    source_dna = first(na.omit(source_dna))
  ) %>%
  ungroup()

# add red list

redlist_cache <- cachem::cache_disk(".cache")
get_redlist <- memoise::memoise(function() {
  redlist <- data.frame()
  page <- 0
  while (TRUE) {
    res <- rl_sp(page, key = "a936c4f78881e79a326e73c4f97f34a6e7d8f9f9e84342bff73c3ceda14992b9")$result
    if (length(res) == 0) {
      break
    }
    redlist <- bind_rows(redlist, res)
    page <- page + 1
  }
  redlist %>%
    filter(is.na(population)) %>%
    select(species = scientific_name, category)
}, cache = redlist_cache)
redlist <- get_redlist() %>%
  mutate(redlist_category = category)

species <- species %>%
  left_join(redlist, by = "species")

# add vernacular

vernaculars <- read.table("supporting_data/vernacularname.txt", sep = "\t", quote = "", header = TRUE) %>%
  mutate(AphiaID = as.numeric(str_replace(taxonID, "urn:lsid:marinespecies.org:taxname:", ""))) %>%
  filter(language %in% c("FRA", "ENG")) %>%
  select(AphiaID, language, vernacularName) %>%
  group_by(AphiaID, language) %>%
  summarize(vernacularName = first(vernacularName)) %>%
  arrange(AphiaID, language) %>%
  group_by(AphiaID) %>%
  summarize(vernacular = paste0(vernacularName, collapse = ", "))

species <- species %>%
  left_join(vernaculars, by = "AphiaID")

# add groups

species <- species %>%
  mutate(
    group = case_when(
      class %in% fish_classes ~ "fish",
      order %in% turtle_orders ~ "turtles",
      class %in% mammal_classes ~ "mammals",
      class %in% bird_classes ~ "birds",
      class %in% amphibia_classes ~ "amphibians",
      phylum %in% mollusc_phyla ~ "molluscs",
      phylum %in% algae_phyla ~ "algae",
      phylum %in% sponge_phyla ~ "sponges",
      phylum %in% cnidaria_phyla ~ "cnidarians",
      phylum %in% unicellular_phyla ~ "unicellular",
      phylum %in% fungi_phyla ~ "fungi",
      phylum %in% worms_phyla ~ "worms",
      phylum %in% bryozoa_phyla ~ "moss animals",
      phylum %in% phoronida_phyla ~ "horseshoe worms",
      phylum %in% echinoderm_phyla ~ "echinoderms",
      class %in% copepod_classes ~ "copepods",
      class %in% crustacean_classes ~ "crustaceans",
      phylum %in% arrowworm_phyla ~ "arrow worms",
      class %in% ascidian_classes ~ "sea squirts"
    )
  )

# fix logical

species <- species %>%
  mutate(across(where(is.logical), ~replace_na(., FALSE)))

# Generate JSON and CSV species lists

for (site_name in sites) {
  
  message(site_name)

  # site species list
  
  site_list <- species %>%
    filter(site == site_name) %>%
    select(-site) %>%
    arrange(group, phylum, class, order, species)
  
  # samples
  
  sample_volumes <- occurrence %>%
    filter(site == site_name) %>%
    distinct(materialSampleID, sampleSize) %>%
    mutate(sampleSize = as.integer(sampleSize))

  sample_species <- occurrence %>%
    filter(site == site_name) %>%
    filter(!remove) %>%
    group_by(materialSampleID) %>%
    summarize(species = n_distinct(na.omit(species)))

  sample_stats <- occurrence %>%
    filter(site == site_name & !remove_reads) %>%
    group_by(locality, materialSampleID) %>%
    summarize(
      reads = sum(organismQuantity),
      asvs = n_distinct(DNA_sequence),
      decimalLongitude = first(decimalLongitude),
      decimalLatitude = first(decimalLatitude)
    ) %>%
    left_join(sample_volumes, by = "materialSampleID") %>%
    left_join(sample_species, by = "materialSampleID")
    
  # sequence stats

  dna_stats <- occurrence %>%
    filter(site == site_name & !remove_reads) %>%
    summarize(reads = sum(organismQuantity), asvs = n_distinct(DNA_sequence))
  
  marker_stats <- occurrence %>%
    filter(site == site_name & !remove_reads) %>%
    group_by(pcr_primer_name_forward) %>%
    summarize(reads = sum(organismQuantity), species = n_distinct(na.omit(species)), asvs = n_distinct(DNA_sequence))

  # red list stats

  cat_obis <- site_list %>%
    filter(source_gbif | source_obis) %>%
    filter(!is.na(category)) %>%
    filter(category %in% threatened_categories) %>%
    group_by(category) %>%
    summarize(obis_species = n_distinct(species))
  
  cat_edna <- site_list %>%
    filter(source_dna) %>%
    filter(!remove) %>%
    filter(category %in% threatened_categories) %>%
    group_by(category) %>%
    summarize(edna_species = n_distinct(species))
  
  redlist_stats <- cat_obis %>%
    left_join(cat_edna, by = "category") %>%
    mutate(fraction = edna_species / obis_species)

  # family marker stats (tree figure)
  
  marker_stats_family <- occurrence %>%
    filter(site == site_name) %>%
    group_by(family, pcr_primer_name_forward) %>%
    summarize(reads = sum(organismQuantity)) %>%
    mutate(pcr_primer_name_forward = factor(pcr_primer_name_forward))
  
  # taxonomic group stats (including database)
  
  group_stats <- site_list %>%
    filter(!remove) %>%
    filter(!is.na(group)) %>%
    group_by(group) %>%
    summarize(n_species = n()) %>%
    data.table::transpose(make.names = TRUE)
  if (nrow(group_stats) == 0) {
    group_stats <- NULL
  }

  # taxonomic group stats (eDNA only)
  
  group_stats_edna <- site_list %>%
    filter(!remove) %>%
    filter(source_dna) %>%
    filter(!is.na(group)) %>%
    group_by(group) %>%
    summarize(n_species = n()) %>%
    data.table::transpose(make.names = TRUE)
  if (nrow(group_stats_edna) == 0) {
    group_stats_edna <- NULL
  }
  
  # source stats
  
  source_stats <- list(
    obis = site_list %>% filter(!remove) %>% filter(source_obis) %>% summarize(n_species = n()) %>% pull(n_species) %>% unbox(),
    gbif = site_list %>% filter(!remove) %>% filter(source_gbif) %>% summarize(n_species = n()) %>% pull(n_species) %>% unbox(),
    db = site_list %>% filter(!remove) %>% filter(source_gbif | source_obis) %>% summarize(n_species = n()) %>% pull(n_species) %>% unbox(),
    edna = site_list %>% filter(!remove) %>% filter(source_dna) %>% summarize(n_species = n()) %>% pull(n_species) %>% unbox(),
    both = site_list %>% filter(!remove) %>% filter(source_dna & (source_obis | source_gbif)) %>% summarize(n_species = n()) %>% pull(n_species) %>% unbox()
  )
  
  # full lists  
  
  json = toJSON(list(
    created = unbox(strftime(as.POSIXlt(Sys.time(), "UTC"), "%Y-%m-%dT%H:%M:%S")),
    species = site_list %>% filter(!remove),
    stats = list(
      samples = sample_stats,
      dna = unbox(dna_stats),
      markers = marker_stats,
      redlist = redlist_stats,
      markers_family = marker_stats_family,
      groups = unbox(group_stats),
      groups_edna = unbox(group_stats_edna),
      source = source_stats
    )
  ), pretty = TRUE)
  write(json, glue("lists_full/json/{site_name}.json"))

  write.csv(site_list %>% filter(!remove), glue("lists_full/csv/{site_name}.csv"), row.names = FALSE, na = "")
  
  # eDNA only lists

  json = toJSON(list(
    created = unbox(strftime(as.POSIXlt(Sys.time(), "UTC"), "%Y-%m-%dT%H:%M:%S")),
    species = site_list %>% filter(!remove) %>% filter(source_dna),
    stats = list(
      samples = sample_stats,
      dna = unbox(dna_stats),
      markers = marker_stats,
      redlist = redlist_stats,
      groups_edna = unbox(group_stats_edna),
      source = source_stats
    )
  ), pretty = TRUE)
  write(json, glue("lists/json/{site_name}.json"))
  
  write.csv(site_list %>% filter(!remove) %>% filter(source_dna), glue("lists/csv/{site_name}.csv"), row.names = FALSE, na = "")

}

# DEBUGGING

# everglades <- species %>%
#   filter(site == "everglades_national_park" & source_dna) %>%
#   arrange(phylum, class, species) %>%
#   View()
