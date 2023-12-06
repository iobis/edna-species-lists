library(jsonlite)
library(glue)
library(purrr)
library(dplyr)
library(stringr)
library(robis)

include_dna <- FALSE
markers <- c("16s", "coi", "mifish", "mimammal", "teleo")
fish_classes <- c("Actinopteri", "Cladistii", "Coelacanthi", "Elasmobranchii", "Holocephali", "Myxini", "Petromyzonti", "Teleostei")
turtle_orders <- c("Testudines")
mammal_classes <- c("Mammalia")
bird_classes <- c("Aves")
mollusc_phyla <- c("Mollusca")
amphibia_classes <- c("Amphibia")
algae_phyla <- c("Chlorophyta", "Haptophyta", "Rhodophyta", "Ochrophyta", "Bacillariophyta")
starfish_phyla <- c("Echinodermata")
sponge_phyla <- c("Porifera")
jelly_phyla <- c("Cnidaria", "Ctenophora")
unicellular_phyla <- c("Cercozoa", "Amoebozoa", "Myzozoa")
fungi_phyla <- c("Ascomycota", "Oomycota")
worms_phyla <- c("Nemertea", "Gnathostomulida", "Annelida")
filter_feeders_phyla <- c("Phoronida", "Bryozoa")
copepod_classes <- c("Copepoda")
crustacean_classes <- c("Malacostraca")

# Read OBIS species lists from https://github.com/iobis/mwhs-obis-species

sites <- fromJSON("https://raw.githubusercontent.com/iobis/mwhs-obis-species/master/lists/sites.json")

obis_species <- map(sites, function(site_name) {
  site_list <- fromJSON(glue("https://raw.githubusercontent.com/iobis/mwhs-obis-species/master/lists/{site_name}.json"))$species
  if (length(site_list) > 0) {
    site_list$site <- site_name
  }
  return(site_list)
}) %>%
  bind_rows()

# Read eDNA species lists from this repository

dna_occurrence <- map(sites, function(site_name) {
  for (marker in markers) {
    occurrence_file <- glue("pipeline_results/{site_name}/{marker}/Occurrence.txt")
    dna_file <- glue("pipeline_results/{site_name}/{marker}/DNADerivedData.txt")
    if (file.exists(occurrence_file)) {
      occurrence <- read.table(occurrence_file, sep = "\t", header = TRUE, na.strings = "")
      dna <- read.table(dna_file, sep = "\t", header = TRUE, na.strings = "")
      occurrence <- occurrence %>%
        left_join(dna, by = "occurrenceID")
      occurrence$site <- site_name
      return(occurrence)
    }
  }  
}) %>%
  bind_rows()

dna_species <- dna_occurrence %>%
  filter(taxonRank == "species") %>%
  rename(species = scientificName) %>%
  mutate(AphiaID = as.numeric(str_extract(scientificNameID, "[0-9]+"))) %>%
  group_by(phylum, class, order, family, genus, species, AphiaID, site) %>%
  summarize(target_gene = paste0(sort(unique(target_gene)), collapse = ";")) %>%
  mutate(
    group = case_when(
      class %in% fish_classes ~ "fish",
      order %in% turtle_orders ~ "turtles",
      class %in% mammal_classes ~ "mammals",
      class %in% bird_class ~ "birds",
      class %in% amphibian_class ~ "amphibians",
      phylum %in% molluscs_phyla ~ "molluscs",
      phylum %in% algae_phyla ~ "algae",
      phylum %in% sponge_phyla ~ "sponges",
      phylum %in% jelly_phyla ~ "jellyfish",
      phylum %in% unicellular ~ "single-cell",
      phylum %in% fungi_phyla ~ "fungi",
      phylum %in% worms_phyla ~ "worms",
      phylum %in% filter_feeders ~ "filter-feeders",
      phylum %in% starfish_phyla ~ "starfish",
      class %in% copepod_class ~ "copepods",
      class %in% crustacean_class ~ "crustaceans"
    )
  )
  # %>%filter(!is.na(group) & species != "Homo sapiens")

dna_species_obis <- taxon(dna_species$AphiaID) %>%
  bind_rows(tibble(category = character())) %>%
  select(AphiaID = taxonID, redlist_category = category)

dna_species <- dna_species %>%
  left_join(dna_species_obis, by = "AphiaID") %>%
  mutate(source_dna = TRUE)

if (include_dna) {
  combined_species <- bind_rows(obis_species, dna_species)
} else {
  combined_species <- obis_species %>% mutate(target_gene = NA, source_dna = NA)
}

species <- combined_species %>%
  group_by(AphiaID, phylum, class, order, family, genus, species, site, group) %>%
  summarize(
    redlist_category = first(na.omit(redlist_category)),
    records = first(na.omit(records)),
    max_year = first(na.omit(max_year)),
    target_gene = first(na.omit(target_gene)),
    source_obis = first(na.omit(obis)),
    source_gbif = first(na.omit(gbif)),
    source_dna = first(na.omit(source_dna))
  ) %>%
  ungroup()
  
# Generate JSON and CSV species lists

for (site_name in sites) {
  site_list <- species %>%
    filter(site == site_name) %>%
    select(-site) %>%
    arrange(group, phylum, class, order, species)

  group_stats <- site_list %>% group_by(group) %>% summarize(n_species = n()) %>% data.table::transpose(make.names = TRUE)
  if (nrow(group_stats) == 0) {
    group_stats <- NULL
  }
  
  site_stats <- list(
    groups = group_stats %>% unbox(),
    redlist = site_list %>% filter(!is.na(redlist_category)) %>% summarize(n_species = n()) %>% pull(n_species) %>% unbox(),
    source = list(
      obis = site_list %>% filter(source_obis) %>% summarize(n_species = n()) %>% pull(n_species) %>% unbox(),
      gbif = site_list %>% filter(source_gbif) %>% summarize(n_species = n()) %>% pull(n_species) %>% unbox(),
      db = site_list %>% filter(source_gbif | source_obis) %>% summarize(n_species = n()) %>% pull(n_species) %>% unbox(),
      edna = site_list %>% filter(source_dna) %>% summarize(n_species = n()) %>% pull(n_species) %>% unbox(),
      both = site_list %>% filter(source_dna & (source_obis | source_gbif)) %>% summarize(n_species = n()) %>% pull(n_species) %>% unbox()
    )
  )
  
  json = toJSON(list(
    created = unbox(strftime(as.POSIXlt(Sys.time(), "UTC"), "%Y-%m-%dT%H:%M:%S")),
    species = site_list,
    stats = site_stats
  ), pretty = TRUE)
  write(json, glue("lists/json/{site_name}.json"))
  write.csv(site_list, glue("lists/csv/{site_name}.csv"), row.names = FALSE, na = "")
}
