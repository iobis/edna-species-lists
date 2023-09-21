library(jsonlite)
library(glue)
library(purrr)
library(dplyr)
library(stringr)

markers <- c("16s", "coi", "mifish", "mimammal", "teleo")
fish_classes <- c("Actinopteri", "Cladistii", "Coelacanthi", "Elasmobranchii", "Holocephali", "Myxini", "Petromyzonti", "Teleostei")
turtle_orders <- c("Testudines")
mammal_classes <- c("Mammalia")

# Read OBIS species lists from https://github.com/iobis/mwhs-obis-species

sites <- fromJSON("https://raw.githubusercontent.com/iobis/mwhs-obis-species/master/lists/sites.json")

obis_species <- map(sites, function(site_name) {
  site_list <- fromJSON(glue("https://raw.githubusercontent.com/iobis/mwhs-obis-species/master/lists/{site_name}.json"))$species
  site_list$site <- site_name
  return(site_list)
}) %>%
  bind_rows() %>%
  mutate(source_obis = TRUE)

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
      order %in% turtle_orders ~ "turtle",
      class %in% mammal_classes ~ "mammal"
    )
  ) %>%
  filter(!is.na(group) & species != "Homo sapiens")

dna_species_obis <- taxon(dna_species$AphiaID) %>%
  bind_rows(tibble(category = character())) %>%
  select(AphiaID = taxonID, redlist_category = category)

dna_species <- dna_species %>%
  left_join(dna_species_obis, by = "AphiaID") %>%
  mutate(source_dna = TRUE)

species <- bind_rows(obis_species, dna_species) %>%
  group_by(AphiaID, phylum, class, order, family, genus, species, site, group) %>%
  summarize(
    redlist_category = first(na.omit(redlist_category)),
    records = first(na.omit(records)),
    max_year = first(na.omit(max_year)),
    target_gene = first(na.omit(target_gene)),
    source_obis = first(na.omit(source_obis)),
    source_dna = first(na.omit(source_dna))
  ) %>%
  ungroup()
  
# Generate JSON and CSV species lists

for (site_name in sites) {
  site_list <- species %>%
    filter(site == site_name) %>%
    select(-site) %>%
    arrange(group, phylum, class, order, species)
  json = toJSON(list(
    created = unbox(strftime(as.POSIXlt(Sys.time(), "UTC"), "%Y-%m-%dT%H:%M:%S")),
    species = site_list
  ), pretty = TRUE)
  write(json, glue("lists/json/{site_name}.json"))
  write.csv(site_list, glue("lists/csv/{site_name}.csv"), row.names = FALSE, na = "")
}
