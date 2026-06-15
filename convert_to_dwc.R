# convert_to_dwc.R
# Converts Tampa Bay seagrass transect data to Darwin Core Archive format
# for OBIS upload.
#
# Outputs:
#   dwc/occurrence.csv
#   dwc/emof.csv        (ExtendedMeasurementOrFact extension)
#
# Requires:
#   install.packages(c("dplyr", "tidyr", "readr", "lubridate", "sf", "worrms"))
#
# Before uploading to OBIS:
#   - Verify all VERIFY comments below
#   - Add EML metadata (required by OBIS — see https://obis.org/manual/dataformat/)
#   - Validate the archive at https://obis.org/manual/processing/

library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(sf)
library(worrms)

# ---------------------------------------------------------------------------
# Configuration — adjust before running
# ---------------------------------------------------------------------------

INSTITUTION_CODE <- "TBEP"
DATASET_NAME     <- "Tampa Bay Seagrass Transect Monitoring"  # VERIFY
COLLECTION_CODE  <- "seagrass"
DATASET_ID       <- ""        # OBIS will assign this UUID after registration
LICENSE          <- "https://creativecommons.org/licenses/by/4.0/"  # VERIFY

# Taxon to attach "No Cover" absence records to.
# "Alismatales" covers all seagrass orders found in Tampa Bay.
ABSENCE_TAXON       <- "Alismatales"
ABSENCE_TAXON_RANK  <- "order"

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------

# trnsct is assumed already loaded in your environment.
# VERIFY: replace `transect_sf` with the name of your sf object and
# VERIFY: replace `Transect` with the column in that sf object that
#         matches the Transect column in trnsct.
transect_locs <- transect_sf |>
  mutate(
    decimalLongitude = st_coordinates(geometry)[, 1],
    decimalLatitude  = st_coordinates(geometry)[, 2]
  ) |>
  st_drop_geometry() |>
  select(Transect, decimalLatitude, decimalLongitude)

dat <- trnsct |>
  left_join(transect_locs, by = "Transect")

# ---------------------------------------------------------------------------
# 2. Resolve taxa to WoRMS
# ---------------------------------------------------------------------------

presence_species <- dat |>
  filter(Species != "No Cover", !is.na(Species), Species != "") |>
  distinct(Species) |>
  pull(Species)

message("Looking up ", length(presence_species), " taxa in WoRMS...")
worms_raw        <- wm_records_names(name = presence_species, fuzzy = FALSE, marine_only = TRUE)
names(worms_raw) <- presence_species

species_lookup <- bind_rows(worms_raw, .id = "Species") |>
  filter(status == "accepted") |>
  group_by(Species) |>
  slice(1) |>
  ungroup() |>
  transmute(
    Species,
    scientificName   = scientificname,
    scientificNameID = paste0("urn:lsid:marinespecies.org:taxname:", AphiaID),
    taxonRank        = tolower(rank),
    kingdom, phylum, class, order, family, genus
  )

unmatched <- setdiff(presence_species, species_lookup$Species)
if (length(unmatched) > 0)
  warning("No accepted WoRMS match for: ", paste(unmatched, collapse = ", "),
          "\nThese rows will be excluded. Consider fuzzy=TRUE or manual lookup.")

# Resolve the absence taxon
message("Looking up absence taxon '", ABSENCE_TAXON, "' in WoRMS...")
absence_worms <- wm_records_names(name = ABSENCE_TAXON, fuzzy = FALSE)[[1]] |>
  filter(status == "accepted") |>
  slice(1)

absence_taxon_row <- tibble(
  scientificName   = absence_worms$scientificname,
  scientificNameID = paste0("urn:lsid:marinespecies.org:taxname:", absence_worms$AphiaID),
  taxonRank        = tolower(absence_worms$rank),
  kingdom          = absence_worms$kingdom,
  phylum           = absence_worms$phylum,
  class            = absence_worms$class,
  order            = absence_worms$order,
  family           = NA_character_,
  genus            = NA_character_
)

# ---------------------------------------------------------------------------
# 3. Shared event/location columns
# ---------------------------------------------------------------------------

# Depth: source values are negative centimetres; DwC expects positive metres
depth_to_m <- function(d) if_else(d != 0, abs(d) / 100, NA_real_)

base_cols <- dat |>
  mutate(
    occurrenceID         = paste(INSTITUTION_CODE, COLLECTION_CODE, ID, sep = ":"),
    eventID              = paste(INSTITUTION_CODE, COLLECTION_CODE, "event",
                                 Transect,
                                 as.Date(ymd_hms(ObservationDate, quiet = TRUE)),
                                 sep = ":"),
    eventDate            = format(ymd_hms(ObservationDate, quiet = TRUE), "%Y-%m-%dT%H:%M:%S"),
    year                 = year(ymd_hms(ObservationDate,  quiet = TRUE)),
    month                = month(ymd_hms(ObservationDate, quiet = TRUE)),
    day                  = day(ymd_hms(ObservationDate,   quiet = TRUE)),
    minimumDepthInMeters = depth_to_m(Depth),
    maximumDepthInMeters = minimumDepthInMeters,
    country              = "United States",
    countryCode          = "US",
    stateProvince        = "Florida",
    waterBody            = BaySegment,
    locality             = paste(BaySegment, Transect),
    locationID           = paste(INSTITUTION_CODE, COLLECTION_CODE, "loc", Transect, sep = ":"),
    samplingProtocol     = "seagrass transect point-intercept survey",
    institutionCode      = INSTITUTION_CODE,
    datasetName          = DATASET_NAME,
    collectionCode       = COLLECTION_CODE,
    datasetID            = DATASET_ID,
    geodeticDatum        = "EPSG:4326",
    license              = LICENSE,
    recordedBy           = Crew
  )

shared_fields <- c(
  "occurrenceID", "eventID",
  "eventDate", "year", "month", "day",
  "decimalLatitude", "decimalLongitude", "geodeticDatum",
  "minimumDepthInMeters", "maximumDepthInMeters",
  "country", "countryCode", "stateProvince", "waterBody", "locality", "locationID",
  "samplingProtocol", "institutionCode", "datasetName", "collectionCode", "datasetID",
  "license", "recordedBy"
)

# ---------------------------------------------------------------------------
# 4. Build occurrence core
# ---------------------------------------------------------------------------

# 4a. Presence rows
presence_occ <- base_cols |>
  filter(Species != "No Cover", !is.na(Species), Species != "") |>
  inner_join(species_lookup, by = "Species") |>     # inner_join drops unmatched taxa
  mutate(
    basisOfRecord    = "HumanObservation",
    occurrenceStatus = "present"
  ) |>
  select(all_of(shared_fields),
         basisOfRecord, occurrenceStatus,
         scientificName, scientificNameID, taxonRank,
         kingdom, phylum, class, order, family, genus)

# 4b. Absence rows ("No Cover" — one record per point, taxon = Alismatales)
absence_occ <- base_cols |>
  filter(Species == "No Cover") |>
  mutate(
    basisOfRecord    = "HumanObservation",
    occurrenceStatus = "absent",
    scientificName   = absence_taxon_row$scientificName,
    scientificNameID = absence_taxon_row$scientificNameID,
    taxonRank        = absence_taxon_row$taxonRank,
    kingdom          = absence_taxon_row$kingdom,
    phylum           = absence_taxon_row$phylum,
    class            = absence_taxon_row$class,
    order            = absence_taxon_row$order,
    family           = absence_taxon_row$family,
    genus            = absence_taxon_row$genus
  ) |>
  select(all_of(shared_fields),
         basisOfRecord, occurrenceStatus,
         scientificName, scientificNameID, taxonRank,
         kingdom, phylum, class, order, family, genus)

occurrence <- bind_rows(presence_occ, absence_occ)

# ---------------------------------------------------------------------------
# 5. Build ExtendedMeasurementOrFact (eMoF) extension
# ---------------------------------------------------------------------------

# Extract the leading integer from SpeciesAbundance, e.g. "2 = 6%-25%" -> 2
extract_cover_code <- function(x) as.integer(sub("^(\\d+).*", "\\1", trimws(x)))

emof_src <- base_cols |>
  filter(Species != "No Cover", !is.na(Species), Species != "") |>
  inner_join(select(species_lookup, Species), by = "Species") |>
  mutate(cover_code = extract_cover_code(SpeciesAbundance))

# Helper: pivot one measurement per occurrence row into eMoF long format
make_emof <- function(df, type, type_id, value_col, unit, unit_id) {
  df |>
    filter(!is.na(.data[[value_col]]), .data[[value_col]] != "") |>
    transmute(
      occurrenceID      = occurrenceID,
      measurementType   = type,
      measurementTypeID = type_id,
      measurementValue  = as.character(.data[[value_col]]),
      measurementUnit   = unit,
      measurementUnitID = unit_id
    )
}

# P01 vocabulary codes — VERIFY these against https://vocab.nerc.ac.uk/
emof <- bind_rows(
  make_emof(emof_src,
    "seagrass percent cover (Braun-Blanquet class)",
    "",           # No standard P01 code for BB class; leave blank or use a local term
    "cover_code", "Braun-Blanquet scale", ""),

  make_emof(emof_src |> filter(BladeLength_Avg > 0),
    "seagrass blade length arithmetic mean",
    "",           # VERIFY P01 code
    "BladeLength_Avg", "mm",
    "http://vocab.nerc.ac.uk/collection/P06/current/UXMM/"),

  make_emof(emof_src |> filter(ShootDensity_Avg > 0),
    "seagrass shoot density arithmetic mean",
    "",           # VERIFY P01 code
    "ShootDensity_Avg", "shoots m-2",
    "http://vocab.nerc.ac.uk/collection/P06/current/UPMS/"),

  make_emof(emof_src |> filter(!is.na(EpiphyteDensity), EpiphyteDensity != ""),
    "epiphyte density (qualitative)",
    "",
    "EpiphyteDensity", "", ""),

  make_emof(emof_src |> filter(!is.na(SedimentType), SedimentType != ""),
    "sediment type",
    "",
    "SedimentType", "", "")
)

# ---------------------------------------------------------------------------
# 6. Write output
# ---------------------------------------------------------------------------

dir.create("dwc", showWarnings = FALSE)
write_csv(occurrence, "dwc/occurrence.csv", na = "")
write_csv(emof,       "dwc/emof.csv",       na = "")

message(
  "\nDone.\n",
  "  occurrence.csv : ", nrow(occurrence), " rows  ",
  "(", nrow(presence_occ), " presence / ", nrow(absence_occ), " absence)\n",
  "  emof.csv       : ", nrow(emof), " rows\n",
  "\nNext steps:\n",
  "  1. Check species_lookup to confirm WoRMS matches are correct\n",
  "  2. Fill in NERC P01 codes for eMoF measurementTypeID fields\n",
  "  3. Register the dataset at https://obis.org and get a datasetID\n",
  "  4. Author EML metadata (required alongside the CSVs)\n",
  "  5. Validate at https://obis.org/manual/processing/ before upload"
)
