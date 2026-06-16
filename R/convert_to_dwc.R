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
library(tbeptools)
library(here)
library(obistools)

transect_sf <- trnpts
# trnsct <- read_transect(raw = T)

# write.csv(trnsct, here('data', 'trnsct.csv'), row.names = FALSE)
trnsct <- read.csv(here('data', 'trnsct.csv'))

# ---------------------------------------------------------------------------
# Configuration — adjust before running
# ---------------------------------------------------------------------------

INSTITUTION_CODE <- "TBEP"
DATASET_NAME     <- "Tampa Bay Interagency Seagrass Monitoring Program"  # VERIFY
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
  filter(Metermark == 0) |>
  mutate(
    decimalLongitude = st_coordinates(geometry)[, 1],
    decimalLatitude  = st_coordinates(geometry)[, 2]
  ) |>
  st_drop_geometry() |>
  select(Transect = TRAN_ID, decimalLatitude, decimalLongitude)

dat <- trnsct |>
  left_join(transect_locs, by = "Transect") |>
  mutate(
    Species = gsub('^.*\\:\\s|\\n$', '', Species),
    Species = gsub('\\s+spp\\.$', '', Species),
    Species = gsub('intestinales$', 'intestinalis', Species),
    # Species = gsub('(^Ulva).*', '\\1', Species),
    # Map informal drift algae field codes to accepted WoRMS taxon names
    Species = case_when(
      trimws(tolower(Species)) == "drift brown" ~ "Phaeophyceae",
      trimws(tolower(Species)) == "drift green" ~ "Chlorophyta",
      trimws(tolower(Species)) == "drift reds"  ~ "Rhodophyta",
      trimws(tolower(Species)) == "drift red" ~ "Rhodophyta",
      Species == 'Lyngbya/Dapis' ~ 'Dapis',
      Species == 'Halodule' ~ 'Halodule wrightii',
      Species == 'Thalassia' ~ 'Thalassia testudinum',
      Species == 'Ruppia' ~ 'Ruppia maritima',
      Species == 'Syringodium' ~ 'Syringodium filiforme',
      Species == 'Ulva fasciata' ~ 'Ulva lactuca',
      TRUE ~ Species
    )
  ) |>
  group_by(IDall) |>
  mutate(transect_date = min(as.Date(ymd_hms(ObservationDate, quiet = TRUE)), na.rm = TRUE)) |>
  ungroup()

# ---------------------------------------------------------------------------
# 2. Resolve taxa to WoRMS
# ---------------------------------------------------------------------------

presence_species <- dat |>
  filter(Species != "No Cover", !is.na(Species), Species != "") |>
  distinct(Species) |>
  pull(Species)

message("Looking up ", length(presence_species), " taxa in WoRMS...")
worms_raw        <- wm_records_names(name = presence_species, fuzzy = TRUE, marine_only = TRUE)
names(worms_raw) <- presence_species

# Pinned AphiaIDs for taxa where WoRMS returns multiple accepted matches.
# Maps the Species name (as it appears in dat) to the correct AphiaID.
aphia_overrides <- c(
  "Chondria"     = 143906L,
  "Digenea"      = 143909L,
  "Polysiphonia" = 143853L,
  "Ulva"         = 144296L,
  "Halophila"    = 144192L
)

species_lookup <- bind_rows(worms_raw, .id = "Species") |>
  filter(status == "accepted") |>
  mutate(override_id = aphia_overrides[Species]) |>
  filter(is.na(override_id) | AphiaID == override_id) |>
  group_by(Species) |>
  slice(1) |>
  ungroup() |>
  select(-override_id) |>
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
    parentEventID        = paste(INSTITUTION_CODE, COLLECTION_CODE, "event",
                                 Transect,
                                 transect_date,
                                 sep = ":"),
    eventID              = paste(INSTITUTION_CODE, COLLECTION_CODE, "event",
                                 Transect,
                                 as.Date(ymd_hms(ObservationDate, quiet = TRUE)),
                                 Site,
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
    locationID           = paste(INSTITUTION_CODE, COLLECTION_CODE, "loc", Transect, Site, sep = ":"),
    samplingProtocol     = "seagrass transect survey",
    institutionCode      = INSTITUTION_CODE,
    datasetName          = DATASET_NAME,
    collectionCode       = COLLECTION_CODE,
    datasetID            = DATASET_ID,
    geodeticDatum        = "EPSG:4326",
    license              = LICENSE,
    recordedBy           = MonitoringAgency,
    eventType            = "Point"
  )

event_fields <- c(
  "eventID", "parentEventID", "eventType",
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
  select(occurrenceID, eventID,
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
  select(occurrenceID, eventID,
         basisOfRecord, occurrenceStatus,
         scientificName, scientificNameID, taxonRank,
         kingdom, phylum, class, order, family, genus)

occurrence <- bind_rows(presence_occ, absence_occ)

# ---------------------------------------------------------------------------
# 5. Build ExtendedMeasurementOrFact (eMoF) extension
# ---------------------------------------------------------------------------

emof_src <- base_cols |>
  filter(Species != "No Cover", !is.na(Species), Species != "") |>
  inner_join(select(species_lookup, Species), by = "Species") |>
  mutate(
    cover_value   = trimws(sub("\\s*=.*$", "", SpeciesAbundance)),
    cover_remarks = ifelse(grepl("=", SpeciesAbundance, fixed = TRUE),
                           trimws(sub("^[^=]+=\\s*", "", SpeciesAbundance)),
                           NA_character_)
  )

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

# P01 vocabulary codes from https://vocab.nerc.ac.uk/collection/P01/current/
#
# cover_code    — Braun-Blanquet ordinal class (0-5), not a true percentage.
#                 PCOV7736 ("Proportion coverage...") is the nearest accepted term
#                 but expects % values; left blank here. Request a BB-specific term
#                 at https://github.com/nvs-vocabs/OBISVocabs/issues if needed.
# BladeLength   — OBSINDLX: "Length of biological entity specified elsewhere" (mm)
# ShootDensity  — SDBIOL02: "Abundance of biological entity specified elsewhere
#                 per unit area of the bed" (n m-2)
# EpiphyteDensity / SedimentType — no accepted P01 code; left blank.
#                 Request terms via https://github.com/nvs-vocabs/OBISVocabs/issues

emof <- bind_rows(
  emof_src |>
    filter(!is.na(cover_value), cover_value != "") |>
    transmute(
      occurrenceID       = occurrenceID,
      measurementType    = "seagrass percent cover (Braun-Blanquet scale)",
      measurementTypeID  = "",
      measurementValue   = cover_value,
      measurementUnit    = "Braun-Blanquet scale",
      measurementUnitID  = "",
      measurementRemarks = cover_remarks
    ),

  make_emof(emof_src |> filter(BladeLength_Avg > 0),
    "seagrass blade length arithmetic mean",
    "http://vocab.nerc.ac.uk/collection/P01/current/OBSINDLX/",
    "BladeLength_Avg", "mm",
    "http://vocab.nerc.ac.uk/collection/P06/current/UXMM/"),

  make_emof(emof_src |> filter(ShootDensity_Avg > 0),
    "seagrass shoot density arithmetic mean",
    "http://vocab.nerc.ac.uk/collection/P01/current/SDBIOL02/",
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
# 6. Build event core
# ---------------------------------------------------------------------------

# Child events: one row per meter mark visit
child_events <- base_cols |>
  distinct(eventID, .keep_all = TRUE) |>
  select(all_of(event_fields))

# Parent events: one row per transect visit; depth varies per point so omit it
parent_events <- base_cols |>
  distinct(parentEventID, .keep_all = TRUE) |>
  transmute(
    eventID              = parentEventID,
    parentEventID        = NA_character_,
    eventType            = "Transect",
    eventDate            = format(transect_date, "%Y-%m-%d"),
    year                 = year(transect_date),
    month                = month(transect_date),
    day                  = day(transect_date),
    decimalLatitude, decimalLongitude, geodeticDatum,
    minimumDepthInMeters = NA_real_,
    maximumDepthInMeters = NA_real_,
    country, countryCode, stateProvince, waterBody, locality,
    locationID           = paste(institutionCode, collectionCode, "loc", Transect, sep = ":"),
    samplingProtocol, institutionCode, datasetName, collectionCode, datasetID, license, recordedBy
  )

event <- bind_rows(parent_events, child_events)

# ---------------------------------------------------------------------------
# 7. Write output
# ---------------------------------------------------------------------------

dir.create("dwc", showWarnings = FALSE)
write_csv(event,      "dwc/event.csv",      na = "")
write_csv(occurrence, "dwc/occurrence.csv", na = "")
write_csv(emof,       "dwc/emof.csv",       na = "")

message(
  "\nDone.\n",
  "  event.csv      : ", nrow(event), " rows  ",
  "(", nrow(parent_events), " transect visits / ", nrow(child_events), " meter mark visits)\n",
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

# qc with obistools

# # check species - note that this only checks those that are ambiguous from worms
# # not those that are have more than one in your actual data
# match_taxa(unique(occurrence$scientificName))

# # check all required fiels are presnt
# chk <- check_fields(event)

# # check coords are in correct region
# plot_map_leaflet(event)

# # check depth of records (within expected range)
# check_depth(event, report = T, depthmargin = 20)

# # check event dates
# check_eventdate(event)

# # check event ids
# check_eventids(event)

# # check extension event ids
# check_extension_eventids(event, emof)
# check_extension_eventids(event, occurrence)

# Hmisc::describe(event)
