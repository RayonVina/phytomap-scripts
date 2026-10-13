# PHYTO-MAP DATA INGESTION ──────────────────────────────────────────────────
# Dataset: Wietkamp-et-al-2020
# Fernando Rayón Viña
# last update: 20261013 ADT
#
# Wietkamp, S. et al. (2020)
# Distribution and abundance of azaspiracid-producing dinophyte species
# and their toxins in North Atlantic and North Sea waters in summer 2018.
# PLOS ONE, 15(6): e0235015. https://doi.org/10.1371/journal.pone.0235015

gc()
rm(list = ls())

# LOAD LIBRARIES ────────────────────────────────────────────────────────────

library(tidyverse)
library(readxl)
library(lubridate)
library(janitor)
library(DBI)
library(RSQLite)
library(this.path)

# LOAD PATHS ────────────────────────────────────────────────────────────────

script_dir <- this.dir()

setwd(script_dir)

fn_xlsx <- file.path(script_dir, "original-data-wietkamp-et-al-2018.xlsx")

fn_ctd <- file.path(script_dir, "HE516-CTD.tab")

fn_s2 <- file.path(script_dir, "s2_microscopy_verified.csv")

# IMPORT DATA ───────────────────────────────────────────────────────────────

station_raw <- read_excel(fn_xlsx, sheet = "station-environment-data") |>
  clean_names() |>
  remove_empty(c("rows", "cols"))

s2_raw <- read_csv(
  fn_s2,
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
) |>
  clean_names()

if (nrow(problems(s2_raw)) > 0 ||
    !all(c("station_id", "microscopy_cells_l") %in% names(s2_raw))) {
  stop("Invalid structure in the verified S2 microscopy table.")
}
if (nrow(s2_raw) != 75 ||
    anyDuplicated(s2_raw$station_id) ||
    !setequal(s2_raw$station_id, paste0("St.", 1:75))) {
  stop("The verified S2 table must contain St.1 to St.75 exactly once.")
}

ctd_raw <- read_delim(
  fn_ctd,
  delim = "\t",
  skip = 103,
  show_col_types = FALSE,
  trim_ws = TRUE
) |>
  clean_names()

# PANGAEA COORDINATES & TIMESTAMPS ──────────────────────────────────────────

pangaea_casts <- ctd_raw |>
  group_by(event) |>
  slice_head(n = 1) |>
  ungroup() |>
  transmute(
    event,
    latitude,
    longitude,
    datetime = date_time
  )

# STATION METADATA ──────────────────────────────────────────────────────────

station_index <- station_raw |>
  mutate(idx = row_number()) |>
  select(
    idx,
    station_id,
    latitude,
    longitude,
    water_column_depth_m,
    year,
    month,
    day
  )

pangaea_casts_ordered <- pangaea_casts |>
  arrange(datetime) |>
  mutate(idx = row_number())

ctd_dcm <- ctd_raw |>
  filter(!is.na(fluorometer_arbitrary_units)) |>
  group_by(event) |>
  slice_max(fluorometer_arbitrary_units, n = 1, with_ties = FALSE) |>
  ungroup() |>
  transmute(event, dcm_depth = depth_water_m)

station_metadata <- station_index |>
  left_join(pangaea_casts_ordered, by = "idx") |>
  left_join(ctd_dcm, by = "event") |>
  mutate(
    latitude = coalesce(latitude.y, latitude.x),
    longitude = coalesce(longitude.y, longitude.x),
    year = year(datetime),
    month = month(datetime),
    day = day(datetime),
    datetime_iso = format(datetime, "%Y-%m-%dT%H:%M:%SZ"),
    depth_label = "wc",
    depth_label_min = pmin(3, 10, dcm_depth),
    depth_label_max = pmax(3, 10, dcm_depth),
    # Round in decimetres to preserve decimal ties-to-even.
    depth = round(
      (round(depth_label_min * 10) + round(depth_label_max * 10)) / 2
    ) / 10
  ) |>
  select(
    station_id,
    event,
    latitude,
    longitude,
    datetime = datetime_iso,
    year,
    month,
    day,
    depth,
    depth_label,
    depth_label_min,
    depth_label_max
  )

# MICROSCOPY COUNT DATA ─────────────────────────────────────────────────────

# Verified transcription of the Microscopy column in S2:
# n.d. = not detected; n.a. = no data available.
s2_raw <- s2_raw |>
  mutate(microscopy_cells_l = str_trim(microscopy_cells_l))

if (anyNA(s2_raw$microscopy_cells_l) ||
    any(!str_detect(s2_raw$microscopy_cells_l, "^(n\\.d\\.|n\\.a\\.|[0-9]+)$")) ||
    sum(s2_raw$microscopy_cells_l == "n.d.") != 10 ||
    sum(s2_raw$microscopy_cells_l == "n.a.") != 30) {
  stop("Unexpected microscopy values or missing-data marks in S2.")
}

micro_long <- s2_raw |>
  filter(microscopy_cells_l != "n.a.") |>
  transmute(
    station_id,
    abundance = parse_double(
      if_else(microscopy_cells_l == "n.d.", "0", microscopy_cells_l)
    ),
    taxon = "amphidomataceae"
  )

build_records <- function(df, method) {
  df |>
    left_join(
      station_metadata |>
        select(
          station_id,
          latitude,
          longitude,
          datetime,
          year,
          month,
          day,
          depth,
          depth_label,
          depth_label_min,
          depth_label_max
        ),
      by = "station_id"
    ) |>
    filter(
      !is.na(latitude),
      !is.na(longitude),
      !is.na(depth)
    ) |>
    mutate(
      taxon = taxon |>
        str_replace_all("_", " ") |>
        str_to_lower() |>
        str_squish(),
      internal_id = paste0(station_id, "-", method),
      source = "wietkamp-2020",
      n_notes = case_when(
        abundance == 0 & method == "microscopy" ~
          paste0(
            "S2 n.d. = not detected; abundance=0 denotes microscopy ",
            "non-detection, not absolute absence; microscopy detection ",
            "limit approximately 20-50 cells/L."
          ),
        TRUE ~ NA_character_
      ),
      s_count_method = "direct count",
      s_enumeration_method = "microscopy",
      s_sampling_method = "niskin bottle (10 L)",
      s_concentration_method = "gravity filtration",
      s_volume_sampled = NA_real_,
      s_microscopy_method = if_else(
        method == "microscopy",
        "light microscope",
        NA_character_
      ),
      s_light_microscope = if_else(
        method == "microscopy",
        "inverted",
        NA_character_
      ),
      s_counting_chamber = if_else(
        method == "microscopy",
        "utermohl",
        NA_character_
      ),
      s_magnification = if_else(
        method == "microscopy",
        "640x-1000x",
        NA_character_
      )
    ) |>
    select(
      internal_id,
      source,
      datetime,
      year,
      month,
      day,
      latitude,
      longitude,
      depth,
      depth_label,
      depth_label_min,
      depth_label_max,
      abundance,
      taxon,
      n_notes,
      s_count_method,
      s_enumeration_method,
      s_sampling_method,
      s_concentration_method,
      s_volume_sampled,
      s_microscopy_method,
      s_light_microscope,
      s_counting_chamber,
      s_magnification
    )
}

# TRANSFORM COUNT DATA ──────────────────────────────────────────────────────

count_data_micro <- build_records(micro_long, "microscopy")

count_data <- count_data_micro

if (nrow(count_data) != 45 ||
    anyDuplicated(count_data$internal_id) ||
    anyNA(count_data$abundance) ||
    sum(count_data$abundance == 0) != 10 ||
    sum(count_data$abundance > 0) != 35) {
  stop("Expected 45 microscopy records: 35 quantified and 10 non-detections.")
}

# PREPARE DATABASE RECORDS ─────────────────────────────────────────────────────
count_data <- count_data |> mutate(across(where(is.factor), as.character))

column_aliases <- c(
  taxon = "t_taxon",
  abundance_l = "abundance",
  abundance_L = "abundance",
  depth_m = "depth",
  datetime_utc = "datetime",
  count_method = "s_count_method",
  taxon_original = "t_taxon",
  taxon_clean = "t_taxon_clean",
  notes = "n_notes"
)

for (old_name in intersect(names(column_aliases), names(count_data))) {
  new_name <- column_aliases[[old_name]]
  if (new_name %in% names(count_data)) {
    old_value <- as.character(count_data[[old_name]])
    new_value <- as.character(count_data[[new_name]])
    conflict <- !is.na(old_value) & !is.na(new_value) & old_value != new_value
    if (any(conflict)) {
      stop(sprintf("Conflicting columns: %s and %s.", old_name, new_name))
    }
    missing_value <- is.na(count_data[[new_name]])
    count_data[[new_name]][missing_value] <- count_data[[old_name]][
      missing_value
    ]
    count_data[[old_name]] <- NULL
  } else {
    names(count_data)[names(count_data) == old_name] <- new_name
  }
}

protocol_fields <- c(
  "s_count_method",
  "s_enumeration_method",
  "s_sampling_method",
  "s_sampling_method_notes",
  "s_concentration_method",
  "s_concentration_method_notes",
  "s_volume_sampled",
  "s_microscopy_method",
  "s_microscopy_notes",
  "s_preservation",
  "s_volume_concentrated",
  "s_volume_counted",
  "s_magnification",
  "s_counting_chamber"
)
count_data <- count_data |>
  mutate(across(any_of(protocol_fields), as.character)) |>
  mutate(across(where(is.factor), as.character))

if ("datetime" %in% names(count_data)) {
  count_data$datetime <- as.character(count_data$datetime)
}
if (anyDuplicated(names(count_data))) {
  stop("Duplicate output column names.")
}
required_fields <- c("source", "t_taxon", "abundance")
if (!all(required_fields %in% names(count_data))) {
  stop("Output must contain source, t_taxon, and abundance.")
}

# EXPORT TO DATABASE ──────────────────────────────────────────────────────────
db_path <- Sys.getenv("PHYTOMAP_DB")

con <- dbConnect(RSQLite::SQLite(), db_path, flags = RSQLite::SQLITE_RW)
tryCatch(
  {
    if (!dbExistsTable(con, "abundance_raw")) {
      stop("The working database must already contain abundance_raw.")
    }
    dbWithTransaction(con, {
      new_fields <- setdiff(
        names(count_data),
        dbListFields(con, "abundance_raw")
      )
      for (field in new_fields) {
        field_type <- dbDataType(con, count_data[[field]])
        dbExecute(
          con,
          paste(
            "ALTER TABLE abundance_raw ADD COLUMN",
            dbQuoteIdentifier(con, field),
            field_type
          )
        )
      }
      dbWriteTable(
        con,
        name = "abundance_raw",
        value = count_data,
        append = TRUE,
        row.names = FALSE
      )
    })
  },
  finally = dbDisconnect(con)
)
rm(con)
