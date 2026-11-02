### BGC-Argo NITRATE & pH Corrections DMQC Pipeline
### Integrated DMQC processing using SAGE v2 per-cycle parameters
### Dataset: Sage_v2DMQC_14Oct2026_jm.csv
### Updated: Distinct mapping for A-mode (ADJUSTED) and D-mode (FINAL)
### Updated: November 2026

# =========================================================================
# 1. ENVIRONMENT CLEANUP & PACKAGE INITIALIZATION
# =========================================================================

cat("\014")
rm(list = ls())

required_pkgs <- c("dplyr", "ggplot2", "lubridate", "gridExtra", "tidyverse", 
                   "ggeffects", "mapdata", "suncalc", "ncdf4", "oce", "terra", 
                   "sf", "stars", "gstat", "rnaturalearth", "rnaturalearthdata", 
                   "viridis", "knitr", "patchwork", "stringr", "interp")

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

options(timeout = 600)
options(HTTPUserAgent = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)")

# Helper function to convert numeric and character representation of fill values to NaN
convert_fill_to_nan <- function(df, fill_vals = c("-99999", "-99999.0", "99999", "99999.0", "-9999")) {
  df %>%
    mutate(
      across(where(is.numeric), ~ ifelse(. %in% as.numeric(fill_vals), NaN, .)),
      across(where(is.character), ~ {
        trimmed <- str_trim(.)
        ifelse(trimmed %in% fill_vals | trimmed %in% c("NA", "NaN"), NA_character_, .)
      })
    )
}

clean_corrupt_file <- function(filepath) {
  if (file.exists(filepath)) {
    if (file.info(filepath)$size == 0) {
      unlink(filepath)
      return(TRUE)
    }
  }
  return(FALSE)
}

download_if_missing <- function(local_path, primary_url, fallback_url = NULL) {
  if (clean_corrupt_file(local_path) || !file.exists(local_path)) {
    try(download.file(primary_url, local_path, quiet = TRUE, mode = "wb", method = "libcurl", cacheOK = TRUE), silent = TRUE)
    clean_corrupt_file(local_path)
  }
  if (!file.exists(local_path) && !is.null(fallback_url)) {
    try(download.file(fallback_url, local_path, quiet = TRUE, mode = "wb", method = "libcurl", cacheOK = TRUE), silent = TRUE)
    clean_corrupt_file(local_path)
  }
  return(file.exists(local_path))
}

# =========================================================================
# 2. CONFIGURATION & DIRECTORY MANAGEMENT
# =========================================================================

dir_main_root <- "C:/Users/Jennifer.McWhorter/Documents/GitHub/DMQC_NITRATE_PH"

dir_functions    <- file.path(dir_main_root, "Functions")
dir_output       <- file.path(dir_main_root, "Output")
dir_plot_base    <- file.path(dir_main_root, "Plots")
dir_script       <- file.path(dir_main_root, "Script")
dir_tables       <- file.path(dir_main_root, "Tables")
dir_profs        <- file.path(dir_main_root, "Tables/profs")

dir_data_root    <- file.path(dir_main_root, "Data")
combo_output_dir <- file.path(dir_output, "COMBO_NITRATE_PH_DMODE")

for (d in c(dir_functions, dir_output, dir_plot_base, dir_script, dir_tables, dir_profs, dir_data_root, combo_output_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

func.sources <- list.files(dir_functions, pattern = "\\.R$", full.names = TRUE)
for (f in func.sources) {
  source(f)
}

# --- IN-MEMORY DICTIONARY MAPPING FOR LAUNCH DATES ---
launch_date_dict <- tibble::tribble(
  ~float_ids, ~launch_date_str,
  4903622,    "2021-10-04 23:36:00",
  2904010,    "2025-05-22 04:32:00",
  2904011,    "2025-05-19 12:15:00",
  4903624,    "2021-09-25 09:51:00",
  4903625,    "2021-09-22 10:53:00",
  4903904,    "2025-06-24 07:00:00",
  6999992,    "2025-06-24 11:41:00",
  7902327,    "2026-05-21 04:14:00",
  3902693,    "2026-05-20 18:58:00",
  1902800,    "2026-05-29 09:33:00",
  7901009,    "2023-06-05 10:30:00",
  4903911,    "2026-05-15 21:52:00", 
  4903910,    "2026-05-20 18:00:00", 
  5907206,    "2026-05-29 09:33:00"
)

# --- LOAD SAGE V2 CALIBRATION DATA ---
calib_csv_path <- file.path(dir_data_root, "Sage_v2_outputs/Sage_v2DMQC_2Nov2026_jm.csv")

if (file.exists(calib_csv_path)) {
  calib_df <- read.csv(calib_csv_path, stringsAsFactors = FALSE) %>%
    convert_fill_to_nan() %>%
    mutate(
      float_ids    = as.numeric(WMOID),
      CYCLE_NUMBER = as.numeric(CYCLE)
    ) %>%
    rename(
      nitrate_gain   = NITRATE_GAIN,
      nitrate_offset = NITRATE_OFFSET,
      nitrate_drift  = NITRATE_DRIFT,
      ph_gain        = PH_GAIN,
      ph_offset      = PH_OFFSET,
      ph_drift       = PH_DRIFT
    ) %>%
    select(-any_of(c("LAUNCH_DATE", "DEPTH_COR_MIN", "DEPTH_COR_MAX"))) %>%
    left_join(launch_date_dict, by = "float_ids") %>%
    mutate(
      launch_date = lubridate::ymd_hms(launch_date_str, tz = "UTC")
    )
} else {
  stop(paste0("Specified calibration file not found at: ", calib_csv_path))
}

# --- LOAD BOTTLE VALIDATION DATA ---
path_bottle_val <- file.path(dir_data_root, "Bottle_vali/Master_Vali_24_July_2025.csv")

if (file.exists(path_bottle_val)) {
  Bottle_Vali <- read.csv(path_bottle_val, stringsAsFactors = FALSE) %>%
    convert_fill_to_nan() %>%
    mutate(
      Float_num    = as.numeric(Float_num),
      CYCLE_NUMBER = as.numeric(CYCLE_NUMBER),
      PRES         = as.numeric(PRES),
      NITRATE_Vali = as.numeric(NITRATE_Vali),
      PH_Vali      = as.numeric(PH_IN_SITU_TOTAL_Vali)
    )
  message("Successfully ingested bottle validation file -> ", path_bottle_val)
} else {
  warning("Bottle validation file not found at path: ", path_bottle_val)
  Bottle_Vali <- data.frame(
    Float_num = numeric(0), CYCLE_NUMBER = numeric(0), 
    PRES = numeric(0), NITRATE_Vali = numeric(0), PH_Vali = numeric(0)
  )
}

float_ids <- unique(calib_df$float_ids)
today_str <- format(Sys.Date(), "%Y-%m-%d")

url_argo_index <- "ftp://ftp.ifremer.fr/ifremer/argo/argo_bio-profile_index.txt"

# Standard BGC-Argo QC Color Palette
qc_color_scale <- c(
  "0" = "gray70", "1" = "#33A02C", "2" = "#B2DF8A", 
  "3" = "#FDBF6F", "4" = "#E31A1C", "5" = "#1F78B4", 
  "8" = "#CAB2D6", "9" = "black"
)

# Multi-mode plot aesthetics matching CHLA/DOXY pipeline:
# - R-mode: Red, dotted line (Raw)
# - A-mode: Brown, dashed line (Adjusted)
# - D-mode: Green, solid line (Final DMQC)
# - Bottle: Black points
mode_levels      <- c("R-mode", "A-mode", "D-mode", "Bottle")
color_mapping    <- c("R-mode" = "#E31A1C", "A-mode" = "#8C510A", "D-mode" = "#33A02C", "Bottle" = "#000000")
linetype_mapping <- c("R-mode" = "dotted",  "A-mode" = "dashed",  "D-mode" = "solid",   "Bottle" = "blank")
linewidth_map    <- c("R-mode" = 0.8,       "A-mode" = 1.2,       "D-mode" = 0.9,       "Bottle" = 0)
shape_mapping    <- c("R-mode" = NA,        "A-mode" = NA,        "D-mode" = NA,        "Bottle" = 16)

# =========================================================================
# 3. GLOBAL DATA INGESTION & PREPARATION
# =========================================================================

index_file_path <- file.path(dir_tables, 'index.txt')
download_if_missing(index_file_path, url_argo_index)

index      <- read.table(index_file_path, skip = 9, sep = ",", stringsAsFactors = FALSE)
files      <- as.character(index[, 1])
ident      <- strsplit(files, "/")
ident      <- matrix(unlist(ident), ncol = 4, byrow = TRUE)
dac        <- ident[, 1]
wod        <- as.character(ident[, 2])
prof_id    <- ident[, 4]
time_index <- index[, 2]

# =========================================================================
# 4. BATCH EXECUTION LOOP (INDIVIDUAL FLOAT PROCESSING)
# =========================================================================

for (woko in float_ids) {
  
  float_id   <- as.numeric(woko)
  active_wmo <- as.character(woko)
  
  dir_date_base   <- file.path(dir_plot_base, today_str)
  dir_float_plots <- file.path(dir_date_base, active_wmo)
  
  dir_nitrate_plots <- file.path(dir_float_plots, "NITRATE")
  dir_ph_plots      <- file.path(dir_float_plots, "PH")
  dir_qc_summary   <- file.path(dir_float_plots, "QC_SUMMARIES")
  
  for (p_dir in c(dir_nitrate_plots, dir_ph_plots, dir_qc_summary)) {
    if (!dir.exists(p_dir)) dir.create(p_dir, recursive = TRUE)
  }
  
  unlink(list.files(dir_qc_summary, pattern = "\\.png$", full.names = TRUE))
  
  cat("\n==========================================\n")
  cat(sprintf("  STARTING PROCESSING FOR FLOAT: %s\n", active_wmo))
  cat("==========================================\n\n")
  
  float_data <- NULL
  crit <- which(wod == active_wmo)
  
  if (length(crit) == 0) {
    warning(paste("No matching index records found on GDAC for WMO:", active_wmo))
    next
  }
  
  # --- STEP 1: DOWNLOAD & EXTRACT PROFILES ---
  for (kiwi in crit) {
    if (substr(prof_id[kiwi], 14, 14) == "D") next 
    if (is.na(time_index[kiwi])) next 
    
    www <- substr(prof_id[kiwi], 3, 13)
    profile_bio <- NULL
    profile_physic <- NULL
    
    meta_path <- file.path(dir_profs, paste0(woko, "_meta.nc"))
    meta_url  <- paste0('ftp://ftp.ifremer.fr/ifremer/argo/dac/', dac[kiwi], '/', woko, '/', woko, '_meta.nc')
    download_if_missing(meta_path, meta_url)
    
    bd_file <- file.path(dir_profs, paste0("BD", www, ".nc"))
    br_file <- file.path(dir_profs, paste0("BR", www, ".nc"))
    d_file  <- file.path(dir_profs, paste0("D", www, ".nc"))
    r_file  <- file.path(dir_profs, paste0("R", www, ".nc"))
    
    if (!file.exists(bd_file) && !file.exists(br_file)) {
      primary_bio_url  <- paste0("ftp://ftp.ifremer.fr/ifremer/argo/dac/", files[kiwi])
      fallback_bio_url <- paste0("ftp://ftp.ifremer.fr/ifremer/argo/dac/", dac[kiwi], "/", wod[kiwi], "/profiles/", prof_id[kiwi])
      download_if_missing(file.path(dir_profs, prof_id[kiwi]), primary_bio_url, fallback_bio_url)
    }
    if (!file.exists(bd_file) && !file.exists(br_file)) next
    
    if (!file.exists(d_file) && !file.exists(r_file)) {
      primary_d_url  <- paste0("ftp://ftp.ifremer.fr/ifremer/argo/dac/", dac[kiwi], "/", wod[kiwi], "/profiles/D", www, ".nc")
      fallback_d_url <- paste0("ftp://ftp.ifremer.fr/ifremer/argo/dac/", dac[kiwi], "/", wod[kiwi], "/D", www, ".nc")
      download_if_missing(d_file, primary_d_url, fallback_d_url)
    }
    if (!file.exists(d_file) && !file.exists(r_file)) {
      primary_r_url  <- paste0("ftp://ftp.ifremer.fr/ifremer/argo/dac/", dac[kiwi], "/", wod[kiwi], "/profiles/R", www, ".nc")
      fallback_r_url <- paste0("ftp://ftp.ifremer.fr/ifremer/argo/dac/", dac[kiwi], "/", wod[kiwi], "/R", www, ".nc")
      download_if_missing(r_file, primary_r_url, fallback_r_url)
    }
    if (!file.exists(d_file) && !file.exists(r_file)) next
    
    active_bio_file <- if (file.exists(bd_file)) bd_file else br_file
    active_phy_file <- if (file.exists(d_file)) d_file else r_file
    
    profile_bio    <- tryCatch(nc_open(active_bio_file, readunlim = FALSE, write = FALSE), error = function(e) NULL)
    profile_physic <- tryCatch(nc_open(active_phy_file, readunlim = FALSE, write = FALSE), error = function(e) NULL)
    
    if (is.null(profile_bio) || is.null(profile_physic)) {
      if (!is.null(profile_bio)) nc_close(profile_bio)
      if (!is.null(profile_physic)) nc_close(profile_physic)
      next
    }
    
    lat_val <- ncvar_get(profile_physic, "LATITUDE")[1]
    lon_val <- ncvar_get(profile_physic, "LONGITUDE")[1]
    if (is.na(lat_val) || is.na(lon_val)) {
      nc_close(profile_bio); nc_close(profile_physic); next
    }
    
    jd <- ncvar_get(profile_physic, "JULD")[1]
    origin <- as.POSIXct("1950-01-01 00:00:00", tz = "UTC") 
    time_val <- origin + jd * 3600 * 24
    if (is.na(time_val)) {
      nc_close(profile_bio); nc_close(profile_physic); next
    }
    
    pres <- as.vector(ncvar_get(profile_physic, "PRES"))
    pres_qc_raw <- ncvar_get(profile_physic, "PRES_QC")
    pres_qc_str <- if (is.matrix(pres_qc_raw)) apply(pres_qc_raw, 2, function(x) paste(x, collapse = "")) else paste(pres_qc_raw, collapse = "")
    pres_qc <- as.numeric(unlist(strsplit(pres_qc_str[1], "")))[1:length(pres)]
    
    temp <- rep(NaN, length(pres)); temp_qc <- rep(NA, length(pres))
    sal  <- rep(NaN, length(pres)); sal_qc  <- rep(NA, length(pres))
    
    if ("TEMP" %in% names(profile_physic$var)) {
      temp <- as.vector(ncvar_get(profile_physic, "TEMP"))
      temp_qc_raw <- ncvar_get(profile_physic, "TEMP_QC")
      temp_qc_str <- if (is.matrix(temp_qc_raw)) apply(temp_qc_raw, 2, function(x) paste(x, collapse = ""))[1] else paste(temp_qc_raw, collapse = "")
      temp_qc <- as.numeric(unlist(strsplit(temp_qc_str, "")))[1:length(pres)]
    }
    
    if ("PSAL" %in% names(profile_physic$var)) {
      sal <- as.vector(ncvar_get(profile_physic, "PSAL"))
      sal_qc_raw <- ncvar_get(profile_physic, "PSAL_QC")
      sal_qc_str <- if (is.matrix(sal_qc_raw)) apply(sal_qc_raw, 2, function(x) paste(x, collapse = ""))[1] else paste(sal_qc_raw, collapse = "")
      sal_qc <- as.numeric(unlist(strsplit(sal_qc_str, "")))[1:length(pres)]
    }
    
    nitrate <- rep(NaN, length(pres)); nitrate_adj <- rep(NaN, length(pres))
    nitrate_qc <- rep(NA, length(pres)); nitrate_adj_qc <- rep(NA, length(pres))
    ph      <- rep(NaN, length(pres)); ph_adj      <- rep(NaN, length(pres))
    ph_qc   <- rep(NA, length(pres)); ph_adj_qc   <- rep(NA, length(pres))
    
    if ("NITRATE" %in% names(profile_bio$var)) {
      nitrate <- as.vector(ncvar_get(profile_bio, "NITRATE"))
      nit_qc_raw <- ncvar_get(profile_bio, "NITRATE_QC")
      nit_qc_str <- if (is.matrix(nit_qc_raw)) apply(nit_qc_raw, 2, function(x) paste(x, collapse = ""))[1] else paste(nit_qc_raw, collapse = "")
      nitrate_qc <- as.numeric(unlist(strsplit(nit_qc_str, "")))[1:length(pres)]
    }
    if ("NITRATE_ADJUSTED" %in% names(profile_bio$var)) {
      nitrate_adj <- as.vector(ncvar_get(profile_bio, "NITRATE_ADJUSTED"))
      nit_adj_qc_raw <- ncvar_get(profile_bio, "NITRATE_ADJUSTED_QC")
      nit_adj_qc_str <- if (is.matrix(nit_adj_qc_raw)) apply(nit_adj_qc_raw, 2, function(x) paste(x, collapse = ""))[1] else paste(nit_adj_qc_raw, collapse = "")
      nitrate_adj_qc <- as.numeric(unlist(strsplit(nit_adj_qc_str, "")))[1:length(pres)]
    }
    if ("PH_IN_SITU_TOTAL" %in% names(profile_bio$var)) {
      ph <- as.vector(ncvar_get(profile_bio, "PH_IN_SITU_TOTAL"))
      ph_qc_raw <- ncvar_get(profile_bio, "PH_IN_SITU_TOTAL_QC")
      ph_qc_str <- if (is.matrix(ph_qc_raw)) apply(ph_qc_raw, 2, function(x) paste(x, collapse = ""))[1] else paste(ph_qc_raw, collapse = "")
      ph_qc <- as.numeric(unlist(strsplit(ph_qc_str, "")))[1:length(pres)]
    }
    if ("PH_IN_SITU_TOTAL_ADJUSTED" %in% names(profile_bio$var)) {
      ph_adj <- as.vector(ncvar_get(profile_bio, "PH_IN_SITU_TOTAL_ADJUSTED"))
      ph_adj_qc_raw <- ncvar_get(profile_bio, "PH_IN_SITU_TOTAL_ADJUSTED_QC")
      ph_adj_qc_str <- if (is.matrix(ph_adj_qc_raw)) apply(ph_adj_qc_raw, 2, function(x) paste(x, collapse = ""))[1] else paste(ph_adj_qc_raw, collapse = "")
      ph_adj_qc <- as.numeric(unlist(strsplit(ph_adj_qc_str, "")))[1:length(pres)]
    }
    
    data_sub <- data.frame(
      PRES = as.numeric(pres), PRES_QC = pres_qc,
      TEMP = as.numeric(temp), TEMP_QC = temp_qc,
      PSAL = as.numeric(sal), PSAL_QC = sal_qc,
      NITRATE = as.numeric(nitrate), NITRATE_QC = nitrate_qc,
      NITRATE_ADJUSTED = as.numeric(nitrate_adj), NITRATE_ADJUSTED_QC = nitrate_adj_qc,
      PH_IN_SITU_TOTAL = as.numeric(ph), PH_IN_SITU_TOTAL_QC = ph_qc,
      PH_IN_SITU_TOTAL_ADJUSTED = as.numeric(ph_adj), PH_IN_SITU_TOTAL_ADJUSTED_QC = ph_adj_qc,
      stringsAsFactors = FALSE
    ) %>%
      convert_fill_to_nan() %>%
      filter(!is.na(PRES)) %>%
      arrange(PRES)
    
    if (nrow(data_sub) > 0) {
      data_sub$TIME <- time_val
      data_sub$LONGITUDE <- lon_val
      data_sub$LATITUDE  <- lat_val
      data_sub$CYCLE_NUMBER <- as.numeric(str_remove(substr(www, 9, 11), "^0+"))
      data_sub$float_num <- substr(www, 1, 7)
      
      float_data <- rbind(float_data, data_sub)
    }
    
    nc_close(profile_bio)
    nc_close(profile_physic)
  }
  
  if (is.null(float_data) || nrow(float_data) == 0) {
    warning(paste("Skipping WMO", active_wmo, "- No valid profiles collected."))
    next
  }
  
  # --- STEP 2: SAGE V2 CALIBRATION & DMQC EVALUATION ---
  df_all <- float_data %>%
    mutate(
      float_ids     = as.numeric(float_num),
      CYCLE_NUMBER  = as.numeric(CYCLE_NUMBER),
      DATE_TIME     = TIME,
      PRES_ADJUSTED = PRES
    ) %>%
    left_join(calib_df, by = c("float_ids", "CYCLE_NUMBER")) %>%
    convert_fill_to_nan()
  
  delta_days <- as.numeric(difftime(df_all$DATE_TIME, df_all$launch_date, units = "days"))
  
  df_all <- df_all %>%
    mutate(
      nitrate_gain_eff = coalesce(nitrate_gain, 1) * (1 + (coalesce(nitrate_drift, 0) / 100) * (delta_days / 365.25)),
      ph_gain_eff      = coalesce(ph_gain, 1) * (1 + (coalesce(ph_drift, 0) / 100) * (delta_days / 365.25)),
      
      # Argo formula: ADJUSTED = (RAW * GAIN) + OFFSET
      NITRATE_CALCULATED_ADJ = ifelse(!is.na(NITRATE), (NITRATE * nitrate_gain_eff) + coalesce(nitrate_offset, 0), NaN),
      PH_CALCULATED_ADJ      = ifelse(!is.na(PH_IN_SITU_TOTAL), (PH_IN_SITU_TOTAL * ph_gain_eff) + coalesce(ph_offset, 0), NaN),
      
      # Distinct A-mode parameter assignments:
      # If raw NetCDF lacked an adjusted profile, populate it with calculated adjusted data so A-mode is complete across cycles
      NITRATE_ADJUSTED = coalesce(NITRATE_ADJUSTED, NITRATE_CALCULATED_ADJ),
      NITRATE_ADJUSTED_QC = coalesce(NITRATE_ADJUSTED_QC, NITRATE_QC, 1),
      
      PH_IN_SITU_TOTAL_ADJUSTED = coalesce(PH_IN_SITU_TOTAL_ADJUSTED, PH_CALCULATED_ADJ),
      PH_IN_SITU_TOTAL_ADJUSTED_QC = coalesce(PH_IN_SITU_TOTAL_ADJUSTED_QC, PH_IN_SITU_TOTAL_QC, 2),
      
      # Distinct D-mode parameter assignments (Final scientific DMQC data)
      NITRATE_FINAL_RAW = coalesce(NITRATE_CALCULATED_ADJ, NITRATE_ADJUSTED, NITRATE),
      PH_FINAL_RAW      = coalesce(PH_CALCULATED_ADJ, PH_IN_SITU_TOTAL_ADJUSTED, PH_IN_SITU_TOTAL),
      
      pres_missing_or_bad = is.na(PRES) | PRES_QC %in% c("9", 9),
      
      NITRATE_FINAL = ifelse(pres_missing_or_bad, NaN, NITRATE_FINAL_RAW),
      PH_FINAL      = ifelse(pres_missing_or_bad, NaN, PH_FINAL_RAW)
    )
  
  # Diagnostic Print
  adj_avail_summary <- df_all %>%
    summarise(
      Total_Rows = n(),
      Raw_Nitrate_Count = sum(!is.na(NITRATE)),
      Adjusted_Nitrate_A_Mode_Count = sum(!is.na(NITRATE_ADJUSTED)),
      Final_Nitrate_D_Mode_Count = sum(!is.na(NITRATE_FINAL)),
      Raw_pH_Count = sum(!is.na(PH_IN_SITU_TOTAL)),
      Adjusted_pH_A_Mode_Count = sum(!is.na(PH_IN_SITU_TOTAL_ADJUSTED)),
      Final_pH_D_Mode_Count = sum(!is.na(PH_FINAL))
    )
  cat("\n==============================================================================\n")
  cat(sprintf("  ADJUSTED PARAMETER AVAILABILITY SUMMARY FOR FLOAT: %s\n", active_wmo))
  cat("==============================================================================\n")
  print(knitr::kable(adj_avail_summary, format = "simple"))
  cat("==============================================================================\n\n")
  
  # --- STEP 2.5: BGC-ARGO COOKBOOK COMPLIANT QC & ERROR EVALUATION ---
  df_all <- df_all %>%
    group_by(CYCLE_NUMBER) %>%
    group_modify(~ {
      df_cyc <- .x
      n_obs  <- nrow(df_cyc)
      
      # 1. Nitrate QC
      nit_v <- df_cyc$NITRATE_FINAL
      is_nit_adjusted <- any(!is.na(df_cyc$NITRATE_CALCULATED_ADJ) | !is.na(df_cyc$NITRATE_ADJUSTED))
      qc_nitrate <- rep(ifelse(is_nit_adjusted, 1, 3), n_obs)
      
      qc_nitrate[!is.na(nit_v) & (nit_v < -0.5 | nit_v > 60.0)] <- 4
      
      valid_nit_idx <- which(!is.na(nit_v))
      if (length(valid_nit_idx) >= 2) {
        max_idx <- valid_nit_idx[which.max(df_cyc$PRES[valid_nit_idx])]
        sec_idx <- valid_nit_idx[order(df_cyc$PRES[valid_nit_idx], decreasing = TRUE)[2]]
        if (abs(nit_v[max_idx] - nit_v[sec_idx]) > 2.0) {
          qc_nitrate[max_idx] <- 3
        }
      }
      
      if ("TEMP_QC" %in% names(df_cyc)) qc_nitrate[df_cyc$TEMP_QC %in% c("4", 4)] <- 4
      if ("PSAL_QC" %in% names(df_cyc)) qc_nitrate[df_cyc$PSAL_QC %in% c("4", 4)] <- 4
      qc_nitrate[df_cyc$PRES_QC %in% c("4", 4)] <- 4
      qc_nitrate[is.na(nit_v)] <- 9
      qc_nitrate[df_cyc$PRES_QC %in% c("9", 9)] <- 9
      df_cyc$NITRATE_FINAL_QC <- qc_nitrate
      
      # 2. pH QC
      ph_v <- df_cyc$PH_FINAL
      is_ph_adjusted <- any(!is.na(df_cyc$PH_CALCULATED_ADJ) | !is.na(df_cyc$PH_IN_SITU_TOTAL_ADJUSTED))
      qc_ph <- rep(ifelse(is_ph_adjusted, 2, 3), n_obs)
      
      qc_ph[!is.na(ph_v) & (ph_v < 7.0 | ph_v > 8.8)] <- 4
      
      valid_ph_idx <- which(!is.na(ph_v))
      if (length(valid_ph_idx) >= 2) {
        max_idx <- valid_ph_idx[which.max(df_cyc$PRES[valid_ph_idx])]
        sec_idx <- valid_ph_idx[order(df_cyc$PRES[valid_ph_idx], decreasing = TRUE)[2]]
        if (abs(ph_v[max_idx] - ph_v[sec_idx]) > 0.05) {
          qc_ph[max_idx] <- 3
        }
      }
      
      if ("TEMP_QC" %in% names(df_cyc)) qc_ph[df_cyc$TEMP_QC %in% c("4", 4)] <- 4
      qc_ph[df_cyc$PRES_QC %in% c("4", 4)] <- 4
      if ("PSAL_QC" %in% names(df_cyc)) {
        psal_bad_idx <- df_cyc$PSAL_QC %in% c("4", 4)
        qc_ph[psal_bad_idx & qc_ph != 4] <- 3
      }
      qc_ph[is.na(ph_v)] <- 9
      qc_ph[df_cyc$PRES_QC %in% c("9", 9)] <- 9
      df_cyc$PH_FINAL_QC <- qc_ph
      
      return(df_cyc)
    }) %>%
    ungroup() %>%
    mutate(
      NITRATE_ADJUSTED_ERROR          = ifelse(NITRATE_FINAL_QC %in% c(1, 2, 3), 1.00, NaN),
      PH_IN_SITU_TOTAL_ADJUSTED_ERROR = ifelse(PH_FINAL_QC %in% c(1, 2, 3), 0.007, NaN)
    ) %>%
    select(-c(NITRATE_FINAL_RAW, PH_FINAL_RAW, pres_missing_or_bad))
  
  # Bottle samples
  df_nit_bottle_master <- Bottle_Vali %>%
    filter(Float_num == float_id & !is.na(NITRATE_Vali) & !is.na(PRES)) %>%
    select(PRES, Nitrate = NITRATE_Vali) %>%
    mutate(Mode = factor("Bottle", levels = mode_levels))
  
  df_ph_bottle_master <- Bottle_Vali %>%
    filter(Float_num == float_id & !is.na(PH_Vali) & !is.na(PRES)) %>%
    select(PRES, pH = PH_Vali) %>%
    mutate(Mode = factor("Bottle", levels = mode_levels))
  
  # ==============================================================================
  # STEP 3: RENDER PER-CYCLE PLOTS (R-mode, A-mode, D-mode & Bottle)
  # - A-mode explicitly plots NITRATE_ADJUSTED and PH_IN_SITU_TOTAL_ADJUSTED
  # - D-mode explicitly plots NITRATE_FINAL and PH_FINAL
  # ==============================================================================
  all_cycles <- sort(unique(df_all$CYCLE_NUMBER[!is.na(df_all$CYCLE_NUMBER)]))
  message(sprintf("Rendering profile plots for WMO %s across %d cycles...", active_wmo, length(all_cycles)))
  
  for (cyc in all_cycles) {
    df_sub <- df_all %>% filter(CYCLE_NUMBER == cyc)
    if (nrow(df_sub) == 0) next
    
    cyc_lat      <- df_sub$LATITUDE[1]
    cyc_lon      <- df_sub$LONGITUDE[1]
    cyc_date_str <- format(df_sub$DATE_TIME[1], "%Y-%m-%d")
    
    # --- 1. NITRATE PROFILE ---
    df_nit_cycle_clean <- df_sub %>%
      select(PRES, NITRATE, NITRATE_ADJUSTED, NITRATE_FINAL) %>%
      pivot_longer(
        cols      = c("NITRATE", "NITRATE_ADJUSTED", "NITRATE_FINAL"),
        names_to  = "Mode_Raw",
        values_to = "Nitrate"
      ) %>%
      mutate(
        Mode_Name = case_when(
          Mode_Raw == "NITRATE"          ~ "R-mode",
          Mode_Raw == "NITRATE_ADJUSTED" ~ "A-mode",
          Mode_Raw == "NITRATE_FINAL"    ~ "D-mode"
        ),
        Mode = factor(Mode_Name, levels = mode_levels)
      ) %>%
      filter(!is.na(Nitrate) & !is.na(PRES))
    
    has_nit_bottle <- (cyc == 1) && (nrow(df_nit_bottle_master) > 0)
    df_nit_combined <- if (has_nit_bottle) bind_rows(df_nit_cycle_clean, df_nit_bottle_master) else df_nit_cycle_clean
    
    # Layer order: D-mode plotted first, A-mode plotted on top, R-mode, Bottle last
    # Ensures A-mode dashed line is fully visible even if D-mode has identical values
    df_nit_combined <- df_nit_combined %>% 
      arrange(factor(Mode, levels = c("D-mode", "A-mode", "R-mode", "Bottle")), PRES)
    
    if (nrow(df_nit_combined) > 0) {
      p_nit_cycle <- ggplot() +
        geom_path(
          data = filter(df_nit_combined, Mode != "Bottle"),
          aes(x = Nitrate, y = PRES, color = Mode, linetype = Mode, linewidth = Mode),
          alpha = 0.95
        ) +
        geom_point(
          data = filter(df_nit_combined, Mode == "Bottle"),
          aes(x = Nitrate, y = PRES, color = Mode, shape = Mode),
          size = 2.5, alpha = 0.95
        ) +
        scale_color_manual(values = color_mapping, limits = mode_levels, drop = FALSE) +
        scale_linetype_manual(values = linetype_mapping, limits = mode_levels, drop = FALSE) +
        scale_linewidth_manual(values = linewidth_map, limits = mode_levels, drop = FALSE) +
        scale_shape_manual(values = shape_mapping, limits = mode_levels, drop = FALSE) +
        scale_y_reverse(limits = c(max(c(df_nit_combined$PRES, 2000), na.rm = TRUE), 0)) +
        scale_x_continuous(position = "top") +
        labs(
          title    = paste0("Float ", float_id, " Cycle ", cyc, " Nitrate Profiles (R-mode, A-mode, D-mode)"),
          subtitle = paste0("Lat: ", round(cyc_lat, 2), "°, Lon: ", round(cyc_lon, 2), "° | Cycle Date: ", cyc_date_str),
          caption  = paste0("Plot generated on: ", today_str),
          x        = expression(Nitrate ~ (mu * "mol/kg")),
          y        = "Pressure (dbar)",
          color    = "Data Source"
        ) +
        theme_bw(base_size = 13) +
        theme(
          aspect.ratio      = 1,
          legend.position   = "right",
          legend.title      = element_text(face = "bold", size = 11),
          legend.text       = element_text(size = 10),
          legend.key.width  = unit(1.4, "cm"),
          legend.key.height = unit(0.5, "cm")
        ) +
        guides(
          color = guide_legend(
            override.aes = list(
              shape     = shape_mapping[mode_levels],
              linetype  = linetype_mapping[mode_levels],
              linewidth = linewidth_map[mode_levels],
              size      = ifelse(mode_levels == "Bottle", 2.5, NA)
            )
          ),
          linetype  = "none",
          shape     = "none",
          linewidth = "none"
        )
      
      nit_file_path <- file.path(dir_nitrate_plots, sprintf("Float_%d_Cycle_%03d_Nitrate_Profile.png", float_id, cyc))
      ggsave(nit_file_path, plot = p_nit_cycle, width = 9.5, height = 7, dpi = 300)
    }
    
    # --- 2. pH PROFILE ---
    df_ph_cycle_clean <- df_sub %>%
      select(PRES, PH_IN_SITU_TOTAL, PH_IN_SITU_TOTAL_ADJUSTED, PH_FINAL) %>%
      pivot_longer(
        cols      = c("PH_IN_SITU_TOTAL", "PH_IN_SITU_TOTAL_ADJUSTED", "PH_FINAL"),
        names_to  = "Mode_Raw",
        values_to = "pH"
      ) %>%
      mutate(
        Mode_Name = case_when(
          Mode_Raw == "PH_IN_SITU_TOTAL"          ~ "R-mode",
          Mode_Raw == "PH_IN_SITU_TOTAL_ADJUSTED" ~ "A-mode",
          Mode_Raw == "PH_FINAL"                  ~ "D-mode"
        ),
        Mode = factor(Mode_Name, levels = mode_levels)
      ) %>%
      filter(!is.na(pH) & !is.na(PRES))
    
    has_ph_bottle <- (cyc == 1) && (nrow(df_ph_bottle_master) > 0)
    df_ph_combined <- if (has_ph_bottle) bind_rows(df_ph_cycle_clean, df_ph_bottle_master) else df_ph_cycle_clean
    
    # Layer order: D-mode plotted first, A-mode plotted on top, R-mode, Bottle last
    df_ph_combined <- df_ph_combined %>% 
      arrange(factor(Mode, levels = c("D-mode", "A-mode", "R-mode", "Bottle")), PRES)
    
    if (nrow(df_ph_combined) > 0) {
      p_ph_cycle <- ggplot() +
        geom_path(
          data = filter(df_ph_combined, Mode != "Bottle"),
          aes(x = pH, y = PRES, color = Mode, linetype = Mode, linewidth = Mode),
          alpha = 0.95
        ) +
        geom_point(
          data = filter(df_ph_combined, Mode == "Bottle"),
          aes(x = pH, y = PRES, color = Mode, shape = Mode),
          size = 2.5, alpha = 0.95
        ) +
        scale_color_manual(values = color_mapping, limits = mode_levels, drop = FALSE) +
        scale_linetype_manual(values = linetype_mapping, limits = mode_levels, drop = FALSE) +
        scale_linewidth_manual(values = linewidth_map, limits = mode_levels, drop = FALSE) +
        scale_shape_manual(values = shape_mapping, limits = mode_levels, drop = FALSE) +
        scale_y_reverse(limits = c(max(c(df_ph_combined$PRES, 2000), na.rm = TRUE), 0)) +
        scale_x_continuous(position = "top") +
        labs(
          title    = paste0("Float ", float_id, " Cycle ", cyc, " pH Profiles (R-mode, A-mode, D-mode)"),
          subtitle = paste0("Lat: ", round(cyc_lat, 2), "°, Lon: ", round(cyc_lon, 2), "° | Cycle Date: ", cyc_date_str),
          caption  = paste0("Plot generated on: ", today_str),
          x        = "pH (In-Situ Total Scale)",
          y        = "Pressure (dbar)",
          color    = "Data Source"
        ) +
        theme_bw(base_size = 13) +
        theme(
          aspect.ratio      = 1,
          legend.position   = "right",
          legend.title      = element_text(face = "bold", size = 11),
          legend.text       = element_text(size = 10),
          legend.key.width  = unit(1.4, "cm"),
          legend.key.height = unit(0.5, "cm")
        ) +
        guides(
          color = guide_legend(
            override.aes = list(
              shape     = shape_mapping[mode_levels],
              linetype  = linetype_mapping[mode_levels],
              linewidth = linewidth_map[mode_levels],
              size      = ifelse(mode_levels == "Bottle", 2.5, NA)
            )
          ),
          linetype  = "none",
          shape     = "none",
          linewidth = "none"
        )
      
      ph_file_path <- file.path(dir_ph_plots, sprintf("Float_%d_Cycle_%03d_pH_Profile.png", float_id, cyc))
      ggsave(ph_file_path, plot = p_ph_cycle, width = 9.5, height = 7, dpi = 300)
    }
  }
  
  # ==============================================================================
  # STEP 4: SECTION CONTOUR & QC SUMMARY PLOTS (ALL CYCLES)
  # - A-mode plots NITRATE_ADJUSTED and PH_IN_SITU_TOTAL_ADJUSTED
  # - D-mode plots NITRATE_FINAL and PH_FINAL
  # ==============================================================================
  df_summary_base <- df_all %>% 
    filter(!is.na(PRES)) %>%
    mutate(
      CYCLE_NUMBER = as.numeric(CYCLE_NUMBER),
      PRES         = as.numeric(PRES)
    )
  
  create_interpolated_grid <- function(df, val_col, grid_nx = 100, grid_ny = 100) {
    clean_df <- df %>% 
      filter(!is.na(!!sym(val_col))) %>% 
      select(x = CYCLE_NUMBER, y = PRES, z = !!sym(val_col)) %>%
      group_by(x, y) %>%
      summarise(z = mean(z, na.rm = TRUE), .groups = "drop")
    
    if (nrow(clean_df) < 5 || length(unique(clean_df$x)) < 2 || length(unique(clean_df$y)) < 2) return(NULL)
    if (diff(range(clean_df$z, na.rm = TRUE)) == 0) return(NULL)
    
    interp_res <- suppressWarnings(tryCatch({
      interp::interp(
        x = clean_df$x, y = clean_df$y, z = clean_df$z,
        xo = seq(min(clean_df$x), max(clean_df$x), length.out = grid_nx),
        yo = seq(min(clean_df$y), max(clean_df$y), length.out = grid_ny),
        linear = TRUE, extrap = FALSE, duplicate = "mean"
      )
    }, error = function(e) return(NULL)))
    
    if (is.null(interp_res) || is.null(interp_res$z) || !is.matrix(interp_res$z)) return(NULL)
    if (nrow(interp_res$z) != grid_nx || ncol(interp_res$z) != grid_ny) return(NULL)
    
    z_vector <- as.vector(interp_res$z)
    grid_df  <- expand.grid(x = interp_res$x, y = interp_res$y)
    if (length(z_vector) != nrow(grid_df) || all(is.na(z_vector))) return(NULL)
    
    grid_df$z <- z_vector
    return(grid_df %>% filter(!is.na(z)))
  }
  
  generate_section_qc_plot <- function(df_base, var_name, qc_var_name, title_prefix, file_suffix, 
                                       col_palette = "viridis", unit_expr = "") {
    if (!var_name %in% names(df_base) || !any(!is.na(df_base[[var_name]]))) return(NULL)
    grid_df <- create_interpolated_grid(df_base, var_name)
    
    df_qc <- df_base %>% filter(!is.na(!!sym(var_name)))
    if (qc_var_name %in% names(df_qc)) {
      df_qc$QC_VAL <- as.character(coalesce(as.character(df_qc[[qc_var_name]]), "1"))
    } else {
      df_qc$QC_VAL <- "1"
    }
    
    if (!is.null(grid_df)) {
      p_data <- ggplot() +
        geom_contour_filled(data = grid_df, aes(x = x, y = y, z = z), bins = 15) +
        scale_y_reverse(limits = c(max(df_base$PRES, na.rm = TRUE), 0)) +
        scale_fill_viridis_d(option = col_palette, name = unit_expr) +
        labs(
          title = paste0("Float ", float_id, " - Interpolated ", title_prefix, " Section Profile"),
          x     = "Cycle Number",
          y     = "Pressure (dbar)"
        ) +
        theme_bw(base_size = 12) + theme(legend.position = "right")
    } else {
      p_data <- ggplot(df_qc, aes(x = CYCLE_NUMBER, y = PRES)) +
        geom_tile(aes(fill = !!sym(var_name)), width = 1.2, height = 10) +
        scale_y_reverse(limits = c(max(df_base$PRES, na.rm = TRUE), 0)) +
        scale_fill_viridis_c(option = col_palette, name = unit_expr) +
        labs(
          title = paste0("Float ", float_id, " - ", title_prefix, " Section Profile"),
          x     = "Cycle Number",
          y     = "Pressure (dbar)"
        ) +
        theme_bw(base_size = 12) + theme(legend.position = "right")
    }
    
    p_qc <- ggplot(df_qc, aes(x = CYCLE_NUMBER, y = PRES, color = QC_VAL)) +
      geom_point(alpha = 0.7, size = 1.2) +
      scale_y_reverse(limits = c(max(df_base$PRES, na.rm = TRUE), 0)) +
      scale_color_manual(values = qc_color_scale, drop = FALSE, name = "QC Flag") +
      labs(title = paste0(title_prefix, " QC Flags Section Map"), x = "Cycle Number", y = "Pressure (dbar)") +
      theme_bw(base_size = 12) + theme(legend.position = "right")
    
    p_combined <- (p_data / p_qc) +
      plot_annotation(
        title    = paste0("BGC-Argo Float ", float_id, " - ", title_prefix, " Section DMQC Summary"),
        subtitle = paste0("Generated on: ", today_str),
        theme    = theme(plot.title = element_text(face = "bold", size = 14))
      )
    
    ggsave(
      file.path(dir_qc_summary, sprintf("Float_%d_%s_Section_Summary.png", float_id, file_suffix)), 
      plot = p_combined, width = 12, height = 8, dpi = 150
    )
  }
  
  # Produce Section Plots for Nitrate & pH (R-mode, A-mode, D-mode)
  generate_section_qc_plot(df_summary_base, "NITRATE", "NITRATE_QC", "NITRATE RAW (R-mode)", "NITRATE_RAW", col_palette = "viridis", unit_expr = expression(Nitrate ~ (mu * "mol/kg")))
  generate_section_qc_plot(df_summary_base, "NITRATE_ADJUSTED", "NITRATE_ADJUSTED_QC", "NITRATE_ADJUSTED (A-mode)", "NITRATE_ADJUSTED", col_palette = "viridis", unit_expr = expression(Nitrate ~ (mu * "mol/kg")))
  generate_section_qc_plot(df_summary_base, "NITRATE_FINAL", "NITRATE_FINAL_QC", "NITRATE_FINAL (D-mode)", "NITRATE_FINAL", col_palette = "viridis", unit_expr = expression(Nitrate ~ (mu * "mol/kg")))
  
  generate_section_qc_plot(df_summary_base, "PH_IN_SITU_TOTAL", "PH_IN_SITU_TOTAL_QC", "PH RAW (R-mode)", "PH_RAW", col_palette = "magma", unit_expr = "pH (In-Situ)")
  generate_section_qc_plot(df_summary_base, "PH_IN_SITU_TOTAL_ADJUSTED", "PH_IN_SITU_TOTAL_ADJUSTED_QC", "PH_ADJUSTED (A-mode)", "PH_ADJUSTED", col_palette = "magma", unit_expr = "pH (In-Situ)")
  generate_section_qc_plot(df_summary_base, "PH_FINAL", "PH_FINAL_QC", "PH_FINAL (D-mode)", "PH_FINAL", col_palette = "magma", unit_expr = "pH (In-Situ)")
  
  # --- STEP 5: EXPORT COMBINED CSV ---
  combo_output_path <- file.path(combo_output_dir, sprintf("NITRATE_PH_%d.csv", float_id))
  
  df_export_prep <- df_all %>%
    mutate(
      float_num = coalesce(as.character(float_ids), as.character(float_num)),
      across(where(is.logical), as.numeric)
    )
  
  df_combo_export <- df_export_prep %>%
    select(-any_of(c(
      "WMOID", "ID", "NITRATE_QC", "PH_IN_SITU_TOTAL_QC", 
      "NITRATE_CALCULATED_ADJ", "PH_CALCULATED_ADJ", 
      "float_ids", "TIME", "launch_date", "LAUNCH_DATE", 
      "launch_date_str"
    ))) %>%
    rename_with(toupper) %>%
    relocate(
      any_of(c(
        "FLOAT_NUM", "CYCLE_NUMBER", "PRES", "TEMP", "PSAL", 
        "NITRATE", "NITRATE_ADJUSTED", "NITRATE_FINAL", "NITRATE_FINAL_QC", "NITRATE_ADJUSTED_ERROR", 
        "PH_IN_SITU_TOTAL", "PH_IN_SITU_TOTAL_ADJUSTED", "PH_FINAL", "PH_FINAL_QC", "PH_IN_SITU_TOTAL_ADJUSTED_ERROR"
      ))
    ) %>%
    convert_fill_to_nan() %>%
    mutate(
      across(where(is.numeric), ~ ifelse(is.na(.), NaN, .))
    )
  
  write.csv(df_combo_export, combo_output_path, row.names = FALSE, na = "NaN")
  message(sprintf("Successfully exported combined CSV -> %s", combo_output_path))
  message(sprintf("=== COMPLETED PROCESSING FOR FLOAT: %d ===\n", float_id))
}