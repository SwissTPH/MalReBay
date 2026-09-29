# Generates the data-entry templates in inst/templates/ (one per data type).
# Run from the package root:  Rscript data-raw/make_templates.R
# Requires openxlsx (not a package dependency -- only needed to rebuild these).
#
# The layout matches what import_data() reads: sheet 1 = late failures,
# sheet 2 = additional Day 0 samples, both as PatientID | Day | Site | alleles.
# The Instructions sheet must stay after them (import_data() reads by position).
# After changing this script, check each template still imports, e.g.
#   import_data("inst/templates/MalReBay_template_msp_glurp.xlsx",
#               system.file("extdata", "makers_details.xlsx", package = "MalReBay"))

library(openxlsx)

out_dir <- "inst/templates"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

font <- "Arial"
hdr_meta   <- createStyle(fontName = font, fontSize = 10, textDecoration = "bold",
                          fgFill = "#1F4E79", fontColour = "#FFFFFF", halign = "center",
                          border = "Bottom")
hdr_allele <- createStyle(fontName = font, fontSize = 10, textDecoration = "bold",
                          fgFill = "#DDEBF7", fontColour = "#1F1F1F", halign = "center",
                          border = "Bottom")
example_st <- createStyle(fontName = font, fontSize = 10, fontColour = "#7F7F7F",
                          textDecoration = "italic", halign = "center")
body_st    <- createStyle(fontName = font, fontSize = 10)
title_st   <- createStyle(fontName = font, fontSize = 14, textDecoration = "bold",
                          fontColour = "#1F4E79")
h2_st      <- createStyle(fontName = font, fontSize = 11, textDecoration = "bold",
                          fontColour = "#1F4E79")
wrap_st    <- createStyle(fontName = font, fontSize = 10, wrapText = TRUE, valign = "top")
bold_st    <- createStyle(fontName = font, fontSize = 10, textDecoration = "bold",
                          wrapText = TRUE, valign = "top")
tbl_hdr_st <- createStyle(fontName = font, fontSize = 10, textDecoration = "bold",
                          fgFill = "#DDEBF7", border = "Bottom", wrapText = TRUE)

# markers: named vector marker -> number of allele columns
allele_cols <- function(markers, sep) {
  unlist(lapply(names(markers), function(m) paste0(m, sep, seq_len(markers[[m]]))))
}

# rows: list of lists(PatientID, Day, Site, alleles = named list col -> value)
make_sheet_df <- function(cols, rows) {
  df <- as.data.frame(matrix(NA, nrow = length(rows), ncol = length(cols) + 3))
  names(df) <- c("PatientID", "Day", "Site", cols)
  for (i in seq_along(rows)) {
    r <- rows[[i]]
    df$PatientID[i] <- r$id; df$Day[i] <- r$day; df$Site[i] <- r$site
    for (col in names(r$alleles)) df[[col]][i] <- r$alleles[[col]]
  }
  df
}

add_data_sheet <- function(wb, sheet, df, notes, day0_only = FALSE) {
  addWorksheet(wb, sheet, gridLines = TRUE)
  writeData(wb, sheet, df, keepNA = FALSE)
  n_col <- ncol(df)
  addStyle(wb, sheet, hdr_meta, rows = 1, cols = 1:3, gridExpand = TRUE)
  addStyle(wb, sheet, hdr_allele, rows = 1, cols = 4:n_col, gridExpand = TRUE)
  if (nrow(df) > 0)
    addStyle(wb, sheet, example_st, rows = 1 + seq_len(nrow(df)), cols = 1:n_col, gridExpand = TRUE)
  addStyle(wb, sheet, body_st, rows = (nrow(df) + 2):500, cols = 1:n_col, gridExpand = TRUE)
  setColWidths(wb, sheet, cols = 1, widths = 18)
  setColWidths(wb, sheet, cols = 2:3, widths = 10)
  setColWidths(wb, sheet, cols = 4:n_col, widths = "auto")
  freezePane(wb, sheet, firstActiveRow = 2, firstActiveCol = 4)

  # Day: whole number >= 0 (and exactly 0 on the Additional sheet)
  if (day0_only) {
    dataValidation(wb, sheet, cols = 2, rows = 2:500, type = "whole", operator = "equal",
                   value = 0, showInputMsg = TRUE, showErrorMsg = TRUE)
  } else {
    dataValidation(wb, sheet, cols = 2, rows = 2:500, type = "whole",
                   operator = "greaterThanOrEqual", value = 0,
                   showInputMsg = TRUE, showErrorMsg = TRUE)
  }

  for (col in names(notes)) {
    writeComment(wb, sheet, col = which(names(df) == col), row = 1,
                 comment = createComment(notes[[col]], author = "MalReBay",
                                         width = 3, height = 5, visible = FALSE))
  }
}

add_instructions <- function(wb, type_label, allele_desc, sep, marker_rows, extra_notes) {
  s <- "Instructions"
  addWorksheet(wb, s, gridLines = FALSE)
  setColWidths(wb, s, cols = 1:3, widths = c(26, 60, 40))
  r <- 1
  put <- function(x, style, col = 1) { writeData(wb, s, x, startCol = col, startRow = r); addStyle(wb, s, style, rows = r, cols = col) }

  put(paste0("MalReBay data template: ", type_label), title_st); r <- r + 2

  put("What to edit", h2_st); r <- r + 1
  steps <- c(
    "1. On 'Late Treatment Failures', delete the grey example rows and enter one row per sample: a Day 0 row and a recurrence row for every patient with a recurrent infection.",
    "2. On 'Additional' (optional), delete the grey example row and enter Day 0 samples from other patients in the same study (e.g. patients without recurrence). They are used only to estimate background allele frequencies. Leave the sheet empty (headers only) if you have none.",
    "3. Rename / add / remove the marker columns to match your panel (see 'Allele columns' below). Both data sheets must have the same marker columns.",
    "4. Keep the sheet order: MalReBay reads the FIRST sheet as late failures and the SECOND as additional data, regardless of their names. Do not insert sheets before them.",
    "5. Run: import_data(filepath = \"your_file.xlsx\", marker_filepath = \"your_marker_details.xlsx\")"
  )
  for (st in steps) { put(st, wrap_st); mergeCells(wb, s, cols = 1:3, rows = r); setRowHeights(wb, s, r, 30); r <- r + 1 }
  r <- r + 1

  put("Columns", h2_st); r <- r + 1
  tbl <- data.frame(
    Column = c("PatientID", "Day", "Site", allele_desc$column),
    Description = c(
      "Patient identifier. Must be identical on the patient's Day 0 and recurrence rows, and unique within a site.",
      "Day of sampling. 0 = baseline (before treatment). Any other number = the recurrence day (e.g. 28, 42). Each patient needs exactly one Day 0 row and one recurrence row. On 'Additional', Day is always 0.",
      "Study site name or code. Patients are analysed per site.",
      allele_desc$description),
    Example = c("PT-001", "0 / 42", "SiteA", allele_desc$example),
    stringsAsFactors = FALSE)
  writeData(wb, s, tbl, startRow = r, headerStyle = tbl_hdr_st)
  addStyle(wb, s, bold_st, rows = r + seq_len(nrow(tbl)), cols = 1, gridExpand = TRUE)
  addStyle(wb, s, wrap_st, rows = r + seq_len(nrow(tbl)), cols = 2:3, gridExpand = TRUE)
  setRowHeights(wb, s, r + seq_len(nrow(tbl)), 45)
  r <- r + nrow(tbl) + 2

  put("Allele columns", h2_st); r <- r + 1
  notes <- c(
    paste0("Name each allele column <marker>", sep, "<n>, e.g. ", marker_rows[1], sep, "1, ",
           marker_rows[1], sep, "2, ... One column per allele observed in a sample; add as many as the highest number of alleles (multiplicity of infection) seen for that marker."),
    "Marker names must match the marker_id column of the marker details file exactly (case-sensitive). Markers not found there are skipped with a message.",
    "Fill alleles from the left (_1 first) and leave the remaining cells empty.",
    "Missing / failed values: leave the cell empty. NA, -, N/A, Failed and 0 are also treated as missing.",
    extra_notes
  )
  for (st in notes) { put(paste0("• ", st), wrap_st); mergeCells(wb, s, cols = 1:3, rows = r); setRowHeights(wb, s, r, 30); r <- r + 1 }
  r <- r + 1

  put("Marker details file", h2_st); r <- r + 1
  md <- c(
    "A separate xlsx with one row per marker and columns: marker_id, markertype, repeatlength, binning_method.",
    "MalReBay ships one with common markers: system.file(\"extdata\", \"makers_details.xlsx\", package = \"MalReBay\"). Copy it and add rows for any markers of yours that are missing."
  )
  for (st in md) { put(paste0("• ", st), wrap_st); mergeCells(wb, s, cols = 1:3, rows = r); setRowHeights(wb, s, r, 30); r <- r + 1 }
  r <- r + 1

  put("Legend", h2_st); r <- r + 1
  legend <- list(list("Dark blue header", "Required sample information columns", hdr_meta),
                 list("Light blue header", "Allele columns: rename/add/remove to match your panel", hdr_allele),
                 list("Grey italic rows", "Examples only: delete before running MalReBay", example_st))
  for (lg in legend) {
    writeData(wb, s, lg[[1]], startCol = 1, startRow = r); addStyle(wb, s, lg[[3]], rows = r, cols = 1)
    writeData(wb, s, lg[[2]], startCol = 2, startRow = r); addStyle(wb, s, wrap_st, rows = r, cols = 2)
    r <- r + 1
  }
  writeData(wb, s, "Hover over a column header on the data sheets for a short note.", startRow = r + 1)
  addStyle(wb, s, wrap_st, rows = r + 1, cols = 1)
}

meta_notes <- list(
  PatientID = "Same ID on the patient's Day 0 and recurrence rows.",
  Day       = "0 = baseline (before treatment); otherwise the day of recurrence, e.g. 28 or 42.",
  Site      = "Study site name or code."
)

build <- function(file, type_label, markers, sep, lf_rows, add_rows, allele_desc, extra_notes) {
  cols <- allele_cols(markers, sep)
  wb <- createWorkbook(creator = "MalReBay")
  modifyBaseFont(wb, fontName = font, fontSize = 10)
  first_allele_note <- setNames(list(paste0("Alleles for marker ", names(markers)[1],
    ": one column per allele (", names(markers)[1], sep, "1, ", names(markers)[1], sep,
    "2, ...). Rename to your markers.")), cols[1])
  add_data_sheet(wb, "Late Treatment Failures", make_sheet_df(cols, lf_rows),
                 c(meta_notes, first_allele_note))
  add_data_sheet(wb, "Additional", make_sheet_df(cols, add_rows),
                 c(meta_notes[c("PatientID", "Site")],
                   list(Day = "Always 0: only baseline samples go on this sheet."),
                   first_allele_note), day0_only = TRUE)
  add_instructions(wb, type_label, allele_desc, sep, names(markers), extra_notes)
  saveWorkbook(wb, file.path(out_dir, file), overwrite = TRUE)
}

al <- function(...) list(...)

# ---- Microsatellites (7 neutral microsatellites) ----------------------------
build("MalReBay_template_microsatellite.xlsx", "microsatellites",
  markers = c("313" = 3, "383" = 3, TA1 = 5, POLYA = 4, PFPK2 = 5, "2490" = 3, TA109 = 5),
  sep = "_",
  lf_rows = list(
    list(id = "EXAMPLE-001", day = 0,  site = "SiteA",
         alleles = al(`313_1` = 248, `383_1` = 133, TA1_1 = 174, POLYA_1 = 153, PFPK2_1 = 159, `2490_1` = 78, TA109_1 = 172)),
    list(id = "EXAMPLE-001", day = 42, site = "SiteA",
         alleles = al(`313_1` = 246, `383_1` = 123, TA1_1 = 177, POLYA_1 = 105, PFPK2_1 = 183, `2490_1` = 81, TA109_1 = 184)),
    list(id = "EXAMPLE-002", day = 0,  site = "SiteA",
         alleles = al(`313_1` = 230, `383_1` = 137, TA1_1 = 162, POLYA_1 = 159, PFPK2_1 = 162, PFPK2_2 = 171, `2490_1` = 81, TA109_1 = 160)),
    list(id = "EXAMPLE-002", day = 28, site = "SiteA",
         alleles = al(`313_1` = 230, `383_1` = 137, TA1_1 = 162, POLYA_1 = 159, PFPK2_1 = 171, `2490_1` = 81, TA109_1 = 160))),
  add_rows = list(
    list(id = "EXAMPLE-101", day = 0, site = "SiteA",
         alleles = al(`313_1` = 242, `383_1` = 149, TA1_1 = 165, POLYA_1 = 150, PFPK2_1 = 180, `2490_1` = 78, TA109_1 = 160))),
  allele_desc = data.frame(
    column = "313_1, 313_2, ...",
    description = "Fragment length (bp) of each allele, as a number. Alleles are binned using the marker's repeat length from the marker details file.",
    example = "248", stringsAsFactors = FALSE),
  extra_notes = "Enter fragment lengths as plain numbers (no 'bp').")

# ---- MSP1 / MSP2 / GLURP ----------------------------------------------------
build("MalReBay_template_msp_glurp.xlsx", "MSP1 / MSP2 / GLURP",
  markers = c(K1 = 4, MAD20 = 3, RO33 = 3, "3D7" = 5, FC27 = 4, glurp = 3),
  sep = "_",
  lf_rows = list(
    list(id = "EXAMPLE-001", day = 0,  site = "SiteA",
         alleles = al(K1_1 = 182, K1_2 = 219, RO33_1 = 150, `3D7_1` = 177, FC27_1 = 399, glurp_1 = 615)),
    list(id = "EXAMPLE-001", day = 42, site = "SiteA",
         alleles = al(MAD20_1 = 187, `3D7_1` = 333, glurp_1 = 941)),
    list(id = "EXAMPLE-002", day = 0,  site = "SiteA",
         alleles = al(K1_1 = 182, K1_2 = 219, MAD20_1 = 183, RO33_1 = 150, `3D7_1` = 389, FC27_1 = 388, glurp_1 = 903)),
    list(id = "EXAMPLE-002", day = 42, site = "SiteA",
         alleles = al(K1_1 = 182, MAD20_1 = 183, RO33_1 = 150, `3D7_1` = 389, FC27_1 = 388, glurp_1 = 903))),
  add_rows = list(
    list(id = "EXAMPLE-101", day = 0, site = "SiteA",
         alleles = al(K1_1 = 200, MAD20_1 = 177, `3D7_1` = 166, FC27_1 = 396, glurp_1 = 988))),
  allele_desc = data.frame(
    column = "K1_1, MAD20_1, RO33_1, 3D7_1, FC27_1, glurp_1, ...",
    description = "Fragment length (bp) of each allele, as a number, per allelic family: MSP1 families K1, MAD20, RO33; MSP2 families 3D7 (IC) and FC27; and GLURP.",
    example = "182", stringsAsFactors = FALSE),
  extra_notes = c(
    "Keep MSP1 and MSP2 split by allelic family as above: MalReBay recognises the families from the column names (K1, MAD20, RO33 / R033; 3D7 or IC; FC27) and combines them into MSP1 and MSP2 calls.",
    "The third marker can be glurp or a single microsatellite; use exactly one."))

# ---- Amplicon sequencing ----------------------------------------------------
build("MalReBay_template_ampseq.xlsx", "amplicon sequencing (AmpSeq)",
  markers = c(cpmp = 5, cpp = 4, amaD3 = 6),
  sep = "_allele_",
  lf_rows = list(
    list(id = "EXAMPLE-001", day = 0,  site = "SiteA",
         alleles = al(cpmp_allele_1 = "cpmp-1", cpp_allele_1 = "cpp-1", amaD3_allele_1 = "ama1-D3-1")),
    list(id = "EXAMPLE-001", day = 19, site = "SiteA",
         alleles = al(cpmp_allele_1 = "cpmp-1", cpp_allele_1 = "cpp-1", amaD3_allele_1 = "ama1-D3-1")),
    list(id = "EXAMPLE-002", day = 0,  site = "SiteA",
         alleles = al(cpmp_allele_1 = "cpmp-24", cpmp_allele_2 = "cpmp-3", cpp_allele_1 = "cpp-14", amaD3_allele_1 = "ama1-D3-12")),
    list(id = "EXAMPLE-002", day = 28, site = "SiteA",
         alleles = al(cpmp_allele_1 = "cpmp-7", cpp_allele_1 = "cpp-2", amaD3_allele_1 = "ama1-D3-4"))),
  add_rows = list(
    list(id = "EXAMPLE-101", day = 0, site = "SiteA",
         alleles = al(cpmp_allele_1 = "cpmp-5", cpp_allele_1 = "cpp-9", amaD3_allele_1 = "ama1-D3-2"))),
  allele_desc = data.frame(
    column = "cpmp_allele_1, cpmp_allele_2, ...",
    description = "Haplotype name/label of each allele (text). Matched exactly between Day 0 and recurrence, so use the same labels across all samples.",
    example = "cpmp-1", stringsAsFactors = FALSE),
  extra_notes = "Use consistent haplotype labels across all samples and both sheets: 'cpmp-1' and 'CPMP-1' are treated as different haplotypes.")

