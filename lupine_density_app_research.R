# =============================================================================
# Two-species range overlap explorer for iNaturalist / GBIF occurrence exports
#
# Run with:  shiny::runApp("lupine_density_app_research.R")
# Data:      put the exported files in ./data (or next to this script).
#            Both GBIF Darwin Core downloads (tab-separated) and iNaturalist
#            "export observations" CSVs are recognised automatically.
# =============================================================================

# ---- Settings ----------------------------------------------------------------
APP_TITLE    <- "Lupine range overlap"
APP_SUBTITLE <- "iNaturalist research-grade observations, North America"
DATA_DIR     <- if (dir.exists("data")) "data" else getwd()
MAX_POINTS   <- 4000   # cap on individually drawn records (map speed)

COL_A    <- "#5E3C99"  # species A (purple)
COL_B    <- "#E66101"  # species B (orange)
COL_BOTH <- "#1A9E8F"  # cells holding both species (teal)

# North America Albers equal-area conic: cells have equal area everywhere
AEA <- "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

# ---- Packages ----------------------------------------------------------------
required_packages <- c("shiny", "bslib", "leaflet", "leaflet.extras", "dplyr",
                       "readr", "ggplot2", "sf", "maps", "DT")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) {
  install.packages(new_packages, repos = "https://cran.r-project.org/")
}
suppressPackageStartupMessages({
  library(shiny); library(bslib); library(leaflet); library(leaflet.extras)
  library(dplyr); library(readr); library(ggplot2); library(sf); library(DT)
})

# ---- Reading and cleaning ----------------------------------------------------
# Take the first non-empty value across several possible column names, so GBIF
# and iNaturalist exports end up in the same schema.
coalesce_cols <- function(df, candidates) {
  out <- rep(NA_character_, nrow(df))
  for (cl in candidates[candidates %in% names(df)]) {
    v <- df[[cl]]
    v[!nzchar(v)] <- NA_character_
    out <- dplyr::coalesce(out, v)
  }
  out
}

read_one <- function(path) {
  header <- readLines(path, n = 1, warn = FALSE)
  if (length(header) == 0L || !nzchar(header[[1]])) {
    stop("Cannot read an empty file or a file without a header: ", basename(path))
  }
  n_tab   <- nchar(gsub("[^\t]", "", header))
  n_comma <- nchar(gsub("[^,]", "", header))
  is_tsv  <- n_tab > n_comma
  df <- read_delim(
    path,
    delim = if (is_tsv) "\t" else ",",
    quote = if (is_tsv) "" else "\"",   # GBIF TSVs contain unescaped quotes
    col_types = cols(.default = col_character()),
    na = c("", "NA"), progress = FALSE, show_col_types = FALSE
  )

  lat <- suppressWarnings(as.numeric(coalesce_cols(df, c("decimalLatitude", "latitude"))))
  lon <- suppressWarnings(as.numeric(coalesce_cols(df, c("decimalLongitude", "longitude"))))
  if (all(is.na(lat)) || all(is.na(lon))) {
    message("Skipping ", basename(path), " (no latitude/longitude columns).")
    return(NULL)
  }

  rec_id   <- coalesce_cols(df, c("gbifID", "id", "occurrenceID"))
  date_chr <- coalesce_cols(df, c("eventDate", "observed_on", "time_observed_at"))
  obs_date <- suppressWarnings(as.Date(substr(date_chr, 1, 10), format = "%Y-%m-%d"))
  yr_col   <- suppressWarnings(as.integer(coalesce_cols(df, "year")))
  mo_col   <- suppressWarnings(as.integer(coalesce_cols(df, "month")))
  quality  <- tolower(coalesce_cols(df, "quality_grade"))
  obscured <- tolower(coalesce_cols(df, "coordinates_obscured")) %in% c("true", "t", "1")
  status   <- toupper(coalesce_cols(df, "occurrenceStatus"))
  sp_raw   <- coalesce_cols(df, c("species", "scientific_name", "scientificName"))
  url      <- coalesce_cols(df, c("url", "references"))

  tibble(
    record_id = rec_id,
    # collapse subspecies/varieties to binomials
    species   = sub("^([A-Z][A-Za-z-]+ [a-z-]+).*$", "\\1", sp_raw),
    obs_date  = obs_date,
    year      = ifelse(!is.na(obs_date), as.integer(format(obs_date, "%Y")), yr_col),
    month     = ifelse(!is.na(obs_date), as.integer(format(obs_date, "%m")), mo_col),
    latitude  = lat,
    longitude = lon,
    unc_km    = suppressWarnings(as.numeric(coalesce_cols(
                  df, c("coordinateUncertaintyInMeters", "positional_accuracy")))) / 1000,
    obscured  = obscured,
    country   = coalesce_cols(df, c("countryCode", "place_country_name")),
    state     = coalesce_cols(df, c("stateProvince", "place_state_name")),
    url       = ifelse(!is.na(url), url,
                       ifelse(!is.na(rec_id) & !is.na(coalesce_cols(df, "gbifID")),
                              paste0("https://www.gbif.org/occurrence/", rec_id), NA_character_)),
    quality   = quality,
    status    = status
  )
}

read_occurrences <- function(root_dir = DATA_DIR) {
  files <- list.files(root_dir, recursive = TRUE, pattern = "\\.(csv|tsv|txt)$",
                      full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) {
    stop("No data files found in '", root_dir, "'. Put your iNaturalist/GBIF exports in a 'data' folder.")
  }
  occ <- bind_rows(lapply(files, read_one))
  if (nrow(occ) == 0) stop("No usable records found in the data files.")

  occ <- occ %>%
    filter(
      is.na(record_id) | !duplicated(record_id),           # same record in two files
      !is.na(species), !is.na(year), year > 0,
      between(latitude, -90, 90), between(longitude, -180, 180),
      !(latitude == 0 & longitude == 0),
      is.na(status) | status != "ABSENT",
      is.na(quality) | quality == "research"
    )

  # Project once so cell assignment is cheap later
  xy <- occ %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(AEA) %>%
    st_coordinates()
  occ$x <- xy[, 1]
  occ$y <- xy[, 2]
  occ %>% select(-status, -quality)
}

occ <- read_occurrences()

species_counts <- sort(table(occ$species), decreasing = TRUE)
if (length(species_counts) < 2) stop("Need at least two species in the data to compare.")
species_choices <- setNames(names(species_counts),
                            paste0(names(species_counts), " (", format(as.integer(species_counts), big.mark = ","), ")"))
all_years  <- range(occ$year)
countries  <- sort(unique(na.omit(occ$country)))

# Country / state outlines for maps
borders <- tryCatch({
  w <- st_as_sf(maps::map("world", regions = c("USA", "Canada", "Mexico"), fill = TRUE, plot = FALSE))
  s <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
  list(world  = st_transform(st_make_valid(w), AEA),
       states = st_transform(st_make_valid(s), AEA))
}, error = function(e) NULL)
borders_ll <- if (!is.null(borders)) {
  boundary_lines <- function(x) {
    st_sf(geometry = st_boundary(st_geometry(st_transform(x, 4326))))
  }
  list(world = boundary_lines(borders$world),
       states = boundary_lines(borders$states))
} else {
  NULL
}

# ---- Helpers -----------------------------------------------------------------
heat_gradient <- function(col) {
  c(colorRampPalette(c("#FFFFFF", col))(4)[2], col, colorRampPalette(c(col, "#000000"))(4)[2])
}
fmt_int <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)

theme_range <- function() {
  theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "top",
          legend.text = element_text(face = "italic"),
          plot.title.position = "plot")
}

# ---- Theme and UI --------------------------------------------------------------
app_theme <- bs_theme(
  version = 5,
  bg = "#F7F6FB", fg = "#1F2133",
  primary = COL_A, secondary = "#6B6F85",
  base_font    = font_collection(font_google("Public Sans", local = FALSE), "system-ui", "sans-serif"),
  heading_font = font_collection(font_google("Newsreader", local = FALSE), "Georgia", "serif")
) %>%
  bs_add_rules("
    .app-title { font-family: 'Newsreader', Georgia, serif; font-size: 1.4rem; font-weight: 600; }
    .app-sub   { margin-left: .8rem; color: #6B6F85; font-size: .9rem; }
    .sp-a .selectize-input { border-left: 6px solid #5E3C99; }
    .sp-b .selectize-input { border-left: 6px solid #E66101; }
    .accordion-button { font-weight: 600; }
    .map-note  { color: #6B6F85; font-size: .85rem; margin: .4rem 0 0; }
    .status    { color: #3A3D52; font-size: .9rem; }
    .methods   { white-space: pre-wrap; font-family: inherit; background: transparent; border: 0; padding: 0; }
  ")

ui <- page_sidebar(
  title = tagList(span(APP_TITLE, class = "app-title"), span(APP_SUBTITLE, class = "app-sub")),
  theme = app_theme,
  sidebar = sidebar(
    width = 340,
    accordion(
      multiple = TRUE, open = c("Species", "Records"),
      accordion_panel(
        "Species",
        div(class = "sp-a", selectizeInput("sp_a", "Species A", choices = species_choices,
                                           selected = names(species_counts)[1])),
        div(class = "sp-b", selectizeInput("sp_b", "Species B", choices = species_choices,
                                           selected = names(species_counts)[2]))
      ),
      accordion_panel(
        "Records",
        sliderInput("years", "Years", min = all_years[1], max = all_years[2],
                    value = all_years, step = 1, sep = ""),
        checkboxGroupInput("months", "Months", choices = setNames(1:12, month.abb),
                           selected = 1:12, inline = TRUE),
        sliderInput("max_unc", "Max coordinate uncertainty (km)", min = 1, max = 100, value = 10, step = 1),
        checkboxInput("keep_unknown", "Keep records with no uncertainty value", TRUE),
        checkboxInput("drop_obscured", "Drop obscured (jittered) locations", TRUE),
        selectizeInput("countries", "Countries (blank = all)", choices = countries, multiple = TRUE)
      ),
      accordion_panel(
        "Grid analysis",
        sliderInput("cell_km", "Grid cell size (km)", min = 10, max = 300, value = 100, step = 1),
        checkboxInput("thin", "One record per species, cell and year", TRUE),
        helpText("Thinning reduces the pull of heavily visited sites on overlap statistics.")
      ),
      accordion_panel(
        "Display",
        sliderInput("heat_radius", "Heat radius (px)", min = 5, max = 40, value = 15, step = 1),
        sliderInput("heat_blur", "Heat blur (px)", min = 5, max = 60, value = 25, step = 1),
        sliderInput("heat_max", "Saturation point", min = 0.05, max = 1, value = 0.6, step = 0.05),
        checkboxInput("equalize", "Give both species equal total weight", TRUE),
        checkboxInput("show_points", "Show individual records", FALSE),
        textOutput("points_note")
      )
    )
  ),

  navset_card_underline(
    id = "tabs",

    nav_panel(
      "Map",
      div(class = "d-flex flex-wrap align-items-center gap-3",
          radioButtons("map_mode", NULL, inline = TRUE,
                       choices = c("Density" = "heat", "Occupied cells" = "cats", "Share ratio" = "ratio")),
          actionButton("fit", "Zoom to records", class = "btn-sm btn-outline-primary"),
          span(class = "status", textOutput("status", inline = TRUE))),
      leafletOutput("map", height = "calc(100vh - 235px)"),
      p(class = "map-note", textOutput("mode_note", inline = TRUE))
    ),

    nav_panel(
      "Overlap",
      layout_columns(
        col_widths = c(5, 7),
        card(card_header("Overlap metrics"), tableOutput("stats_table")),
        card(card_header("Records per cell"), plotOutput("cell_scatter", height = 340))
      ),
      layout_columns(
        col_widths = c(5, 7),
        card(
          card_header("Do the distributions differ?"),
          selectInput("n_perm", "Permutations", choices = c(199, 499, 999), selected = 499),
          actionButton("run_perm", "Run permutation test", class = "btn-primary"),
          verbatimTextOutput("perm_result"),
          helpText("Species labels are shuffled among the pooled records, which holds observer effort constant. ",
                   "A low P means the species are more spatially distinct than shared sampling effort would produce.")
        ),
        card(card_header("Methods text"), tags$pre(class = "methods", textOutput("methods")))
      )
    ),

    nav_panel(
      "Phenology and trends",
      layout_columns(
        col_widths = c(6, 6, 12),
        card(card_header("Flowering-season timing (share of each species' records)"), plotOutput("plot_month", height = 300)),
        card(card_header("Records per year"), plotOutput("plot_year", height = 300)),
        card(card_header("Latitudinal distribution"), plotOutput("plot_lat", height = 300))
      )
    ),

    nav_panel(
      "Figure",
      layout_sidebar(
        sidebar = sidebar(
          width = 240, position = "right",
          helpText("Static equal-area map for slides or manuscripts. It uses the grid and the cell display chosen on the Map tab ",
                   "(density view falls back to occupied cells)."),
          downloadButton("dl_png", "Download PNG (300 dpi)"),
          br(), br(),
          downloadButton("dl_pdf", "Download PDF")
        ),
        plotOutput("static_map", height = 640)
      )
    ),

    nav_panel(
      "Data",
      div(class = "d-flex gap-2 mb-2",
          downloadButton("dl_records", "Filtered records (CSV)"),
          downloadButton("dl_cells", "Grid cells (CSV)"),
          downloadButton("dl_stats", "Overlap metrics (CSV)")),
      DTOutput("records_table")
    )
  )
)

# ---- Server ------------------------------------------------------------------
server <- function(input, output, session) {

  sp_cols <- reactive(c(setNames(COL_A, input$sp_a), setNames(COL_B, input$sp_b)))
  cat_cols <- reactive({
    setNames(c(COL_A, COL_BOTH, COL_B),
             c(paste("Only", input$sp_a), "Both", paste("Only", input$sp_b)))
  })

  # Records after filters, with grid cell assignment
  recs <- reactive({
    validate(need(input$sp_a != input$sp_b, "Choose two different species."))
    df <- occ %>%
      filter(species %in% c(input$sp_a, input$sp_b),
             year >= input$years[1], year <= input$years[2],
             is.na(month) | month %in% as.integer(input$months))

    unc_ok <- (!is.na(df$unc_km) & df$unc_km <= input$max_unc) |
              (is.na(df$unc_km) & input$keep_unknown)
    df <- df[unc_ok, ]
    if (input$drop_obscured) df <- df[!df$obscured, ]
    if (length(input$countries)) df <- df[df$country %in% input$countries, ]
    validate(need(nrow(df) > 0, "No records match the current filters."))

    s <- input$cell_km * 1000
    df %>% mutate(grp = if_else(species == input$sp_a, "A", "B"),
                  cx = floor(x / s), cy = floor(y / s),
                  cell_id = paste(cx, cy))
  })

  # Records used for grid statistics (optionally thinned)
  recs_thin <- reactive({
    df <- recs()
    if (input$thin) df <- distinct(df, species, cell_id, year, .keep_all = TRUE)
    df
  })

  # One row per occupied cell
  cells <- reactive({
    cl <- recs_thin() %>%
      group_by(cx, cy) %>%
      summarise(a = sum(grp == "A"), b = sum(grp == "B"), .groups = "drop")
    validate(need(sum(cl$a) > 0 && sum(cl$b) > 0, "Both species need at least one record."))
    lv <- c(paste("Only", input$sp_a), "Both", paste("Only", input$sp_b))
    cl %>% mutate(
      pa = a / sum(a), pb = b / sum(b),
      category = factor(case_when(a > 0 & b > 0 ~ "Both",
                                  a > 0 ~ lv[1], TRUE ~ lv[3]), levels = lv),
      log2_ratio = if_else(a > 0 & b > 0, log2(pa / pb), NA_real_)
    )
  })

  cells_aea <- reactive({
    cl <- cells(); s <- input$cell_km * 1000
    polys <- lapply(seq_len(nrow(cl)), function(i) {
      x0 <- cl$cx[i] * s; y0 <- cl$cy[i] * s
      st_polygon(list(matrix(c(x0, y0, x0 + s, y0, x0 + s, y0 + s, x0, y0 + s, x0, y0),
                             ncol = 2, byrow = TRUE)))
    })
    st_sf(cl, geometry = st_sfc(polys, crs = AEA))
  })
  cells_ll <- reactive(st_transform(cells_aea(), 4326))

  # ---- Overlap statistics ----
  stats_list <- reactive({
    cl <- cells(); df <- recs_thin()
    a <- cl$a; b <- cl$b; both <- a > 0 & b > 0
    da <- df[df$grp == "A", ]; db <- df[df$grp == "B", ]
    list(
      n_a = sum(a), n_b = sum(b),
      cells_a = sum(a > 0), cells_b = sum(b > 0),
      shared = sum(both), occupied = nrow(cl),
      jaccard = sum(both) / nrow(cl),
      a_in_b = sum(both) / sum(a > 0),
      b_in_a = sum(both) / sum(b > 0),
      schoener = 1 - 0.5 * sum(abs(cl$pa - cl$pb)),
      rho = tryCatch(suppressWarnings(cor(a, b, method = "spearman")), error = function(e) NA_real_),
      centroid_km = sqrt((mean(da$x) - mean(db$x))^2 + (mean(da$y) - mean(db$y))^2) / 1000,
      lat_a = median(da$latitude), lat_b = median(db$latitude)
    )
  })

  stats_df <- reactive({
    s <- stats_list(); a <- input$sp_a; b <- input$sp_b
    data.frame(
      Metric = c(paste("Records in analysis:", a), paste("Records in analysis:", b),
                 paste("Occupied cells:", a), paste("Occupied cells:", b),
                 "Cells holding both species",
                 "Jaccard index of occupied cells",
                 paste0("Share of ", a, " cells also holding ", b),
                 paste0("Share of ", b, " cells also holding ", a),
                 "Schoener's D (0 = no overlap, 1 = identical)",
                 "Spearman correlation of cell counts",
                 "Distance between record centroids (km)",
                 paste("Median latitude:", a), paste("Median latitude:", b)),
      Value = c(fmt_int(s$n_a), fmt_int(s$n_b), fmt_int(s$cells_a), fmt_int(s$cells_b), fmt_int(s$shared),
                sprintf("%.2f", s$jaccard), sprintf("%.2f", s$a_in_b), sprintf("%.2f", s$b_in_a),
                sprintf("%.2f", s$schoener), sprintf("%.2f", s$rho), fmt_int(s$centroid_km),
                sprintf("%.1f", s$lat_a), sprintf("%.1f", s$lat_b)),
      stringsAsFactors = FALSE
    )
  })
  output$stats_table <- renderTable(stats_df(), striped = TRUE, spacing = "s", width = "100%")

  perm <- eventReactive(input$run_perm, {
    df <- recs_thin()
    idx <- match(df$cell_id, unique(df$cell_id))
    K <- max(idx); lab <- df$grp == "A"
    D <- function(l) {
      a <- tabulate(idx[l], K); b <- tabulate(idx[!l], K)
      1 - 0.5 * sum(abs(a / sum(a) - b / sum(b)))
    }
    n <- as.integer(input$n_perm)
    obs <- D(lab)
    null <- numeric(n)
    withProgress(message = "Running permutations", value = 0, {
      for (i in seq_len(n)) {
        null[i] <- D(sample(lab))
        if (i %% 25 == 0) incProgress(25 / n)
      }
    })
    list(obs = obs, null = null, n = n, p = (1 + sum(null <= obs)) / (n + 1))
  })
  output$perm_result <- renderText({
    r <- perm()
    q <- quantile(r$null, c(0.025, 0.975))
    sprintf("Observed D = %.3f\nNull mean = %.3f (95%% range %.3f to %.3f)\nP = %.3f (%d permutations, one-sided)",
            r$obs, mean(r$null), q[1], q[2], r$p, r$n)
  })

  output$methods <- renderText({
    s <- stats_list(); df <- recs()
    unc <- paste0(" with coordinate uncertainty of ", input$max_unc, " km or less",
                  if (input$keep_unknown) " (or unreported)" else "")
    paste0(
      "Research-grade iNaturalist observations of ", input$sp_a, " (n = ", fmt_int(sum(df$grp == "A")),
      ") and ", input$sp_b, " (n = ", fmt_int(sum(df$grp == "B")), ") from ", input$years[1], " to ", input$years[2],
      " were retained", unc, if (input$drop_obscured) ", excluding obscured locations" else "",
      ". Records were assigned to ", input$cell_km, "-km equal-area grid cells (Albers conic)",
      if (input$thin) " and thinned to one record per species, cell and year" else "",
      ". The species shared ", fmt_int(s$shared), " of ", fmt_int(s$occupied), " occupied cells (Jaccard = ",
      sprintf("%.2f", s$jaccard), "), and overlap in the distribution of records across cells was Schoener's D = ",
      sprintf("%.2f", s$schoener), "."
    )
  })

  # ---- Map ----
  output$map <- renderLeaflet({
    leaflet(options = leafletOptions(minZoom = 2, preferCanvas = TRUE)) %>%
      addProviderTiles(providers$OpenStreetMap.Mapnik, group = "Light") %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
      setView(-98, 45, 3) %>%
      addScaleBar("bottomleft")
  })

  observeEvent(list(input$fit, input$sp_a, input$sp_b), {
    d <- tryCatch(recs(), error = function(e) NULL)
    if (is.null(d) || nrow(d) == 0) return()
    q <- function(v, p) unname(quantile(v, p))
    leafletProxy("map") %>%
      fitBounds(q(d$longitude, 0.01), q(d$latitude, 0.01), q(d$longitude, 0.99), q(d$latitude, 0.99))
  })

  observe({
    d <- recs(); cl <- cells_ll(); mode <- input$map_mode
    a <- input$sp_a; b <- input$sp_b
    px <- leafletProxy("map") %>%
      clearHeatmap() %>% clearGroup(c("Grid cells", "Records", "Boundaries")) %>% removeControl("legend")
    overlays <- character()

    if (mode == "heat") {
      n_a <- sum(d$grp == "A"); n_b <- sum(d$grp == "B")
      for (g in c("A", "B")) {
        dd <- d[d$grp == g, ]
        if (!nrow(dd)) next
        w <- if (input$equalize) min(n_a, n_b) / nrow(dd) else 1
        px <- addHeatmap(px, data = dd, lng = ~longitude, lat = ~latitude,
                         intensity = rep(w, nrow(dd)),
                         radius = input$heat_radius, blur = input$heat_blur, max = input$heat_max,
                         minOpacity = 0.25, gradient = heat_gradient(if (g == "A") COL_A else COL_B),
                         group = if (g == "A") a else b)
      }
      overlays <- c(a, b)
      px <- addLegend(px, position = "bottomright", colors = c(COL_A, COL_B), labels = c(a, b),
                      opacity = 1, title = "Density", layerId = "legend")
    } else {
      popup <- sprintf("<i>%s</i>: %s records<br><i>%s</i>: %s records",
                       a, fmt_int(cl$a), b, fmt_int(cl$b))
      if (mode == "cats") {
        pal <- colorFactor(c(COL_A, COL_BOTH, COL_B), levels = levels(cl$category))
        fill <- pal(cl$category); op <- 0.75
        px <- addLegend(px, position = "bottomright", pal = pal, values = cl$category,
                        opacity = 1, title = paste0(input$cell_km, "-km cells"), layerId = "legend")
      } else {
        m <- suppressWarnings(max(abs(cl$log2_ratio), na.rm = TRUE))
        lim <- if (is.finite(m)) max(1, ceiling(m)) else 1
        pal <- colorNumeric(c(COL_B, "#F7F7F7", COL_A), domain = c(-lim, lim), na.color = "#CFCFE0")
        fill <- pal(cl$log2_ratio); op <- ifelse(is.na(cl$log2_ratio), 0.35, 0.8)
        px <- addLegend(px, position = "bottomright", pal = pal, values = c(-lim, lim), opacity = 1,
                        title = paste0("log2 share<br>", a, " / ", b), layerId = "legend")
      }
      px <- addPolygons(px, data = cl, fillColor = fill, fillOpacity = op,
                        color = "#FFFFFF", weight = 0.6, opacity = 0.6, popup = popup, group = "Grid cells",
                        highlightOptions = highlightOptions(weight = 2, color = "#1F2133", bringToFront = TRUE))
      overlays <- "Grid cells"
    }

    if (input$show_points) {
      dd <- d
      if (nrow(dd) > MAX_POINTS) dd <- dd[sample.int(nrow(dd), MAX_POINTS), ]
      link <- ifelse(is.na(dd$url), "",
                     sprintf("<br><a href='%s' target='_blank' rel='noopener'>Open record</a>", dd$url))
      px <- addCircleMarkers(px, data = dd, lng = ~longitude, lat = ~latitude,
                             radius = 3.5, stroke = TRUE, color = "#FFFFFF", weight = 0.8,
                             fillColor = ifelse(dd$grp == "A", COL_A, COL_B), fillOpacity = 0.85,
                             popup = sprintf("<i>%s</i><br>%s%s", dd$species,
                                             ifelse(is.na(dd$obs_date), dd$year, format(dd$obs_date)), link),
                             group = "Records")
      overlays <- c(overlays, "Records")
    }

    if (!is.null(borders_ll)) {
      px <- px %>%
        addPolylines(data = borders_ll$world, color = "#4C4C4C", weight = 1.2,
                     opacity = 0.9, group = "Boundaries") %>%
        addPolylines(data = borders_ll$states, color = "#686868", weight = 0.6,
                     opacity = 0.75, group = "Boundaries")
      overlays <- c(overlays, "Boundaries")
    }

    addLayersControl(px, baseGroups = c("Light", "Satellite"), overlayGroups = overlays,
                     options = layersControlOptions(collapsed = TRUE))
  })

  output$status <- renderText({
    d <- recs()
    sprintf("%s: %s records.  %s: %s records.",
            input$sp_a, fmt_int(sum(d$grp == "A")), input$sp_b, fmt_int(sum(d$grp == "B")))
  })
  output$points_note <- renderText({
    n <- nrow(recs())
    if (input$show_points && n > MAX_POINTS)
      sprintf("Drawing a random %s of %s records.", fmt_int(MAX_POINTS), fmt_int(n)) else ""
  })
  output$mode_note <- renderText({
    switch(input$map_mode,
           heat  = "Kernel density of records. Each species has its own color; overlap appears where both are dense.",
           cats  = paste0("Equal-area ", input$cell_km, "-km cells occupied by one or both species (after thinning, if selected)."),
           ratio = paste0("Cells holding both species, colored by relative share of records: purple means a larger share of ",
                          input$sp_a, ", orange a larger share of ", input$sp_b, ". Grey cells hold only one species."))
  })

  # ---- Phenology and trends ----
  output$plot_month <- renderPlot({
    d <- recs() %>% filter(!is.na(month)) %>% count(species, month) %>%
      group_by(species) %>% mutate(share = n / sum(n)) %>% ungroup()
    ggplot(d, aes(month, share, colour = species)) +
      geom_line(linewidth = 1.1) + geom_point(size = 2) +
      scale_x_continuous(breaks = 1:12, labels = month.abb) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      scale_colour_manual(values = sp_cols(), name = NULL) +
      labs(x = NULL, y = "Share of records") + theme_range()
  })
  output$plot_year <- renderPlot({
    d <- recs() %>% count(species, year)
    ggplot(d, aes(factor(year), n, fill = species)) +
      geom_col(position = "dodge", width = 0.75) +
      scale_fill_manual(values = sp_cols(), name = NULL) +
      labs(x = NULL, y = "Records", caption = "The most recent year may be incomplete.") + theme_range()
  })
  output$plot_lat <- renderPlot({
    ggplot(recs(), aes(latitude, fill = species, colour = species)) +
      geom_density(alpha = 0.35, linewidth = 0.8) +
      scale_fill_manual(values = sp_cols(), name = NULL) +
      scale_colour_manual(values = sp_cols(), name = NULL) +
      labs(x = "Latitude (degrees N)", y = "Density") + theme_range()
  })

  output$cell_scatter <- renderPlot({
    cl <- cells()
    ggplot(cl, aes(log10(1 + a), log10(1 + b), colour = category)) +
      geom_abline(slope = 1, intercept = 0, colour = "grey70", linetype = "dashed") +
      geom_jitter(width = 0.04, height = 0.04, alpha = 0.7, size = 2) +
      scale_colour_manual(values = cat_cols(), name = NULL, drop = FALSE) +
      labs(x = paste0("log10(1 + records): ", input$sp_a),
           y = paste0("log10(1 + records): ", input$sp_b)) +
      theme_range() + theme(legend.text = element_text(face = "plain"))
  })

  # ---- Static figure ----
  static_plot <- reactive({
    cl <- cells_aea(); s <- input$cell_km * 1000
    bb <- st_bbox(cl); pad <- 2 * s
    p <- ggplot()
    if (!is.null(borders)) p <- p + geom_sf(data = borders$world, fill = "#EDECF4", colour = NA)
    if (input$map_mode == "ratio") {
      m <- suppressWarnings(max(abs(cl$log2_ratio), na.rm = TRUE))
      lim <- if (is.finite(m)) max(1, ceiling(m)) else 1
      p <- p + geom_sf(data = cl, aes(fill = log2_ratio), colour = NA) +
        scale_fill_gradient2(low = COL_B, mid = "#F7F7F7", high = COL_A, midpoint = 0,
                             limits = c(-lim, lim), na.value = "#CFCFE0",
                             name = paste0("log2 share\n", input$sp_a, " / ", input$sp_b))
    } else {
      p <- p + geom_sf(data = cl, aes(fill = category), colour = NA) +
        scale_fill_manual(values = cat_cols(), drop = FALSE, name = NULL)
    }
    if (!is.null(borders)) {
      p <- p + geom_sf(data = borders$states, fill = NA, colour = "#A9A8BC", linewidth = 0.15) +
        geom_sf(data = borders$world, fill = NA, colour = "#7E7D93", linewidth = 0.25)
    }
    p + coord_sf(crs = AEA, default_crs = AEA,
                 xlim = c(bb[["xmin"]] - pad, bb[["xmax"]] + pad),
                 ylim = c(bb[["ymin"]] - pad, bb[["ymax"]] + pad), expand = FALSE) +
      labs(caption = paste0(input$cell_km, "-km equal-area grid. iNaturalist research-grade observations, ",
                            input$years[1], "-", input$years[2], ".")) +
      theme_void(base_size = 13) +
      theme(legend.position = "bottom", legend.text = element_text(face = "italic"),
            plot.background = element_rect(fill = "white", colour = NA),
            plot.caption = element_text(colour = "#6B6F85", hjust = 0, margin = margin(t = 8)),
            plot.margin = margin(10, 10, 10, 10))
  })
  output$static_map <- renderPlot(static_plot(), res = 110)

  # ---- Data tab and downloads ----
  records_out <- reactive({
    recs() %>% select(species, obs_date, year, month, latitude, longitude, unc_km, obscured,
                      country, state, url, record_id)
  })
  output$records_table <- renderDT(
    datatable(records_out(), rownames = FALSE, filter = "top",
              options = list(pageLength = 15, scrollX = TRUE))
  )
  output$dl_records <- downloadHandler(
    filename = function() paste0("filtered_records_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(records_out(), file, row.names = FALSE, na = ""))
  output$dl_cells <- downloadHandler(
    filename = function() paste0("grid_cells_", input$cell_km, "km_", Sys.Date(), ".csv"),
    content  = function(file) {
      cl <- cells(); s <- input$cell_km * 1000
      ctr <- st_transform(st_as_sf(data.frame(x = (cl$cx + 0.5) * s, y = (cl$cy + 0.5) * s),
                                   coords = c("x", "y"), crs = AEA), 4326)
      out <- cbind(cell_center_lon = st_coordinates(ctr)[, 1],
                   cell_center_lat = st_coordinates(ctr)[, 2],
                   setNames(as.data.frame(cl[, c("a", "b", "pa", "pb", "category", "log2_ratio")]),
                            c(paste0("n_", input$sp_a), paste0("n_", input$sp_b),
                              "share_a", "share_b", "category", "log2_share_ratio_a_over_b")))
      write.csv(out, file, row.names = FALSE, na = "")
    })
  output$dl_stats <- downloadHandler(
    filename = function() paste0("overlap_metrics_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(stats_df(), file, row.names = FALSE))
  output$dl_png <- downloadHandler(
    filename = function() paste0("range_overlap_", Sys.Date(), ".png"),
    content  = function(file) ggsave(file, static_plot(), width = 10, height = 7, dpi = 300, bg = "white"))
  output$dl_pdf <- downloadHandler(
    filename = function() paste0("range_overlap_", Sys.Date(), ".pdf"),
    content  = function(file) ggsave(file, static_plot(), width = 10, height = 7, device = "pdf"))
}

shinyApp(ui = ui, server = server)
