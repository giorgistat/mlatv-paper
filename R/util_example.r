# util_example.r — shared helpers for the analysis scripts in test/.
# Source from the repository root with: source("test/util_example.r")

# Collapse observations to unique spatial locations, with optional subsampling.
#
# Expects a data frame that already carries raw-metre coordinates (long_m,
# lat_m, used as the location key), centred-metre coordinates (x_m, y_m), and
# the model columns age, y_1, y_2. Optionally keep only n_sub unique locations
# and cap the individuals per location, then return the arrays the spatial
# fitters consume.
#
#   dat          data frame with long_m, lat_m, x_m, y_m, age, y_1, y_2
#   n_sub        unique locations to keep, or NULL to keep all (full data)
#   max_per_loc  max individuals per unique location; applied whenever n_sub is
#                set (matches the original scripts' behaviour)
#   verbose      print the one-line summaries when TRUE
#
# Returns a list: dat, loc, loc_unique, n_unique, n_obs, age, y_1, y_2.
example_prepare_locations <- function(dat, n_sub = NULL, max_per_loc = NULL,
                                      verbose = TRUE) {
  loc_key     <- paste(dat$long_m, dat$lat_m, sep = "_")
  unique_keys <- unique(loc_key)
  if (verbose)
    cat(sprintf("Unique locations after filter: %d\n", length(unique_keys)))

  # optional subsample of unique locations
  if (!is.null(n_sub) && n_sub < length(unique_keys)) {
    chosen_keys <- sample(unique_keys, n_sub)
    dat <- dat[loc_key %in% chosen_keys, ]
    if (verbose) cat(sprintf("Subsampled to %d unique locs\n", n_sub))
  }

  # optional cap on individuals per location (any non-NULL n_sub triggers it)
  if (!is.null(n_sub)) {
    loc_key2 <- paste(dat$long_m, dat$lat_m, sep = "_")
    dat <- do.call(rbind, lapply(split(dat, loc_key2), function(g) {
      g[seq_len(min(nrow(g), max_per_loc)), ]
    }))
  }

  loc_key_final <- paste(dat$long_m, dat$lat_m, sep = "_")
  loc_uniq_keys <- unique(loc_key_final)
  first_idx     <- match(loc_uniq_keys, loc_key_final)

  loc_unique <- as.matrix(dat[first_idx, c("x_m", "y_m")])
  loc        <- as.matrix(dat[, c("x_m", "y_m")])
  colnames(loc_unique) <- colnames(loc) <- c("x", "y")

  out <- list(
    dat        = dat,
    loc        = loc,
    loc_unique = loc_unique,
    n_unique   = length(loc_uniq_keys),
    n_obs      = nrow(dat),
    age        = dat$age,
    y_1        = dat$y_1,
    y_2        = dat$y_2
  )
  if (verbose)
    cat(sprintf("n_unique: %d  n_obs: %d  age range: [%.1f, %.1f]\n",
                out$n_unique, out$n_obs, min(out$age), max(out$age)))
  out
}

# Apply the KISUMU HDSS geographic filter: drop the southern arm and the
# north-west outlier cluster. Thresholds are dataset-specific; the defaults are
# the values shared across the example scripts.
#
#   dat              data frame with columns latitude, longitude
#   drop_south_lat   drop rows south of this latitude
#   drop_nw_lon_max  \ together drop the north-west cluster: rows west of
#   drop_nw_lat_min  / drop_nw_lon_max and north of drop_nw_lat_min
#   verbose          print the removed/retained counts when TRUE
#
# Returns the filtered data frame.
example_geographic_filter <- function(dat,
                                      drop_south_lat  = -0.62,
                                      drop_nw_lon_max =  34.80,
                                      drop_nw_lat_min = -0.46,
                                      verbose = TRUE) {
  drop_flag <- dat$latitude < drop_south_lat |
               (dat$longitude < drop_nw_lon_max & dat$latitude > drop_nw_lat_min)
  if (verbose)
    cat(sprintf("Geographic filter: removed %d, retained %d\n",
                sum(drop_flag), sum(!drop_flag)))
  dat[!drop_flag, ]
}

# Build a converter from centred-metre coordinates (x_m, y_m) back to lon/lat,
# by regressing raw lon/lat on raw metres in `dat`. The spatial examples use it
# to draw meshes and fitted fields on a geographic axis.
#
#   dat          data frame with long_m, lat_m, longitude, latitude
#   site_center  named vector c(long_m =, lat_m =) used to centre the metres
#
# Returns a function: xy_centred (2-col matrix) -> data.frame(longitude, latitude).
example_lonlat_converter <- function(dat, site_center) {
  longitude_model <- stats::lm(longitude ~ long_m, data = dat)
  latitude_model  <- stats::lm(latitude  ~ lat_m,  data = dat)
  function(xy_centred) {
    long_m <- xy_centred[, 1L] + site_center["long_m"]
    lat_m  <- xy_centred[, 2L] + site_center["lat_m"]
    data.frame(
      longitude = as.numeric(stats::predict(longitude_model,
                             newdata = data.frame(long_m = long_m))),
      latitude  = as.numeric(stats::predict(latitude_model,
                             newdata = data.frame(lat_m = lat_m)))
    )
  }
}

# Save a diagnostic map of the mesh over the observation sites, on a lon/lat
# axis: all sites faint, the unique fit locations darker, mesh edges overlaid.
# Shared across the spatial examples (example 20-23).
#
#   mesh              fmesher mesh
#   dat               data frame with long_m, lat_m, longitude, latitude
#   loc_unique        2-col matrix of unique fit locations (centred metres)
#   site_center       named vector c(long_m =, lat_m =) used to centre metres
#   n_unique, n_obs   counts shown in the title
#   out_path          PDF path to write
#   width, height     device size in inches
#
# Returns the ggplot invisibly.
example_save_mesh_plot <- function(mesh, dat, loc_unique, site_center,
                                   n_unique, n_obs, out_path,
                                   width = 8, height = 5.5) {
  to_lonlat <- example_lonlat_converter(dat, site_center)

  # mesh triangle edges -> unique undirected edges -> lon/lat segments
  tri <- mesh$graph$tv
  ep  <- do.call(rbind, lapply(seq_len(nrow(tri)), function(i)
           rbind(tri[i, c(1L, 2L)], tri[i, c(2L, 3L)], tri[i, c(3L, 1L)])))
  ep  <- t(apply(ep, 1L, sort))
  ep  <- unique(data.frame(a = ep[, 1L], b = ep[, 2L]))
  nodes <- to_lonlat(mesh$loc[, 1:2, drop = FALSE])
  edges_df <- data.frame(
    longitude     = nodes$longitude[ep$a],
    latitude      = nodes$latitude [ep$a],
    longitude_end = nodes$longitude[ep$b],
    latitude_end  = nodes$latitude [ep$b]
  )

  # all sites (faint) + convex hull, and the unique fit locations (darker)
  sites_all_df  <- data.frame(longitude = dat$longitude, latitude = dat$latitude)
  sites_used_df <- to_lonlat(loc_unique)
  hull    <- grDevices::chull(sites_all_df$longitude, sites_all_df$latitude)
  hull_df <- sites_all_df[c(hull, hull[1L]), , drop = FALSE]

  mesh_plot <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = hull_df, ggplot2::aes(x = longitude, y = latitude),
      fill = "grey96", color = "grey55", linewidth = 0.35) +
    ggplot2::geom_segment(
      data = edges_df,
      ggplot2::aes(x = longitude, y = latitude,
                   xend = longitude_end, yend = latitude_end),
      linewidth = 0.18, alpha = 0.55, color = "grey35") +
    ggplot2::geom_point(
      data = sites_all_df, ggplot2::aes(x = longitude, y = latitude),
      color = "grey20", alpha = 0.035, size = 0.06) +
    ggplot2::geom_point(
      data = sites_used_df, ggplot2::aes(x = longitude, y = latitude),
      color = "black", alpha = 0.08, size = 0.08) +
    ggplot2::coord_equal(
      xlim = range(c(edges_df$longitude, edges_df$longitude_end)),
      ylim = range(c(edges_df$latitude,  edges_df$latitude_end)),
      expand = FALSE) +
    ggplot2::labs(
      title = sprintf("Mesh: %d nodes  |  %d unique locs  |  %d obs",
                      mesh$n, n_unique, n_obs),
      x = "longitude", y = "latitude") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  ggplot2::ggsave(out_path, mesh_plot, width = width, height = height)
  cat(sprintf("Saved: %s\n\n", out_path))
  invisible(mesh_plot)
}
