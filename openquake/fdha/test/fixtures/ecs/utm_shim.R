# =============================================================================
# DEV-ONLY UTM shim for the ECS oracle (NOT shipped, NOT runtime).
#
# Newer sp::spTransform delegates CRS transforms to sf, which needs a heavy
# system geo stack. The ECS reference impl only uses sp for WGS84 longlat<->UTM
# (same-datum, no datum shift), so we replace longlat2UTM / UTM2longlat with a
# self-contained Transverse Mercator (Snyder 1987, USGS PP-1395) good to ~mm in
# a 6-degree zone. Source this AFTER ecs_functions.R to override those two
# functions. Engine parity: both this shim and the Python side (pyproj) use a
# plain "+proj=utm +zone=N +datum=WGS84" definition (NO +south, FN=0), so a
# southern-hemisphere northing is negative in BOTH -- consistency is what the
# GC2 gate checks. Validated against pyproj to < 1e-3 m on the Calingiri points.
# =============================================================================

.WGS84_A  <- 6378137.0
.WGS84_F  <- 1.0 / 298.257223563
.UTM_K0   <- 0.9996
.UTM_FE   <- 500000.0
.UTM_FN   <- 0.0   # no +south, matching "+proj=utm +zone=N +datum=WGS84"

.utm_lon0_deg <- function(zone) zone * 6 - 183  # central meridian of the zone

# forward: WGS84 (lon,lat deg) -> UTM (x=easting, y=northing) m
.ll2utm_xy <- function(lon, lat, zone) {
  a <- .WGS84_A; f <- .WGS84_F; e2 <- f * (2 - f); ep2 <- e2 / (1 - e2); k0 <- .UTM_K0
  d2r <- pi / 180
  phi <- lat * d2r
  lam <- lon * d2r
  lam0 <- .utm_lon0_deg(zone) * d2r
  N <- a / sqrt(1 - e2 * sin(phi)^2)
  T <- tan(phi)^2
  C <- ep2 * cos(phi)^2
  A <- (lam - lam0) * cos(phi)
  M <- a * ((1 - e2/4 - 3*e2^2/64 - 5*e2^3/256) * phi
            - (3*e2/8 + 3*e2^2/32 + 45*e2^3/1024) * sin(2*phi)
            + (15*e2^2/256 + 45*e2^3/1024) * sin(4*phi)
            - (35*e2^3/3072) * sin(6*phi))
  x <- .UTM_FE + k0 * N * (A + (1 - T + C) * A^3 / 6
                           + (5 - 18*T + T^2 + 72*C - 58*ep2) * A^5 / 120)
  y <- .UTM_FN + k0 * (M + N * tan(phi) * (A^2/2
                       + (5 - T + 9*C + 4*C^2) * A^4 / 24
                       + (61 - 58*T + T^2 + 600*C - 330*ep2) * A^6 / 720))
  list(x = x, y = y)
}

# inverse: UTM (x,y) m -> WGS84 (lon,lat deg)
.utm2ll <- function(x, y, zone) {
  a <- .WGS84_A; f <- .WGS84_F; e2 <- f * (2 - f); ep2 <- e2 / (1 - e2); k0 <- .UTM_K0
  r2d <- 180 / pi
  e1 <- (1 - sqrt(1 - e2)) / (1 + sqrt(1 - e2))
  M <- (y - .UTM_FN) / k0
  mu <- M / (a * (1 - e2/4 - 3*e2^2/64 - 5*e2^3/256))
  phi1 <- mu + (3*e1/2 - 27*e1^3/32) * sin(2*mu) +
          (21*e1^2/16 - 55*e1^4/32) * sin(4*mu) +
          (151*e1^3/96) * sin(6*mu) + (1097*e1^4/512) * sin(8*mu)
  C1 <- ep2 * cos(phi1)^2
  T1 <- tan(phi1)^2
  N1 <- a / sqrt(1 - e2 * sin(phi1)^2)
  R1 <- a * (1 - e2) / (1 - e2 * sin(phi1)^2)^1.5
  D <- (x - .UTM_FE) / (N1 * k0)
  phi <- phi1 - (N1 * tan(phi1) / R1) * (D^2/2
         - (5 + 3*T1 + 10*C1 - 4*C1^2 - 9*ep2) * D^4 / 24
         + (61 + 90*T1 + 298*C1 + 45*T1^2 - 252*ep2 - 3*C1^2) * D^6 / 720)
  lam0 <- .utm_lon0_deg(zone) * pi / 180
  lam <- lam0 + (D - (1 + 2*T1 + C1) * D^3 / 6
         + (5 - 2*C1 + 28*T1 - 3*C1^2 + 8*ep2 + 24*T1^2) * D^5 / 120) / cos(phi1)
  list(lon = lam * r2d, lat = phi * r2d)
}

# --- overrides matching ecs_functions.R signatures ---------------------------
longlat2UTMzone <- function(long, lat) {
  as.numeric((floor((long + 180) / 6) %% 60) + 1)
}

longlat2UTM <- function(long, lat, utm_zone = NA) {
  mean_long <- mean(long)
  if (is.na(utm_zone)) utm_zone <- longlat2UTMzone(mean_long, mean(lat))
  xy <- .ll2utm_xy(long, lat, utm_zone)
  df <- data.frame(x = xy$x, y = xy$y)
  list(df, as.numeric(utm_zone))
}

UTM2longlat <- function(x_utm, y_utm, utm_zone) {
  ll <- .utm2ll(x_utm, y_utm, utm_zone)
  data.frame(Longitude = ll$lon, Latitude = ll$lat)
}

message("utm_shim.R loaded: longlat2UTM / UTM2longlat overridden (pure-R TM).")
