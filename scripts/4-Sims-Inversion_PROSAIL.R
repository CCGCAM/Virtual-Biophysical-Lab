# ----------------------------------------------------------------------------
# PROSAIL model and Inversions
# ----------------------------------------------------------------------------
# OPTIONAL: clear environment (comment out if you don't like this behaviour)
rm(list = ls())

library(SCOPEinR)
library(ToolsRTM)

if (!require("dplyr")) { install.packages("dplyr"); require("dplyr") }  ###
if (!require("readr")) { install.packages("readr"); require("readr") }  ###
if (!require("ggplot2")) { install.packages("ggplot2"); require("ggplot2") }  ###
if (!require("plotly")) { install.packages("plotly"); require("plotly") }  ###
if (!require("tidyr")) { install.packages("tidyr"); require("tidyr") }  ###
if (!require("fmsb")) { install.packages("fmsb"); require("fmsb") }  ###
if (!require("doParallel")) { install.packages("doParallel"); require("doParallel") }  ###
if (!require("parallel")) { install.packages("parallel"); require("parallel") }  ###
if (!require("caret")) { install.packages("caret"); require("caret") }  ###
if (!require("randomForest")) { install.packages("randomForest"); require("randomForest") }  ###
if (!require("pls")) { install.packages("pls"); require("pls") }  ###
if (!require("e1071")) { install.packages("e1071"); require("e1071") }  ###
if (!require("gbm")) { install.packages("gbm"); require("gbm") }  ###


source('codes/compute_stats.R')


# ----------------------------------------------------------------------------
# 1. Create the LUT
# ----------------------------------------------------------------------------
# Here we assume you already generated a LUT (Leaf/Canopy parameters).
# The goal of this block is to:
#   • build a soil reflectance spectrum
#   • simulate canopy reflectance (BRF) using PROSAIL (FourSAIL)
#   • extract per-wavelength mean and standard deviation across all runs
# ----------------------------------------------------------------------------


#Get default PROSAIL input list
inputs <- ToolsRTM::inputsPROSAIL

# (Everything else stays as in defaults; you can tweak similarly if needed)

# 2.2 Generate a LUT with 2000 samples
set.seed(1234)
n.samples = 2000
LUT <- as.data.frame(ToolsRTM::getLUT(inputs = inputs, nLUT = n.samples, setseed = 1234))

# Uniform pigment distributions
LUT$Car = stats::runif(n.samples,min = 0,max=25)
LUT$Anth = stats::runif(n.samples,min = 0,max=5)

LUT$EWT <- stats::rnorm(n.samples, mean = 0.02, sd = 0.001)



# Choose traits to inspect (students can modify this vector)
traits_to_plot <- c("Cab", "Car", "Anth", "LAI", 'EWT','N' )

# Keep only traits that actually exist in the LUT
traits_to_plot <- traits_to_plot[traits_to_plot %in% names(lut)]

if (length(traits_to_plot) == 0) {
  stop("None of the selected traits are present in inputLUT.csv")
}

# Long format for ggplot
lut_long <- LUT |>
  select(all_of(traits_to_plot)) |>
  pivot_longer(cols = everything(),
               names_to = "Trait",
               values_to = "Value")

# Histogram (barplot) per trait
ggplot(lut_long, aes(x = Value)) +
  geom_histogram(fill = "grey70", color = "black", bins = 20) +
  facet_wrap(~ Trait, scales = "free") +
  theme_bw() +
  xlab("Value") +
  ylab("Frequency") +
  ggtitle("Histogram of LUT input traits (SCOPE simulations)")



# Summarise min–max per parameter
lut_ranges <- data.frame(
  Parameter = names(LUT),
  Min = round(sapply(LUT, min, na.rm = TRUE),3),
  Max = round(sapply(LUT, max, na.rm = TRUE),3)
)

# Display in an interactive datatable
DT::datatable(
  lut_ranges,
  caption = "Ranges of sampled input parameters in the LUT",
  options = list(pageLength = 10)
)


# ---------------------------------------------------------
# 2. Build soil reflectance profile (dry–wet mixture)
# ---------------------------------------------------------
# SCOPE/PROSAIL requires a soil reflectance curve (400–2500 nm).
# ToolsRTM already provides two standard soil spectra:
#   - dry soil
#   - wet soil
# Here we mix them linearly using parameter psoil (0 = wet, 1 = dry).
# This allows controlling soil moisture optically.

wl_grid <- c(400:2500)   # wavelength grid required by FourSAIL (nm)

rsoil.dry <- ToolsRTM::dataSpec_PDB$dry_soil    # dry soil reflectance (400–2500 nm)
rsoil.wet <- ToolsRTM::dataSpec_PDB$wet_soil    # wet soil reflectance (lower reflectance)

psoil <- 0.5   # 50% dry + 50% wet

# Linear interpolation between dry and wet soil reflectance
rsoil.default <- psoil * rsoil.dry + (1 - psoil) * rsoil.wet

# Build a tidy dataframe for ggplot
df.soil <- data.frame(
  Wavelength = wl_grid,
  Reflectance = rsoil.default
)


# ---------------------------------------------------------
# 3. Plot soil reflectance spectrum
# ---------------------------------------------------------
# Students should notice:
#   • Dry soil is brighter (higher reflectance)
#   • Wet soil absorbs more in NIR/SWIR
#   • Soil optical properties affect canopy reflectance when LAI is low

ggplot(df.soil, aes(x = Wavelength, y = Reflectance)) +
  geom_line(color = "brown4", size = 1) +
  theme_bw() +
  labs(
    title = "Soil Reflectance Spectrum",
    x = "Wavelength (nm)",
    y = "Reflectance"
  ) +
  ylim(0, max(df.soil$Reflectance) * 1.1)



# ---------------------------------------------------------
# 4. Simulate canopy reflectance using PROSAIL (FourSAIL)
# ---------------------------------------------------------
# For each row in the LUT:
#   1) Run PROSAIL-D through ToolsRTM::foursail()
#   2) Combine all directional + diffuse components into BRF
#   3) Return a tibble with run index, wavelength, and BRF value
# ----------------------------------------------------------------------------

sims <- lapply(seq_len(nrow(LUT)), function(i) {

  # FourSAIL simulation for one LUT configuration
  out <- suppressMessages(
    ToolsRTM::foursail(
      inputLUT     = LUT[i, ],
      rsoil        = rsoil.default,
      LeafModel    = "PROSPECT-D",    # biochemical leaf model
      spectrum.all = TRUE               # output full spectral components
    )
  )

  # Compute BRF (Bidirectional Reflectance Factor)
  # BRF combines directional/diffuse reflectance according to solar zenith angle (tts)
  brf <- ToolsRTM::Compute_BRF(
    rdot = out$rdot,                    # diffuse → observer
    rsot = out$rsot,                    # direct sun → observer
    tts  = LUT[i, "tts"],               # solar zenith angle
    data.light = ToolsRTM::dataSpec_PDB
  )

  # Return tidy tibble for this simulation
  tibble(
    run = i,
    wl  = wl_grid,
    rho = as.numeric(brf)
  )
})

# ---------------------------------------------------------
# 5. Combine all runs + compute spectral statistics
# ---------------------------------------------------------
# After simulating all canopy reflectance curves:
#   • merge into a single table
#   • compute mean + standard deviation per wavelength
# ---------------------------------------------------------

brf_runs <- dplyr::bind_rows(sims)

spec_stats <- brf_runs %>%
  group_by(wl) %>%
  summarise(
    mean_rho = mean(rho, na.rm = TRUE),
    sd_rho   = sd(rho, na.rm = TRUE),
    .groups  = "drop"
  )

# ----------------------------------------------------------------------------
# 6.Get LUT and fet indices
# ----------------------------------------------------------------------------

brf_wide <- bind_rows(sims) %>%
  mutate(wl = paste0("R.", wl)) %>%        # column names like R.400, R.405, ...
  tidyr::pivot_wider(
    id_cols     = run,                     # one row per simulation
    names_from  = wl,                      # each wavelength becomes a column
    values_from = rho
  ) %>%
  arrange(run)

# Remove 'run' if LUT has no explicit ID column
refl_mat <- brf_wide %>%
  dplyr::select(-run)

# Check dimensions: should match nrow(LUT)
stopifnot(nrow(LUT) == nrow(refl_mat))

# Build final LUT + reflectance table
LUT_refl <- cbind(LUT, refl_mat)
# Identificar columnas espectrales
spec_cols <- grep("^R\\.", names(LUT_refl), value = TRUE)

df.indices <- ToolsRTM::getIndices(LUT_refl[,spec_cols],  pattern.rfl      = "R.",
  spectral.domain  = "VNIR")
df.indices <- df.indices[, colSums(is.na(df.indices)) == 0]
#nzv <- caret::nearZeroVar(df.indices)
#df.indices <- df.indices[, -nzv]
head(df.indices)
# ----------------------------------------------------------------------------
# 7. Visualise all PROSAIL reflectance simulations + mean ± SD envelope
# ----------------------------------------------------------------------------
# This plot shows:
#   • All individual BRF spectra (one per LUT row) in faint grey lines
#   • The mean reflectance spectrum across all simulations (blue line)
#   • A shaded band showing ±1 standard deviation (blue ribbon)
# ----------------------------------------------------------------------------
p <- ggplot() +
  geom_line(data = brf_runs, aes(x = wl, y = rho, group = run), alpha = 0.12) +
  geom_ribbon(data = spec_stats,
              aes(x = wl, ymin = mean_rho - sd_rho, ymax = mean_rho + sd_rho),
              inherit.aes = FALSE, alpha = 0.2, fill = "blue") +
  geom_line(data = spec_stats, aes(x = wl, y = mean_rho),
            inherit.aes = FALSE, color = "blue") +
  labs(title = "PROSAIL (custom LUT): mean ± 1 sd (BRDF)",
       x = "Wavelength (nm)", y = "Reflectance (BRDF)") +
  theme_bw()
print(p)



# ----------------------------------------------------------------------------
# 7. Define target trait and predictor variables
# ----------------------------------------------------------------------------

target_trait <- "Cab"   # <-- choose the trait you want to retrieve

#Target wavelengths every 5 nm from 400 to 800
target_wavelengths <- paste('R.',seq(400, 800, by = 5),sep='')


# Combine all predictors you want to use
predictors <- c(target_wavelengths)

# Build the modelling dataset: 1 column for Cab + predictors
dataset <- cbind(
  Cab = LUT_refl[[target_trait]],
  LUT_refl[, predictors],
  df.indices
)

# Remove any duplicated columns (can happen if names overlap)
dataset <- dataset[, !duplicated(names(dataset))]

head(dataset)



df_spec <- dataset %>%
  select(starts_with("R.")) %>%       # keep only spectral columns
  mutate(sim = 1:n()) %>%             # simulation index
  pivot_longer(cols = starts_with("R."),
               names_to = "Wavelength",
               values_to = "Reflectance")

# Convert "R.400" → 400
df_spec$Wavelength <- as.numeric(gsub("R\\.", "", df_spec$Wavelength))

ggplot(df_spec, aes(x = Wavelength, y = Reflectance, group = sim)) +
  geom_line(alpha = 0.25, color = "darkblue") +
  theme_bw() +
  xlab("Wavelength (nm)") +
  ylab("Reflectance") +
  ggtitle("Reflectance spectra for all simulations")



# ----------------------------------------------------------------------------
# 8. Optional: remove highly collinear predictors using VIF
# ----------------------------------------------------------------------------

set.seed(42)
rows.vif <- sample(nrow(dataset), min(200, nrow(dataset)))  # subset for speed

vif_keep <- ToolsRTM::getVIF(
  dataset[rows.vif, -1],   # all predictors (exclude first column = Cab)
  thresh = 10              # VIF threshold
)
print(vif_keep)
selected_keep <-  c("MCARI","SRPI","DNCabxc",
                    'Chlred.edge','PSSRb','NDVI',
                    "IRECI", "PRI515","PRIn","PRI_CI","BF1","BF5",
                    "CR.red.nir.2", "CR.red.nir.4", "CR.red.nir.5")



# Keep only predictors selected by VIF
dataset_vif <- cbind(Cab = dataset$Cab, dataset[, vif_keep])




# ----------------------------------------------------------------------------
# 9. Inversion using GB &RF
# ----------------------------------------------------------------------------



model.rf=get.inversion(data = dataset_vif, depVar = "Cab", inputs = vif_keep,
                        algorithm  = "RF", n.samples = 500,seed = 123)

rf_fit <- model.rf$model   # caret::train object




model.gb=get.inversion(data = dataset_vif, depVar = "Cab", inputs = vif_keep, algorithm  = "GB",
                       n.samples = 300,seed = 123)

gb_fit <- model.gb$model   # caret::train object






# ----------------------------------------------------------------------------
# 10. Predictions in a real scenario
# ----------------------------------------------------------------------------


rfl.sensor <- read.csv("Tables/FieldData/reflectance_plot.csv")
head(rfl.sensor)


field.maize <- read.csv("Tables/FieldData/Maize.csv")
field.potato <- read.csv("Tables/FieldData/Potato.csv")
field.data <- rbind(field.maize,field.potato)
colnames(field.data) <- c('PlotNum', colnames(field.data)[2:10])
head(field.data)





# rfl.sensor has: Wavelengths, Reflectance, Plot, PlotNum, Crop, Date
# Step 1: round wavelengths to integer nm and build column names like "R.400"
rfl.sensor2 <- rfl.sensor %>%
  mutate(
    wl_round = round(Wavelengths),             # 398.017 → 398, etc.
    band_name = paste0("R.", wl_round)
  )

# Step 2: pivot to wide format: one row per Plot/Date, one column per wavelength
rfl.wide <- rfl.sensor2 %>%
  select(Plot, PlotNum, Crop, Date, band_name, Reflectance) %>%
  group_by(Plot, PlotNum, Crop, Date, band_name) %>%
  summarise(Reflectance = mean(Reflectance, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = band_name,
    values_from = Reflectance
  )

# Identify reflectance columns ("R.xxx")
rfl_cols <- grep("^R\\.", names(rfl.wide), value = TRUE)

# Extract numeric wavelength (400, 402, 405, ...)
wl_numeric <- as.numeric(gsub("R\\.", "", rfl_cols))

# Order columns by wavelength
sorted_cols <- rfl_cols[order(wl_numeric)]

# Create a reordered dataset
rfl.wide.sorted <- rfl.wide[, c("Plot", "PlotNum", "Crop", "Date", sorted_cols)]

# Check first 12 columns after metadata
head(rfl.wide.sorted[, 1:12])


# Unir reflectancias y datos de campo por PlotNum + Crop
field.spec <- field.data %>%
  left_join(rfl.wide.sorted,
            by = c("PlotNum", "Crop"),
            suffix = c(".spec", ".field"))
field.spec <- field.spec %>%
  mutate(
    Date.spec  = as.Date(Date.spec,  format = "%d/%m/%Y"),
    Date.field = as.Date(as.character(Date.field), format = "%Y%m%d")
  )# Identificar columnas espectrales
spec_cols <- grep("^R\\.", names(field.spec), value = TRUE)

head(field.spec)

df_field_spec <- field.spec %>%
  select(Plot = PlotNum, Crop, N_treatment, Date.spec, all_of(spec_cols)) %>%
  pivot_longer(cols = all_of(spec_cols),
               names_to  = "Wavelength",
               values_to = "Reflectance")

# Convertir "R.400" → 400
df_field_spec$Wavelength <- as.numeric(gsub("R\\.", "", df_field_spec$Wavelength))

# Asegurar orden lógico de niveles de N
df_field_spec$N_treatment <- factor(df_field_spec$N_treatment,
                                    levels = c("N0", "N80", "N120"))
df_field_spec_clean <- df_field_spec %>%
  filter(!is.na(N_treatment))

ggplot(df_field_spec,
       aes(x = Wavelength,
           y = Reflectance,
           group = Plot,
           color = N_treatment)) +
  geom_line(alpha = 0.5, linewidth = 0.7) +
  facet_wrap(~ Crop, scales = "free_y") +
  theme_bw() +
  xlab("Wavelength (nm)") +
  ylab("Reflectance") +
  ggtitle("Field Reflectance Spectra by Nitrogen Treatment") +
  scale_color_brewer(palette = "Set1",
                     name = "N treatment")

# ----------------------------------------------------------------------------
# 2.Get indices
# ----------------------------------------------------------------------------


df.indices.field <-ToolsRTM::getIndices(field.spec[,spec_cols], pattern.rfl = 'R.', spectral.domain = 'VNIR')
df.indices.field <- df.indices.field[, colSums(is.na(df.indices.field)) == 1]
head(df.indices.field)


# Build the modelling dataset: 1 column for Cab + predictors
dataset.field <- cbind(
  field.spec, df.indices.field
)
# Subset predictors
dataset.field <- na.omit(dataset.field)
dim(dataset.field)


# Extract the trained RF model
rf_fit <- model.rf$model

# Make sure all needed predictors exist in the field indices
all(vif_keep %in% colnames(dataset.field))
# If FALSE, check name mismatches / missing indices


dim(dataset.field)
# Predict Cab for each Plot/Date
dataset.field$Cab_pred.rf <- predict(rf_fit, dataset.field[, vif_keep])
dataset.field$Cab_pred.gb <- predict(gb_fit, dataset.field[, vif_keep])




ggplot(dataset.field, aes(x = CR.red.nir, y = Chl)) +
  geom_point(size = 3, alpha = 0.7, color = "darkgreen") +
  geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  labs(title = "CR.red.nir vs Chl",
       x = "CR.red.nir (NIR / Red)",
       y = "Chlorophyll (Chl)")


ggplot(dataset.field, aes(x = Cab_pred.rf, y = NBI)) +
  geom_point(size = 3, alpha = 0.7, color = "darkgreen") +
  geom_smooth(method = "lm", color = "black") +
  theme_bw() +
  labs(title = "CR.red.nir vs Chl",
       x = "Predicted Cab using RF",
       y = "measured Chlorophyll (Chl)")


ggplot(dataset.field, aes(x = Cab_pred.rf,y = Chl,color = N_treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  facet_wrap(~ Crop, scales = "free_y") +
  labs(
    x = "Predicted Cab (RF inversion)",
    y = "Measured Chlorophyll (Chl)",
    color = "N treatment"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

ggplot(dataset.field, aes(x = Cab_pred.rf,y = Chl,color = N_treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", color = "black") +
  facet_wrap(~ Crop, scales = "free_y") +
  theme_bw() +
  labs(
       x = "Predicted Cab using GB",
       y = "measured Chlorophyll (Chl)")

