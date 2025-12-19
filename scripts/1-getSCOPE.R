
# ============================================================================
# Running SCOPE v2.1 with SCOPEinR + ToolsRTM
# Soil Canopy Observation, Photochemistry and Energy balance model
# Example script for students
# Author: Carlos Camino (adapted for teaching)
# ============================================================================

rm(list = ls())   # Clear workspace

# ----------------------------------------------------------------------------
# 0. Load or install required packages
# ----------------------------------------------------------------------------

pkgs <- c(
  "parallel", "doParallel",
  "ToolsRTM", "SCOPEinR",
  "dplyr", "tidyr", "ggplot2"
)

for (p in pkgs) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}

# ------------------------------------------------------------------------------
# Step 1: Install & Load ToolsRTM from GitLab
# ------------------------------------------------------------------------------

if (!requireNamespace("ToolsRTM", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  # Public GitLab repo (no token needed). Set upgrade="never" for reproducibility.
  remotes::install_gitlab("caminoccg/toolsrtm", upgrade = "never")
}
library(ToolsRTM)

cat("\n ToolsRTM is ready: ", as.character(packageVersion("ToolsRTM")), "\n", sep = "")
#Check the version is 0.65

# ------------------------------------------------------------------------------
# Step 2: Install & Load SCOPE from GitLab
# ------------------------------------------------------------------------------
if (!requireNamespace("SCOPEinR", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  # Public GitLab repo (no token needed). Set upgrade="never" for reproducibility.
  remotes::install_gitlab("caminoccg/scopeinR", upgrade = "never")
}
library(SCOPEinR)

cat("\n SCOPEinR is ready: ", as.character(packageVersion("SCOPEinR")), "\n", sep = "")

#Check the version is 0.48


# ----------------------------------------------------------------------------
# 1. Load SCOPE model options
#    These options control radiative transfer, energy balance, and model numerics.
# ----------------------------------------------------------------------------

table.with.opts <- read.table(
  "Tables/inputs/setoptions.csv",
  header = TRUE, sep = ","
)

# ----------------------------------------------------------------------------
# 2. Build the LUT (Look-Up Table) with model input parameters
#    The LUT contains one row per SCOPE simulation.
# ----------------------------------------------------------------------------

N.Samples <- 100
file.LUT  <- "LUT"   # choose "default" or "LUT"

start_time <- Sys.time()

# ----------------------------------------------------------------------------
# 2A. Option A: Use predefined LUT (one simulation)
# ----------------------------------------------------------------------------
if (file.LUT == "default") {

  Table.LUT <- read.table(
    "Tables/inputs/LUT_input.csv",
    header = TRUE, sep = ","
  )

  db.sim <- SCOPEinR::get.SCOPE(
    LUT           = Table.LUT[1, ],
    options.SCOPE = table.with.opts,
    optipar       = SCOPEinR::optipar2021.Pro.CX,
    leaf.model    = "fluspect-CX",
    canopy.model  = "fourSAIL",
    get.outputs   = "ALL",
    get.plots     = FALSE
  )

  # ----------------------------------------------------------------------------
  # 2B. Option B: Generate a random LUT with N.Samples simulations
  # ----------------------------------------------------------------------------
} else {

  inputLUT <- read.table(
    "Tables/inputs/inputs_SCOPE.csv",
    header = TRUE, sep = ","
  )

  Table.LUT <- getLUT.SCOPE(
    inputLUT = inputLUT,
    nLUT     = N.Samples
  )

  db.sim <- get.SCOPE(
    LUT           = Table.LUT,
    n.LUT         = N.Samples,
    options.SCOPE = table.with.opts,
    optipar       = SCOPEinR::optipar2021.Pro.CX,
    leaf.model    = "fluspect-CX",
    canopy.model  = "fourSAIL",
    get.outputs   = "ALL",
    get.plots     = FALSE
  )
}

end_time <- Sys.time()

cat("\nExecution time:", end_time - start_time, "\n")

# ----------------------------------------------------------------------------
# 3. Extract and save SCOPE outputs
#    This saves reflectance, fluorescence, radiance, fluxes, and derived variables.
# ----------------------------------------------------------------------------

get.SCOPE.outputs(
  data.sim        = db.sim,
  N.sims          = N.Samples,
  LUT             = Table.LUT,
  get.outputs = 'ALL',
  path.out        = "outs/",
  get.more.inputs = c("refl", "lidf", "LIDFb", "Ft_Fo", "rdo"),
  get.plots       = TRUE
)
