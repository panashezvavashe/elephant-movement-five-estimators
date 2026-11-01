# Data

The elephant GPS locations are **not** distributed with this repository because of poaching risk. They are held by the Zimbabwe Parks and Wildlife Management Authority (ZPWMA) and are available on reasonable request to ZPWMA (Scientific Services Unit, Marongora) and the corresponding author, subject to a data-sharing agreement.

- `raw/`: put the Africa Wildlife Tracking (AWT) exports here, one file per collar, named `Elephant <ID>.csv`.
- `processed/`: `code/01_screening/build_screened_v6.R` writes `all_elephants_screened_v6.csv` here. Every later script reads that file.

Both folders are excluded by `.gitignore`, so location data cannot be committed by accident.
