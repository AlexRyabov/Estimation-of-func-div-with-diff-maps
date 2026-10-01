# DM4FuncDiv: diffusion maps for functional diversity

This repository contains the data and MATLAB source code for the analysis of the
Lithuanian phytoplankton monitoring data (south-eastern Baltic Sea, 1993–2015)
presented in

> Ryabov A. et al. (2022). Estimation of functional diversity and species traits
> from ecological monitoring data. *Proceedings of the National Academy of
> Sciences* 119(43), e2118156119.
> https://doi.org/10.1073/pnas.2118156119

The method uses a diffusion map to estimate species traits ("i-traits")
from how often species occur together in monitoring samples. Functional
diversity is then computed from these traits.
It works in four steps:

1. A similarity matrix between species is computed from the co-occurrence
   of their biomass in the samples (Spearman correlation).
2. The similarity matrix is sparsified (k nearest neighbours) and the
   eigenvectors of its random-walk Laplacian are computed. The eigenvectors
   scaled by 1/eigenvalue are the species i-traits.
3. Functional diversity of each sample is the Rao quadratic entropy, with
   the diffusion distance between species used as the trait distance.
4. The i-traits are compared with the biomass-weighted mean environmental
   conditions of each species, and trends of functional diversity are
   analysed for each station and season.

## Repository structure

```
DM4FuncDiv/
├── README.md
├── data/
│   ├── LT_phytoplankton_Cleaned_Plankton.mat           phytoplankton biomass (table tPlAbnd)
│   ├── Chem_Hidrology_TempSal_Avrg_1-10.mat            temperature and salinity, 1-10 m average (table tTempSalDataAvgd)
│   ├── Chem_Hidrology_Chem_Avrg_1-10m.mat              nutrients and pH, 1-10 m average (table tChemDataAvgd)
│   └── LT_phytoplankton_Cleaned_StationCoordinates.xlsx station coordinates
├── source/
│   ├── DM_CompareSpec.m         main script: run this
│   ├── DM_Calculations.m        per-sample relative biomass and diversity indexes
│   ├── DM_Adds_Calculations.m   figures of the paper
│   ├── f_DM_DMit.m              diffusion map (similarity, Laplacian, eigenvectors)
│   └── f_DiversityMetrics.m     functional diversity (Rao), richness, Simpson ESN
└── functions/                   general plotting and helper functions
```

## Usage

1. Open MATLAB and run `source/DM_CompareSpec.m`. The script can be started from any
   folder, because all paths are resolved relative to the script.
2. The parameters are set at the top of the script:
   - `DataID`: which stations to use. `'60vis'` (stations visited at least
     60 times) was used in the paper. Other options are `'20vis'`, `'Regs'`,
     `'2004on'` and `'All'`.
   - `PlotsID`: extra plots, e.g. `{'VisStations'}`.
   - `SaveFigures`: if `true`, the figures `Fig1.png`, `Fig3.png` and
     `FigSI_EnvITraits.png` are saved to the folder `figs/`, which is created
     automatically.

## Requirements

- MATLAB R2021a or newer. The code was developed with R2021. Some functions
  need at least R2020a: `boxchart`, `exportgraphics`, `tiledlayout` and
  `arguments` blocks.
- Statistics and Machine Learning Toolbox. It provides `grpstats`, `fitlm`,
  `fitnlm`, `corr`, `pdist`, `boxplot`, `quantile` and `nanmean`.
- Internet connection for the map basemaps used by `geobubble`.

The sign of eigenvectors is arbitrary. The script flips the sign of the 2nd and
3rd eigenvectors so that the figures have the same orientation as with
MATLAB R2021, which was used for the paper.

## Data

The data come from the Lithuanian phytoplankton and hydrochemical monitoring
in the south-eastern Baltic Sea, 1993–2015, provided by I. Olenina. See the
paper for details.

- `LT_phytoplankton_Cleaned_Plankton.mat`: species biomass per sample.
  Rows without a species name and "total" rows were removed from the original
  data.
- `Chem_Hidrology_*.mat`: temperature, salinity, pH and nutrient
  concentrations averaged over 1–10 m depth.
- `LT_phytoplankton_Cleaned_StationCoordinates.xlsx`: coordinates of the
  stations (WGS 84).

Each phytoplankton sample is matched with the environmental measurement made
at the same station within ±7 days.

## Citation

If you use this code or data, please cite the paper above.
