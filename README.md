
# JetStream_WarmingHole

This repository contains the Matlab code and analysis scripts used in the manuscript:

**A Wavier Polar Jet Stream Contributed to the Mid-20th Century Winter Warming Hole in the United States**  
*AGU Advances*

**Note:** Raw reanalysis data files are not included due to size restrictions. These are publicly available in online databases.


## Repository explanation

```
/Data/        # Processes reanalysis and other datasets used in figure generation
/Functions/   # Matlab functions used in analysis and figure generation (written for this study, if not otherwise attributed)
/SOMs/        # Contains .mat file with master SOM

/MatFig_WavinessClimatology.m   # Plots LWA and MCI waviness climatologies over US (Fig 1 in paper)
/MatFig_InsideOutside.m         # Plots US temperature inside/outside WH with map (Fig 2 in paper)
/MatFig_SOMcomposite.m          # Plots composite height (Fig 3) and winds (Fig 4) for SOM classes
/MatFig_Waviness.m              # Plots LWA, MCI, SOM Troughing Index, and U-wind in Warming Hole area (Fig 5 in paper)
/MatFig_Impact_Periods.m        # Plots warming hole impact across the 2 periods (Fig 6 in paper)
/MatFig_Impact_Time.m           # Plots annual warming hole impact across time (Fig 7 in paper)
/MatFig_GlobalWaviness.m        # Plots LWA and MCI waviness across Northern Hemisphere (Fig 8 in paper)

/SUPP_BerkeleyComparison.m       # Plots temperature and explores bias-correction for Berkeley Earth datasets (Supplementary Figs 1-2)
/SUPP_Impact_Time_Livneh.m       # temperature decomposition using Livneh instead of Berkeley Earth (Supplementary Fig 3)
/SUPP_WavinessIdentification.m   # Scatterplot of LWA and MCI for each SOM class (Supplementary Fig 4)
/SUPP_WHPhases.m                 # Warming Hole temperature with phases of Warming Hole distinctly marked (Supplementary Fig 5)
/SUPP_CEDS.m                     # Plots US SO2 emissions (Supplementary Fig 6)
/SUPP_WavinessPrecedent.m        # Timeseries of LWA and MCI over 20th century (Supplementary Fig 7)

** NEED RAW DATA BELOW **
/ReadRean_InsideOutside.m    # Reads temperature data from various datasets
/ReadRean_LWA.m              # Reads 500 hPa geopotential height data and calculates LWA from various datasets
/ReadRean_MCI.m              # Reads 300 hPa winds data and calculates MCI from various datasets
/ReadRean_LWAclimatology.m   # Calculates LWA climatology from NCEP/NCAR data
/ReadRean_MCIclimatology.m   # Calculates MCI climatology from NCEP/NCAR data
/ReadRean_ZonalWinds.m       # Reads 300 hPa winds data from various datasets
/ReadUSHCN.m                 # Reads USHCN temperature data
```
## Citation

If you use this code, please cite:

> Chalif, J.I., et al. (2025). *A Wavier Polar Jet Stream Contributed to the Mid-20th Century Winter Warming Hole in the United States*. AGU Advances. DOI: [to be added]

## Contact
Contact Jacob Chalif (jacob.i.chalif.gr@dartmouth.edu) with any questions.

---
