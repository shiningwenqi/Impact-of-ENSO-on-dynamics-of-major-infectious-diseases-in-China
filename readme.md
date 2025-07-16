# Impact of El Niño Southern Oscillation on dynamics of major infectious diseases in China

Core code and data for "Impact of El Niño Southern Oscillation on dynamics of major infectious diseases in China".

## System Requirements

#### Hardware Requirements

The package requires only a standard computer with enough RAM to support the operations defined by a user. For minimal performance, this will be a computer with about 2 GB of RAM. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB
CPU: 4+ cores, 3.3+ GHz/core

#### Software requirements and Installation guide

 All code were performed using R (v 4.3.2) or python (v 3.8.18). The code is  tested on *Windows* operating systems. Our code  mainly depends on the following

* **python dependencies**: pandas, numpy, matplotlib, os, cartopy, geopandas, descartes,  Affine, seaborn, and xarray. You can **install them from PyPi** , for example:

  ```
  pip3 install pandas
  ```

  It takes about an hour and a half to install all the above python dependencies.
* **R package dependencies**: ggplot2, margins, sandwich, stargazer, lfe, dplyr, lemon,  lspline, fixest, cowplot, texreg and magrittr. You can install them with the following sentence:

  ```
  install.packages(c('ggplot2', 'margins', 'sandwich', 'stargazer', 'lfe,', 'dplyr', 'lemon', 'lspline', 'fixest', 'cowplot', 'texreg', 'magrittr'))
  ```

  which will install in about 3 minutes on a machine.

## Organization

The repository is organized into **Scripts/**,  and **Data/** folders.

**Data/**: This folder includes intermediate and processed summary data that enable replication of most the figures and numbers cited in the text. However, disease incidence/case data in the panel (log_rate column) is simulated dataset used to test the code, as it not publicly available and are protected by data privacy laws. Some of the files for the climate model projections, raw observational data, and future risk estimates are quite large, so they are not provided here. But the download URL has been provided.

**Scripts/**: Code required to reproduce the findings of our work is included in this folder. Take scarlet fever as an example.

## Data

Raw climate datasets not available in this repository are publicly available as follows, details can also see in Supplementary Table 5:

* The **HadISST** sea surface temperature data are available from HadISST v1.1, https://www.metoffice.gov.uk/hadobs/hadisst/
* The **Berkeley Earth** surface temperature data are available https://berkeleyearth.org/data/
* The **precipitation data** are available https://crudata.uea.ac.uk/cru/data/hrg/
* **CMIP6 temperature, precipitation, and SST** data are generally available from the https://esgf-node.llnl.gov/projects/cmip6/;
* **Historical Population** are available from Gridded
  Population of the World, https://sedac.ciesin.columbia.edu/data/collection/gpw-v4 ; **Projection population** are available from NASA Socioeconomic Data and Applications Center:https://sedac.ciesin.columbia.edu/data/set/popdynamics-1-km-downscaled-pop-base-year-projection-ssp-2000-2100-rev01
* **Climate zoning data** are available from Data Center for Resources and Environmental Sciences, Chinese Academy of Science (RESDC: http://www.resdc.cn)
* **COVID-19: Stringency Index** are available from Oxford
  Coronavirus Government Response Tracker, https://ourworldindata.org/covid-stringency-index.
* **Data/Cumulative_effect/** folder: This file contains the results of the  cumulative effects of ENSO on all infectious diseases in our study.
* **Data/Future/Future risk/** folder: This file contains the results of the projected cumulative case numbers from 2030s to 2092s for infectious diseases under different climate change scenarios (SSP1-2.6, SSP2-4.5, SSP3-7.0, and SSP5-8.5) and the results of  the impact of  global warming on the predicted  infectious dieases cases in China. But future estimates are quite large, so we only provided results of Scarlet fever.

## Scripts and expected output and run time

* The calculation of E-index, C-index  please  refer to code of  "Increased variability of eastern Pacific  El Niño under greenhouse warming"
* The calculation of teleconnection please  refer to code of "Persistent effect of El Niño on global economic growth"
* **DLM.r**: performs the main distributed lag models analysis, which are performed  with 1000 bootstraps. The model result are saved in "Data/Cumulative_effect/","Data/significance/" and "Data/Coefficients/".  It can take about one hour to run in full.
* **Calculate_cumulative_effect.ipynb** and **Cumulative_effect.ipynb**: calculates and plot current cumulative effect of ENSO on infectious dieaseas. The current cumulative effect of ENSO on all infectious dieaseas in our study are shown in "Data/Cumulative_effect/". It can take about one minutes to run.
* **Future cumulative cases.ipynb:** predictes future 10-year cumulative cases from 2030s to 2080s for selected infectious diseases under different
  climate change scenarios (SSP1-2.6, SSP2-4.5, SSP3-7.0, and SSP5-8.5). The ipynb file records the expected    run time.
* **ENSO_change_risk.ipynb:** calculates the impacts of warming-driven ENSO changes on infectious dieases.  It can take hours  to run in full.
