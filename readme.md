# Impact of El Niño Southern Oscillation on dynamics of major infectious diseases in China

Core code and data for "Impact of El Niño Southern Oscillation on dynamics of major infectious diseases in China"

### Organization

The repository is organized into **Scripts/**,  and **Data/** folders.

**Data/**: This folder includes intermediate and processed summary data that enable replication of most the figures and numbers cited in the text. Some of the files for the climate model projections, raw observational data, and future risk estimates are quite large, so they are not provided here. But the download URL has been provided.

**Scripts/**: Code required to reproduce the findings of our work is included in this folder. Take scarlet fever as an example

### Data

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
* **Data/Future/Future risk/** folder: This file contains the results of the projected cumulative case numbers from 2030s to 2092s for infectious diseases under different climate change scenarios (SSP1-2.6, SSP2-4.5, SSP3-7.0, and SSP5-8.5) and the results of  the impact of  global warming on the predicted  infectious dieases cases in China. But future estimates are quite large, so we only provided results of Scarlet fever in figshare (Link please see Data availability in our manuscript).

### Scripts

* The calculation of E-index, C-index  please  refer to code of  "Increased variability of eastern Pacific  El Niño under greenhouse warming"
* The calculation of teleconnection please  refer to code of "Persistent effect of El Niño on global economic growth"
* **DLM.r**: performs the main distributed lag models analysis, which are performed  with 1000 bootstraps. The model result are saved in "Data/Cumulative_effect/","Data/significance/" and "Data/Coefficients/"
* **Calculate_cumulative_effect.ipynb** and **Cumulative_effect.ipynb**: calculates and plot current cumulative effect of ENSO on infectious dieaseas. The current cumulative effect of ENSO on all infectious dieaseas in our study are shown in "Data/Cumulative_effect/".
* **Future cumulative cases.ipynb:** predictes future 10-year cumulative cases from 2030s to 2080s for selected infectious diseases under different
  climate change scenarios (SSP1-2.6, SSP2-4.5, SSP3-7.0, and SSP5-8.5).
* **ENSO_change_risk.ipynb:** calculates the impacts of warming-driven ENSO changes on infectious dieases.
