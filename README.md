# TraitsWinnersLosers
This repository contains data and code for the manuscript 'Plant traits poorly predict winner and loser shrub species in a warming tundra biome'

## Description
This repo contains data and code for the manuscript "Plant traits poorly predict winner and loser shrub species in a warming tundra biome".

## Authors
Mariana García Criado, Isla H. Myers-Smith, Anne D. Bjorkman, Signe Normand, Anne Blach-Overgaard, Haydn J.D. Thomas, Anu Eskelinen, Konsta Happonen, Juha M. Alatalo, Alba Anadon-Rosell, Isabelle Aubin, Mariska te Beest, Katlyn R. Betway-May, Daan Blok, Allan Buras, Bruno E.L. Cerabolini, Katherine Christie, J. Hans C. Cornelissen, Bruce C. Forbes, Esther R. Frei, Paul Grogan, Luise Hermanutz, Robert D. Hollister, James Hudson, Maitane Iturrate-Garcia, Elina Kaarlejärvi, Michael Kleyer, Laurent J. Lamarque, Jonas J. Lembrechts, Esther Lévesque, Miska Luoto, Petr Macek, Jeremy L. May, Janet S. Prevéy, Gabriela Schaepman-Strub, Serge N. Sheremetiev, Laura Siegwart Collier, Nadejda A. Soudzilovskaia, Andrew Trant, Susanna E. Venn and Anna-Maria Virkkala

Contact: Mariana García Criado, mariana.garcia.criado@gmail.com

## Data use guidelines
(this needs updating) Data are publicly available using a Creative Commons Attribution 4.0 International copyright (CC BY 4.0; see license.txt). Data are fully public but should be appropriately referenced by citing the paper. Although not mandatory, we additionally suggest that data users contact and collaborate with data contributors if this dataset will contribute a substantial proportion of observations used in a particular paper or analysis. DOI for this dataset is xxx

## Data availability & access
The data and code for this manuscript will be mantained at this GitHub repository (https://github.com/marianagarciacriado/TraitsWinnersLosers). 

## Citation
(needs updating)

## Data
Trait data are available at https://www.try-db.org/TryWeb/Home.php (TRY) and https://tundratraitteam.github.io/ (TTT). Cover change over time data will be published at https://github.com/annebj/ITEX30_VegComp. A previous version of this dataset can be accessed at http://polardata.ca/, CCIN Reference Number 10786. Species range data are available as a summarised dataset, together with the rest of input data necessary to reproduce figures and analyses, in the `data` folder.

Scripts 0_cover_change_over_time and 0_TRY_TTT_cleaning are included in this repository to show the process of generating input files, but the original raw files are not provided in this repository. When a raw data file is not included in this repository, this is specified in the script as "this file is not available in the repo as it contains raw data". All summarised data files are provided in this repository as input files and enable the reproducibility of the figures and analyses of this manuscript.

## Scripts
All the analyses undertaken for this manuscript are split between multiple R scripts within the `scripts`folder.
They can be followed in a sequential order (i.e., 1 to 8), with both 0_ scripts showing the process of generating summarised input files.

## Figures
The figures generated in R are stored in the `figures` folder.

## Model outputs
Full model outputs for statistical analyses are stored in the `models` folder.

## Software requirements
(needs updating) R version 3.4.3 or greater.

(needs updating) R packages: `tidyverse, proj4, scales, ggalt, ggplot2, dplyr, cowplot, MCMCglmm, ggbpubr, raster, rgdal, rasterVis, sp, tidyr, readr, broom, stargazer`
