# TraitsWinnersLosers

## Description
This repository contains data and code for the manuscript "Plant traits poorly predict winner and loser shrub species in a warming tundra biome".

## Authors
Mariana García Criado, Isla H. Myers-Smith, Anne D. Bjorkman, Signe Normand, Anne Blach-Overgaard, Haydn J.D. Thomas, Anu Eskelinen, Konsta Happonen, Juha M. Alatalo, Alba Anadon-Rosell, Isabelle Aubin, Mariska te Beest, Katlyn R. Betway-May, Daan Blok, Allan Buras, Bruno E.L. Cerabolini, Katherine Christie, J. Hans C. Cornelissen, Bruce C. Forbes, Esther R. Frei, Paul Grogan, Luise Hermanutz, Robert D. Hollister, James Hudson, Maitane Iturrate-Garcia, Elina Kaarlejärvi, Michael Kleyer, Laurent J. Lamarque, Jonas J. Lembrechts, Esther Lévesque, Miska Luoto, Petr Macek, Jeremy L. May, Janet S. Prevéy, Gabriela Schaepman-Strub, Serge N. Sheremetiev, Laura Siegwart Collier, Nadejda A. Soudzilovskaia, Andrew Trant, Susanna E. Venn and Anna-Maria Virkkala

Contact: Mariana García Criado, mariana.garcia.criado@gmail.com

## Data use guidelines
Data output from this manuscript are publicly available using a Creative Commons Attribution 4.0 International copyright (CC BY 4.0; see license.txt). Data are fully public but should be appropriately referenced by citing the paper. Although not mandatory, we additionally suggest that data users contact and collaborate with data contributors if this dataset will contribute a substantial proportion of observations used in a particular paper or analysis. DOI for this dataset is xxx

## Data availability & access
The data and code for this manuscript will be mantained at this GitHub repository (https://github.com/marianagarciacriado/TraitsWinnersLosers). 

## Citation
García Criado, M., Myers-Smith, I.H., Bjorkman, A.D., Normand, S., Blach-Overgaard, A., Thomas, H.J.D., Eskelinen, A., Happonen, K., Alatalo, J.M., Anadon-Rosell, A., Aubin, I., te Beest, M., Betway-May, K.R., Blok, D., Buras, A., Cerabolini, B.E.L., Christie, K., Cornelissen, J.H.C., Forbes, B.C., Frei, E.R., Grogan, P., Hermanutz, L., Hollister, R.D., Hudson, J., Iturrate-Garcia, M., Kaarlejärvi, E., Kleyer, M., Lamarque, L.J., Lembrechts, J.J., Lévesque, E., Luoto, M., Macek, P., May, J.L., Prevéy, J.S., Schaepman-Strub, G., Sheremetiev, S.N., Siegwart Collier, L., Soudzilovskaia, N.A., Trant, A., Venn, S.E., and Virkkala, A-M. (2023). Plant traits poorly predict winner and loser shrub species in a warming tundra biome. <i>Nature Communications</i>.

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
R version 4.2.0. or greater.

R packages: `AMR, brms, broom, cowplot, data.table, dplyr, ggbpubr, ggplot2, ggOceanMaps, ggrepel, modelr, plyr, reshape2, rjags, rstan, RVAideMemoire, R2jags, stargazer, Taxonstand, tidyverse, vegan`

