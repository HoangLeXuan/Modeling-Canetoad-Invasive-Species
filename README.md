# README for Modeling Invasive Species Dynamics and Distribution Repository
This repository contains the research and code used in the study "Modeling Invasive Species Dynamics and Distribution to Guide Management Strategies: A Case Study on Cane Toads in Australia" by[ Hoang Le](https://github.com/HoangLeXuan). The research aims to comprehensively assess potential cane toad distribution in Australia using statistical and machine learning models to inform management strategies and combat the spread of this invasive species. 

### About the Author
[Hoang. X. Le](https://github.com/HoangLeXuan)
[le_h2@Denison.edu](le_h2@Denison.edu)

### Research Overview
The study employs several modeling approaches to predict the presence and spread of Cane Toads in Australia, including Generalized Linear Models (GLM), BIOCLIM, Domain, Support Vector Machine (SVM), MaxEnt, and Random Forest.

The repository includes the final paper that provides a detailed account of the methodology, the data analyzed, model evaluations, and the results, including a geospatial prediction map. The findings are significant for biodiversity conservation and the development of targeted management strategies to mitigate the impact of Cane Toads.

### Purpose
Fun of learning

### Data
data.zip and input.zip contain data about the species occurrence, and Raster Layers are >100mb. Currently, I am figuring out a way to push those files without extensive problems.
For Terrestrial Ecoregions of the World data, refer to the package's documentation and download the full documentation at [here](https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world)

### Code

1. AbsenceBGData.Rmd: R Markdown file pulling and analyzing absence and background data of invasive species in various locations.
2. DataPipeline_ctsel.Rmd: R Markdown script creating a data pipeline for Cane Toad (Rhinella marina) under the ctsel dataset.
3. DifferentModelMethods.Rmd: An R Markdown document exploring different algorithms and machine learning models to study the dynamics and spread of Cane Toads in Australia.
4. ModelFitting.Rmd: R Markdown file where models are fitted to the data to predict the presence and spread of Cane Toads.
5. Figures: Folder containing visual outputs from the analyses and models.
6. DA401 - Autobiography - Hoang Le.pdf: Autobiographical document detailing the study's context and the researcher's background.
7. README.md: Markdown document (this file) describing the repository's content and purpose.
8. LICENSE.txt: The license under which this repository's content is distributed. We work under a BSD-3-Clause license.
9. .Rproj, .Rhistory, .gitignore: Files associated with the R project configuration, history, and version control.

### How to Use This Repository
Researchers and practitioners in invasive species management and ecological modeling can refer to the R Markdown files for methodology and code for implementing similar models in their work. The Figures folder contains visual aids that are used in the final paper of the research and can be used for presentations or further analysis.

### Citation
To cite this repository in your work, please use the following:

Le, H. X. (2024). Modeling the Dynamics of Invasive Species with Distribution Models to Inform Management Strategies: A Case Study on Cane Toads. GitHub.
https://github.com/HoangLeXuan/modeling-invasive-species

### Acknowledgements
The author thanks the contributors to the WorldClim database and The Global Biodiversity Information Facility (GBIF) and all referenced works that have made this study possible. The author would love to extend thanks to Dr. Supp for their guidance and Lena Le at [here](https://github.com/lenaledenison) for their extended reviews of the research, without their the research wouldn't be completed. 

### Contact
For any inquiries regarding the study or collaborations, please contact [le_h2@denison.edu] provided in the paper.
