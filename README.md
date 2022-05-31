# ploidy_project

The programs in this directory were used to complete the study entitled "Global biodiversity data suggests allopolyploid plants do not occupy larger ranges or harsher conditions compared to their progenitors", which is in the process of being published. 

ploidy_project_loop.R
<br> This program calculates the environmental and geographical range and overlap, centroid, maximum latitude, and climatic variables (temp_max, temp_min, temp_breadth, precip_max, precip_min, precip_seasonality). This loop only works for species with 1 or 2 parents since I calculated the species with 3 parents separately. 

sdm_code.R
<br> This program creates maxent models using the species occurrence data and then saves them to my local computer. It also calculates the suitability measures, B1 and B2. 

chi-square.R
<br> This program performs chi square tests for geographical and environmental ranges, B1, B2, and the climatic variables. This is done with data including all of the species and again with only diploid parents.
