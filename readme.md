# Automated Calculation & Visualization of % Inhibition 
## This method is applicable to 96-well plate derived raw data having triplicate samples


This pipeline is based on my usual routine of analyzing 96-well based data from antimicrobial assay.
The script is written in Python3 using different packages namely Python3, NumPy, matplotlib, pandas, os, math, and re.
The script is designed to be automated, such that user will only input few details of their filenames and the location (96-well) of their growth and media controls. This also allows users to have visualized heatmaps, bargraphs, summarized mean % inhibition, and list of active samples having high antimicrobial activity. 

##USAGE: python inhibition.py <metadata.csv> <filename1> <filename2> <filename...n> <media control wells> <growth control wells>

##USAGE example: python inhibition.py metadata.csv plate1.csv plate2.csv plate3.csv plate4.csv H10,H11,H12 H7,H8,H9

##Files needed
- [inhibition.py](https://github.com/JakiPeppero/Automated-Calculation-of-Bacterial-Inhibition---based-on-microdilution-assay-on-96-well-plate/blob/main/inhibition.py)
- metadata.csv - example is the followingg but users can opt to add more details if needed [metadata.csv](https://github.com/JakiPeppero/Automated-Calculation-of-Bacterial-Inhibition---based-on-microdilution-assay-on-96-well-plate/blob/main/metadata.csv)
- Rawdata - examples are
  - [plate1.csv](https://github.com/JakiPeppero/Automated-Calculation-of-Bacterial-Inhibition---based-on-microdilution-assay-on-96-well-plate/blob/main/plate1.csv) 
  - [plate2.csv](https://github.com/JakiPeppero/Automated-Calculation-of-Bacterial-Inhibition---based-on-microdilution-assay-on-96-well-plate/blob/main/plate2.csv) 
  - [plate3.csv](https://github.com/JakiPeppero/Automated-Calculation-of-Bacterial-Inhibition---based-on-microdilution-assay-on-96-well-plate/blob/main/plate3.csv) 
  - [plate4.csv](https://github.com/JakiPeppero/Automated-Calculation-of-Bacterial-Inhibition---based-on-microdilution-assay-on-96-well-plate/blob/main/plate4.csv) 
  - NOTE: Plate 4 has errors on them (overflow, missing/ blank data, negative values) so this can be included as a test for errors in future runs
  
##Expected Output files
- combined_summary
- active_hits.csv - summarized the active samples having high antimicrobial activity. Users can change the filtering for "active hits" but the default setting here are those having > 60% inhibition activity
- error_output.csv - only if errors are detected 
- QC_facet_bars.png - bar graph
- FACET_heatmaps.png - tiled heatmap per plate, showing all % inhibition per well
- COMBINED_heatmap.png - one heatmap summarizing all the %inhibition for all the plates having same sample_codes

Note: This project aside from helping my bioassay computation, is also a part of the prerequisite of BIO_539.