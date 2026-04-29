# Automated Calculation & Visualization of % Bacterial Inhibition Assay
## This method is applicable to 96-well plate derived raw data having triplicate samples


This pipeline is based on my usual routine of analyzing 96-well based data from antimicrobial assay.
The script is written in Python3 using different packages namely Python3, NumPy, matplotlib, pandas, os, math, and re.
The script is designed to be automated, such that user will only input few details of their filenames and the location of their growth and media controls (based on the 96-well format). This also allows users to have visualized heatmaps, bargraphs, summarized mean % inhibition, and list of active samples having high antimicrobial activity. 

## USAGE:  
python inhibition.py metadata.csv rawfile1.csv rawfile2.csv ... rawfileN.csv media_control_wells growth_control_wells

## USAGE example: 
python inhibition.py metadata.csv plate1.csv plate2.csv plate3.csv plate4.csv H10,H11,H12 H7,H8,H9

## FILES needed
- [inhibition.py](https://github.com/JakiPeppero/Automated-Calculation-of-Bacterial-Inhibition---based-on-microdilution-assay-on-96-well-plate/blob/main/inhibition.py)
- metadata.csv - metadata example is included but users can opt to add more columns (i.e. test pathogens) and revise some codes if needed [metadata.csv](https://github.com/JakiPeppero/Automated-Calculation-of-Bacterial-Inhibition---based-on-microdilution-assay-on-96-well-plate/blob/main/metadata.csv)
- Rawdata - users can add as many as files as they want
  - [plate1.csv](https://github.com/JakiPeppero/Automated-Calculation-of-Bacterial-Inhibition---based-on-microdilution-assay-on-96-well-plate/blob/main/plate1.csv) 
  - [plate2.csv](https://github.com/JakiPeppero/Automated-Calculation-of-Bacterial-Inhibition---based-on-microdilution-assay-on-96-well-plate/blob/main/plate2.csv) 
  - [plate3.csv](https://github.com/JakiPeppero/Automated-Calculation-of-Bacterial-Inhibition---based-on-microdilution-assay-on-96-well-plate/blob/main/plate3.csv) 
  - [plate4.csv](https://github.com/JakiPeppero/Automated-Calculation-of-Bacterial-Inhibition---based-on-microdilution-assay-on-96-well-plate/blob/main/plate4.csv) NOTE: Plate 4 has errors on them (overflow, missing/ blank data, negative values) so this can be included as a test for errors in future runs
  
## Expected Output files
- combined_summary - shows the computed mean % inhibition with standard error
- active_hits.csv - summarized the active samples having high antimicrobial activity. Users can change the filtering for "active hits" but the default setting here are those having > 60% inhibition activity
- error_output.csv - only if errors are detected (see plate 4 above)
- QC_facet_bars.png - tiled bar graph per plate 
- FACET_heatmaps.png - tiled heatmap per plate, showing all % inhibition per well
- COMBINED_heatmap.png - one heatmap summarizing all the %inhibition for all the plates having same sample_codes

## Other potential use
- This code can also be used to compute % bacterial growth proliferation. You just need to change the computation part in inhibition.py 
- This is also applicable to other biological assays in 96-well format such as mammalian cell culture assays, protein quantification, etc. Few tweaks can be made to retrofit to each users prerequisite.

Note: This project aside from helping my bioassay computation, is also a part of the prerequisite of BIO_539. Thanks Dr. Schwartz!