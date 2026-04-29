#!/usr/bin/env python3
# usage: python inhibition.py <metadata.csv> <filename1> <filename2> <filename...n> <media control wells> <growth control wells>
# usage example: python inhibition.py metadata.csv plate1.csv plate2.csv plate3.csv plate4.csv H10,H11,H12 H7,H8,H9

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import re


# ERROR STORAGE
ALL_ERRORS = [] #stores possible errors under list
                # possible errors to be detected in  raw data: overflow readings,
                # missing values, negative values


# LOADING OF DATA VIA ARGUMENTS
if len(sys.argv) < 5:
                    #0-python script i.e. script.py
                    #1-filename1
                    #2-filename2
                    #3-metadata
                    #4-media control
                    #5-growth control
    sys.exit("Usage: python y12.py metadata.csv plate1.csv plate2.csv ... MEDIA_WELLS GROWTH_WELLS")

#This is similar to $1 in Shell scripting where we can incorporate the filename to our phython command in the terminal.
#But in this case, we added not just the filenames, but the metadata, and the wells for growth and media control
#The following commands identifies which argument is the metadata, media blanks, and growth
metadata_file = sys.argv[1]
media_arg = sys.argv[-2]
growth_arg = sys.argv[-1]

#This easily identifies our files so that even if we add more files on the middle, we dont have to edit the script. 
#Hence the file_list are all from index 2 until before the last two argument.
file_list = sys.argv[2:-2]

if len(file_list) == 0:
    sys.exit("ERROR: No plate files provided.") #just a loading error check if no files are provided


#LOADING METADATA
#-------------------------------------------------
#Metadata is needed to incorporate the sample codes in the setup. It also prevents hardcoding for users which may
#have different sample codes and sample concentrations(dose)
metadata = pd.read_csv(metadata_file)
metadata["row"] = metadata["row"].str.strip() #str.strip() cleans up whitespaces in row
metadata["col"] = metadata["col"].astype(int) #str.strip() cleans up whitespaces in column


#Converts the well location into usable string coordinates (a1,a2,a3...h12) --> [("a,1"), ("a,2"), ("a,3")]
def parse_wells(s):
    return [(w[0], int(w[1:])) for w in s.split(",")]

#get the media_blank and growth values based on the parsed values
MEDIA_BLANK = parse_wells(media_arg) 
GROWTH_CONTROL = parse_wells(growth_arg)


#LOADING OF RAWDATA FILES
#-------------------------------------------------


def load_plate(file_path):

    plate_name = os.path.splitext(os.path.basename(file_path))[0]

    # Read everything as string and store as dataframe
    df = pd.read_csv(file_path, index_col=0, dtype=str)

    df.index = df.index.str.strip()  #str.strip() cleans up whitespaces in row/index
    df.columns = [str(c).strip() for c in df.columns]#str.strip() cleans up whitespaces in columns

    # Keep only 96-well structure columns
    df = df[[c for c in df.columns if c.isdigit()]] #c.isdigit keeps only numeric values (NOT letters)
    df.columns = [int(c) for c in df.columns]

    df = df.dropna(how='all') #remove empty rows

    cleaned = pd.DataFrame(index=df.index, columns=df.columns)
    
#CLEANUP of DATAFARME
#-------------------------------------------------
# This section checks if my raw data values are okay for processing. This reflects real-life cases where 
# the raw data might contain possible errors. For example, the most common issues in reading plates would
# be a sudden "OVERFLOW"  when the reader's capacity has been maxed out. Possible human error could also occur 
# happen in acquiring the data such as missing, negative, or letter-based values which are deemed erroneous.
# This section allows raw data values to be converted into a float ONLY IF our values are valid.
#I asked chatgpt for a script to log these potential errors, AND if observed, tabulate them into error_output.csv 
# This allows to easily trace the location of wells having errors, so that users can re-evaluate the raw files.
    
    for r in df.index:
        for c in df.columns: #Go through every column and row

            val = df.loc[r, c] #get the values for each well [r, c] and check for the ff:

            #1. MISSING VALUES
            if pd.isna(val) or str(val).strip() == "": #== "" empty/ blank cells
                ALL_ERRORS.append({                    #if present, log/ append the error 
                    "plate": plate_name,               #log the plate number/name
                    "sample_code": f"{r}{c}",          #log the sample_code
                    "error": "missing_value"           #and write missing_value as the type of error  
                })
                cleaned.loc[r, c] = np.nan
                continue

            val_str = str(val).strip()

            #2. NUMERIC WITH COMMAS
            if re.match(r'^-?\d{1,3}(,\d{3})*(\.\d+)?$', val_str): #this detects 1,234 and replace it with 1234. 

                num = float(val_str.replace(',', ''))

                if num < 0:                                        #This also logs negative values as error
                    ALL_ERRORS.append({                            #Same thing above, which logs the information about the errors for traceability 
                        "plate": plate_name,
                        "sample_code": f"{r}{c}",
                        "error": "neg_value_raw"                   #log/append the error as neg_value_raw in our error_output.csv
                    })
                    cleaned.loc[r, c] = np.nan
                    continue

                cleaned.loc[r, c] = num
                continue


            #3. INVALID STRING (OVERFLOW, abc, etc.)
            ALL_ERRORS.append({
                "plate": plate_name,
                "sample_code": f"{r}{c}",
                "error": val_str
            })

            cleaned.loc[r, c] = np.nan #any error will be removed from calculations and will be considered NaN automatically

    # Final numeric conversion safety net
    cleaned = cleaned.astype(float) #cleaned dataset contains all the valid numbers stored as float

    return cleaned


# COMPUTATION OF OUR MEDIA (STERILITY) and GROWTH (negative) CONTROLS
#-------------------------------------------------
#Here we used the numpy to have faster math operations and to avoid extensive scripts (if loop is used)
#First we get the mean values of the media and growth for each plate. The values are stored as numpy arrays here
def compute_controls(plate):
    media_vals = np.array([plate.loc[r, c] for r, c in MEDIA_BLANK])
    growth_vals = np.array([plate.loc[r, c] for r, c in GROWTH_CONTROL])

    return np.mean(media_vals), np.mean(growth_vals) #this computes the average/mean for both controls. 
                                                          


# COMPUTES FOR % BACTERIAL INHIBITION
#-------------------------------------------------
#This formula here can also be changed if other analysis will be needed (i.e. %bacterial growth proliferation which is just the inverse of the % inhibition)

def compute_inhibition(plate, media_avg, growth_avg):
    values = plate.values #creates numpy values from each well the samples
    inhibition = (1 - ((values - media_avg) / (growth_avg - media_avg))) * 100 #formula for getting the %inhibition on EACH WELL

    return pd.DataFrame(inhibition, index=plate.index, columns=plate.columns)


# NOTE: This % inhibition values are computed per well. Later we want to summarize them per sample 
# as mean_inhibition with se (standard error), which will be summarized later in the output file "combined_summary.csv" 
# Our intended output file contains the ff info: plate,sample_code,dose,mean_inhibition,se
# As you can see, some of the details included here (i.e. sample_code, dose) must be taken from the metadata.csv.
# Hence, we have to prepare our dataset into a compatible format for merging

# CONVERTING PANDAS DATAFRAME (containing Numpy array) INTO LONG TABLE FORMAT - suggested by chatgpt 
#-------------------------------------------------
# This step is similar to data preparation in Rscript prior to joining metadata to our dataframe.
# We need this step to be able to merge data from Metadata.csv later.
# This converts numpy backed dataframe into long table format which is compatible later during the "df.merge(metadata, on=["row", "col"])

def plate_to_long(inhibition_plate):
    return (
        inhibition_plate
        .reset_index()
        .melt(id_vars="index", var_name="col", value_name="inhibition")
        .rename(columns={"index": "row"})
    )


# GROUPS AND SUMMARIZE THE % INHIBITION VALUES
#-------------------------------------------------
#Note: Our % inhibition values above are computed per well. But since our setup is organized in triplicate values, we want to
#obtain the mean % inhibition as well as the standard error to reveal the variation in our sample. This accounts for possible pipetting errors, etc.

def compute_summary(long_df, plate_name):
    df = long_df.merge(metadata, on=["row", "col"]) #this now merge the data "long_df"" which contains the inhibition values, 
                                                    # and "metadata" which contains the row, col, sample_code, and dose
    df["plate"] = plate_name                        #plate stored as plate_name

    df = df.dropna(subset=["inhibition"])           #remove the row if NaN is observed

    grouped = df.groupby(["plate", "sample_code", "dose"])["inhibition"]    #This creates our sample grouping 
    
                        #Meaning a1, a2, and a3 from Plate1 having same dose is considered as triplicates having same condition
                        #this is important since we want to get the mean %inhibition per triplicate of the samplee
    return grouped.agg(
        mean_inhibition="mean",
        se=lambda x: x.std(ddof=1) / np.sqrt(len(x)) #computes for standard error
    ).reset_index()


# PLOTTING OF HEATMAPS - this is mostly via chatgpt since we only had like 15mins discussion about matplotlib.pyplot
#-------------------------------------------------
#Heatmaps are very important to easily detect trends accross our plates. This identifies the active sample easily accross the samples. 
#At first my plan was to export my data and run the export file into Rscript using ggplot function.
#But since it was discussed that matplotlib.pyplot can also do the same thing with ggplot, 
#I used matplotlib.pyplot package to have just 1 whole python script for everything.

def plot_heatmap_facet(all_inhibition, base_dir): #heatmap_facet to have a tiled heatmaps for all the plates analyzed

    n = len(all_inhibition) # the values to be plotted would be the inhibition values in our plates
    cols = math.ceil(math.sqrt(n))
    rows = math.ceil(n / cols)

    vmin = min(mat.values.min() for _, mat in all_inhibition) #sets the color scaling to the minimum values found in our calculation
    vmax = max(mat.values.max() for _, mat in all_inhibition) #sets the color scaling for max values

    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows))
    axes = np.array(axes).reshape(-1)

    cmap = plt.cm.viridis #theme - viridis

    for i, (plate_name, mat) in enumerate(all_inhibition): #this loops through each plate data

        ax = axes[i]

        mat = mat.sort_index() #sorts our rows (a-h) so that they are not shuffled
        mat = mat.reindex(sorted(mat.columns), axis=1) #sorts the columns (1 to 12)

        data = mat.values
        masked = np.ma.masked_invalid(data) #this handles the missing/invalid values. 
                                            #i.e. in Plate4.csv there are deliberate error values which should leave a blank on the heatmap

        im = ax.imshow(
            masked,
            cmap=cmap,  #colors
            vmin=vmin,
            vmax=vmax,
            aspect="equal",   #square wells
            interpolation="none"
        )
        
  # the following sets the labels and ticks on heatmap
        ax.set_xlabel(plate_name, labelpad=10)
        ax.set_xticks(np.arange(mat.shape[1]))
        ax.set_yticks(np.arange(mat.shape[0]))
        ax.set_xticklabels(mat.columns)
        ax.set_yticklabels(mat.index)
        ax.tick_params(axis="x", rotation=0, labelsize=6)
        ax.tick_params(axis="y", rotation=0, labelsize=6)

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

  # cbar function creates a shared colorbar for all the plates
    # reserve space on right for colorbar
    fig.subplots_adjust(top=0.85)
    cbar_ax = fig.add_axes([0.25, 0.90, 0.5, 0.02])  # top horizontal bar
    cbar = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("% Inhibition")
  
  # output file
    out = os.path.join(base_dir, "FACET_heatmaps.png")
    plt.savefig(out, dpi=100, bbox_inches="tight")
    plt.close()

    print(f"Saved: {out}")


# HEATMAP-but creates a single heatmap containing mean %inhibition and grouped sample_code
# This would output a heatmap but will group similar sample_codes together. Like in the example raw files, all of the plates have same sample_codes but having different pathogen.
# This easily determines the active hits accross different pathogens.

def plot_combined_heatmap(combined_summary, base_dir): #this time, it gets the values from the combined_summary.csv where the mean_inhibition is found

    df = combined_summary.copy()

    #dont include media and growth in the plot
    df = df[~df["sample_code"].isin(["media", "growth"])] 

    # combine plate + dose
    df["plate_dose"] = df["plate"] + "_D" + df["dose"].astype(str) 

    # pivot into one matrix
    pivot = df.pivot(
        index="sample_code",
        columns="plate_dose",
        values="mean_inhibition"
    )

    # sort rows to ascending order (Y axis) to avoid 1-1, 1,10, and 1,2 to appear on the y-axis
    pivot = pivot.reindex(sorted(pivot.index, key=natural_sort_key))

    # sort columns 
    pivot = pivot.reindex(
        sorted(
            pivot.columns,
            key=lambda x: (x.split("_D")[0], int(x.split("_D")[1]))
        ),
        axis=1
    )

    data = pivot.values
    
    #show blank spaces for missing values/ NaN
    masked = np.ma.masked_invalid(data)

    # setting of figure size (prevents overlap)
    fig, ax = plt.subplots(
        figsize=(max(10, len(pivot.columns) * 0.6),
                 max(6, len(pivot.index) * 0.4))
    )

    cmap = plt.cm.viridis #theme/ colormap
    cmap.set_bad(color="white")

    im = ax.imshow(
        masked,
        cmap=cmap,
        aspect="auto"
    )


    # Set out the Axes
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_yticks(np.arange(len(pivot.index)))

    ax.set_xticklabels(pivot.columns, rotation=90, fontsize=7)
    ax.set_yticklabels(pivot.index, fontsize=7)

    # eep Y axis on LEFT (default, explicit for clarity)
    ax.yaxis.tick_left()
    ax.yaxis.set_label_position("left")

    # djust the spacing so labels don't crowd heatmap
    ax.tick_params(axis="y", pad=8)

    ax.set_title("Combined Heatmap (Mean % Inhibition)", pad=20)

    # eave space on right side for colorbar
    fig.subplots_adjust(right=0.85)

    # colorbar axis on far right
    cbar_ax = fig.add_axes([0.88, 0.2, 0.02, 0.6])  # [left, bottom, width, height]

    cbar = fig.colorbar(im, cax=cbar_ax, orientation="vertical")
    cbar.set_label("Mean % Inhibition")

    # Heatmap Layout
    plt.tight_layout(rect=[0, 0, 0.85, 1])

    out = os.path.join(base_dir, "COMBINED_heatmap.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Saved: {out}")



# PLOTTING OF BAR PLOTS
#-------------------------------------------------
#Same thing with the heatmap, most of the commands here were derived from chatgpt.

def plot_bars(df, base_dir):

    df = df.copy()

    COLOR_MAP = {                       #This specifically distinguish the samples type (via dose) from the positive control
        ("media", None): "#808080",
        ("positive_control", 5): "#e41a1c",
        ("positive_control", 50): "#ff7f00",
        ("sample", 5): "#377eb8",
        ("sample", 50): "#4daf4a"
    }

    #creates a facet layout for the bar graphs and automatically change it depending on the number of plates read
    plates = df["plate"].unique()
    cols = math.ceil(math.sqrt(len(plates)))
    rows = math.ceil(len(plates) / cols)

    #creates unique subplot axes. Sometimes the % inhibition varies from negative values to high values.
    fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 4 * rows)) 
    axes = np.array(axes).reshape(-1)

    for i, plate in enumerate(plates): #loops through our plate data

        ax = axes[i]
        sub = df[df["plate"] == plate]

        #defines the x-axis which will be the different sample_code
        samples = sorted(sub["sample_code"].unique(), key=natural_sort_key)
        #this creates a clustered bar graph where sample1 would have two bars representing different doses
        doses = sorted(sub["dose"].unique())

        x = np.arange(len(samples))
        width = 0.8 / len(doses)

        for j, dose in enumerate(doses): #then loops through our dose

            dsub = sub[sub["dose"] == dose] #so that each dose would have a set of bar graphs

            #variables to be plotted
            means = []
            errors = []
            colors = []

            for s in samples: #now we loop through each sample 

                row = dsub[dsub["sample_code"] == s] #and extract our sample_code

                if len(row) > 0:
                    means.append(row["mean_inhibition"].iloc[0]) #if mean_inhibition exists, 
                    errors.append(row["se"].iloc[0]) 

                    #assign the following colors per sample type
                    if s == "media":
                        colors.append(COLOR_MAP[("media", None)])
                    elif s == "positive_control":
                        colors.append(COLOR_MAP[("positive_control", dose)])
                    else:
                        colors.append(COLOR_MAP[("sample", dose)])
                        
              #for missing data, it creates invisible bars (NaN)
                else:
                    means.append(np.nan)
                    errors.append(0)
                    colors.append("black")
            
        #bar axis formatting
            ax.bar(x + j * width, means, width, yerr=errors, capsize=3, color=colors)

        ax.set_title(f"Plate {plate}")
        ax.set_xticks(x + width * (len(doses) / 2))
        ax.set_xticklabels(samples, rotation=90, fontsize=6)
        ax.set_ylabel("% Inhibition")
    
    #loop through the plate data and remove any empty subplots
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    #output parameters. You can increase resolution if needed
    out = os.path.join(base_dir, "QC_facet_bars.png")
    plt.tight_layout()
    plt.savefig(out, dpi=100, bbox_inches="tight")
    plt.close()

    print(f"Saved: {out}")
    
def natural_sort_key(s):
    return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', str(s))]


# OUTPUT SUMMARY FILE
#-------------------------------------------------

all_summary = [] # this contains the list for the summary per plate
all_inhibition = [] #creates list for the inhibition which contains plate_name, inhibition_matrix

for FILE_PATH in file_list: #go through each plate

    plate_name = os.path.splitext(os.path.basename(FILE_PATH))[0] #extract the plate name (i.e. /data/plate1.csv will retrieve just the plate1.csv)
    

    plate = load_plate(FILE_PATH) #loads the folder path where the file is located

    media_avg, growth_avg = compute_controls(plate) #then compute for the controls per plate

    inhibition_plate = compute_inhibition(plate, media_avg, growth_avg) #proceed with computing the %inhibition per well

    all_inhibition.append((plate_name, inhibition_plate)) #store the inhibition data into inhibition_plate. This will be used for heatmap

    long_df = plate_to_long(inhibition_plate) #converts to compatible long format (i.e well a1 to row a and column 1)

    summary = compute_summary(long_df, plate_name) #now let;s merge the metadata with the compute_summary

    all_summary.append(summary) #and store it as summary

BASE_DIR = os.path.dirname(file_list[0]) #setting our output directory

combined_summary = pd.concat(all_summary, ignore_index=True)

combined_summary.to_csv(os.path.join(BASE_DIR, "combined_summary.csv"), index=False)

#This filters out samples having >60% inhibition. this is helpful to have a summarized list of samples which we can consider as active antimicrobial agent.
active_hits = combined_summary[combined_summary["mean_inhibition"] > 60] 
active_hits.to_csv(os.path.join(BASE_DIR, "active_hits.csv"), index=False)

print("CSV files saved.")


# ERROR OUTPUT FILE
#-------------------------------------------------
#this is to print out any detected errors especially in our raw data
if ALL_ERRORS:  #if errors were detected
    error_df = pd.DataFrame(ALL_ERRORS)
    error_df.to_csv(os.path.join(BASE_DIR, "error_output.csv"), index=False) #log them into the error_output.csv
    print(f" !!!!ERRORS FOUND!!! {len(ALL_ERRORS)} (saved to error_output.csv)") #and caution ERRORS FOUND
else:
    print("No data errors found.") #otherwise, no error_output would be created.
    
#NOTE: Raw data plate1.csv, plate2.csv, and plate3.csv have no errors.
# Plate4.csv has deliberate errors on it so if this will be included in the pipeline, error_output.csv should be present


# PLOT our heatmap and bar graph
#-------------------------------------------------
plot_heatmap_facet(all_inhibition, BASE_DIR)
plot_combined_heatmap(combined_summary, BASE_DIR)
plot_bars(combined_summary, BASE_DIR)

print("Done.")
