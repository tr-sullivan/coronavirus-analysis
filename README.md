# Torus Biosystems Coronavirus Panel Analysis Pipeline
## Section I. System Requirements
Software is designed to be run on MATLAB on a Windows computer. All MATLAB versions from 2020b to 2022b have been used without error. The Image Processing and Computer Vision toolbox is required for the spot-finding function used during masking. 

## Section II. Installation Guide
Aside from MATLAB and the required toolbox, no other installation steps are necessary. Ensure all included code is on the MATLAB path. 

## Section III. Demo
A set of experimental data used in the manuscript is provided for demo purposes. The structure of the included code and data is organized as follows: 

The main folder (Analysis Pipeline) includes 3 subfolders. 
- Sample data contains several subfolders, each of which contains 60 images acquired from RT-caPCR experimental runs as described in the manuscript. Images were taken every 30 seconds and are saved in the .raw format. 
- MainAnalysis contains all functions needed to run the main code. To run analysis for each experimental dataset, follow these instructions:
-- Open the AnalyzeImages_grid.m script in the MATLAB editor and navigate into the directory containing the images for the run to be analyzed. 
-- Run the AnalyzeImages_grid.m script.
-- Move the cursor over the top-left, top-right, bottom-left, and bottom-right corner spots of the array when prompted, and click at each location (these instructions are also displayed in MATLAB's command window).
-- Figures will be generated showing the 'raw', 'baseline-subtracted', and 'normalized' data for the run. Color coding on these plots specifies if the target was detected (red) to be amplifying or not based on the criteria described in the manuscript. 
-- Images are also generated showing the 'subtracted' behavior of the run, i.e., the final image minus the second image, which gives a quick readout for what was found to be present. 
-- Other generated files include .jpg versions of each image for quick viewing; MATLAB files storing the images, generated masks, assay information, and data; and Microsoft Excel files including the data and calculated threshold times (a measure of target abundance, as in standard qPCR). 

- PositiveNegativeOverlays and SNPCalling contains two other scripts that are used to generate the two-plot overlays shown in the manuscript for each dataset and to automatically determine the dominant SNV at each SARS-CoV-2 mutation site based on the criteria specified in the manuscript. 
-- GenerateFluorescenceOverlays_Respiratory.m is designed to be run in the PositiveNegativeOverlays and SNPCalling folder and will generate two-plot summary figures for each analyzed dataset. The left-hand plot shows fluorescence traces for all targets expected to be amplifying based on the known sample input and the right-hand plot shows the same for targets expected to be absent. Note that some SNV-sensitive probes will show a low degree of off-target activity in the presence of similar sequences. Traces are color-coded based on probe (red for robust viral-detection; blue for SNV-sensitive discriminating) and sample type (gray for negative-control probes; green for extraction and amplification controls). 
-- SNPCalling_Loop.m is also designed to be run in the PositiveNegativeOverlays and SNPCalling folder and will iterate through all analyzed datasets and will generate predicted 'winners' for each of the five SNV clusters contained on the panel. Each individual dataset's prediction will be reported in SNP Calls.txt in the dataset's folder, and all predictions are written to SNP Calls_all.txt in the PositiveNegativeOverlays and SNPCalling folder. Note that predictions are generated for all datasets regardless of whether SARS-CoV-2 was found to be present, but these should be ignored unless the latter criterion is satisfied. 

- Each dataset can be analyzed in <5 minutes. Downloading the images is the largest overall contribution to the runtime. Post-processing for SNV detection and generation of positive-negative plots is quick. All analyses can be run on ThinkPad laptops. 

## Section IV. Instructions for Use
See above section describing how to run demo data. 
