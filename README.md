# SRS-FISH drugs

This repository is adapted from the codebase used to produce the results in the paper "The Parkinson’s drug entacapone disrupts gut microbiome homeostasis via iron sequestration".

## Requirements
The code in this repo can be implemented with the following software
* CellProfiler 4.2.1
* Matlab 2023a

## Files
This repo should contain the following files:
* CD_ratio_V10_FISH masked_CHnoFISH masked_Bac Ratio.cpproj - single cell statistical analysis on FISH-masked-SRS data with CH masks excluding FISH masks
* V2.1_Drug_distribution_FISH_masked_drug_intensity.cpproj - single cell statistical analysis on FISH-masked drug intensity data
* V2.1_Drug_distribution_FISH_masked_image_save.cpproj - FISH masked drug distribution image processing for visualization
* Drug_Intensity_BoxPlot_v3.m - statistical check (box plot) on single cell drug intensity
* WidefieldPointscanRegistration_V2.m - widefield and point scanning image registration

## Contact
If you have any questions, please contact xwge@bu.edu.
