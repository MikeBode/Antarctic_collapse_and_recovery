This repository contains the data and MATLAB code needed to reproduce the
analysis and figures in the manuscript: 

"Collapse and recovery of decision-making in Antarctic governance"
Bode M, Botelho LL, Yaryab P, Liggett D.

Running one script rebuilds everything from the raw data.


WHAT IS IN THE REPOSITORY

RawData/
  MeasureCorpus.csv          Every Measure, Recommendation, Decision and
                             Resolution adopted under the Antarctic Treaty,
                             with the year adopted and whether it is in force.
  ATCM_IPWP.xlsx             Working and Information Papers submitted to each
                             Antarctic Treaty Consultative Meeting (ATCM),
                             with the Parties that submitted them.
  Country_joined.xlsx        Treaty Parties and the year each became a
                             Consultative Party.
  ATCM_roman_numerals.xlsx   List of ATCMs with their number, year and host.

Code
  RunAllAnalyses.m                Main script. Runs everything below in order.
  ExtractMeasureData.m            Counts Measures adopted and entered into
                                  force at each ATCM.
  ExtractNetworkData.m            Builds networks of which Consultative
                                  Parties co-authored Working Papers, over
                                  windows of neighbouring meetings.
  CreateModelTimeseries.m         Runs the consensus model on the networks
                                  over a grid of parameter values (slow).
  FitModel.m                      Finds the parameter values that best match
                                  the observed data, with a bootstrap.
  Figure1_MeasuresAndNetworks.m   Draws Figure 1.
  Figure2_ModelledObserved.m      Draws Figure 2.
  AddNetworkToFigure.m            Draws a single network in Figure 1.
  PrintRunConfig.m                Prints the settings each step is using.
  MakeTIFF.m                      Saves a figure as a TIFF.


REQUIREMENTS

  MATLAB (developed with R2025b)
  Statistics and Machine Learning Toolbox
  Curve Fitting Toolbox

The figures also use slanCM, a colormap function by Zhaoxu Liu, from the
MATLAB File Exchange ("200 colormaps"):
  https://www.mathworks.com/matlabcentral/fileexchange/120088-200-colormaps
Put slanCM.m and slanCM_Data.mat in the repository folder.


HOW TO RUN IT

1. Open MATLAB and change directory to the repository folder.
2. Type:  RunAllAnalyses

The script creates the folders CreatedData/ and Figures/ if they do not exist.

The model simulation (step 3) takes around 9 hours at full resolution. If
CreatedData/Timeseries/pk_flat_store.mat is already there, you can skip it
by setting RunMode = 'figures' at the top of RunAllAnalyses.m. The script
checks that the stored file was made with the same settings before using
it, and the rest then runs in a minute or two.

All the model settings (parameter grid, number of Monte Carlo runs, random
seed, spans searched, fit years) are set at the top of RunAllAnalyses.m.


WHAT IT PRODUCES

Figures/
  Figure1_MeasuresAndNetworks.tif   Figure 1
  Figure2_ModelledObserved.tif      Figure 2

CreatedData/
  MeasureData.mat              Measures adopted and in force per meeting
  AveragedNetworks/            Co-authorship networks, one file per span
  Timeseries/pk_flat_store.mat Model output for every parameter set
  BestFit.mat                  Best-fitting parameters, SSR for every
                               parameter set, and bootstrap results
  Figure1_MovingAverage.mat    Moving average line used in both figures

When it finishes, the script prints the best-fitting parameters and warns
if the best fit is on the edge of the parameter grid.


LICENCE

See License file. 

Note that: slanCM.m and slanCM_Data.mat are by Zhaoxu Liu and are not 
covered by this licence. They are distributed under their own licence 
from the MATLAB File Exchange.

CONTACT
Michael Bode
michael.bode@qut.edu.au
