% Spikes Toolbox.
%   Version 2010.01.13
%  
%   Toolbox to load and analyze electrode data in the Carandini lab.
%  
%   To run interactively, type "ExptPick". To work in batch, read below and
%   consider starting from examples in SpikesExamples.m. 
%  
%   The toolbox is based on four data structures:
%   Protocol      the parameters of an experiment
%   Experiment    the raw data for an experiment
%   Channel       the data for one channel, with sorted spikes
%   Unit          the data for one discriminated neuron
%  
%   Data are saved in Animals, which contain Series, which contain
%   Experiments, which contain Stimuli.
%   
%   A global DIRS needs to be declared with fields
%   DIRS.data     The directory where the raw data are kept
%   DIRS.spikes   The directory where the analyzed data are kept
%  
%  % --------------------------- Functions ----------------------------------
%  
%   SpikesExamples            - Examples of batch scripts to analyze data
%  
%%
%  % GRAPHICAL USER INTERFACES
%   ExptPick                  - (formerly SPIKES) defines a few environment variables and launches the GUI
%   ExptPickFcn               - (formerly SPIKES) defines a few environment variables and launches the GUI
%   DirDialog                 - Open directory dialog window (for ExptPick)
%   FigPicker                 - Matlab-generated figure for the FigPicker GUI
%   FigPicker_callbacks       - switchyard for the FigPicker GUI
%   FigSpikeSorter            - Matlab-generated figure for the FigSpikeSorter GUI
%   FigSpikeSorter_callbacks  - switchyard for the FigSpikeSorter figure
%   FigTuning                 - Matlab-generated figure for the FigSpikeSorter GUI
%   FigTuning_callbacks       - switchyard for the FigTuning figure
%   BringToTop                - Brings to the top the objects with the same UserData as the callback object (in FigSpikeSorter)
%   ProtocolInspect           - MATLAB-generated figure for ProtocolInspect GUI
%   ProtocolInspect           - shows a table of the parameters values in a protocol
%%
%  
%  % EXPERIMENT - BASIC
%   ExptInitialize            - Initializes an Experiment
%   ExptLoad                  - Loads or updates an Experiment
%   ExptLoadMichigan          - Loads or updates a Michigan experiment
%   ExptSaveAnalysis          - Creates Channels and Units and saves them
%   ExptReadLogFile           - Reads the log file for an Experiment
%%
%  
%  % EXPERIMENT - UTILITIES
%   ExptEstimateNoise         - Estimates the noise as described by Blanche et al
%   ExptFilter                - Filters the traces if filter coefficients are provided
%   ExptTemplateMatch         - Work in progress by Matteo, does not create any Units, just inspection by eye
%%
%  
%  % PROTOCOL
%   ProtocolInitialize        - Initializes a Protocol data structure
%   ProtocolLoad              - Loads a Protocol data structure
%   ProtocolSave              - Saves a Protocol in the appropriate directory
%   ProtocolTimingIsBuggy     - Checks whether a Protocol might have timing problems
%   ScreenLogLoad             - Loads the screenlog of an experiment
%%
%  
%   ProtocolAbsentMatFiles    - See which files are present and gives a report
%   ProtocolGetGraphInfo      - Gets info on what to plot for a given protocol
%   ProtocolGetGratingInfo    - Gives information about gratings in a stimulus
%%
%  
%  % UNIT - BASIC FUNCTIONS
%   UnitCreate                - Creates an empty Unit data structure
%   UnitLoad                  - Loads a Unit data structure
%   UnitMake                  - Creates a Unit from an Expt and Chan data structure
%   UnitSave                  - saves a Unit data structure
%%
%   
%  % UNIT - ANALYSIS FUNCTIONS
%   UnitGetRates              - A Unit's firing rates
%   UnitGetFilename           - The file name of a file for a Unit data structure
%   UnitGetCycles             - Cycle averages of a Unit's firing rates
%   UnitGetDC                 - Mean of a Unit's firing rates 
%   UnitGetF1                 - First harmonic of a Unit's firing rates
%   UnitGetHarm               - Frequency component of a Unit's firing rates
%   UnitGetStimFreqs          - UnitGetStimFreqs returns the frequency of the stimulus
%%
%  
%  % GRAPHICS
%   ProtocolGetGraphInfo      - Gets info on what to plot for a given protocol
%   UnitPlotCycles            - Plots cycle hystograms
%   UnitPlotFreq              - Plots power spectrum of a Unit's spike trains
%   UnitPlotHistos            - Plots the histograms of a Unit's firing rates
%   UnitPlotRasters           - plots the rasters of a Unit's spike trains
%   UnitPlotTuning            - plots the tuning curves for a Unit
%%
%  
%  % UTILITIES
%   SetDefaultDirs            - Sets DIRS variable
%   SetDefaultDirs_Znetgear1  - Dummy SetDefaultDirs replacement, in case primary server is down
%   ListAnimals               - Lists all the animals in a data directory
%   ListSeries                - Lists all the series that were run for a given animal
%   ListExpts                 - Lists all the expts that were run for a given series
%   ProtocolFind              - Finds experiments with specified active parameters
%   ProtocolGetGratingInfo    - Gives information about gratings in a stimulus
%   XFileLoad                 - Loads an x file
%   LoadMichiganFPs           - Loads the field potentials for Michigan probe recordings
%%
%  
%  % channel - MOSTLY USED BY THE SPIKE SORTER
%   ChanLoad                  - Loads a saved Channel
%   ChanImport                - Imports Channels for an Experiment
%   ChanInitialize            - Initializes a Channel
%   ChanGetFileName           - Gives the file name of a Channel
%   ChanFindCandidates        - Finds candidate spikes by thresholding
%   ChanFindPrototypes        - Finds prototype spikes using unsupervised learning.
%   ChanSpikeClassify         - Classifies spikes for a Channel
%%
%  
%  % INTERFACE WITH VISBOX TOOLBOX
%   ShowStimulusFrame         - Shows a frame of the visual stimulus
%  
%  %
% 
% 
%   For processing spike files
%  
%   CCGBinContents - finds spike pairs corresponding to a certail CCG bin
%   CSD - current source density analysis
%   ClusterFind - locates a cluster in a .clu.n file from a cluster in a.clu file
%   ClusteringErrors - Evaluate false positive and negative %ages of a clustering
%   Feature - a matlab version of sfeature
%   FiringRate - computes firing rate from .clu and .res files
%   ISIGrams - collect and plot a load of ISI histograms for several cells  
%   LoadClu - loads a .clu file
%   LoadFet - loads a .fet file
%   LoadPar - loads a .par parameter file
%   LoadPar1 - loads a .par.n parameter file
%   LoadSpk - loads a .spk file
%   LoadSpkPartial - loads a subset of the spikes from a .spk file (to save memory)
%   LoadData - loads in .res .clu and .spk file
%   LoadSingleChan - loads a single channel from a .dat file
%   MahalView - Plots mahalanobis distance of fale positives and negatives
%   Peak2Peak - Detects peak to peak amplitudes out of a spike array.
%   PointCorrel - plots a cross correlogram
%   PointCorrelM - cross correlograms of multiple spike trains
%   ResampleRealignSpikes - resamples a bunch of spikes and realigns on their peaks
%   SaveClu - saves a .clu file
%   SaveFet - saves a .fet file
%   SaveSpk - saves a .spk file
%   SmoothFiring - turns a vector of spike times into an estimated firing rate function
%   SpikePCA - calculate PCA coefficients given a spike file AND PC waveforms
%   SpikeExtract - gets spikes from a .dat file using intracellular spike info
%   SpikeExtractRealign - gets spikes from a .dat file using intracellular spike info, resamples and realigns
