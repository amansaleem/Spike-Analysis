function expt = ExptInitialize
% ExptInitialize initializes an expt data structure
%
% expt = ExptInitialize
%
% part of Spikes
%
% MC 09 2000 created
% MC 12 2009 added ChanDescriptions

expt.animal           = '';
expt.iseries          = NaN;
expt.iexp             = NaN;
expt.nrepeats         = 0;
expt.nstim            = NaN;
expt.stimdurs         = [];
expt.nchans           = NaN;
expt.data             = [];
expt.timestamp        = NaN;
expt.channames        = [];
expt.unitspervolt     = NaN;
expt.samplerate       = NaN;
expt.threshsigns      = [];
expt.threshvalues     = [];
expt.ChanDescriptions = {};


