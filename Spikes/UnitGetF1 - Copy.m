function [F1mean, F1err] = UnitGetF1( varargin )
% UnitGetF1 returns the first harmonic firing rate response 
%
% [F1mean, F1std] = UnitGetF1( unit, protocol ) 
% returns the mean and the standard deviation of the F1 responses (computed at frequencies fs)
%
% The F1 is reported PEAK to PEAK
%
% [F1mean, F1sem] = UnitGetF1( unit, protocol, 'std' ) 
% does the same
%
% [F1mean, F1sem] = UnitGetF1( unit, protocol, 'sem' ) 
% returns the mean and the standard error of the mean of F1 the responses
%
% [F1mean, F1sem] = UnitGetF1( unit, protocol, ..., 'broad' ) 
% looks around for the biggest response near the indicated frequencies. Useful for those
% case when you are not entirely sure of the stimulus frequency.
%
% OBSOLETE: REPLACED BY UnitGetHarm(1,...)
%
% part of Spikes
%
% 2000-11 MC
% 2001-07 MC moved all the code to UnitGetHarm

disp('---------> Calls to UnitGetF1 are obsolete. Replaced by UnitGetHarm(1,...)');

[F1mean, F1err] = UnitGetHarm(1, varargin{:});
