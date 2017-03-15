function chanfilename = ChanGetFileName(spikedir,animal,iseries,iexp,DataType)
% ChanGetFileName Gives the file name of a channel file
%
% chanfilename = ChanGetFileName(spikedir,animal,iseries,iexp)
%
% chanfilename = ChanGetFileName(spikedir,animal,iseries,iexp,DataType)
% lets you specify the type of data, which can be 'Michigan', or
% 'Multispike' (DEFAULT).
%
% 2000 Matteo Carandini created
% 2007-09 MC added 5th argument

if nargin<5
    DataType = 'Multispike';
end

switch DataType
    case 'Multispike'
        chanfilename = ['discrim_pars_expt_' num2str(iexp) '.mat'];
    case 'Michigan'
        chanfilename = ['discrim_pars_Michigan_expt_' num2str(iexp) '.mat'];
    case 'CerebusTraces'
        chanfilename = ['DiscrimParsCerebusTraces_' num2str(iexp) '.mat'];
    otherwise
        error('Do not know this DataType');
end

chanfilename = fullfile(spikedir,animal,num2str(iseries),chanfilename);
