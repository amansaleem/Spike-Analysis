function ExptPlotLFPHarm( expt, channame, harm )
% plots a harmonic component of the LFP
%
% ExptPlotLFPHarm( expt, harm )
%
% Example:
% expt = ExptLoad(PICK);
% ExptPlotLFPHarm( expt, 6, 2 );
%
% 2010-03 MC

% harm = 2;
% channame = 6;

if isempty(expt)
    error('expt is empty');
end

p = ProtocolLoad(expt);

ichan = (expt.channames == channame);
fs = expt.samplerate;

for istim = 1:p.nstim
    
    StimFreq = p.pfilefreqs(istim);
    
    for irepeat = 1:expt.nrepeats
        
        yy = expt.data{ichan}{istim,irepeat};
        yy = double(yy);
        yy = yy - mean(yy);
        
         tt = (1:length(yy))/fs;
         
%          ff = freq(length(yy), length(yy)/fs);
%          figure; stem( ff, abs(fft(yy))); set(gca,'xlim',[0 10]);
         
        rr(istim,irepeat) = sum( yy .*  exp(2*pi*1i*tt*harm*StimFreq))/length(tt);
         
    end
end


figure; plot( mean(abs(rr),2) )
