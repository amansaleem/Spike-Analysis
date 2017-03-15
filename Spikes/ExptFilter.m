function expt = ExptFilter(expt,chan)
% ExptFilter filters the traces if filter coefficients are provided
%
% expt = ExptFilter(expt,chan);
% (uses chan.filtercoeffs)
%
% 2009-02 MC extracted from ChanFindCandidates
% 2010-03 MC removed mean before filtering
% 2010-06 MC moved the test for class type downward
% 2010-06 MC added a test for empty vv

ichan = chan.iexptchan;

if ~isempty(chan.filtercoeffs)
    b = chan.filtercoeffs{1};
    a = chan.filtercoeffs{2};
    fprintf(1,'Filtering...');   
    for istim = 1:expt.nstim
        for irepeat = 1:expt.nrepeats
            ExptNumType = class(expt.data{ichan}{istim,irepeat});
            % typically int16 (for TDT data) or single (for multispike data)
            vv = single(expt.data{ichan}{istim,irepeat});
            if ~isempty(vv) % MC added this 2010-06-14
                vv = vv-mean(vv); % MC added this 2010-03-16
                % ramp up the first 2 ms of data
                n = round(2* expt.samplerate/1000);
                vv(1:n) = vv(1:n).* (1:n)/n;
                vv = filter(b,a,vv);
            end
            expt.data{ichan}{istim,irepeat} = cast(vv,ExptNumType);
            % figure; plot(vv-mean(vv)); hold on;
            % plot(expt.data{ichan}{istim,irepeat},'r');
        end
    end
    fprintf(1,'done\n');
end