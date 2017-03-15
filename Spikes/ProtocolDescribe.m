function ppp = ProtocolDescribe(p)
% Describes the active and passive parameters of a protocol
%
% ppp = ProtocolDescribe(p)
%
% 2010-05 MC taken out from FigPick_callbacks

% p = PICK.protocol;

if isempty(p)
    ppp = '';
    return
end

nstim = p.nstim;
nonblanks = setdiff(1:nstim,p.blankstims);

if isempty(nonblanks)
    % they are all blank
    ppp = {'All stimuli are blank'};
    return
end

ppp = cell(p.npars,1);
for ipar = 1:p.npars
    ppp{ipar} = [ p.pardefs{ipar} ' -- ' p.parnames{ipar} ];
    if ismember(ipar, [p.activepars{:}]);
        parvalue = NaN;
        ppp{ipar} = [ppp{ipar} ' is ACTIVE'];
    else
        parvalue = unique(p.pars(ipar,nonblanks));
        if length(parvalue)>1,
            parvalue = parvalue(1);
            disp('Warning: inactive parameter has more than one value!');
        end
        ppp{ipar} = [ppp{ipar} ' = ' num2str(parvalue)];
    end
    if p.adapt.flag
        adaptvalues = unique(p.pars(ipar,nstim+(1:2) ));
        if adaptvalues~=parvalue
            ppp{ipar} = [ppp{ipar} ' --- ADAPT: ' mat2str(adaptvalues) ];
        end
    end
    % disp(ppp{ipar})
end