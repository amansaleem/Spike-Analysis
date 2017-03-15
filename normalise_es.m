function [es, bins]         = normalise_es(es_input, c, pGrid, vrGrid, type)

es = es_input;

[bins.P_bins es.traj] = normalise1var(es.traj, pGrid,c);
[bins.V_bins es.trajspeed] = normalise1var(es.trajspeed, vrGrid,c);
[bins.R_bins es.ballspeed] = normalise1var(es.ballspeed, vrGrid,c);

if strcmp(type,'OL')
    for iangle = 1:numAngles
        [bins.projBins{iangle} es.projSpd{iangle}] = normalise1var(es.projSpd{iangle}, vrGrid,c);
    end
end

end