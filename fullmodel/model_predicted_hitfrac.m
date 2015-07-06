function [hitfrac simdata] = model_predicted_hitfrac(rawdata, likey)
simdata = rawdata;
for i = 1:numel(simdata), 
    if simdata(i).pokedR==1, 
        simdata(i).pokedR = likey(i); 
    else
        simdata(i).pokedR = 1-likey(i); 
    end; 
end;

hits = zeros(numel(simdata), 1);
samples = zeros(numel(simdata), 1);
nleftbups = zeros(numel(simdata), 1);
nrightbups = zeros(numel(simdata), 1);
for i = 1:numel(simdata),
    samples(i)    = rawdata(i).T;
    nleftbups(i)  = sum(rawdata(i).leftbups < samples(i));
    nrightbups(i) = sum(rawdata(i).rightbups < samples(i));
    dbups         = nrightbups(i) - nleftbups(i);
    
    if dbups > 0,
        hits(i) = simdata(i).pokedR;
    elseif dbups < 0,
        hits(i) = 1 - simdata(i).pokedR;
    else % if there's no right answer, flip a coin
        hits(i) = 0.5;
    end;
end;

hitfrac = mean(hits);