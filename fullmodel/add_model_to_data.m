
function simdata=add_model_to_data(rawdata,likey)

if numel(likey)~=numel(rawdata)
    error('# trials in data and model must be equal')
end

    simdata = rawdata;
    for tx = 1:numel(simdata)
        if rawdata(tx).hit==1
            simdata(tx).probHit = likey(tx);
        else
            simdata(tx).probHit = (1-likey(tx));
        end
        
        if rawdata(tx).pokedR==1
            simdata(tx).probR = likey(tx);
        else
            simdata(tx).probR = (1-likey(tx));
        end
        
        
    end