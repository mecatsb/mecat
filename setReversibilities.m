function model = setReversibilities(model, revMap)    
   
    reactions = keys(revMap);  
    
    for k = 1:numel(reactions)      
        r = reactions(k);
        IndexC = strfind(model.rxns, r{:});
        Index = not(cellfun('isempty', IndexC)); 
        model.rev(Index) = revMap(r{:});
    end
end