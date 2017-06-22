function reactions = findReactionWithSameEductProduct(reaction_map)
    
   reaction_ids = keys(reaction_map);
   reactions =[];
   
   for ind = 1:numel(reaction_ids)
        id = reaction_ids(ind);
        data = reaction_map(id{:});
        equation = data.equation;
        
        [metaboliteList,stoich,~] = parseRxnFormula(equation);
        
        educts = metaboliteList(stoich < 0);
        products =metaboliteList(stoich > 0);
        
        if numel(unique(metaboliteList)) < numel(metaboliteList)
           reactions = [reactions; id]
        end
        
   end    
end