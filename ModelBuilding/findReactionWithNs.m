function reactions = findReactionWithNs(reaction_map)
    
   reaction_ids = keys(reaction_map);
   reactions = {};
   
   for ind = 1:numel(reaction_ids)
        id = reaction_ids(ind);
        data = reaction_map(id{:});
        equation = data.equation;
        
        [metaboliteList,~,~] = parseRxnFormula(equation);
        
        metabolites = regexp(metaboliteList, 'C\d{5}|G\d{5}', 'match');
        num_metab = numel(metabolites);
        
        metabolites(cellfun('isempty',metabolites)) = [];
        
        if num_metab > numel(metabolites)
%             fprintf(fid, '%s\n', data.entry); 
            reactions = [reactions; {data.entry}];
        end
                
   end    
end