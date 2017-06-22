function start_reactions = findStartReactions(usable_metabolite_ids, metabolites, reactions, reaction_map)

    usable_metabolites = metabolites(usable_metabolite_ids);
    usable_metabolites = regexprep(usable_metabolites,'\[c\]','');

    num_reactions = numel(reactions);
    start_reactions = zeros(num_reactions,1);

    for id = 1:num_reactions
        r = reactions{id};
        
        if ~isempty(strfind(r, '-'))
            r = regexprep(r, '-', '');
            reaction_data = reaction_map(r);
            educts = reaction_data.products;
        else
            reaction_data = reaction_map(r);
            educts = reaction_data.educts;      
        end

        start_reactions(id) = all(ismember(educts, usable_metabolites));            
    end
    l = logical(start_reactions);
    start_reactions = find(l);
    
end