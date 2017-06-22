function reaction_lists = setReactionLists(reaction_lists, start_reactions)

    for k = 1:numel(reaction_lists)
            reaction_list = reaction_lists{k};
            
            if isempty(reaction_list)
                continue;
            end
            
            reaction_ids = reaction_list{1,1};
            reacs = ismember(reaction_ids, start_reactions);
            if any(reacs)   
                reaction_lists{k}{1,1} = reaction_ids(reacs);
            else
                reaction_lists{k} = [];
            end       
    end   
end