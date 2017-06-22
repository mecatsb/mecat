function [reaction_ids, metabolite_reaction] = findReaction(reactands, reaction_map, reactions)
    
    reaction_ids = {};
    metabolite_reaction = containers.Map(reactands, zeros(numel(reactands),1));
    
    for k = 1:numel(reactions)
       ind = reactions(k);
       reaction_data = reaction_map(ind{:});
       educts = reaction_data.educts;
       products = reaction_data.products;
       
       if all(ismember(educts, reactands))
            reaction_ids = [reaction_ids; ind];  
            
            for m = 1:numel(educts)
            	metabolite_reaction(educts{m}) = metabolite_reaction(educts{m})+1;    
            end
            
        elseif all(ismember(products, reactands))
            reaction_ids = [reaction_ids; strcat('-',ind)];
            
            for m = 1:numel(products)
            	metabolite_reaction(products{m}) = metabolite_reaction(products{m})+1;     
            end
       end    
    end  
end