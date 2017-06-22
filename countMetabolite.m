function cofac_balance = countMetabolite(metabolites, reaclist, reactions, reaction_map)

    reactand_map = containers.Map();
    num_reactions = numel(reaclist);
    num_metabolites = numel(metabolites);
    
    cofac_balance = zeros(num_metabolites,1);

    for id = 1:num_reactions
        if reaclist(id) == 0
            continue;
        end
        
        r = reactions(reaclist(id));
        [~,educts , products, stoich_e, stoich_p] = getReactionPartners(r, reaction_map);

        num_educts = numel(educts);
        num_products = numel(products);

        for num = 1:num_educts
           if isKey(reactand_map, educts{num}) 
               reactand_map(educts{num}) = reactand_map(educts{num}) - str2double(stoich_e{num});
           else
               reactand_map(educts{num}) = - str2double(stoich_e{num});
           end
        end
        
        for num = 1:num_products
               if isKey(reactand_map, products{num})
                   reactand_map(products{num}) = reactand_map(products{num}) + str2double(stoich_p{num});
               else
                   reactand_map(products{num}) = str2double(stoich_p{num});
               end
        end 
    end
    
    for k = 1:num_metabolites   
        if isKey(reactand_map,metabolites{k})
            cofac_balance(k) = reactand_map(metabolites{k});
        else
            cofac_balance(k) = 0;
        end
    end
    
end