function number_cofac = countCofactors(cofactors, metabolites, reaclist, reactions, reaction_map)
   
    kegg_ids = metabolites(cofactors);
    kegg_ids = regexprep(kegg_ids,'\[c\]','');
    cm  = countMetabolite(kegg_ids, reaclist, reactions, reaction_map);
    
    %count cofactors that are used
    number_cofac = abs(sum(cm(cm < 0)));
end