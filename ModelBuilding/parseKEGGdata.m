 function parseKEGGdata(compound_db, reaction_db, rpair_db, enzyme_db)
    compound_map = containers.Map();
    compound_ids = keys(compound_db);

    for ci = 1:numel(compound_ids)
            abbrev = compound_ids(ci);
            abbrev = abbrev{:};
            data = compound_db(abbrev);
            [~, compound] = parseCOMPOUND(data);
            compound_map(abbrev) = compound;
    end

    save(['parsed_KEGG_COMPOUND_' date '.mat'], 'compound_map')

    reaction_map = containers.Map();
    reaction_ids = keys(reaction_db);

    for ri = 1:numel(reaction_ids)
            abbrev = reaction_ids(ri);
            abbrev = abbrev{:};
            data = reaction_db(abbrev);
            [~, parsed] = parseREACTION(data);
            reaction_map(abbrev) = parsed;
    end  
    
    save(['parsed_KEGG_REACTION_' date '.mat'], 'reaction_map')

    rpair_map = containers.Map();
    rpair_ids = keys(rpair_db);

    for ri = 1:numel(rpair_ids)
            abbrev = rpair_ids(ri);
            abbrev = abbrev{:};
            data = rpair_db(abbrev);
            [~, parsed] = parseRPAIR(data);
            rpair_map(abbrev) = parsed;
    end
    
    save(['parsed_KEGG_RPAIR_' date '.mat'], 'rpair_map')

    enzyme_map = containers.Map();
    enzyme_ids = keys(enzyme_db);

    for ri = 1:numel(enzyme_ids)
            abbrev = enzyme_ids(ri);
            abbrev = abbrev{:};
            data = enzyme_db(abbrev);
            parsed = parseENZYME(data);
            
            e = parsed.entry; 
            e = e{:};
            enzyme_map(e) = parsed;
    end
     save('parsed_KEGG_ENZYME_mat', 'enzyme_map')
    save(['parsed_KEGG_ENZYME_' date '.mat'], 'enzyme_map')
end