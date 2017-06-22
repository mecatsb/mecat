function writeOverallPathwayBalance(reaclist, reactions, reaction_map, thermodynamics_map, outfile)

    reaclist = reaclist{:};
    reaclist(reaclist==0) = [];
    reactand_map = containers.Map();
    fid = fopen(outfile, 'w+');
    
    num_reactions = numel(reaclist);
        
    pathway_balance.reactions = cell(num_reactions,1);
        
    for id = 1:num_reactions
        r = reactions(reaclist(id));
        pathway_balance.reactions(id,1) = r;

        if ~isempty(strfind(r, '-'))
            r = regexprep(r, '-', '');
            reaction_data = reaction_map(r{:});
            products = reaction_data.educts;
            educts = reaction_data.products;
            stoich_p = reaction_data.stoich_e;
            stoich_e = reaction_data.stoich_p;
        else     

           reaction_data = reaction_map(r);
            
            if isKey(thermodynamics_map, r)                   
                tdata = thermodynamics_map(r);
                rev = tdata.rev; 
                
                if (rev == -1)
                    products = reaction_data.educts;
                    educts = reaction_data.products;
                    stoich_p = reaction_data.stoich_e;
                    stoich_e = reaction_data.stoich_p; 
                else
    
                    educts = reaction_data.educts;
                    products = reaction_data.products;
                    stoich_e = reaction_data.stoich_e;
                    stoich_p = reaction_data.stoich_p;  
                end
            else
                educts = reaction_data.educts;
                products = reaction_data.products;
                stoich_e = reaction_data.stoich_e;
                stoich_p = reaction_data.stoich_p; 
            end  
        end

        fprintf(fid,'\n%s %s\n', 'reaction', reaction_data.entry);
        num_educts = numel(educts);
        num_products = numel(products);

        fprintf(fid,'\n%s\n', 'substrates:');       
                
        for num = 1:num_educts
        fprintf(fid, '\t%s %s\n', stoich_e{num}, educts{num});  
           if isKey(reactand_map, educts{num}) 
               reactand_map(educts{num}) = reactand_map(educts{num}) - str2double(stoich_e{num}); 
           else
               reactand_map(educts{num}) = - str2double(stoich_e{num});
           end
        end

        fprintf(fid,'\n%s\n', 'products:');
        for num = 1:num_products
           fprintf(fid, '\t%s %s\n', stoich_p{num}, products{num});  
           if isKey(reactand_map, products{num})
               reactand_map(products{num}) = reactand_map(products{num}) + str2double(stoich_p{num});  
           else
               reactand_map(products{num}) = str2double(stoich_p{num});
           end
        end  
    end

    fprintf(fid,'\n%s\n', 'overall:');
    reactand_keys = keys(reactand_map);
    pathway_balance.overall_stoichiometry = cell(numel(reactand_keys),1);

    for k = 1:numel(reactand_keys)
        p = reactand_keys(k);
        fprintf(fid, '%s\t%d\n', p{:}, reactand_map(p{:}));
        pathway_balance.overall_stoichiometry{k,1} = {p{:}, reactand_map(p{:})};
    end

    fprintf(fid,'\n%s\n', 'not balanced:');
    for k = 1:numel(reactand_keys)
        p = reactand_keys(k);
        if reactand_map(p{:}) ~= 0      
        fprintf(fid, '%s\t%d\n', p{:}, reactand_map(p{:}));
        end

    end
        num_reactions = numel(reaclist);
        for id = 1:num_reactions
            r = reactions(reaclist(id));
        end

    fclose(fid);

end