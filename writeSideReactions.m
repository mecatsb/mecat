function writeSideReactions(reaclist, reactions, position, reac_names, reaction_map, compound_map, fid) 

    reactand_map = containers.Map();
    reactant_name_map = containers.Map();
    reactand_molweight_map = containers.Map();
    reactand_formula_map = containers.Map();
    
    reaclist = reaclist(find(reaclist));
    
    fprintf(fid, '%s%d\n', 'position: ', position);
    
    for t = 1:length(reaclist) 
        r = strjoin(reactions(reaclist(t)),', ');    
        fprintf(fid, '\n%s%s%s\n', r, ' : ', reac_names{reaclist(t)});
             
        r = reactions(reaclist(t));
        r = r{:};
        if ~isempty(strfind(r, '-'))
            r = regexprep(r, '-', '');
            reaction_data = reaction_map(r);
            products = reaction_data.educts;
            educts = reaction_data.products;
            stoich_p = reaction_data.stoich_e;
            stoich_e = reaction_data.stoich_p;
        else
           reaction_data = reaction_map(r);
                  
            educts = reaction_data.educts;
            products = reaction_data.products;
            stoich_e = reaction_data.stoich_e;
            stoich_p = reaction_data.stoich_p;                  
        end
        num_educts = numel(educts);
        num_products = numel(products);

        fprintf(fid,'\n\t%s\n', 'substrates:');       
        ed = cell(num_educts,1);
        stoich = zeros(num_educts,1);
        for num = 1:num_educts
           data = compound_map(educts{num});
           reactant_name_map(educts{num}) = data.names{1};
           
           if isfield(data, 'molweight')
                reactand_molweight_map(educts{num}) = data.molweight;
           else
                reactand_molweight_map(educts{num}) = 'n/a';
           end
           
           if isfield(data, 'formula')
                reactand_formula_map(educts{num}) = data.formula;
           else
                reactand_formula_map(educts{num}) = 'n/a';
           end
           
           if isKey(reactand_map, educts{num}) 
               reactand_map(educts{num}) = reactand_map(educts{num}) - str2double(stoich_e{num}); 
           else
               reactand_map(educts{num}) = - str2double(stoich_e{num});
           end
           ed{num,1} = educts{num};
           stoich(num,1) = str2double(stoich_e{num});          
        end
        
        [sorted_stoich, sorted_indices] = sort(stoich);
        sorted_ed = ed(sorted_indices);
        for num = 1:num_educts       
            fprintf(fid, '\t\t%d\t%s\t%s\n', sorted_stoich(num), sorted_ed{num}, reactant_name_map(educts{num}));
        end

        fprintf(fid,'\n\t%s\n', 'products:');
        prod = cell(num_products,1);
        stoich = zeros(num_products,1);
        for num = 1:num_products
               data = compound_map(products{num});
               reactant_name_map(products{num}) = data.names{1};
               
               if isfield(data, 'molweight')
                    reactand_molweight_map(products{num}) = data.molweight;
               else
                    reactand_molweight_map(products{num}) = 'n/a';
               end
               
               if isfield(data, 'formula')
                    reactand_formula_map(products{num}) = data.formula;
               else
                    reactand_formula_map(products{num}) = 'n/a';
               end

           if isKey(reactand_map, products{num})
               reactand_map(products{num}) = reactand_map(products{num}) + str2double(stoich_p{num});  
           else
               reactand_map(products{num}) = str2double(stoich_p{num});
           end
           prod{num,1} = products{num};
           stoich(num,1) = str2double(stoich_p{num});
        end  
        
        [sorted_stoich, sorted_indices] = sort(stoich);
        sorted_prod = prod(sorted_indices);
        
        for num = 1:num_products        
            fprintf(fid, '\t\t%d\t%s\t%s\n', sorted_stoich(num), sorted_prod{num}, reactant_name_map(products{num}));
        end
    end
 	fprintf(fid,'\n%s\n\n', '#######################################################################################################');

end