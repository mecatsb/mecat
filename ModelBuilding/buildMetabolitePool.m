function buildMetabolitePool(model, cea_file, compound_file, compound_name_file, compound_map, min, max)

    fid = fopen(compound_file);
    compounds = textscan(fid, '%s', 'delimiter', '\n');
    fclose(fid);
    compounds = compounds{:};
    
    fid = fopen(compound_name_file);
    compound_names = textscan(fid, '%s', 'delimiter', '\n');
    fclose(fid);
    compound_names = compound_names{:};

    [selected_compounds, ~] = findCompoundMolWeightModel(compound_map, compound_file, min, max);   
    [arc_compounds, wo_arcs] = findCompoundsWithArcs(cea_file, compound_file);
    
    [generic_compounds, ~] = findGenericCompound(compound_map, compound_file);

    wo_carbons = findCompoundsWithXCarbonAtoms(model, 0, 1);
    wo_carbons = regexprep(wo_carbons,'\[c\]','');
    
    start_candidates = setdiff(selected_compounds, wo_carbons);     
    start_metabolites = intersect(start_candidates, arc_compounds);
       
    metabolite_pool   = intersect([wo_arcs; wo_carbons], selected_compounds);  
    metabolite_pool = setdiff(metabolite_pool, start_metabolites); %nicht unbedingt nötig
          
    sm_ind = readMetabolites(start_metabolites, compounds, true);
    pool_ind = readMetabolites(metabolite_pool, compounds, true);
    generic_ind = readMetabolites(generic_compounds, compounds, true);
    
    wo_arcs_ind = readMetabolites(wo_arcs, compounds, true);
    arcs_ind = readMetabolites(arc_compounds, compounds, true);
    
    not_pool_indices = setdiff((1:numel(compounds))', pool_ind);    
    not_pool_wo_generic_indices = setdiff(not_pool_indices, generic_ind);
    not_pool_wo_arcs = intersect(not_pool_wo_generic_indices, wo_arcs_ind);
    not_pool_arcs = intersect(not_pool_wo_generic_indices, arcs_ind);

    fid = fopen('startpoints.txt', 'w');
        fprintf(fid, '%d\n', sm_ind);
    fclose(fid);
    
    fid = fopen('excluded_metabolites.txt', 'w');
        fprintf(fid, '%d\n', pool_ind);
    fclose(fid);  
    
    fid = fopen('generic_compounds.txt', 'w');
        fprintf(fid, '%d\n', generic_ind);
    fclose(fid);
       
    writeList('generic_names.txt', generic_ind, compounds, compound_names);
    writeList('startpoints_names.txt', sm_ind, compounds, compound_names);
    writeList('excluded_metabolites_names.txt', pool_ind, compounds, compound_names);
    writeList('not_pool_names.txt', not_pool_indices, compounds, compound_names);
    writeList('not_pool_wo_generic_indices_names.txt', not_pool_wo_generic_indices, compounds, compound_names);
    writeList('not_pool_wo_arcs_names.txt', not_pool_wo_arcs, compounds, compound_names);
    writeList('not_pool_arcs_names.txt', not_pool_arcs, compounds, compound_names);
   
end

function writeList(filename, indices, compounds, compound_names)
    
    fid = fopen(filename', 'w');   
    for ind = 1:numel(indices)
        index = indices(ind);
        id = compounds{index};
        name = compound_names{index};
        fprintf(fid, '%s\t%s\n', id, name);
    end
    fclose(fid);
end