function writeCytoscapeDataTable(model, path_reactionlist, path_metabolitelist, active_reactionlist, reactions, metabolites, filename)
    
    values = cell(numel(reactions),1);
    val2 = cell(numel(metabolites),1);
    
    
    for k = 1:numel(active_reactionlist)      
        r = reactions(active_reactionlist(k));
    values = node_colors(values, r , model.rxns, 'Green');
    end
    
    path_reactionlist = path_reactionlist(find(path_reactionlist));
    path_metabolitelist = path_metabolitelist(find(path_metabolitelist));
    
    for k = 1:numel(path_reactionlist)      
         r = reactions(path_reactionlist(k));     
        values = node_colors(values, r , model.rxns, 'Red');
        
    end
    
    val2 = node_colors(val2, metabolites(path_metabolitelist(1)) , model.mets, 'Yellow');     
    val2 = node_colors(val2, metabolites(path_metabolitelist(numel(path_metabolitelist))) , model.mets, 'Magenta');
    
    for l = 2:numel(path_metabolitelist)-1
        m = metabolites(path_metabolitelist(l));
        val2 = node_colors(val2, m , model.mets, 'Blue');
    end

    fs = '%s\t%s\n';
    fid = fopen([filename '.txt'], 'A');
    fprintf(fid,fs, 'ID', 'node colors');
    for k = 1:numel(model.rxns)  
         fprintf(fid,fs, [model.rxnNames{k}], values{k});
    end
    
    for k = 1:numel(model.mets)  
        fprintf(fid,fs, [model.metNames{k}], val2{k});
    end
    fclose(fid);
end

function values = node_colors(values, id, model_list,  color)

        IndexC = strfind(model_list, id{:});
        Index = not(cellfun('isempty', IndexC)); 
        values{Index,1} = color;
end