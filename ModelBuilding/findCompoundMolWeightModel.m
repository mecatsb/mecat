function [selected_compounds, compound_data] = findCompoundMolWeightModel(compound_map, compound_file, min, max)

    compounds = findCompoundMolWeight(compound_map, min, max);
    compounds = compounds(:,1);   
    [generic, ~] = findGenericCompound(compound_map, compound_file);
    
    fid = fopen(compound_file);
    compound_names = textscan(fid, '%s', 'delimiter', '\n');
    fclose(fid);
    compound_names = compound_names{:};
    compound_names = regexprep(compound_names, '\[c\]', '');
    compounds_in_model = intersect(compounds, compound_names); 
    
    % keine generischen Metabolite
    not_generic = setdiff(compounds_in_model, generic);     
    compound_data = cell(numel(not_generic),4);
    
    for ind = 1:numel(not_generic)
        data = compound_map(not_generic{ind});
        
        entry = data.entry;
        name = '';
        molweight = '';
        cas = '';

        if isfield(data, 'names')
            name = data.names(1);
        end

        if isfield(data, 'molweight')
            molweight = data.molweight;
        end

        if isfield(data, 'cas')
            cas = data.cas;
        end

        compound_data(ind,1) = entry;
        compound_data(ind,2) = name;
        compound_data{ind,3} = str2double(molweight);
        compound_data(ind,4) = {cas};       
    end

    compound_data = sortrows(compound_data,3);
    
    selected_compounds = compound_data(:,1);
  
end