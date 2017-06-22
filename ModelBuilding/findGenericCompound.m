%generic_compound zurückgeben, wenn man alle in KEGG will, nicht nur die
%im Model
function [generic_model, generic_kegg] = findGenericCompound(compound_map, compound_file)

    fid = fopen(compound_file);
    compound_names = textscan(fid, '%s', 'delimiter', '\n');
    fclose(fid);
    compound_names = compound_names{:};

    compound_names = regexprep(compound_names, '\[c\]', '');

    generic_kegg = [];
    compound_ids = keys(compound_map);
    
   for ind = 1:numel(compound_ids)
        id = compound_ids(ind);
        data = compound_map(id{:});
        entry = data.entry;
        entry = entry{:};
        
        if isfield(data, 'comment')
            comment = data.comment;       
            if (~isempty(strfind(comment, 'generic')) || ~isempty(strfind(comment, 'Generic')))
                generic_kegg = [generic_kegg; entry];
                continue;
            end      
        end
        
        if isfield(data, 'formula')
            formula = data.formula;
            if ~isempty(strfind(formula, 'R'))
                generic_kegg = [generic_kegg; entry];
            end
        end     
   end 
   
   generic_model = intersect(compound_names, generic_kegg);
   
end
