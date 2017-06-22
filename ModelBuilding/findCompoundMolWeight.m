function compounds = findCompoundMolWeight(compound_map, min, max)
   compound_ids = keys(compound_map);
   compounds = cell(numel(compound_ids),2);
   for ind = 1:numel(compound_ids)
        id = compound_ids(ind);
        data = compound_map(id{:});
        entry = data.entry;
        entry = entry{:};       
        mol_weight = 0;
        
        if isfield(data, 'molweight')
            mol_weight = data.molweight;
            mol_weight = str2double(mol_weight);            
        end
        
       if mol_weight > min && mol_weight < max         
            compounds{ind,1} = entry;
            compounds{ind,2} = mol_weight;
       end
   end
   
   compounds = compounds.';
   compounds = reshape(compounds(~cellfun(@isempty,compounds)),2,[])';   
   c = sortrows(compounds,2);
   compounds = c;
end