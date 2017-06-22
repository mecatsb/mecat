function compounds = findCompoundsWithXCarbonAtoms(model, min, max)
    metabolites = model.mets;
    formulas = model.metFormulas;
    compounds = {};
    
    % capture number of C atoms
%     num_carbons = regexp(formulas, 'C(\d*)', 'tokens');
    num_carbons = regexp(formulas, 'C((?![a-z])\d*)', 'tokens');
    
    % get indices of metabolites without carbons
    empty = cellfun('isempty', num_carbons);
        
    % get indices of metabolites with at least one carbon
    not_empty_indices = find(~cellfun('isempty', num_carbons));   

    for k = not_empty_indices' 
        n = num_carbons{k}; 
            
        % more than one block of C atoms (mostly for polymers)
        % ignore polymers
        if numel(n) > 1
            continue
        end
        
        n = [n{:}];
        
        if isempty(n{:})
            n = 1;            
        else
           n = n{:};
           n = str2num(n);
        end
           
        if (n < min) || (n > max)
            continue
        end
          
        compounds{end+1} = metabolites{k};
%          end
    end 

    % write metabolites without C atoms
    if min == 0      
        compounds = [compounds'; metabolites(empty)];
    end
end