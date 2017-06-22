function startind = readMetabolites(file, compounds, suffix)
     metabolite_ids = file;  
    startind = zeros(length(metabolite_ids),1);
    
    if (suffix == true)
        suffix = '[c]';
    else
        suffix = '';
    end
    
    for i = 1 : length(metabolite_ids)
        metabolite = [metabolite_ids{i} suffix];
        compTmp = find(ismember(compounds, metabolite));
        if isempty(compTmp)
             disp(['Metabolite ' metabolite ' not found. Check compounds.txt'])
            continue;
        end    
        startind(i) = compTmp;
    end
    
    startind = startind(find(startind));
    
end