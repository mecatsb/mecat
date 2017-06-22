function startind = readStartMetabolites(file, compounds, suffix)

%     fid = fopen(file, 'r');
%     % KEGG ID Name
%     data = textscan(fid, '%s%s%d', 'delimiter', '\t', 'whitespace', '');
%     fclose(fid);   
%     metabolite_ids = data{1,1}; 

       metabolite_ids = file;  

    fid = fopen(compounds);
    compounds = textscan(fid, '%s', 'delimiter', '\n');
    fclose(fid);
    compounds = compounds{:};

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