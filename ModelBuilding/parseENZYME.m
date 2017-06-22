function enzyme = parseENZYME(entry)

    fieldkeys = { 'ENTRY', 'NAME', 'CLASS', 'SYSNAME', 'REACTION', ...
        'ALL_REAC', 'SUBSTRATE', 'PRODUCT', 'COMMENT', 'HISTORY', 'PATHWAY'...
        'ORTHOLOGY', 'GENES', 'REFERENCE', 'DBLINKS'};
    
    pattern = strjoin(fieldkeys, '|');
    [data, keys] = regexp(entry, pattern, 'split', 'match');
    data(cellfun(@isempty,data)) = [];

    if ~numel(data) == numel(keys)
        error('Parsing error');
    end

    gene_pattern = '^[A-Z]{3,}(?=(:|\s))';
    
    for k = 1:numel(keys)
        switch keys{k}
           case 'ENTRY'
               e = regexp(data{k}, '\d[.]\d+[.]\d+[.]\d+', 'match');
               enzyme.entry = unique(e);
           case 'NAME'
              names = strsplit(data{k},';');   
              enzyme.names = strtrim(names);
            case 'GENES'                
                genes = strtrim(strsplit(data{k}, '\n'));
                genes(strcmp('',genes)) = [];
                organism = regexp(genes, gene_pattern, 'match');
                gene = regexprep(genes, [gene_pattern ':'], '');
                gene = strtrim(gene);
                gene = regexp(gene,'\s*', 'split'); 
                enzyme.organism = organism;
                enzyme.genes = gene;
                
           otherwise
             continue;
        end
        
    end
end