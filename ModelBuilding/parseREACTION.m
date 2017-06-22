function [equation, reaction] = parseREACTION(entry)

    fieldkeys = {'ENTRY', 'NAME', 'DEFINITION', 'EQUATION', 'REMARK', 'COMMENT'...
    'RPAIR', 'ENZYME', 'PATHWAY', 'MODULE', 'ORTHOLOGY', 'LINKDB'};

    pattern = strjoin(fieldkeys, '|');
    [data, keys] = regexp(entry, pattern, 'split', 'match');
    data(cellfun(@isempty,data)) = [];

    if ~numel(data) == numel(keys)
        error('Parsing error');
    end

    for k = 1:numel(keys)
        switch keys{k}
           case 'ENTRY'
               e = regexp(data{k}, 'R\d{5}', 'match');
               reaction.entry = e{:};
           case 'NAME'
              reaction.names = strtrim(strsplit(data{k},';'));   
              
%             this is different than the rpair parsing before!
            case 'RPAIR'
               reaction.rpairs = regexp(data{k}, '(RP\d{5})\s*(C\d{5}\_C\d{5})\s*(\w*)','tokens');

            case 'EQUATION'
                equation = strtrim(data{k});
                equation = strtrim(regexprep(equation, '///',''));
                reaction.equation = equation;
                
            case 'ENZYME'
                reaction.enzyme = regexp(data{k}, '\d[.]\d+[.]\d+[.]\d+', 'match');  
             
            case 'COMMENT'
                reaction.comment = strtrim(data{k});
                
           otherwise
             continue;
        end
        
    end 
% parse equation
[split, ~] = splitReaction({equation});

educts = strtrim(regexp(split{1},'(C\d{5}|G\d{5})','match'));
products = strtrim(regexp(split{2},'(C\d{5}|G\d{5})','match'));
educts = [educts{:}];
products = [products{:}];

stoich_e = split{1};
stoich_p = split{2};

parfor k = 1:length(stoich_e)
    f = strfind(stoich_e{k}, '(')
    if (~isempty(f))
        stoich_e{k} = strtrim(regexp(stoich_e{k}, '+(?=\))', 'split'));  
    else
         stoich_e{k} = strtrim(regexp(stoich_e{k}, '+', 'split'));
    end
    stoich_e{k} = strtrim(regexprep(stoich_e{k},'(C\d{5}|G\d{5})',''));
end

parfor k = 1:length(stoich_p)
    f = strfind(stoich_p{k}, '(')
    if (~isempty(f))
        stoich_p{k} = strtrim(regexp(stoich_p{k}, '+(?=\))', 'split'));  
    else
         stoich_p{k} = strtrim(regexp(stoich_p{k}, '+', 'split'));
    end
    stoich_p{k} = strtrim(regexprep(stoich_p{k},'(C\d{5}|G\d{5})',''));
end

stoich_e = [stoich_e{:}];
stoich_p = [stoich_p{:}];

parfor k = 1:length(stoich_e)
    if(isempty(stoich_e{k}))
        stoich_e{k} = '1';
    end
end

parfor k = 1:length(stoich_p)
    if(isempty(stoich_p{k}))
        stoich_p{k} = '1';
    end
end

reaction.educts = educts;
reaction.products = products;
reaction.stoich_e = stoich_e;
reaction.stoich_p = stoich_p;
end