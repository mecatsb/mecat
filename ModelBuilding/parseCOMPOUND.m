function  [entry, compound] = parseCOMPOUND(entry)

    fieldkeys = {'ENTRY', 'NAME', 'FORMULA', 'EXACT_MASS', 'MOL_WEIGHT', ...
        'SEQUENCE', 'REMARK', 'COMMENT', 'REACTION', 'PATHWAY', 'BRITE' ...
        'ENZYME', 'PATHWAY', 'OTHER DBS', 'LINKDB', 'KCF DATA', 'MODULE', ...
        'STRUCTURE', 'ATOM', 'BOND', 'DBLINKS'};

    pattern = strjoin(fieldkeys, '|');
    % tmp = textscan(entry, '%s', 'delimiter', '\n', 'whitespace', '');
    % tmp = tmp{:};

    [data, keys] = regexp(entry, pattern, 'split', 'match');
    data(cellfun(@isempty,data)) = [];

    if ~numel(data) == numel(keys)
        error('Parsing error');
    end

    for k = 1:numel(keys)
        switch keys{k}         
           case 'ENTRY'
                compound.entry = regexp(data{k}, 'C\d{5}', 'match');
                
           case 'NAME'
              compound.names = strtrim(strsplit(data{k},';'));   
              
            case 'FORMULA'
                compound.formula = strtrim(data{k});
                
            case 'REACTION'
               reactions = regexp(entry, 'R\d{5}', 'match');
               compound.reactions = reactions;
            
            case 'COMMENT'
                compound.comment = strtrim(data{k});
                
            case 'MOL_WEIGHT'
                compound.molweight = strtrim(data{k});
              
            case 'DBLINKS'
                  tmp = regexp(data{k}, 'CAS:\s*\d*-\d*-\d*', 'match');
                  if ~isempty(tmp)
                    compound.cas = regexp(tmp{:}, '\d*-\d*-\d*', 'match');
                  end                               
           otherwise
             continue;
        end   
    end
    
    entry = compound.entry;

% % find ENTRY
% e = strfind(tmp, 'ENTRY');
% ind = find(~cellfun('isempty', e), 1);
% entry = tmp(ind);
% % entry = strtrim(regexprep(entry,'ENTRY',''));
% % entry = strtrim(regexprep(entry,'Compound',''));
% entry = regexp(entry, 'C\d{5}', 'match');
% entry = entry{:};
% entry = entry{:};
% compound.entry = entry;
% 
% % line = strfind(tmp{1}, 'NAME');
% % ind = find(~cellfun('isempty', line), 1);
% % name = tmp{1}(ind);
% % name = strtrim(regexprep(name,'NAME',''));
% % name = strtrim(regexprep(name,';',''));
% % name = name{:};   
% % compound.name = name;
% 
% 
% entry_start = strfind(tmp, 'NAME');
% starti = find(~cellfun('isempty', entry_start), 1)
% entry_end = regexp(tmp, pattern);
% 
% 
% entry_end = strfind(tmp{1}, 'FORMULA')
% endi = find(~cellfun('isempty', entry_end), 1);
% 
% if (~isempty(starti) && ~isempty(endi))
%          name = tmp(starti:endi-1)
%          name{1} = regexprep(name{1},'^NAME','');
%          name = strtrim(regexprep(name,';',''));
%          name_lines = strtrim(name);
%          name = strtrim(name_lines);
%          compound.name = name;
% else
%     compound.name = {''};
% end
% 
% % % find comment(reads currently only the first line of comments)
% % line = strfind(tmp{1}, 'COMMENT');
% % ind = find(~cellfun('isempty', line), 1);
% % 
% % if ~isempty(ind)
% %     comment = tmp{1}(ind);
% %     comment = strtrim(regexprep(comment,'COMMENT',''));
% %     comment  = comment{:};   
% %     compound.comment = comment;
% % else
% %     compound.comment = '';
% % end
% 
% entry_start = strfind(tmp{1}, 'COMMENT');
% starti = find(~cellfun('isempty', entry_start), 1);
% entry_end = strfind(tmp{1}, 'REACTION');
% 
% endi = find(~cellfun('isempty', entry_end), 1);
% 
% if (~isempty(starti) && ~isempty(endi))
%          comment = tmp(starti:endi-1);
%          comment{1} = regexprep(comment{1},'^COMMENT','');
%          comment = strtrim(regexprep(comment,';',''));
%          comment_lines = strtrim(comment);
%          comment = strtrim(comment_lines);
%          compound.comment = comment;
% else
%     compound.comment = {''};
% end
% 
% % find formula
% line = strfind(tmp{1}, 'FORMULA');
% ind = find(~cellfun('isempty', line), 1);
% if ~isempty(ind)
%     formula = tmp{1}(ind);
%     formula = strtrim(regexprep(formula,'FORMULA',''));
%     formula  = formula{:};   
%     compound.formula = formula;
% else
%     compound.formula = '';
% end
% 
% % find other links
% entry_start = strfind(tmp, 'DBLINKS');
% starti = find(~cellfun('isempty', entry_start), 1);
% entry_end = strfind(tmp, 'ATOM');
% endi = find(~cellfun('isempty', entry_end), 1);
%   if (~isempty(starti) && ~isempty(endi))
%         links = tmp(starti:endi-1);
%         links{1} = regexprep(links{1},'^DBLINKS','');
%         links_lines = strtrim(links);
%         pubchem_ind = strfind(links_lines, 'PubChem');
%         pubchem_ind = find(~cellfun('isempty', pubchem_ind), 1);
%         pubchem = links_lines(pubchem_ind);
%         pubchem = strtrim(regexprep(pubchem,'PubChem:',''));
%         compound.pubchem = pubchem;
%         
%         chebi_ind = strfind(links_lines, 'ChEBI');
%         chebi_ind = find(~cellfun('isempty', chebi_ind), 1);
%         chebi = links_lines(chebi_ind);
%         chebi = strtrim(regexprep(chebi,'ChEBI:',''));
%         compound.chebi = chebi;
%           
%   else
%     compound.pubchem = '';
%     compound.chebi = '';   
%   end
end
