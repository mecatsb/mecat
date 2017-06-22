function model = KEGGtoCOBRAModel(reaction_map, compound_map, exclude_reaction)
% KEGGtoCOBRAModel create a COBRA model from KEGG data
% reaction_map: map with parsed KEGG data (KEGG Reaction id -> parsed data
% struct
% filename: output name of model
% (optional) exclude_reaction_file: list of KEGG reactions ids to exclude

%     if (nargin < 3)
%         exclude_reaction = [];
%     else
%         exclude_reaction = unique(importdata(exclude_reaction_file));
%     end

    reaction_ids = keys(reaction_map);
    abbreviations = cell(numel(reaction_ids), 1);
    names = cell(numel(reaction_ids), 1) ;
    equations = cell(numel(reaction_ids), 1);
    
    fid = fopen('exclude_reactions_new.txt', 'w+');
    
    for ri = 1:numel(reaction_ids)
        abbrev = reaction_ids(ri);
        abbrev = abbrev{:};
        data = reaction_map(abbrev);
        add = 1;
        
        % filter reactions containing a glycan
        if KEGGReactioncontainsGlycan(data.equation)   
           add = 0;
           fprintf(fid, '%s\n', abbrev);
        end
        
%         %exclude reactions containing generic metabolites
%         reactionpartners = [data.educts data.products];
%         if ismember(reactionpartners, exclude_compound)
%             add = 0;
%             fprintf(fid, '%s\n', abbrev);
%         end
        
        if ismember(abbrev, exclude_reaction)
            add = 0;
            fprintf(fid, '%s\n', abbrev);
        end
        
        if add
            abbreviations{ri} = abbrev;
            
            if isfield(data, 'names')    
                n = data.names;
                names{ri} = n{1};
            else
                names{ri} = ' ';
            end
            
            equations{ri} = data.equation;
        end       
    end
    fclose(fid);

    nonempty_abbreviations = find(~cellfun(@isempty,abbreviations));
    
    new_abbreviations = cell(numel(nonempty_abbreviations),1);
    new_names = cell(numel(nonempty_abbreviations),1);
    new_equations = cell(numel(nonempty_abbreviations),1);
     
    for k = 1:numel(nonempty_abbreviations)
        new_abbreviations(k) = abbreviations(nonempty_abbreviations(k));
        new_names(k) = names(nonempty_abbreviations(k));
        new_equations(k) = equations(nonempty_abbreviations(k));
    end   
    model = createModel(new_abbreviations, new_names, new_equations);
     
     % search for invalid metabolite names (due to incorrect equation parsing in
     % COBRA function)
      
     wellformed = regexp(model.mets, '(^(C\d{5}\[c\]))|^(G\d{5}\[c\])', 'match');
     illformed = cellfun(@isempty,wellformed);
     del_metab = model.mets(illformed);     
     
     model = removeMetabolites(model, del_metab); 
        
%      f = fopen('duplicates', 'w+');
%      
%      % find duplicate metabolites     
%      [metName,metCnt] = countUnique(model.mets);
%       metInd = find(metCnt > 1);
%         if ~isempty(metInd)
%             for i = 1:length(metInd)
%                 thisMetName = metName{metInd(i)};
%                 fprintf(f, '%s\n',thisMetName);
%             end
%         end
%        fclose(f); 
        
    % set compound information   
    metabolite_names = model.metNames;
    new_metabolite_names = cell(numel(metabolite_names),1);
    metabolite_formulas = model.metFormulas;
    new_metabolite_formulas = cell(numel(metabolite_formulas),1);
       
    for ind = 1:numel(metabolite_names)
        index = regexp(metabolite_names{ind}, 'C\d{5}', 'match');
        index = index{:};
        data = compound_map(index);
        
        if isfield(data, 'names')    
            names = data.names;
            new_metabolite_names{ind} = names;
        else
            new_metabolite_names{ind} = '';
        end
        
        
        new_metabolite_names{ind} = names{1};
        
        if isfield(data, 'formula')    
            formula = data.formula;
            new_metabolite_formulas{ind} = formula;
        else
            new_metabolite_formulas{ind} = '';
        end
        
    end
    
    model.metNames = new_metabolite_names;
    model.metFormulas = new_metabolite_formulas;
        
end