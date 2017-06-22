function model = KEGGtoSBML(reactionlist, reactions, reaction_map, compound_map, revMap)

    reactionlist = reactionlist(find(reactionlist));
    reactionlist = unique(reactionlist);
    
    reaction_ids = cell(numel(reactionlist,1));
    reaction_names = cell(numel(reactionlist,1));
    reaction_equations = cell(numel(reactionlist,1));
    
    del_token  = '(-->|<==>|<=>)';
    
    for ind = 1:numel(reactionlist)
        old_id = reactions(reactionlist(ind));
        old_id = old_id{:};
        id = regexprep(old_id,'-','');
        
        data = reaction_map(id);
        reaction_names{ind} = '';
        if (isfield(data, 'names'))          
            reaction_names{ind} = data.names{1};
        end
        
        equation = data.equation;
        
        if ((isKey(revMap, id) && revMap(id) == -1) || ~strcmp(old_id,id))
           del = regexp(equation,del_token,'match');
           split = regexp(equation, del{:}, 'split');
           split = strtrim(split);         
           equation = strcat(split{1,2}, {' '}, '=>',  {' '}, split{1,1});
           equation = equation{:};
%            revMap(id) = 1;
           reaction_names{ind} = strcat('- ', reaction_names{ind});          
        end   
        
        if (isKey(revMap, id) && revMap(id)== 0 || ~isKey(revMap, id))
            revMap(id) = 1; 
        end
        
        reaction_ids{ind} = old_id;
        reaction_equations{ind} = equation;
    end
    
    model = createModel(reaction_ids, reaction_names, reaction_equations); 
    
    for ind = 1:numel(model.metNames)
        id = model.metNames{ind};
        id = regexp(id, 'C\d{5}', 'match');
        id = id{:};
%         id = regexprep(id, '\[.\]', '');
        
        data = compound_map(id);
        name = '';
        if (isfield(data, 'names'))          
            name = data.names{1};
        end
        
        model.metNames{ind} = name;     
    end
    
%     ########################################

    % reversibilities in network are 0: reversible, 1: irreversible from left to
    % right and -1: irreversible from right to left
    r = model.rev;
    if ~isempty(find(r == -1, 1))
        reversed_reactions = find(r==-1);   
    else
        reversed_reactions = [];
    end
    
    for ind = 1:numel(reversed_reactions)
        model.S(:, reversed_reactions(ind)) = -model.S(:, reversed_reactions(ind));
        revMap(reversed_reactions(ind)) = 1; 
    end
    
        model = setReversibilities(model, revMap);
end