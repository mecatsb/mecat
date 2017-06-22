%  function model = readDataFromSBML(sbml_file, del_comp, reactionList, metaboliteList, equivalence_list, reaction_id_map, rpair_map)
function readDataFromSBML(model)

% sbml_file: SBML model of network
% reactionList: list of reaction abbreviations
% metaboliteList: list of metabolite abbreviations

%     if (nargin < 2)
%         del_comp = false;
%     end

    % read SBML Model
%     model = readCbModel(sbml_file);

    % read stoichiometric matrix from SBML file
    S = model.S;
    
    % read metabolites from SBML file
    metabolites = model.mets;
    metabolite_names = model.metNames;
    
%     if del_comp  
%             token = '[[pce]\]';
%             metabolites =  strtrim(regexprep(metabolites, token, ''));
%     end
    
    % and save them to file
    fid = fopen('compounds.txt','w+');
    cnid = fopen('compound_names.txt','w+');
    
    % metabolite abbreviation -> metabolite index
    metabolite_to_index_map = containers.Map();
    
    %don't use parfor!!!
    for ind = 1:numel(metabolites)
        fprintf(fid, char(metabolites(ind))); 
        fprintf(fid, '\n');
        
        fprintf(cnid, char(metabolite_names(ind))); 
        fprintf(cnid, '\n');
        metabolite_to_index_map(char(metabolites(ind))) = ind;
    end
    fclose(fid);

    % read reactions from SBML file
    reactions = model.rxns; 
    reaction_names = model.rxnNames;
    
    % also add reversible reactions
    r = model.rev;
    
    % reversibilities in network are 0: reversible, 1: irreversible from left to
    % right and -1: irreversible from right to left

    if ~isempty(find(r == -1, 1))
        rev = find(r == 0);
        reversed_reactions = find(r==-1);
        reversible = reactions(rev);
        reversible_names = reaction_names(rev);     
    else
        r = logical(r);
        reversible = reactions(r);
        rev = find(r);
        reversed_reactions = [];
        reversible_names = reaction_names(r);
    end
 
    % reaction index -> reaction abbreviation
    r_map = containers.Map('KeyType','uint32','ValueType','any');
%     r_names_map = containers.Map('KeyType','uint32','ValueType','any');
    
%   r_abbrev_map = containers.Map();
     
    % and save them to file
    fid = fopen('reactions.txt','w+'); 
    rnid = fopen('reaction_names.txt','w+'); 
    for ind = 1:numel(reactions)
        fprintf(fid, char(reactions(ind))); 
        fprintf(fid, '\n'); 
        
        fprintf(rnid, char(reaction_names(ind))); 
        fprintf(rnid, '\n'); 
        
        r_map(ind) = char(reactions(ind));
%         r_names_map = char(reaction_names(ind));
%         r_abbrev_map(char(reactions(ind))) = ind;
    end
    
    % also save reversible reactions
    for ind = 1:numel(reversible)
        fprintf(fid, strcat('-',char(reversible(ind)))); 
        fprintf(fid, '\n');
        
        fprintf(rnid, strcat('- ',char(reversible_names(ind)))); 
        fprintf(rnid, '\n');
        % reaction index -> reaction abbreviation
        r_map(numel(reactions)+ind) = strcat('-',char(reversible(ind)));
    end
    fclose(fid);
    fclose(cnid);
    
    % add reversible reactions to stoichiometric matrix
    S = full(S);
    fid = fopen('Rev.txt', 'w+');
    for ind = 1:numel(rev)
        fprintf(fid, '%d\t%d\n', [numel(reactions)+ind'; rev(ind)']); 
        S(:, numel(reactions)+ind) = -S(:, rev(ind));
    end 
    
    % revert reactions that are irreversible from right to left
    for ind = 1:numel(reversed_reactions)
        S(:, reversed_reactions(ind)) = -S(:, reversed_reactions(ind));
    end
    
    S = sparse(S);
    save('S.mat', 'S')
    
    ext = regexp(metabolites, '\[e\]','match');
    external = metabolites(~cellfun('isempty', ext));
    
    id = fopen('ext.txt','w+');
    for k = 1:numel(external)
        fprintf(id, '%d\n', metabolite_to_index_map(char(external(k))));         
    end
    fclose(id); 
    
%     % metabolite abbreviation -> metabolite KEGG ID
%     metab_abbrev_to_kid_map = containers.Map(metabolite_data{1,1}, metabolite_data{1,6}); 

%     reaction_metabolite_map = containers.Map();
%     reaction_ids = keys(reaction_id_map);
%     for k = reaction_ids
%         reaction = reaction_id_map(k{:})        
%         reactants.educts = reaction.educts;
%         reactants.products = reaction.products;            
%         reaction_metabolite_map(k{:}) = reactants;
%     end
    
%    createBIGGReactionMap(reactionList, 'Rev.txt', r_abbrev_map);
       
% %     %create CEAs
%       equival_map = readEquivalenceList(equivalence_list)   
%       buildCEAFromBIGG(reactionList, metaboliteList, reaction_id_map, rpair_map, metabolite_to_index_map, r_map, equival_map, del_comp);
end