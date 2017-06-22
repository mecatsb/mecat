function buildModelFromKEGG(reaction_map, compound_map, rpair_map, thermodynamics_raw_data, exclude_reactions)
%   buildModelFromKEGG build a SBML model from KEGG data
%   reaction_map: parsed data from KEGG REACTION
%   compound_map:  parsed data from KEGG COMPOUND
%   rpair_map: parsed data from KEGG RPAIR
%   exclude_reactions: list of KEGG REACTION Ids that should be excluded

    if (nargin < 5)
        exclude_reaction = [];
    else
        exclude_reaction = unique(importdata(exclude_reactions));
    end
    n_reactions = findReactionWithNs(reaction_map);
    generic_reactions = findGenericReaction(reaction_map);
    illformed_reactions =  findReactionWithSameEductProduct(reaction_map);
    ex = [];
    
    disp('build model')
    model = KEGGtoCOBRAModel(reaction_map, compound_map, [exclude_reaction; n_reactions; generic_reactions; illformed_reactions;ex]);

    disp('update reaction reversibilities')
    readThermodynamicsData(thermodynamics_raw_data, -15.0, 'reversibilities.txt');
    
    for k = 1:numel(model.rev)
        model.rev(k) = 0;
    end
    %%
    
    fid = fopen('reversibilities.txt', 'r');
    data = textscan(fid, '%s%d', 'delimiter', '\t', 'whitespace', '');
    fclose(fid);  
    
    revMap = containers.Map(data{1,1}, data{1,2});   
    save('reversibilities_map', 'revMap');

    %%
    model = setReversibilities(model, revMap);
    
    disp('write model to file')
    writeCbToSBML(model,['KEGG_' date])
   
    disp('extract model data')
    readDataFromSBML(model);
    
    disp('create arcs')
    ArcsFromKEGG('compounds.txt', 'reactions.txt', 'reversibilities.txt', reaction_map, rpair_map, ['RxR_' date() '.txt']);
    
    buildMetabolitePool(model,['RxR_' date() '.txt'], 'compounds.txt', 'compound_names.txt', compound_map, 0, 300);

    disp('read cofactors')
    
    %%
    fid = fopen('possible_terminal_metabolites.txt', 'r');
    data = textscan(fid, '%s%s%d', 'delimiter', '\t', 'whitespace', '');
    fclose(fid);   
    metabolite_ids = data{1,1};     
    startind = readStartMetabolites(metabolite_ids, 'compounds.txt', true);    
    fid = fopen('terminal_metabolites.txt', 'w');
    fprintf(fid, '%d\n', startind);
    fclose(fid); 
    
    %%
    fid = fopen('possible_cofactors.txt', 'r');
    data = textscan(fid, '%s%s%d', 'delimiter', '\t', 'whitespace', '');
    fclose(fid); 
    metabolite_ids = data{1,1}; 
    startind = readStartMetabolites(metabolite_ids, 'compounds.txt', true);    
    fid = fopen('cofactors.txt', 'w');
    fprintf(fid, '%d\n', startind);
    fclose(fid); 
    
    %%
    
    fid = fopen('reactions_eco.txt', 'r');
    data = textscan(fid, '%s', 'delimiter', '\t', 'whitespace', '');
    fclose(fid); 
    
    reaction_ids = data{1,1}; 
    startind = readStartMetabolites(reaction_ids, 'reactions.txt', false);    
    fid = fopen('reactions_eco_indices.txt', 'w');
    fprintf(fid, '%d\n', startind);
    fclose(fid);

end