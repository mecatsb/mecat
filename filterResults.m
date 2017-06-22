function results = filterResults(input_file, results)

    input_data = parseInput(input_file);
    old_folder = cd(input_data.data_path);
    excluded_indices = unique(importdata(input_data.excluded_metabolites));
    cofactors = importdata(input_data.cofactors);
    start_metabolites = importdata(input_data.start_metabolites);
    reactions = importdata(input_data.reactions);
    compounds = importdata(input_data.compounds);
    
    cd(old_folder);
    old_folder = cd(input_data.kegg_data_path);
    
    reaction_map = load(input_data.reaction_map);
    reaction_map = reaction_map.reaction_map;   
    
    cd(old_folder)
    mkdir(input_data.result_path);
    cd(input_data.result_path);   
    
    metabolite_pool = unique([cofactors; start_metabolites; excluded_indices]);
    start_reactions = findStartReactions(metabolite_pool, compounds, reactions, reaction_map);
    
    
    maxz = max(cellfun('length',results.z));
    zs = zeros(size(results.z,1),maxz);
    
    for k = 1:numel(results.z)
        zs(k,1:numel(results.z{k})) = results.z{k}';
    end
    
    [~,ia,~] = unique(zs, 'rows');   
    double =  setdiff(1:size(results.z,1),ia);
   
    results.metabolite_lists(double) = {[]};
    results.lengths(double) = 0;
    results.x(double) = {[]};
    results.v(double) = {[]};
    results.active_reactions(double) = {[]};
    results.number_active_reactions(double) = 0;
    results.reaction_lists(double) = {[]};
    results.z(double) = {[]};
    results.side_reactions(double) = {[]};
     
    results.reaction_lists = setReactionLists(results.reaction_lists, start_reactions);

    empty_ind= cellfun(@isempty,results.reaction_lists);
    results.metabolite_lists(empty_ind) = {[]};
    results.lengths(empty_ind) = 0;
    results.x(empty_ind) = {[]};
    results.v(empty_ind) = {[]};
    results.active_reactions(empty_ind) = {[]};
    results.number_active_reactions(empty_ind) = 0;
    results.z(empty_ind) = {[]};
    results.side_reactions(empty_ind) = {[]};

    save('results_filtered', 'results'); 
    
end