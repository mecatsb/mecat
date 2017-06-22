function rankResults(input_file, results)
     
    % read options and import data
    input_data = parseInput(input_file);
    old_folder = cd(input_data.data_path);

    reactions = importdata(input_data.reactions);    
    compounds = importdata(input_data.compounds);
    cofactors = importdata(input_data.cofactors);
    excluded_metabolites = unique(importdata(input_data.excluded_metabolites));  
    basis_metabolites = importdata(input_data.basis_metabolites);   
    reversible = importdata(input_data.reversible_reactions);
    start_reversible_reactions = min(reversible(:,1));
    reversible_map = containers.Map(reversible(:,1), reversible(:,2));
    reversible_map2 = containers.Map(1:start_reversible_reactions-1, 1:start_reversible_reactions-1);
    reversible_map = [reversible_map; reversible_map2];   
    non_enzymatic_reactions = importdata(input_data.non_enzymatic_reactions);
    thermodynamics_map = load(input_data.thermodynamics_map);
    thermodynamics_map = thermodynamics_map.thermodynamics_map;  
    reversibilities_map = load(input_data.reversibilities_map);
    reversibilities_map = reversibilities_map.revMap;   
    fid = fopen(input_data.reaction_names);
    reaction_names = textscan(fid, '%s', 'delimiter', '\n');
    fclose(fid);
    
    reaction_names = reaction_names{:};
    cd(old_folder);
    old_folder = cd(input_data.kegg_data_path);
    
    reaction_map = load(input_data.reaction_map);
    reaction_map = reaction_map.reaction_map;      
    compound_map = load(input_data.compound_map);
    compound_map = compound_map.compound_map;  
    
    enzyme_map = load(input_data.enzyme_map);
    enzyme_map = enzyme_map.enzyme_map;
   
    % create folder for result files
    cd(old_folder)
    mkdir(input_data.result_path);
    cd(input_data.result_path);
    num = clock;
    folder_name = [input_data.result_target '_' num2str(num(3)) '.' num2str(num(2)) '.' num2str(num(1)) '_' num2str(num(4)) '.' num2str(num(5))];
    mkdir(folder_name)
    cd(folder_name)    
    
    mp_with_basis_metabolites = unique([cofactors; basis_metabolites; excluded_metabolites]);
    start_reactions = findStartReactions(mp_with_basis_metabolites, compounds, reactions, reaction_map); 
    
    max_length = int8(floor(max(results.lengths)));

    all_pathways = [];
    metalist_all = cell(max_length,1);   
    metalist_indices = [];
    
    side_reactions = cell(max_length,1);
    metab_reac = cell(max_length,1);
    side_reaction_indices = [];
    
    v_all =  cell(max_length,1);   
    v_indices = [];

    % go over pathways, expand and concatenate them
    for k = 1:numel(results.lengths)        
        metalist_all{k} = results.metabolite_lists(k);
        v_all{k} = results.v(k);
        reactionlist_all =  results.reaction_lists(k);
        side_reactions{k} = results.side_reactions(k);
        metab_reac{k} = results.metabolite_reaction(k);

        if isempty(reactionlist_all{:})
            continue;
        end
        
         reac_pathways_all = expandPathway(reactionlist_all);
%          reac_pathways_all = reactionlist_all{:}';
        if size(reac_pathways_all,2) < max_length
            reac_pathways_all(:,max_length) = 0;
        end         

        dimensions_all_pathways = size(all_pathways,1);
        dimensions_reac_pathways_all = size(reac_pathways_all,1);

        all_pathways = [all_pathways; reac_pathways_all];   

        metalist_indices(dimensions_all_pathways+1:dimensions_all_pathways+dimensions_reac_pathways_all,1) = k;
        v_indices(dimensions_all_pathways+1:dimensions_all_pathways+dimensions_reac_pathways_all,1) = k;
        side_reaction_indices(dimensions_all_pathways+1:dimensions_all_pathways+dimensions_reac_pathways_all,1) = k;
        
    end 
    
    clear reac_pathways_all dimensions_all_pathways dimensions_reac_pathways_all fid
    
    for s = 1:size(all_pathways,1)
        reactionlist = all_pathways(s, :);
        test = reactionlist;
        test(test==0) = [];
        
        test = arrayfun(@(x) reversible_map(x), test);
        
        if numel(unique(test)) < numel(test)
         all_pathways(s, :) = 0;
        end        
    end
    
   % size depends on number of ranking functions
    pathway_score = zeros(size(all_pathways,1),6);
    scores = cell(size(all_pathways,1),1);
    scores_names = cell(size(all_pathways,1),1);
    deltaGs = cell(size(all_pathways,1),1);
    
    % define organism
     reaction_gene_map =  getEnyzmeInformation (all_pathways, reactions, reaction_map, non_enzymatic_reactions, 'eco', enzyme_map);

    for s = 1:size(all_pathways,1)
        reactionlist = all_pathways(s, :);
               
        if reactionlist == zeros(1,max_length)
            continue;
        end
        
        v = v_all(v_indices(s));
        v = v{:};
        pathway_score(s,1) = numel(find(v{:}));
        pathway_score(s,2) = startsWithbasisMetabolite(reactionlist, start_reactions);
        [pathway_score(s,5), pathway_score(s,4), pathway_score(s,3), deltaGs{s}] = rankThermodynamics(v{:}, thermodynamics_map, reactions);                           
        pathway_score(s,6) = findHeterologousEnzymes(v{:}, reaction_gene_map);
        pathway_score(s,7) = countCofactors(cofactors, compounds, v{:}, reactions, reaction_map); 

        score.num_active = pathway_score(s,1);
        score.basis = pathway_score(s,2);
        score.dG = pathway_score(s,5);
        score.sumdG = pathway_score(s,4);
        score.wodG = pathway_score(s,3);
        score.numhet = pathway_score(s,6);
        score.numcofac = pathway_score(s,7);        
              
        names.num_active = 'number of active reactions: ';
        names.dG = 'dG: ';
        names.sumdG = 'sum (dG + |dG|):';
        names.wodG = 'reactions w/o dG: ';
        names.numhet = 'number of heterologeous enzymes: ';
        names.numcofac = 'number of cofactors: ';
        names.basis = 'starts with basic: ';
        
        scores{s,1} = score;
        scores_names{s,1} = names;
        
    end 
    
    [score, sorted_indices] = sortrows(pathway_score);
    
    sorted_pathways = all_pathways(sorted_indices,:);
    sorted_metalist = metalist_indices(sorted_indices);
    sorted_side_reactions = side_reaction_indices(sorted_indices);    
    
    sorted_v = v_indices(sorted_indices);
    sorted_dG = deltaGs(sorted_indices);
    sorted_scores = scores(sorted_indices);
    sorted_scores_names = scores_names(sorted_indices);
    
    fid = fopen(input_data.result_target, 'A');
    rid = fopen([input_data.result_target '_reactions'], 'A');
    sid = fopen([input_data.result_target '_side'], 'A');
    
    offset = 0;

     for s = 1:numel(sorted_indices) 
        reactionlist = sorted_pathways(s, :);   
        position = s-offset;
        if reactionlist == zeros(1,max_length)
           offset = offset+1;
            continue;
        end
        metalist = metalist_all{sorted_metalist(s)};
        sidelist = side_reactions{sorted_side_reactions(s)};
        v = v_all(sorted_v(s));
        v = v{:}; 
        dG = sorted_dG(s);
        dG = dG{:};
        
        if isempty(sidelist)
            sreactions = [];
        else  
            sreactions = setdiff(sidelist{:},v{:});
        end
        
         numside = numel(sreactions);
         scores_for_output = [sorted_scores{s},sorted_scores_names{s}];
         writePathways(v, reactionlist, reactions, reaction_names, reaction_map, compound_map, thermodynamics_map, fid, rid, position, scores_for_output, numside, dG);
         plotThermodynamics(dG', reactions(reactionlist(find(reactionlist))), ['thermodynamics' num2str(position)]);
         writeOverallPathwayBalance(v, reactions, reaction_map,thermodynamics_map, [input_data.result_target '_balance' int2str(s)]);
         model = KEGGtoSBML(v{:}, reactions, reaction_map, compound_map, reversibilities_map);   
         writeCbToSBML(model, ['result' num2str(s)]);
         writeCytoscapeDataTable(model, reactionlist, metalist{:}, v{:}, reactions, compounds, ['datatable' num2str(position)])        
         writeSideReactions(sreactions, reactions, position, reaction_names, reaction_map, compound_map, sid);
         allreactions = unique([v{:}; sidelist{:}]);
         sidemodel = KEGGtoSBML(allreactions, reactions, reaction_map, compound_map, reversibilities_map);   
         writeCbToSBML(sidemodel, ['side' num2str(s)]);
         
         clear reactionlist metalist v dG;
    end
    fclose(fid);
    fclose(rid);
    fclose(sid);   
    cd(old_folder)
end