function results = findSideReactions(input_file, results)

    input_data = parseInput(input_file);
    old_folder = cd(input_data.data_path);
    reactions = importdata(input_data.reactions);
    reactions_in_host = importdata(input_data.reactions_in_host);
    
    cd(old_folder);
    old_folder = cd(input_data.kegg_data_path);
    
    reaction_map = load(input_data.reaction_map);
    reaction_map = reaction_map.reaction_map;   
    
    metabolite_map = load(input_data.compound_map);
    metabolite_map = metabolite_map.compound_map;   
    
    cd(old_folder)
    mkdir(input_data.result_path);
    cd(input_data.result_path);   
        
    for s = 1:size(results.reaction_lists,1)

        reaclist = results.reaction_lists{s};
        if isempty(reaclist)
            results.side_reactions{s} = [];
            results.metabolite_reaction{s} = [];
            continue;
        end
         
        metab = {};
 
        num_reactions = numel(reaclist);
        reacnames = cell(num_reactions,1);

        for id = 1:num_reactions
            r = reactions(reaclist{id});
            reacnames(id) = r;

          if ~isempty(strfind(r, '-'))
                r = regexprep(r, '-', '');
                reaction_data = reaction_map(r{:});
                products = reaction_data.educts;
                educts = reaction_data.products;
            else
                reaction_data = reaction_map(r);
                educts = reaction_data.educts;
                products = reaction_data.products;         
          end

            educt_ids = values(metabolite_map, educts);
            product_ids = values(metabolite_map, products);
            metab = [metab educt_ids product_ids];
         end

        metab_entries1 = cell(1,numel(metab));

        for m = 1:numel(metab)
            met = metab{m};
            metab_entries1{m} = met.entry{:};
        end
  
        mets = unique(metab_entries1);      
        [reaction_ids, metabolite_reaction] = findReaction(mets, reaction_map, reactions_in_host);        
        reaction_ids2 = zeros(numel(reaction_ids),1);        
    
        for i = 1 : length(reaction_ids)
            tmp = find(ismember(reactions, reaction_ids{i}));
            if isempty(tmp)
                   continue;
            end    
            reaction_ids2(i) = tmp;
        end
    
        reaction_ids2 = reaction_ids2(find(reaction_ids2));
  
       	results.side_reactions{s} = reaction_ids2;
        results.metabolite_reaction{s} = metabolite_reaction;
        clear metab reaction_ids2;
    end
    save('results_side_reactions.mat','results')
end