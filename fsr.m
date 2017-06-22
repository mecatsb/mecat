function fsr(input_file)

    % read input file
    input_data = parseInput(input_file);
    
    % load KEGG data
    old_folder = cd(input_data.kegg_data_path);   
    cd(old_folder)
    old_folder = cd(input_data.data_path);
    
    % load model data
    generic_metabolites = unique(importdata(input_data.generic_metabolites));
    excluded_metabolites = unique(importdata(input_data.excluded_metabolites));
    start_metabolites = importdata(input_data.start_metabolites);
    cofactors = importdata(input_data.cofactors);
    compounds = importdata(input_data.compounds);
    reactions = importdata(input_data.reactions);
    basis_metabolites   = importdata(input_data.basis_metabolites);  
    reversible_reactions = importdata(input_data.reversible_reactions);
    S = load(input_data.stoichiometric_matrix);
    S = S.S;
    
    % adjust arcs
    d = importdata(input_data.arcs);
    d2 = modifyArcs(d, reversible_reactions, [compounds(generic_metabolites); compounds(cofactors); compounds(excluded_metabolites)], compounds);           
    
    % eliminate duplicate arcs and connect them to their reactions
    [C, ~, ic] = unique(d2(:, 1:2), 'rows');
    rpair_reaction_map = cell(numel(C(:,1)),1);
    reaction_rpair_map = cell(numel(reactions),1);
    
    for k = 1:numel(ic)
        rpair_reaction_map{ic(k)} = [rpair_reaction_map{ic(k)} d2(k,3)];
        reaction_rpair_map{d2(k,3)} = [reaction_rpair_map{d2(k,3)} ic(k)];
    end
    
    for k = 1:numel(reaction_rpair_map(:,1))
        reaction_rpair_map{k} = unique(reaction_rpair_map{k});
    end
    
    % metabolite pools
    excluded_indices = zeros(numel(excluded_metabolites),1);
    for ind = 1:numel(excluded_metabolites)
        excluded_indices(ind) = find(ismember(compounds,compounds(excluded_metabolites(ind))));
    end

    metabolite_pool = unique([cofactors; start_metabolites; excluded_metabolites]);
    compound_indices = 1:numel(compounds);    
    external = setdiff(compound_indices, metabolite_pool);
    target = input_data.target_metabolite;

    cd(old_folder)
    mkdir(input_data.result_path)
    copyfile(input_file,input_data.result_path)
    cd(input_data.result_path)
    num = clock;
    folder_name = [num2str(num(3)) '.' num2str(num(2)) '.' num2str(num(1)) '_' num2str(num(4)) '.' num2str(num(5))];
    mkdir(folder_name)
    copyfile(input_file,folder_name)
    delete(input_file);
    cd(folder_name)
    
    % check if target is in model
    T = find(ismember(compounds,target));
    if isempty(T)
        disp(['Target compound  ' target ' not found. Check your target input file']);
        return;
    end 

    % set up basic MILP
    cplex = buildBasicProblem(C, rpair_reaction_map,S, reversible_reactions, external, metabolite_pool, basis_metabolites, T);

    save('parameters.mat','T','rpair_reaction_map', 'C', 'cplex', 'S', 'reversible_reactions');
    
    % solve MILP
    findSynthesisRoutesAll(C, rpair_reaction_map, S, T, cplex);

    clear cplex;
    cd(old_folder)
end
