function input_data = parseInput(input_file)

    input = fileread(input_file);

    fieldkeys = {'DATA_PATH', 'COFACTORS', 'START_METABOLITES', ...
        'EXTERNAL_METABOLITES', 'MEDIUM', 'COMPOUNDS', 'REACTIONS', ...
        'COMPOUND_NAMES', 'REACTION_NAMES', 'ARCS', 'STOICHIOMETRIC_MATRIX', ...
        'TARGET_METABOLITE', 'NUMBER_PATHWAYS', 'REVERSIBLE_REACTIONS',...
        'KEGG_DATA_PATH', 'REACTION_MAP', 'COMPOUND_MAP', 'ENZYME_MAP'...
        'ORGANISM_MAP', 'RESULT_PATH' 'RESULT_TARGET' 'MAX_PATHWAY_LENGTH' ...
        'MODEL', 'EXCLUDED_COMPOUNDS' 'NON_ENZYMATIC_REACTIONS' 'THERMODYNAMICS_MAP'...
        'SBML_MODEL' 'REVERSIBILITIES_MAP' 'BASIS_METABOLITES' 'GENERIC_METABOLITES' 'EXCLUDED_METABOLITES' 'R_IN_HOST'};
    
    pattern = strjoin(fieldkeys, '|');
    [data, keys] = regexp(input, pattern, 'split', 'match');
    data(cellfun(@isempty,data)) = [];
       
     if ~numel(data) == numel(keys)
        error('Parsing error');
    end

    for k = 1:numel(keys)
        switch keys{k}         
            case 'DATA_PATH'
                input_data.data_path = strtrim(regexprep(data{k}, 'DATA_PATH', ''));
            case 'RESULT_PATH'
                input_data.result_path = strtrim(regexprep(data{k}, 'RESULT_PATH', ''));
            case 'COFACTORS'
                input_data.cofactors = strtrim(regexprep(data{k}, 'COFACTORS', ''));
            case 'GENERIC_METABOLITES'
                input_data.generic_metabolites = strtrim(regexprep(data{k}, 'GENERIC_METABOLITES', ''));
            case 'EXCLUDED_METABOLITES'
                input_data.excluded_metabolites = strtrim(regexprep(data{k}, 'EXCLUDED_METABOLITES', ''));
            case 'START_METABOLITES' 
                input_data.start_metabolites = strtrim(regexprep(data{k}, 'START_METABOLITES', ''));
            case 'EXTERNAL_METABOLITES'
                input_data.external_metabolites = strtrim(regexprep(data{k}, 'EXTERNAL_METABOLITES', ''));
            case 'BASIS_METABOLITES'
                input_data.basis_metabolites = strtrim(regexprep(data{k}, 'TERMINAL_METABOLITES', ''));
            case 'MEDIUM'
                input_data.medium = strtrim(regexprep(data{k}, 'MEDIUM', ''));
            case 'COMPOUNDS'
                input_data.compounds = strtrim(regexprep(data{k}, 'COMPOUNDS', ''));
            case 'REACTIONS'
                input_data.reactions = strtrim(regexprep(data{k}, 'REACTIONS', ''));
            case 'COMPOUND_NAMES'
                input_data.compound_names = strtrim(regexprep(data{k}, 'COMPOUND_NAMES', ''));
            case 'REACTION_NAMES'
                input_data.reaction_names = strtrim(regexprep(data{k}, 'REACTION_NAMES', ''));
            case 'ARCS'
                   input_data.arcs = strtrim(regexprep(data{k}, 'ARCS', ''));
            case 'STOICHIOMETRIC_MATRIX'
                input_data.stoichiometric_matrix = strtrim(regexprep(data{k}, 'STOICHIOMETRIC_MATRIX', ''));
            case 'TARGET_METABOLITE'
                input_data.target_metabolite = strtrim(regexprep(data{k}, 'TARGET_METABOLITE', ''));
            case 'NUMBER_PATHWAYS'
                input_data.number_pathways = str2double(strtrim(regexprep(data{k}, 'NUMBER_PATHWAYS', '')));
            case 'REVERSIBLE_REACTIONS'
                input_data.reversible_reactions = strtrim(regexprep(data{k}, 'REVERSIBLE_REACTIONS', ''));
            case 'KEGG_DATA_PATH'
                input_data.kegg_data_path = strtrim(regexprep(data{k}, 'KEGG_DATA_PATH', ''));
            case 'REACTION_MAP'
                input_data.reaction_map = strtrim(regexprep(data{k}, 'REACTION_MAP', ''));
            case 'COMPOUND_MAP'
                input_data.compound_map = strtrim(regexprep(data{k}, 'COMPOUND_MAP', ''));
            case 'ENZYME_MAP'
                input_data.enzyme_map = strtrim(regexprep(data{k}, 'ENZYME_MAP', ''));
            case 'ORGANISM_MAP'
                input_data.organism_map = strtrim(regexprep(data{k}, 'ORGANISM_MAP', ''));
            case 'RESULT_TARGET'
                input_data.result_target = strtrim(regexprep(data{k}, 'RESULT_TARGET', ''));
            case 'MAX_PATHWAY_LENGTH'
                input_data.max_pathway_length = str2double(strtrim(regexprep(data{k}, 'MAX_PATHWAY_LENGTH', '')));
            case 'MODEL' 
                input_data.model = strtrim(regexprep(data{k}, 'MODEL', ''));
            case 'NON_ENZYMATIC_REACTIONS'
                input_data.non_enzymatic_reactions = strtrim(regexprep(data{k}, 'NON_ENZYMATIC_REACTIONS', ''));
            case 'THERMODYNAMICS_MAP'
                 input_data.thermodynamics_map = strtrim(regexprep(data{k}, 'THERMODYNAMICS_MAP', ''));
            case 'SBML_MODEL' 
                input_data.sbml_model = strtrim(regexprep(data{k}, 'SBML_MODEL', ''));
            case 'REVERSIBILITIES_MAP'
                 input_data.reversibilities_map = strtrim(regexprep(data{k}, 'REVERSIBILITIES_MAP', ''));
            case 'R_IN_HOST'
                input_data.reactions_in_host = strtrim(regexprep(data{k}, 'R_IN_HOST', ''));     
                 
                 
            otherwise              
             continue;
        end   
    end
end
