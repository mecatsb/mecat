function [r, substrates , products, stoich_s, stoich_p] = getReactionPartners(r, reaction_map)
    f = strfind(r, '-');
    if ~isempty(f{:})
        r = regexprep(r, '-', '');
        reaction_data = reaction_map(r{:});
        products = reaction_data.educts;
        substrates = reaction_data.products;
        stoich_p = reaction_data.stoich_e;
        stoich_s = reaction_data.stoich_p;
    else
        reaction_data = reaction_map(r{:});
        substrates = reaction_data.educts;
        products = reaction_data.products;
        stoich_s = reaction_data.stoich_e;
        stoich_p = reaction_data.stoich_p ;           
    end
end