function [arc_compounds, wo_arcs] = findCompoundsWithArcs(cea_file, compound_file)

    cea = importdata(cea_file);
    compounds = unique([cea(:,1); cea(:,2)]);
    modelcompounds = regexprep(importdata(compound_file), '\[c\]', '');
    compounds(isnan(compounds)) = [];
    arc_compounds = modelcompounds(compounds);  
    wo_arcs = setdiff(modelcompounds, arc_compounds);
end