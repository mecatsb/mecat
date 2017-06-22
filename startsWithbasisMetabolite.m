function res = startsWithbasisMetabolite(reactionlist, start_reactions)
    reacs = intersect(reactionlist(1), start_reactions);

    if isempty(reacs)
        res = 1;
    else
        res = 0;
    end;  
  
end