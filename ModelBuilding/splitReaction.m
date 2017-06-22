function [split, comp] = splitReaction(reaction)
    % remove compartment
    token = '\<\[[pce]\] : *';
    comp = regexp(reaction, '\<\[[pce]\]', 'match');
    reaction = regexprep(reaction,token,'');

    del_token  = '(-->|<==>|<=>)';
    del = regexp(reaction,del_token,'match');
    
    for j = numel(del)
        d = del{j,1};
        split = regexp(reaction, d, 'split');
        split = strtrim(split);
        split = split{:};
                
        for k = 1:numel(split)
%                  split(k) = strtrim(regexp(split(k), '+(?=\))', 'split')); 
                    split(k) = strtrim(regexp(split(k), '+(?!\w*\))', 'split')); 
        end
    end
end