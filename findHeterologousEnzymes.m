function num_het = findHeterologousEnzymes(reactionlist, reaction_genes)

    max_length = size(reactionlist,2);
    reactionlist(reactionlist==0) = [];
    het = 0;
    num_het = 0;
    num_reactions = numel(reactionlist);
    
    for k = 1:num_reactions
        t = reactionlist(k);
        data = reaction_genes(t);
        
        if (isempty(find(not(data.het), 1)))            
            if (~isfield(data, 'organisms'))               
                het = het+ num_reactions;   
                num_het = num_het+1;
            else
                het = het+max_length;
                num_het = num_het+1;
            end    
        end
    end
end