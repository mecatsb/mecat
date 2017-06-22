function reaction_gene_map = getEnyzmeInformation (reactionlist, reactions, reaction_map, non_enzymatic_reactions, host, enzyme_map)

    reactionids = unique(reactionlist);
    reactionids(reactionids==0) = [];
    
    reaction_gene_map = containers.Map('KeyType','double','ValueType','any');
    organismc = 0;
    
    for t = 1:numel(reactionids)
          
         rid = reactions(reactionids(t));
         rid = regexp(rid, 'R\d{5}', 'match');
         rid = rid{:};
         reac = reaction_map(rid{:});
         
         %enzymes for reaction
         if isfield(reac, 'enzyme')
            data.ec = reac.enzyme;
            data.het = [];
            data.organisms = cell(numel(data.ec),1);

            if (isempty(data.ec))
                data.het(1) = 1;
            end
            
            for ne = 1: numel(data.ec)
                 enzyme_data = enzyme_map(data.ec{ne});   

                 % organisms associated to this enzyme
                 if (isfield(enzyme_data, 'organism'))
                      organismc = organismc+1;
                      data.organisms{ne} = enzyme_data.organism;

                    % look if host organism is in organism list
                    indices = find(cellfun(@(x) strcmpi(x,host), enzyme_data.organism), 1);
                    if (~isempty(indices))
                        data.het(ne) = 0;
                    else
                        data.het(ne) = 1;
                    end
                 else                  
                     % no information: count as heterolog
                     data.het(ne) = 1;
                 end
            end   
        
         % no enzyme given
         else
             if isempty(find(non_enzymatic_reactions == reactionids(t), 1))
                data.het(1) = 1;
                data.ec = {''};
             else
               data.het(1) = 0;
               data.ec = {''};
             end             
         end
         reaction_gene_map(reactionids(t)) = data;   
    end    
end