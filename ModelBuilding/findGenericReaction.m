% function findGenericReaction(reaction_map, outfile)
function generic_reactions = findGenericReaction(reaction_map)
   reaction_ids = keys(reaction_map);
%    fid = fopen(outfile, 'w+');    
    generic_reactions = {};
       
   for ind = 1:numel(reaction_ids)
        id = reaction_ids(ind);
        data = reaction_map(id{:});
        
        if isfield(data, 'comment')
            comment = data.comment;
            names = '';
            if isfield(data, 'names')
                names = data.names;
                names = names{:};
            end
            if ~isempty(strfind(comment, 'generic')) || ~isempty(strfind(comment, 'incomplete')) || ~isempty(strfind(comment, 'general')) || ~isempty(strfind(comment, 'general'))
%                 fprintf(fid, '%s\t%s', data.entry, names);
%                 fprintf(fid, '\t%s\n', comment);
                
                generic_reactions = [generic_reactions; {data.entry}];
            end
        end
   end  
%     fclose(fid);
    
end