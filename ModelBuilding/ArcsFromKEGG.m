function ArcsFromKEGG(compound_file, reaction_file, revfile, reaction_map, rpair_map, outfile)

%    only consider main RPAIRS (can be changed)
%     pattern = '(main)';
    pattern = '(trans)|(main)';
    compounds = importdata(compound_file);
    reactions = importdata(reaction_file);
    
    fid = fopen(revfile, 'r');
    rev_data = textscan(fid, '%s%d', 'delimiter', '\t', 'whitespace', '');
    fclose(fid);

    if ~isempty(rev_data{1,1})
        revMap = containers.Map(rev_data{1,1}, rev_data{1,2});
    else
        revMap = containers.Map();
    end

    fid = fopen(outfile,'w+');
    
    for ri = 1:numel(reactions)
        abbrev = reactions(ri);     
        abbrev = abbrev{:};
        rev = 0;

         if(isKey(reaction_map,abbrev))            
            data = reaction_map(abbrev);
             if isfield(data, 'rpairs')

                rpair_data = data.rpairs;
                rpair_types = cell(1,numel(rpair_data));
                rpairs = cell(1,numel(rpair_data));
                
                for nrp = 1:numel(rpair_data)
                   rp = rpair_data{nrp};
                   rpair_types{nrp} = rp{1,3};
                   rpairs{nrp} =  rp{1,1};
                end
                
                educts = data.educts;               
                
                type = regexp(rpair_types, pattern);
            
            if(isKey(revMap, abbrev))
                rev = revMap(abbrev);
            end

            for ind = 1:numel(type)
                if ~isempty(type{ind})
            % search RPAIR in map
                rpair_id = rpairs(ind);
                    if isKey(rpair_map, rpair_id)
                        rpair = rpair_map(char(rpair_id));
                        cea_compounds = rpair.compounds;

                        index = strfind(compounds, cea_compounds{1});
                        c_index_1 = find(not(cellfun('isempty', index)));

                        index = strfind(compounds, cea_compounds{2});
                        c_index_2 = find(not(cellfun('isempty', index)));
                        
                        if (rev == 0 || rev == 1)   
                            if (find(strcmp(cea_compounds{1}, educts)))
                                fprintf(fid, '%d\t%d\t', c_index_1, c_index_2);
                            else
                                fprintf(fid, '%d\t%d\t', c_index_2, c_index_1);
                            end
                        else
                            if (find(strcmp(cea_compounds{1}, educts)))
                                fprintf(fid, '%d\t%d\t', c_index_2, c_index_1);
                            else
                                fprintf(fid, '%d\t%d\t', c_index_1, c_index_2);
                            end
                        end
                       
                        fprintf(fid,'%d\t',ri);
                        fprintf(fid,'%d',1);
                        fprintf(fid,'\n');    
%                         disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
                    end
                end
            end   
        else
             disp(['Reaction: ' abbrev ' does not contain any RPAIRS']);
            end
        end
    end
    
    fclose(fid);
end