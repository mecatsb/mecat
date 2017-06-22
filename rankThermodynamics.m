function [deltaG, sumdG, num_wo_data, dGs] = rankThermodynamics(reactionlist, thermodynamics_map, reactions)

    reactionlist(reactionlist==0) = [];
    deltaG = 0;
    sumdG = 0;
    num_wo_data = 0;
    isrev = 0;
    dGs = zeros(numel(reactionlist),1);
     
    for k = 1:numel(reactionlist)
        id = reactions{reactionlist(k)};
        
        if strfind(id, '-')   
            id = regexprep(id, '-', '');
            isrev = 1;
        end
        
        if isKey(thermodynamics_map, id)        
            
            tdata = thermodynamics_map(id);
            rev = tdata.rev; 
            dGu = tdata.dG0u;
            dG = tdata.dGm;
            sdG = 0;
            
            uncertain = false;
            
            if (dGu > 1000)
                uncertain = true;
                dG = 0;
            end
            
            if (uncertain == false)               
                sdG = tdata.dGm + abs(tdata.dGm);
            end
            
            if (isnan(dG) )
                num_wo_data = num_wo_data+1;
                continue;
            end
            
            if (isrev && (rev == 0))
                 dG = -dG;
                 sdG = -tdata.dGm + abs(tdata.dGm);
            end
            
            if (uncertain)
                num_wo_data = num_wo_data+1;
            else
                deltaG = deltaG+dG;
                sumdG = sumdG + sdG;
            end
           
            dGs(k) = dG;
        else
            num_wo_data = num_wo_data+1;
        end       
    end
end
