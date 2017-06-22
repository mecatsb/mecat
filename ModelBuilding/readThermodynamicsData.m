function thermodynamics_map = readThermodynamicsData(tdata_file, threshold, outfile)

    fid = fopen(tdata_file, 'r');
    % KEGG dG0' sigma(dG0) ph mM T Note
    tdata = textscan(fid, '%s%f%f%f', 'delimiter', '\t', 'whitespace', '', 'headerlines',1);
    fclose(fid);
    
    dG0 = tdata(1,2);
    dG0 = dG0{1,1};
    dG0u = tdata(1,3);
    dG0u = dG0u{1,1};   
    dGm = tdata(1,4);
    dGm = dGm{1,1}; 
    abbrev = tdata(1,1);
    abbrev = abbrev{1,1};

    thermodynamics_map = containers.Map();    
    
    fid = fopen(outfile, 'w+');
    count = 0;

     for k = 1:length(dGm)
        val = dGm(k);
        dgu = dG0u(k);
        a = abbrev{k};
        
%        if (dgu > 5*abs(val))
%        if(dgu > 392000)
%        if (abs(val) > abs(threshold) && dgu > 5*abs(val))
       if (abs(val) > 30 && dgu > 5*abs(val))
            direction = 1;
            count = count+1;
        elseif (val > -abs(threshold) && val < abs(threshold)) 
            direction = 0;
        elseif val < threshold
            direction = 1;
        elseif val > abs(threshold)
            direction = -1;
        else 
            direction = 1; %Reaktionen ohne dG sind jetzt irreversibel statt reversibel (war: NaN)
        end 

        fprintf(fid, '%s\t%d\n', a, direction);
        
        if direction == -1
            data.dG0 = - dG0(k);
            data.dGm = -dGm(k);
        else
            data.dG0 = dG0(k);
            data.dGm = dGm(k);
        end
        
        data.dG0u = dG0u(k);
        data.rev = direction;      
        
        thermodynamics_map(abbrev{k}) = data;
     end
     
     fclose(fid);
     
     save('thermodynamics_map', 'thermodynamics_map')
     count
end