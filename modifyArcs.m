function d = modifyArcs(d, reversible, cofac, metabolites)
    
    %remove arcs with NaNs
    [row, ~] = find(isnan(d));
    d(row,:) = [];   

    if ~isempty(reversible > 0)
        % reaction in normal direction
        rev1 = reversible(:,2);
        % reaction in reverse direction
        rev2 = reversible(:,1);              
        
        final = length(d(:,1));
        for i=1:length(reversible(:,1))

            % find reversible reactions in arcs
            reactions = find(d(:,3) == rev1(i));

            % add arcs for reverse reaction
            new = (final + 1):(final + length(reactions));
            final = final + length(reactions);
            d(new,1) = d(reactions,2);
            d(new,2) = d(reactions,1);
            d(new,4) = d(reactions,4);
            d(new,3) = rev2(i);
        end
    end

    % validate cofactors given in list
    for i=1:length(cofac)
        cofactorsTmp = find(ismember(metabolites,cofac(i)));
        if isempty(cofactorsTmp)
            c = cofac(i);
            error('Cofactor ''%s'' not found.''\n', c{:});
        end
    
        for k = 1:length(cofactorsTmp)     
            cofac_cea1= d(:,1)==cofactorsTmp(k);
            d(cofac_cea1,:)=[];
            cofac_cea2= d(:,2)==cofactorsTmp(k);
            d(cofac_cea2,:)=[];        
        end
    end
end