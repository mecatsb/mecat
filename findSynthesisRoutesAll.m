function findSynthesisRoutesAll(d, rpair_reaction_map, stoichiometric_matrix, beta, cplex) 
     [numMetabolites, numReactions] = size(stoichiometric_matrix);
     numArcs = size(d,1); 

     % variable that determines which arc is active
     u = 1:numArcs;

     % reaction variables
     z = numArcs+1:numArcs+numReactions;
     
	% variable with value of flux for each reaction
     v = numArcs + numReactions+1:numArcs+2*numReactions;
     
     s=[];

    cplex.Model.sense = 'minimize';
    
    fprintf('Solve Model \n');
    n = 1;
    while true
        cplex.solve(); 
        % (MIP only) Solution is integer infeasible
        if cplex.Solution.status == 103 
            
           cplex.refineConflict
            disp(['Already calculated: ' num2str(n-1)])  
            conflict = cplex.Conflict;
            save('conflict', 'conflict');
            return      
        end

        x = cplex.Solution.x;
        save(['solution' int2str(n)], 'x');
        if isequal(x,xold)
            save(['cplex_same_solution_' int2str(n)], 'cplex');
            return
        end
        us = find(x(u)>0.001);
        zs = find(x(z)>0.001);

        num_rpairs_in_path = numel(us);
        num_metabolites_in_path = numel(us)+1;

        metalist = [zeros(num_metabolites_in_path-1,1);beta];
        reaclist = cell(num_rpairs_in_path,1);
        reaclist_all = cell(num_rpairs_in_path,1);

        for ind = num_rpairs_in_path:-1:1
          % find starting nodes of last arc
          vact = find(d(us,2) == metalist(ind+1));

          if isempty(vact)
             disp('Path is not connected');
             cplex.writeModel(['model' int2str(n) '.lp']);
             save(['cplex' int2str(n)], 'cplex');
              break;
          end

          % append starting node to metabolite list
          metalist(ind) =  d(us(vact),1);

          % reactions with flux       
          reaction_list = intersect(rpair_reaction_map{us(vact)},zs);

          % append coresponding reaction to list
          reaclist{ind,1} = reaction_list;
          reaclist_all{ind,1} = rpair_reaction_map{us(vact)};       
        end

        alpha = metalist(1);
        
        if(~any(x(s)>0.5))
           cplex.addCols(0,spalloc(length(cplex.Model.A(:,1)),1,0),0,1,'B');
           s=numArcs+2*numReactions+numMetabolites+1:numArcs+2*numReactions+numMetabolites+1+length(s);
           curSol=s(end);
           tmp = spalloc(3,length(cplex.Model.A(1,:)),3*num_rpairs_in_path+numArcs);
           tmp(1,u(x(u)>0.5)) = 1;         
           arcsto = d(:,2) == alpha;
           tmp(1,u(arcsto)) = -1;
           tmp(1,curSol) = -1;
           lb=zeros(3,1);
           up=ones(3,1);
           up(1) = length(find(x(u)>0.5))-1;
           lb(1) = -Inf;
           tmp(2,u(x(u)<0.5)) =1;
           tmp(2,curSol) = numArcs;
           up(2)= numArcs;
           lb(2)=-Inf;
           tmp(3,u(x(u)>0.5)) =1;
           tmp(3,curSol) = - num_rpairs_in_path;
           up(3)=Inf;
           lb(3)=0;
           cplex.addRows(lb,tmp,up);
           clear tmp;
        else
           curSol=s(x(s)>0.5);
        end
        
        % expand 
        sizeVec = cellfun('prodofsize',reaclist);
        indices = fliplr(arrayfun(@(n) {1:n}, sizeVec));
        [indices{:}] = ndgrid(indices{:});
        combMat = cellfun(@(c,i) {reshape(c(i(:)), [], 1)}, reaclist, fliplr(indices));
        
        combMat = [combMat{:}];
        
        if (size(combMat,1)) > 1
            disp('mehr als eine Zeile')
        end
        
        tmp = spalloc(size(combMat,1),length(cplex.Model.A(1,:)),numel(us));
        indices = sub2ind(size(tmp),repmat((1:size(combMat,1))',1,size(combMat,2)),z(combMat));
        tmp(indices)=1;
        tmp(:,curSol)=numReactions;
        
        %calculate number of unique reactions on path
        lens = cellfun(@(row) length(unique(row)),num2cell(combMat,2));
        %lens = lens.';
        %cplex.addRows(-Inf, tmp,lens-1+numReactions);
        cplex.addRows(ones(numel(lens),1)*-Inf, tmp,lens-1+numReactions);
        
        clear tmp;
        
        n = n+1;
    end

end

       
