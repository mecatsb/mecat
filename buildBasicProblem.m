function cplex = buildBasicProblem(d, rpair_reaction_map, stoichiometric_matrix, reversible_reactions, E, Em, terminal, beta)   
    fid = fopen('problem_stats','w+');

    [numMetabolites, numReactions] = size(stoichiometric_matrix);
    numArcs = size(d,1);
    numReversible = size(reversible_reactions,1);
    
    fprintf(fid, '%s\t%d\n', 'number arcs:', numArcs);

    metabolite_indices = 1:numMetabolites;
    
    fprintf(fid, '%s\t%d\n', 'number metabolites:', numMetabolites);
    
    M = 1000;
    
    educt_to_reaction = stoichiometric_matrix(1:numMetabolites,:) < 0;
    
    %split stoichiometric matrix
    stoichiometric_matrix(Em,:) = 0;
    Sint = stoichiometric_matrix;
    Sext = stoichiometric_matrix(E,:);
    Sbeta = stoichiometric_matrix(beta,:);
    Sbeta(Sbeta>0) = 0; 
    
    Sint([E Em'],:) = [];
    numExt = length(E);
    numInt = numMetabolites - numExt - length(Em);
    
    fprintf(fid, '%s\t%d\n', 'number external:', numExt);
    fprintf(fid, '%s\t%d\n', 'number internal:', numInt);

    int = setdiff(metabolite_indices,[E Em']);

    numVars = numArcs+2*numReactions+numMetabolites;
    numConstraints = 4*numMetabolites+numArcs+numReversible+2*numReactions+numInt+numExt; 
    
    fprintf(fid, '%s\t%d\n', 'number vars:', numVars);
    fprintf(fid, '%s\t%d\n', 'number constraints:', numConstraints);
    
    % variable that determines which arc is active in CFP
    u = 1:numArcs;
    
    % reaction variables
    z = numArcs+1:numArcs+numReactions;
    
    fprintf(fid, '%s\t%d-%d\n', 'reaction variables :', numArcs+1, numArcs+numReactions);
    
    % variable with value of flux for each reaction
    v = numArcs + numReactions+1:numArcs+2*numReactions;
    
%     metabolite variables, 1: metabolite is part of the path, 0 otherwise
     m = numArcs+2*numReactions+1:numArcs+2*numReactions+numMetabolites;
   
    fprintf(fid, '%s\t%d-%d\n', 'reaction fluxes :', numArcs + numReactions+1, numArcs+2*numReactions);
    
    % set type for each variable
    binary = 'B';
    continuous = 'C';
    ctype(u)= binary(1,ones(1,numArcs));
    ctype(z) = binary(1,ones(1,numReactions));
    ctype(v) = continuous(1,ones(1,numReactions));
    ctype(m) = binary(1,ones(1,numMetabolites));
     
    % model
    cplex = Cplex('');
    cplex.Model.sense = 'minimize';
    cplex.DisplayFunc = [];
        
    % objective function: minimize sum(u)
    obj = zeros(numVars,1);
    obj(u,1) = 10000;
    obj(z,1) = 10000.0/(numReactions+1);
    % lower and upper bounds on variables
    lb = zeros(numVars,1);
    ub = inf*ones(numVars,1);
    
    cplex.addCols(obj, [], lb, ub, ctype);

    % constraint matrix 
    A = spalloc(numConstraints, numVars, 5*numArcs+4*numReversible+4*numReactions);
    
    % left and right hand side of constraints
    lhs = zeros(numConstraints,1);
    rhs = zeros(numConstraints,1);

    fprintf('Constraint 3 + 4 + 5 \n');
    for c = metabolite_indices
        arcsfrom = d(:,1) == c;
        arcsto = d(:,2) == c;
        
        % total number of input arcs <= 1 
        %constraint (5)
        A(c,u(arcsfrom)) = 1;
        
        A(c + numMetabolites, u(arcsfrom)) = 1;
        A(c + numMetabolites, u(arcsto)) = -1;
    end
    
    rhs(metabolite_indices) = 1;
    rhs(numMetabolites+1:numMetabolites+numMetabolites) = 1;

    for im = int
        rhs(im+numMetabolites) = 0;
    end
 
    for s = terminal'
        ends = d(:,2) == s;
        A(s,:) = 0;
        A(s, ends) = 1;
        rhs(s) = 0;
    end
      
    % constraint (1)
    endb = d(:,2)== beta;
  
    A(beta,:)=0; 
    A(beta,endb)=1; 

    lhs(beta,1) = 1; 
    rhs(beta,1) = 1;    
    
    % constraint (2)
    startb = d(:,1)== beta;     
    
    A(numMetabolites+beta,:)=0; 
    A(numMetabolites+beta,startb)=1;

    lhs(beta+numMetabolites,1)=0; 
    rhs(beta+numMetabolites,1)=0;
   
    % stoichiometric constraints
    A(2*numMetabolites+numArcs+numReversible+2*numReactions+1:2*numMetabolites+numArcs+numReversible+2*numReactions+numInt,v) = Sint;
    
    A(2*numMetabolites+numArcs+numReversible+2*numReactions+numInt+1:2*numMetabolites+numArcs+numReversible+2*numReactions+numInt+numExt,v) = Sext;
    rhs(2*numMetabolites+numArcs+numReversible+2*numReactions+numInt+1:2*numMetabolites+numArcs+numReversible+2*numReactions+numInt+numExt,1) = inf;
    
    %   beta can only be produced
    A(2*numMetabolites+numArcs+numReversible+2*numReactions+numInt+numExt+1,v) = Sbeta;
    lhs(2*numMetabolites+numArcs+numReversible+2*numReactions+numInt+numExt+1,1) = 0;
    rhs(2*numMetabolites+numArcs+numReversible+2*numReactions+numInt+numExt+1,1) = inf;
 
    Atmp = spalloc(numReactions, numVars, 2*numReactions);
    Atmp(:,z) = eye(numReactions)*M;
    Atmp(:,v) = -eye(numReactions);
    A(2*numMetabolites+numArcs+numReversible+1:2*numMetabolites+numArcs+numReversible+numReactions,:) = Atmp;
    rhs(2*numMetabolites+numArcs+numReversible+1:2*numMetabolites+numArcs+numReversible+numReactions,1) = inf;      
    clear Atmp
    
    Atmp = spalloc(numReactions, numVars, 2*numReactions);
    Atmp(:,z) = -eye(numReactions);
    Atmp(:,v) = eye(numReactions);
    A(2*numMetabolites+numArcs+numReversible+numReactions+1:2*numMetabolites+numArcs+numReversible+2*numReactions,:) = Atmp;
    rhs(2*numMetabolites+numArcs+numReversible+numReactions+1:2*numMetabolites+numArcs+numReversible+2*numReactions,1) = inf;
   
    clear Atmp
   
    % if network contains reversible reactions
    if numReversible > 0
        Atmp = spalloc(numReversible, numReactions, 2*numReversible);
        indices = sub2ind(size(Atmp),1:numReversible, reversible_reactions(:,2)');
        Atmp(indices) = 1;      
        indices = sub2ind(size(Atmp),1:numReversible, reversible_reactions(:,1)');
        Atmp(indices) = 1;
        A(2*numMetabolites+numArcs+1:2*numMetabolites+numArcs+numReversible,z) = Atmp;
        rhs(2*numMetabolites+numArcs+1:2*numMetabolites+numArcs+numReversible,1) = 1;

        clear Atmp
    end
    
    Atmp = spalloc(numArcs, numVars, 2*numArcs);
    Atmp(:,u) = -eye(numArcs);

    for k = 1:numel(rpair_reaction_map)  
       reaction_list = rpair_reaction_map{k};
       for rind = 1:numel(reaction_list)
            rl = reaction_list(rind);   
            idx = sub2ind(size(Atmp), k, z(rl)); 
            Atmp(idx) = 1;
       end
    end    
    
    A(2*numMetabolites+1:2*numMetabolites+numArcs,:) = Atmp;    
    rhs(2*numMetabolites+1:2*numMetabolites+numArcs,1) = inf;    
    
    % link u_ij's to m_i's
    for c = metabolite_indices
        arcsfrom = d(:,1) == c;
   
        A(2*numMetabolites+numArcs+numReversible+2*numReactions+numInt+numExt+c,u(arcsfrom)) = zeros(1,numel(find(arcsfrom)));
        A(2*numMetabolites+numArcs+numReversible+2*numReactions+numInt+numExt+c,m(c)) = 0; 
       
        rwithc = educt_to_reaction(c,:);
              
        A(3*numMetabolites+numArcs+numReversible+2*numReactions+numInt+numExt+c,z(rwithc)) = 0; 
        A(3*numMetabolites+numArcs+numReversible+2*numReactions+numInt+numExt+c, m(c)) = 0; 
        
        rhs(3*numMetabolites+numArcs+numReversible+2*numReactions+numInt+numExt+c) = 0;        
    end
    
    %beta is always on the path
    A(2*numMetabolites+numArcs+numReversible+2*numReactions+numInt+numExt+beta,:) = 0;
    A(2*numMetabolites+numArcs+numReversible+2*numReactions+numInt+numExt+beta,m(beta)) = 1;
    lhs(2*numMetabolites+numArcs+numReversible+2*numReactions+numInt+numExt+beta)=1;
    rhs(2*numMetabolites+numArcs+numReversible+2*numReactions+numInt+numExt+beta)=1;

    clear Atmp;
   
    cplex.Model.A=A;
    cplex.Model.lhs=lhs;
    cplex.Model.rhs=rhs;
    
    clear A lhs rhs
    
    fclose(fid);
end
