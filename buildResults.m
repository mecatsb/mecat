function result = buildResults(path, parameter_file)
% path: path the result mat-files representing one pathway solution each
% parameter file: mat file containing variables created by the fsr call
%   C: arcs 
%   S: stoichiometric matrix
%   rpair_reaction_map
%   T: target index

    old_folder = cd(path);
    load(parameter_file);
    numArcs = size(C,1); 
    [~, numReactions] = size(S);
      
   	pathways = dir('solution*.mat');
    number_pathways = numel(dir('solution*.mat'));
      
    % data structures for results
    result.metabolite_lists = cell(number_pathways,1);
    result.reaction_lists = cell(number_pathways,1);
    result.lengths = zeros(number_pathways,1);
    result.x = cell(number_pathways,1);
    result.v = cell(number_pathways,1);
    result.active_reactions = cell(number_pathways,1);
    result.number_active_reactions = zeros(number_pathways,1);
   
    u = 1:numArcs;
    z = numArcs+1:numArcs+numReactions;
    v = numArcs + numReactions+1:numArcs+2*numReactions;
    
    for n = 1:number_pathways
        
       load(['solution'  int2str(n) '.mat']);
       us = find(x(u)>0.001);
       zs = find(x(z)>0.001);
       
       num_rpairs_in_path = numel(us);
       num_metabolites_in_path = numel(us)+1;
       
       metalist = [zeros(num_metabolites_in_path-1,1);T];
       reaclist = cell(num_rpairs_in_path,1);
     
       for ind = num_rpairs_in_path:-1:1
          % find starting nodes of last arc
          vact = find(C(us,2) == metalist(ind+1));

          if isempty(vact)
             disp('Path is not connected');
%               break;
          else

          % append starting node to metabolite list          
          metalist(ind) =  C(us(vact),1);
          
          % reactions with flux        
          reaction_list = intersect(rpair_reaction_map{us(vact)},zs);

          % append coresponding reaction to list
          reaclist{ind,1} = reaction_list;
          end
       end 
       
       vs = find(x(v));       
       result.metabolite_lists{n} = metalist;
       result.reaction_lists{n} = reaclist;
       result.lengths(n) = numel(reaclist);
       result.x{n} = x;
       result.v{n} = vs;
       result.z{n} = zs;
       result.active_reactions{n} = num2cell(vs,1);
       result.number_active_reactions(n) = numel(vs);

    end
    result.z = result.z';
    save('results', 'result');
    cd(old_folder);
end