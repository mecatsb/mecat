function combMat = expandPathway(reactionlist)
    sizeVec = cellfun('prodofsize', [reactionlist{:}]);
    indices = fliplr(arrayfun(@(n) {1:n}, sizeVec));
    [indices{:}] = ndgrid(indices{:});
    combMat = cellfun(@(c,i) {reshape(c(i(:)), [], 1)}, [reactionlist{:}], fliplr(indices));
    combMat = [combMat{:}];
end