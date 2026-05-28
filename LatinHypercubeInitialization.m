% ---------------- Initialization Helper Functions ----------------

function pop = LatinHypercubeInitialization(n, dim, low, up)
    tempLHS = lhsdesign(n, dim);
    pop = bsxfun(@plus, low, bsxfun(@times, tempLHS, (up - low)));
end

