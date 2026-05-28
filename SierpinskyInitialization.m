function pop = SierpinskyInitialization(n, dim, low, up)
    midpoint = (low + up)/2;
    pop = bsxfun(@plus, midpoint, bsxfun(@times, (rand(n,dim) - 0.5)*0.25, (up - low)));
    pop = max(min(pop, up), low);
end