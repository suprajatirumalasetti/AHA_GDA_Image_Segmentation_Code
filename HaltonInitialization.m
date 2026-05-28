function pop = HaltonInitialization(n, dim, low, up)
    p = haltonset(dim);
    tempHalton = net(p, n);
    pop = bsxfun(@plus, low, bsxfun(@times, tempHalton, (up - low)));
end

