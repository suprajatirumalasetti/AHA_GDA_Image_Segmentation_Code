function pop = SobolInitialization(n, dim, low, up)
    p = sobolset(dim);
    tempSobol = net(p, n);
    pop = bsxfun(@plus, low, bsxfun(@times, tempSobol, (up - low)));
end

