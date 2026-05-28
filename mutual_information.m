
function mi = mutual_information(h)
    px = sum(h,2); py = sum(h,1); mi = 0;
    for i = 1:256
        for j = 1:256
            if h(i,j) > 0
                mi = mi + h(i,j)*log2(h(i,j)/(px(i)*py(j) + eps));
            end
        end
    end
end