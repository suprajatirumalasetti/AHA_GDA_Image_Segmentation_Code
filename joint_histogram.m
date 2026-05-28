function h = joint_histogram(A, B)
    h = zeros(256,256);
    for i = 1:numel(A)
        h(A(i)+1, B(i)+1) = h(A(i)+1, B(i)+1) + 1;
    end
    h = h / sum(h(:));
end

