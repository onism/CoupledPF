function A = createA(alpha,dim)
    A = zeros(dim,dim);
    for i = 0 : dim-1
        for j = 0 : dim-1
            A(i+1,j+1) = alpha^(1 + abs(i-j));
        end
    end
end