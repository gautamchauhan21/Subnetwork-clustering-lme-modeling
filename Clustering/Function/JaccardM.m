function M = JaccardM(E)

Ne = size(E,2);
M = zeros(Ne,Ne);
parfor i = 1:Ne
    for j = 1:Ne
        A = E(:,i);
        B = E(:,j);
        M(i,j) = sum(A.*B)/(sum(A+B) - sum(A.*B));
    end
end
M(isnan(M)) = 0;
