function M = CosineM(E)

Ne = size(E,2);
M = zeros(Ne,Ne);
for i = 1:Ne
    for j = 1:Ne
        A = double(E(:,i));
        B = double(E(:,j));
        xy   = dot(A,B);
        nx   = norm(A);
        ny   = norm(B);
        nxny = nx*ny;
        Cs   = xy/nxny; 
        M(i,j) = Cs;
    end
end
M(isnan(M)) = 0;
