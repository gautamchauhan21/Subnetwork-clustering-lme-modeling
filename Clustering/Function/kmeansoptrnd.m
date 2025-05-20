function sCl = kmeansoptrnd(E,N,NCl,type)

[Np,Ne] = size(E);

%% Randomization
if strcmp(type, 'var')
    ERnd = zeros(Np,Ne);
    for i = 1:Ne
        ERnd(:,i) = E(randperm(Np),i);
    end

    M = CovarM(ERnd);
else
    M = E;
end
%% k-means
parfor k = 1:N
    IDX = kmeans(M,NCl);
    s = silh(M,IDX);
    IDX0(k,:) = IDX;
    S(k) = median(s);
end
% keep best silhouette
[~,ClOK] = max(S);
IDX = IDX0(ClOK,:);
s = silh(M,IDX);
sCl = zeros(1,NCl);
for i = 1:NCl
    sCl(i) = median(s(IDX==i));
end
sCl = max(sCl);
end

%if sCl==1, keyboard; end
