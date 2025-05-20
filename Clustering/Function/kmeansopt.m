function [IDXs,sCl,M,S] = kmeansopt(E,N,type)

%E p parameters (cells) by N Events
%N number of trials per cluster number
%p p-value for cluster validation

Ne = size(E,2);

%% Covariance matrix
if strcmp(type,'var')
    M = CovarM(E);
else
    M = E;
end

%% k-means loop
p = 1;
P = N*18;
parfor k = 1:N*18
    NCl = floor((k-1)/N) + 2;
    IDX = kmeans(M,NCl); %Kmeans on distance of covariance metric
    s = silh(M,IDX);
    IDX0(k,:) = IDX;
    S(k) = mean(s);
    
end
%Best clustering for each cluster number
IDX1 = zeros(18,Ne);
for i = 1:18
    tmp = IDX0((i-1)*N+(1:N),:);
    [~,idx] = max(S((i-1)*N+(1:N)));
    IDX1(i,:) = tmp(idx,:);
end
%% keep best silhouette
[~,ClOK] = max(S);
NCl = floor((ClOK-1)/N) + 2;
IDX = IDX0(ClOK,:);
s = silh(M,IDX);
sCl = zeros(1,NCl);
for i = 1:NCl
    sCl(i) = median(s(IDX==i));
end

[sCl,xCl] = sort(sCl,'descend');
IDXs = zeros(1,Ne);

for i = 1:NCl
    IDXs = IDXs + (IDX == xCl(i))*i;
end
end
