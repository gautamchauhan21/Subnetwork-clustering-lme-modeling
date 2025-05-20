function s = silh(M,IDX)

NCl = max(IDX);
NRace = size(M,1);
a = zeros(NCl,NRace);
for j = 1:NCl
    for i = 1:NRace
        Events = IDX == j;
        Events(i) = 0;
        a(j,i) = 1 - mean(M(i,Events));
    end
end
b = zeros(1,NRace);
s = zeros(1,NRace);
for i = 1:NRace
    b(i) = min(a(~ismember(1:NCl,IDX(i)),i));
    s(i) = (b(i)-a(IDX(i),i))/max(b(i),a(IDX(i),i));
end