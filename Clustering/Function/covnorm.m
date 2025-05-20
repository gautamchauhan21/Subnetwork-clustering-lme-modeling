function C=covnorm(A,B,win)

N=length(A);
C=xcov(A,B,win)/std(A)/std(B)/N;