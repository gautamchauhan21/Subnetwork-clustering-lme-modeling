function cmapjet = fcn_cmapjet(len)
if nargin == 0
    len = 256;
end
cmapjet = jet(7); cmapjet(4,:) = 1; cmapjet = interp1(linspace(0,1,size(cmapjet,1)),cmapjet,linspace(0,1,len));