function r = tetrac(M)
%tetrac computes pairwise tetrachoric correlation coefficients.
%   r = tetrac(M) returns a vector of pairwise correlation coefficients
%   for an M-by-N data array M. Rows and columns of M correspond to
%   observations and variables, respectively.

%   Source: https://www.statology.org/tetrachoric-correlation/
%           https://link.springer.com/chapter/10.1007/978-3-319-93177-7_2
            
m_col = size(M,2);  % Number of columns in events matrix

r = NaN(size(M,2)); % Declaring dimension of correlation matrix
p = pi();

for col = 1:m_col
    for row = 1:m_col
        
        if row == col % Diagonal elements (correlation with itself)
        
            r(row,col) = 1;
        
        elseif isnan(r(row,col)) == 1 % Checking if element is not calculated (proceed if NaN)            
            
            % z = crosstab(M(:,col),M(:,row)); % fails when all 1 or 0
            
            a = sum(M(:,col)==0 & M(:,row)==0);
            b = sum(M(:,col)==0 & M(:,row)==1);
            c = sum(M(:,col)==1 & M(:,row)==0);
            d = sum(M(:,col)==1 & M(:,row)==1);
            
            ad = a * d; % diagonal
            bc = b * c; % off-diagonal
            
            % 1 if the off-diagonals are 0 and −1 if the diagonals are 0
            
            if bc == 0 % Off-diagonal
                tr = 1;
            else 
                % when diagonal is 0, the formula returns -1
                tr = cos(p/(1 + sqrt(ad/bc)));
            end
            
            % symmetric matrix
            r(row,col) = tr;  
            r(col,row) = tr;    
        
        end
    end           
end    

end
