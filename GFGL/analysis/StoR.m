function [ R ] = StoR( S )
%STOR Converts covariance matrix to correlation matrix

P=length(S);

for i=1:P
    for j=1:P
    R(i,j)=S(i,j)./(sqrt(S(i,i))*sqrt(S(j,j)));    
    end
end


end

