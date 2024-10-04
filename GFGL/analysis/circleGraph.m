function [ xy ] = circleGraph( A,r )
%CIRCLEGRAPH Creates xy position matrix for given adjacency matrix A

[P,~]=size(A);
xy=zeros([P,2]);

theta=2*pi/P;
for i=1:P
xy(i,1)=r*cos(theta*i);
xy(i,2)=r*sin(theta*i);
end

end

