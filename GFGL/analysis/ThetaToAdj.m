function [ A] = ThetaToAdj( Theta,thresh )
%THETATOADJ Summary of this function goes here
%   Converts Theta estimate to graph adjacency matrix


[P,~,T]=size(Theta);
A=zeros([P,P,T]);
for t=1:T
    for i=1:P-1
        for j=i+1:P
            if abs(Theta(i,j,t))>thresh
                A(i,j,t)=1;
                A(j,i,t)=1;
            else
                A(i,j,t)=0;
                A(j,i,t)=0;
            end
        end
    end
end


end

