function [ cpPos ] = cpCluster( Theta,thresh )
%CPCLUSTER outputs histogram of changepoint positions

[P,~,T]=size(Theta);
cpPos=zeros([T,1]);

for t=2:T
   for i=1:P-1
       for j=i+1:P
       if(abs(Theta(i,j,t)-Theta(i,j,t-1))>thresh)
           cpPos(t)=cpPos(t)+1;
       end
       end
   end
end
end

