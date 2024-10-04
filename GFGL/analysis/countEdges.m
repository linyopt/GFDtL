function [ nedge,edges ] = countEdges( G,thresh )
%counts number of edges in adjacency matrix G, also returns edge positions
% Thresh for filtering
p=size(G,2);
edges=[];
nedge=0;

for i=1:p-1
    for j=i+1:p
        if(abs(G(i,j))>thresh)
        edges=[edges;[i,j]];
        nedge=nedge+1; 
        end
    end
end

end

