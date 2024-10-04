function [  ] = plotGraph( Z,r,leg )
%PLOTDIFFGRAPH Plots difference between two precision matrices as a graph

% - Z1 = first precision matrx
% - Z2 = second precision matrix
% - r = radius of graph
[P,~]=size(Z);

[ xy ] = circleGraph( Z,r );
% wgPlot(Z,xy);
% Find degree of node
% Set vertex weight to be proportional to the degree of node
d=ones([1,P]);
for i=1:P
    for j=1:P
        if(Z(i,j)~=0 && i~=j)
          d(i)=d(i)+1;   
        end
   end
end
 wgPlotlw(Z,xy,'vertexWeight',30*d,'vertexScale',1000,'legend',leg);
% wgPlot(Z,xy);

end

