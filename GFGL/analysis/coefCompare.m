function []=coefCompare(edges,actual,estimated)
% plots a comparison between estimated and actual matrix components, is 
% useful for monitoring performance of lasso estimators and assessing
% shrinkage effect.
% coefCompare(edges,actual,estimated)
%
% edges =[[1,2];[1,4];[1,1]] a vector of index pairs specifying location
% in matrix
% actual: pxp simulated covariance/precision matrix
% estimated: pxp estimated covaraince/precision matrix

% Get end bits for trimming due to windowing effects of estimation.
delta=length(actual)-length(estimated);
if(delta ~=0)
trm=ceil(delta/2)
else
    trm=1;
end

p=size(actual,1);   % total number of variables
nedge=size(edges,1);    % number of edges to plot

cc=hsv(nedge);

for i=1:nedge

    act=squeeze(actual(edges(i,1),edges(i,2),:));
    stdact=sqrt(var(act));
    act2=(act);%/stdact;
    plot(act2(trm:(end-trm)),'color',cc(i,:));
    hold on
    
    est=squeeze(estimated(edges(i,1),edges(i,2),:));
    stdest=sqrt(var(est));
    est2=(est);%/stdest;
    plot(est2,'--','color',cc(i,:));
    hold on
%     str=['Estimated std:',num2str(stdest,'% 10.2f'),'     Actual std:',num2str(stdact,'% 10.2f'),'    Actual range:',num2str(max(act)-min(act)),'     Estimated range',num2str(max(est)-min(est))];
%     text(10, max([act2;est2]), str, 'Color', cc(i,:));
    hold on
end

end