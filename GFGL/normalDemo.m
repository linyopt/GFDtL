function []=normalDemo(lambda1G,lambda2G,lambda1I,lambda2I)
%% hdDemo
% This is a wrapper function for publishing the example normalExample(). It
% generates a pdf file which displays the results of the analysis and
% discusses the various outputs. Alternatively, one may simply run the
% function normalExample() directly.

% Alex Gibberd 7/4/2015

if(nargin<1)
    lambda1G=0.35;
end
if(nargin<2)
    lambda2G=20;
end
if(nargin<3)
    lambda1I=0.25;
end
if(nargin<4)
    lambda2I=5;
end

assignin('base','lambda1G',lambda1G);
assignin('base','lambda2G',lambda2G);
assignin('base','lambda1I',lambda1I);
assignin('base','lambda2I',lambda2I);

funpar=['lambda1G',char(10),...
    'lambda2G',char(10), 'lambda1I',char(10), 'lambda2I',char(10),...
    'normalExample(lambda1G,lambda2G,lambda1I,lambda2I)',char(10)];

options=struct('codeToEvaluate',funpar,'format','html');
ref=publish('normalExample',options);
web(ref);   % View output..
display('Note: If the figures dont update, try clicking the refresh button in the MATLAB browser.');
end