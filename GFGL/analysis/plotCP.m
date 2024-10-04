function [] = plotCP(cp1,cp2,data,offset)
%% Plots change points against supplied data vector
% cp is vector of changepoints
% data is a nxp data stream for plotting
% due to windowing n< n_original so use t_new for plotting and time
% adjustment
%
% plotCP(cp,data,t_new)

% Alex Gibberd 2015
%figure(1)

plot(data);
hold on
l=4.*max(sqrt(var(data)));    % limit = +-2sigma
for i=1:size(cp1,1)
    x=cp1(i)+offset;  % Adjust for windowing
    line([x x],[-l -l/2],'Color','Blue','LineWidth',2);
    hold on
end

for i=1:size(cp2,1)
    x=cp2(i)+offset;  % Adjust for windowing
    line([x x],[l/2 l],'Color','Red','LineWidth',2);
    hold on
end
hold off

end

