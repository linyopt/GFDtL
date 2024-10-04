function [ inv0 ] = geninv( p,m,scale )
% Generates random graphical precision matrix, used for
% spawning simulations
%
% [ inv0 ] = geninv0( p,m,scale )
% p = number of variables (dimension of network)
% m = number of edges (off diagonal)
%
% Alex Gibberd 1 May 2013

inv0=[];

if(m>round((p^2-p)/2)) % Check that number of edges isnt too much
    disp('Error: You specified too many edges (m>p(p-1)/2)');
    return
end
% Perform final check to see if is actually positive semi-definite
if(~isempty(inv0))
    [~,E]=chol(inv0);   % Check for cholesky decomposition
        if(E)
             [inv0]=genpsd(p,m,scale);
        else
            return
        end
else        %First run, generate first attempt 
         [inv0]=genpsd(p,m,scale);
end

    %% Create attempts at Positive Semi-Definite Matrix
    function [inv0]=genpsd(p,m,scale)
    %% initiate inv0
% Use diagonal with maginitude equivalent to edge range to ensure positive
% semi-definite...

inv0=0.5*eye([p,p]);


%% Random selection of indicies for new entries

Tran = randomGraph(p,m);
[si,sj]=find(triu(Tran));
list=[si,sj]';

%% Actually add edges
for i=1:m
    %new=0.5-rand(1);    % new weight to add to edge
    new=scale*2*(rand(1)-0.5);  % Scale now goes here 9/10/2014
    
    % This constraint basically limits the minimum value of the edge to
    % scale/2
    while(abs(new)<scale/2)
       new=scale*2*(rand(1)-0.5);   % Keep picking until accepted
    end
    inv0(list(1,i),list(2,i))=-new;
    inv0(list(2,i),list(1,i))=-new;    % make symmetric
    
    %% Try to make positive semi-definite
    inv0(list(1,i),list(1,i))=inv0(list(1,i),list(1,i))+abs(new);
    inv0(list(2,i),list(2,i))=inv0(list(2,i),list(2,i))+abs(new);
end

        sigma=inv(inv0);  % calculate inverse of new sigma
        
        % Create vector of inverse variances 
        v=diag(sigma);
        v=1./sqrt(v);
        V=diag(v);
        % Divide by variance associated with edge updates
        sigma=V*sigma*V';
        
        inv0=inv(sigma);
        
        % Prune for really small values in the inverse (numerical??)
        inv0(find(abs(inv0)<1e-12))=0;
        
                
% Rescale matrix
%inv0=scale*inv0;

function Tran = randomGraph(p,n)
Tran  = tril(ones(p,p),-1)';
idx = find(Tran>0);
len = numel(idx);
eidx = randsample(len,n);
Tran=zeros(p,p);
Tran(idx(eidx))=1;
Tran = Tran+Tran';
end

    end

end
