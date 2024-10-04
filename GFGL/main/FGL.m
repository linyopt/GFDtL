function [ Theta,Z,S,fobj ] = FGL( y,lambda1,lambda2,gamma )
%ADMM to solve Fused Graphical lasso
% This penalises edge variation independently via the FLSA
% Alex Gibberd - UCL

% This is public version...for acompanying JCGS paper
% 5/4/2015

%% Some initial parameter setting
[P,T]=size(y');    % Get size of data array 


%% Calculate covariance dynamic structure

interval=1;  % how much to move the window along each step
maxsub=100;
maxit=1000;

tolDual=1e-4;   % Tolerance for convergence of ADMM
tolPrime=1e-4;

% for testing..
i=0;    % indexes for time-step reference S(:,:,i)..etc
S=[];   % initiate tracking for pre-regularised covariance matrix

window_width=30;

% Select window for kernel estimate

weight=hanning(window_width);

weight=sqrt(weight);    % Take square root to avoid double weighting when calculating (wY-wmu)(wY-wmu)'

% Generate new time axis for plotting (takes into account window width subtraction)
% Ts=ceil(window_width/2);
% Te=T-ceil(window_width/2);

Ts=1;Te=T;

data=y;

%% calculate S
% index ts so not to get confused with scope of t in subproblem
for ts=Ts:interval:Te
%     i=i+1;
% % Set bounds for window
% start=ts-floor(window_width/2)+1;
% fin=ts+ceil(window_width/2);
% 
% Y=[];
% mut=[];
% 
%     for j=1:P
%         Y=[Y, weight(:).*data(start:fin,j)];    % the set of samples to be used for covariance estimation in an NxP
%                             %  matrix where N is the number of samples and P is the number of variables
%     end
% 
%     % Generate empirical estimate of covariance matrix
%     % Use either sample mean or regularization 
%     St=zeros(P);
%     for k=1:window_width
%         St=St+(Y(k,:))'*(Y(k,:));
%     end
% 
%     % Normalize covariance matrix
%     St=St./sum(weight.^2);
%     
% %     % Rescale covariance matrix ( Dont think this works..)
% %     
% %     V=diag(1./sqrt(diag(St)));
% %     St=V*St*V';
% 
%     S=cat(3,S,St);  % Construct Tensor for S
    
%     %% Dont use kernel (added 30/10/2014)
    S=zeros([P,P,T]);
    for ts=1:T
    S(:,:,ts)=y(ts,:)'*y(ts,:)/2;
    end
    
end


% Construct new time indexing T_s:T_w->1:T' , where T'=T-window_width
Tp=T-window_width+1;
[~,~,Tp]=size(S);   % Automatically length find from S

Z=zeros([P,P,Tp]);    % Initialise dual and aux variables
U=Z;

Theta=zeros([P,P,Tp]);
for t=1:Tp
Theta(:,:,t)=eye(P);        % Intiialise theta
end

fobj=[0,findF(Theta,S,lambda1,lambda2)];    % Initiate f;

% Initialize convergence criteria
difPrime=tolPrime+1;
difDual=tolDual+1;
conv=0;
%gamma=gammaMin;

aditer=2;

%% MAIN ADMM LOOP
% while ((abs(fobj(aditer)-fobj(aditer-1))>tol || aditer==2 || aditer<minit ) && aditer<maxit)
% while(conv~=1)
while(difPrime>tolPrime && difDual > tolDual && aditer<maxit)
    % Solve step 1 through eigendecomposition
    for t=1:Tp
                  %[V,Sd]=eig(squeeze(S(:,:,t)-gamma*(Z(:,:,t)-U(:,:,t)) ),'nobalance');
                  [V,Sd]=eig(squeeze(S(:,:,t)-gamma*(Z(:,:,t)-U(:,:,t)) ));

        E=[];
        for r=1:P
            sr=Sd(r,r); % Get rth eigenvalue
              thetar=(-sr+sqrt(sr^2+4*gamma))/(2*gamma); % Solve quadratic for eigenvalues 
            E=[E,thetar];
        end
        E=diag(E);
        Theta(:,:,t)=V*E*V';    % Reconstruct theta
        
    end
    
    
    Zold=Z; % For convergence checking
     Target=U+Theta;
     % Options for fused lasso signal approximator
         % If we need to do smoothing as well
         if(lambda2>0)
            opts=struct('fusedPenalty',lambda2/gamma);
            % Solve for each edge seperately...
             for i=1:P-1
                 for j=i+1:P
                     [Z(i,j,:), ~, ~]=fusedLeastR(eye(Tp), squeeze(Target(i,j,:)), lambda1/gamma, opts);
                    Z(j,i,:)=Z(i,j,:);  % Make symmetric
                 end
             end
         else
             Z=wthresh(Target,'s',lambda1/gamma);
         end
      
      % Update diagonals
      for t=1:Tp
          for i=1:P
         Z(i,i,t)=Theta(i,i,t)+U(i,i,t); 
          end
      end
      
    % Update dual variable
    Uold=U;
     U=Uold+(Theta-Z); % Yuan version
     

%         U=U+(Theta-Z)

    % Check for psd (DEBUG)
    % Check dual and primal feasbility
    difDual=0;
    difPrime=0;
    for t=1:Tp
        % Dual feasability
        difDual=difDual+norm(Z(:,:,t)-Zold(:,:,t))^2;
        
        %Primal feasability
        difPrime=difPrime+norm(Theta(:,:,t)-Z(:,:,t))^2;
        
        % Positive semi-definite Z
        [~,psd]=chol(Z(:,:,t));   % Check for cholesky decomposition
        if(psd~=0)
            display('outpuyt NOT PSD');
        end
    end

    % Update "old" estimates
    Zold=Z;
    
    % Calculate objective
    fnew=findF(Theta,S,lambda1,lambda2);

    fobj=[fobj,fnew];
%     % Check convergence!!
%     if (difPrime < tolPrime && difDual < tolDual)
%        conv=1; 
%     end
    
    % Update gamma according to Boyd strategy
    %gamma=min([1.1*gamma,gammaMax]); % This destroys convergence theory...
    aditer=aditer+1;

end

end
%% END OF MAIN LOOP
    
%% Construct matrix rep
function Thetam=constructMatrix(Theta,P)

[~,~,T]=size(Theta);

 Thetam=zeros([P*(P-1)/2,T]);
% Use full matrix size rather than making use of symmetry (for debugging)
%Thetam=zeros([P*P,T]);
for t=1:T
    k=1;
    % Puts all col/rows in matrix
%    for i=1:P  % iterate over rows
%        for j=1:P % iterate over cols
%            Thetam(k,t)=Theta(i,j,t);
%            k=k+1;
%        end
%    end

    % Puts only edges in matrix
    for i=1:P-1
        for j=i+1:P
            Thetam(k,t)=Theta(i,j,t);
            k=k+1;
        end
    end
end
end

%% Reverse operation to tensor
function Theta=constructTensor(Thetam,Theta0,P)
[~,T]=size(Thetam);

% Theta=Theta0;
% Only changes offdiagonal entries

for t=1:T
    k=1;
   for i=1:P-1  % iterate over rows
       for j=i+1:P % iterate over cols
           Theta(i,j,t)=Thetam(k,t);
           Theta(j,i,t)=Thetam(k,t);
           k=k+1;
       end
   end
end

end

function f=findF(Theta,S,lambda1,lambda2)
[~,P,T]=size(Theta);

mask=ones([P,P]);

mask=mask-eye(P);  % For masking whole problem ie calculating objective

        f=-log(det(squeeze(Theta(:,:,1))))+trace(squeeze(S(:,:,1))*squeeze(Theta(:,:,1)))...
            + lambda1*norm(mask.*squeeze(Theta(:,:,1)),1);
    for t=2:T
        f=f-log(det(squeeze(Theta(:,:,t))))+trace(squeeze(S(:,:,t))*squeeze(Theta(:,:,t)))...
            + lambda1*norm(mask.*squeeze(Theta(:,:,t)),1)...
            + lambda2*norm(mask.*squeeze(Theta(:,:,t)-Theta(:,:,t-1)));
    end
end

%% Reconstructs signal from gflasso Bleakley
function Yest=reconstruct(res,Y,T)

beta=res.value{1,1};    % Value of beta at jump point
active=res.jump{1,1};   % Active jump points
K=length(active);   % Number of changepoints
dg=size(beta,2);
% 
% if(option.weights)
%     X=ones([T,T-1]);
% else
X=zeros([T,T-1]);

for i=1:T
    for j=1:T-1
        if(i>j)
        X(i,j)=sqrt(T/(j*(T-j)));   % For boundary effects
        else
            X(i,j)=0;
        end
    end
end
% end

Beta=zeros([T-1,dg]);
for k=1:K
    Beta(active(k),:)=beta(k,:);
end

offset=ones([1,T])*(Y-X*Beta)./T;

% Center X
Yest=ones([T,1])*offset + X*Beta;

end
