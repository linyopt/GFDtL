function [ sigmainv,sigma ] = genPeriodGraph( P,M,scale,K,T )
%GENPERIODGRAPH Summary of this function goes here
%   Generates precision matrices accoridng to a
%   periodic graph sturcture with:
%   P - dimension
%   M - Number of edges
%   scale - size of off-diagonal entries
%   K - number of changepoints
%   T - total time-series length
%
%   JCGS version
%   Alex Gibberd 1/4/2015 - UCL Department of Statistical Science

% Generate initial graph
[ inv0 ] = geninv( P,M,scale );
sigma0=inv(inv0);  % calculate inverse of new sigma
% Generate alternate graph
inv1=addrmp(inv0,M,M,scale);
sigma1=inv(inv1);

sigmainv=zeros([P,P,T]);
sigma=sigmainv;

if(mod(T,K+1)~=0)
   display('Time length not divisible by number of Changepoints');
   return
else
    Tphase=T/(K+1);
end
% Find breakpoints
bpc=floor(linspace(1,T,K+2));
bpc(end)=T+1;   % To correct for spacing to be inclusive

% Loop over changepoints alternating graph structure

for k=1:K+1
    if(mod(k,2)==0)
        for t=bpc(k):bpc(k+1)-1
            sigmainv(:,:,t)=inv0;
            sigma(:,:,t)=sigma0;
        end
    else
        for t=bpc(k):bpc(k+1)-1
            sigmainv(:,:,t)=inv1;
            sigma(:,:,t)=sigma1;
        end
    end
end

end

