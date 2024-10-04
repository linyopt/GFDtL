function [ y ] = simulateData( sigma )
%SIMULATEDATA Produces a realisation for a given MV Gaussian noise process

[P,P,T]=size(sigma)
y=[];
for t=1:T
y=[y;randn([1,P])*squeeze(sigma(:,:,t))];
end

end

