function y=GMVP(sigma)

% Global Minimum Variance Portfolio (uncstr)

N=size(sigma,2);
I=ones(N,1);
C=I'*inv(sigma)*I;
y=(inv(sigma)*I)/C;   