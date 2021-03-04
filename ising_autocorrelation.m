function [C,x] = ising_autocorrelation(spin,displ)

% [C,x] = ising_autocorrelation(spin,displ)
%
% ising_autocorrelation computes some of the autocorrelation
% function of an Ising lattice along vertical and horizontal directions
% 		spin        2D matrix of +/- 1's
%       displ       disply results if displ>0, Default is with
%		C 			autocorrelation function
%		x			corresponding lag

%   TDdW 3/2019

if nargin<2, displ = 1; end

N = size(spin,1);

y = spin;
y = y - mean(y(:));
y = y / std(y(:));

maxlag = ceil(N/5);
x = 0:maxlag;
nx = maxlag+1;
C = zeros(maxlag,2*N);

for i=1:nx
    lag = x(i);
    w = 1:nx-lag;
    nw = length(w);
    for k = 1:N
        C(i,k)   = spin(w,k)'*spin(w+lag,k)/nw;
        C(i,k+N) = spin(k,w)*spin(k,w+lag)'/nw;
    end
end
Cave = mean(C,2);
Cstd = std(C')';

if displ
    clf
    h = plot(x,Cave,'-',x,Cave+Cstd/2,'--',x,Cave-Cstd/2,'--');
    set(h(2),'Color',get(h(1),'Color'));
    set(h(3),'Color',get(h(1),'Color'));
    set(h(1),'Linewidth',1.4);
    grid on
    xlabel('lag')
    ylabel('autocorrelation fctn')
    title(['N = ',int2str(N)])
    axis([0 x(end) -1 1])
end
    
