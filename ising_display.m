function ising_display(spin)

% ising_display(spin)
%
% ising_display displays a 2D lattice of spins
% spin is a 2D matrix of +/- 1's.

%   TDdW 3/2019


image((spin+1)*128);
N = size(spin,1);           % nr of spins per dimension
Ms = sum(sum(spin))/(N*N);  % magnetisation per spin

set(gca,'YTickLabel',[],'XTickLabel',[]);
axis square
colormap bone
title(sprintf('N = %d   Ms = %0.2f', N, Ms));
drawnow

