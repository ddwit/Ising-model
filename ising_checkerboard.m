function [ind1,ind2] = ising_checkerboard(N)

% [ind1,ind2] = ising_checkerboard(N)
%
% For a 2D Ising lattice, ising_checkerboard computes the indices of the
% spins that are distributed on checkerboard (i.e. one out of two). This
% function is used by ising_metropolis to avoid flipping simultaneously a 
% spin and on of its nearest neighbours.
%   N       number of spins per dimension
%   ind1    indices of elements that are located on a checkerboard [array]
%   ind2    idem, but on complementary checkerboard [array]

%   T. Dudok de Wit, 3/2019

parity = mod(1:N,2);         % parity
r = repmat(parity,  [N 1]);  % expand the parities to rows
c = repmat(parity', [1 N]);  % expand the parities to columns
M = xor(r,c);                % apply the xor operator

ind1 = find(M>0);
ind2 = find(M<1);