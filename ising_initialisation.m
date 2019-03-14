function spin = ising_initialisation(N,p)

% spin = ising_initialisation(N,p)
%
%   ising_initialisation initialises a 2D Ising lattice by
%   returning a configuration with N spins along each 
%   dimension and a proportion p (0≤p≤1) of them
%   pointing upwards. 
%		spin  	N x N matrix of +/- 1's
%   	p       probability of having spins pointing upwards
%				must be 0≤p≤1. Default is p=0.5 

%   T. Dudok de Wit 3/2019



if nargin<2,	p = 0.5;	end

if p<0 | p>1, error('** 0 ≤ p ≤ 1  is required **'), end

spin = sign(p - rand(N,N));
