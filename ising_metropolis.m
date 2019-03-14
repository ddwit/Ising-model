function [spin,Es,Ms] = ising_metropolis(spin,J,niter,displ)

% [spin,Es,Ms] = ising_metropolis(spin,J,niter,displ)
% 
%   ising_metropolis runs the Metropolis-Hastings algorithm on a
%   configuration of spins with a coupling coefficient J at a
%   temperature kT. spin is a matrix of +/- 1's.
%       spin    [N,N]  matrix of initial spins, as obtained from
%               ising_initialisation
%       niter   number of iterations, ideally niter N*N*2^8
%       J       interaction strength
%       displ   display iterations every displ'th time step if displ>0. 
%               Default displ = 1
%
%       spin    [N,N] matrix of final spin state
%       Es      [niter,1] time series of energy per spin
%       Ms      [niter,1] time series of magnetisation per spin

%   inspired by algorithm by Tobias Fricke
%   T. Dudok de Wit, 3/2019


if nargin<4, displ = 1; end
if nargin<3, niter = 500; end
if displ, colormap(bone); end


kT = 1;                 % k * T
flip_fraction = 0.1;    % fractional amount of spins that may flip
                        % at each iteration; should be << 1
                        % larger implies more frequent updates of lattice

N = size(spin,1);
N2 = N*N;

Ms = zeros(niter,1);
Es = zeros(niter,1);

% generate indices for selecting four nearest neighbours
index = reshape(1:N*N,N,N);
indexR = circshift(index,[ 0  1]);  % for selecting neighbours on the right
indexL = circshift(index,[ 0 -1]);
indexT = circshift(index,[ 1  0]);
indexB = circshift(index,[-1  0]);

% flip only those spins that are located on a checkerboard. This is to avoid
% flipping simultaneously a spin and some of its four nearest neigbours
[index1,index2] = ising_checkerboard(N);

% iterate to burn-in
for i=1:niter,
      
    % switch checkerboard at each iteration
    if mod(i,2)==1
        index = index1;
    else
        index = index2;
    end
    Nindex = length(index);
        
        
    neighbours = spin(indexR) + spin(indexL) + spin(indexT) + spin(indexB);

    % change in energy when flipping a spin
    deltaE = 2*J*spin.*neighbours;
    
    % transition probability (only on checkerboard)
    prob = exp(-deltaE(index)/kT);
    
    % build mask with transitions that will occur. This is the most
    % time-consuming operation
    flip = (rand(Nindex,1) < prob).*(rand(Nindex,1) < flip_fraction)*(-2) + 1;
    
    % flip spins that are located on checkerboard
    spin(index) = spin(index).*flip;
    
    % sum up the variables of interest
    Ms(i) =  sum(sum(spin))/N2;
    Es(i) = -sum(sum(deltaE))/2/N2; 
        
    % if desired, display the current state of the system every displ time
    % steps
    if displ & mod(i,displ)==0
        clf
        image((spin+1)*128);
        title(sprintf('# %d   J=%0.2f  Ms=%0.2f   Es=%0.2f',i, J, Ms(i), Es(i)));
        set(gca,'YTickLabel',[],'XTickLabel',[]);
        axis square
        drawnow
    end
end