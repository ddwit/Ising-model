function [E,Es,M,Ms,chi,chis] = ising_stats(spin,J,displ)

% [E,Es,M,Ms,chi,chis] = ising_stats(spin,J,displ)
%
% ising_stats computes some of the main statistical
% quantities associated with a 2D Ising lattice
% 		spin        N x N matrix of +/- 1's
%       J           dimensionless coupling strength
%       displ       display results if displ>0, Default is with display
%       E           total energy
%       Es          energy per spin
%		M 			total magnetisation
%		Ms 			magnetisation per spin
%		chi			total magnetic susceptibility
%		chis	    magnetic susceptibility per spin

%   T. Dudok de Wit, 3/2019

if nargin<3, displ = 1; end

N = size(spin,1);
N2 = N*N;


%%%%% magnetisation and susceptibility

M = sum(sum(spin));
Ms = M/N2;
z = spin - Ms;
chi = sum(sum(z.*z));
chis = chi/N2;


%%%%% now compute the energy

neighbors = circshift(spin, [ 0  1]) + ...
            circshift(spin, [ 0 -1]) + ...
            circshift(spin, [ 1  0]) + ...
            circshift(spin, [-1  0]);

% change in energy of flipping a spin
dE = 2 * J * (spin .* neighbors);
E = -sum(sum(dE))/2;
Es = E/N2;
    
    
    
%%%%% With the image toolbox, uncomment the latter to compute connectivity

% [L,Nconnect] = bwlabel(spin == 1, 4); 




%%%%% display results

if displ
    fprintf('\n')
    disp(['Lattice size                         N    = ',int2str(N)])
    disp(['Total energy                         E    = ',num2str(E)])
    disp(['Energy per spin                      Es   = ',num2str(Es)])
    disp(['Total magnetisation                  M    = ',num2str(M)])
    disp(['Magnetisation per spin               Ms   = ',num2str(Ms)])
    disp(['Motal magnetic susceptibility        chi  = ',num2str(chi)])
    disp(['Magnetic susceptibility per spin     chis = ',num2str(chis)])
    fprintf('\n')
end
