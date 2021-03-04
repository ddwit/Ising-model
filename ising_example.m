%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ising_example
%
% example script to simulate a 2D Ising model using the
% Metropolis-Hastings algorithm
%
% T. Dudok de Wit, 3/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start with simple example

displ = 1;      % turn display on
N = 50;         % number of spins per dimension
J = 0.4;        % dimensionless coupling strength
niter = 200;    % number of iterations

% initialise a given configuration with on average as
% many ups and downs
spin0 = ising_initialisation(N,0.5);

% visualise the initial state
ising_display(spin0)

% start iterating the Metropolis Hastings algorithm
[spin,Es,Ms] = ising_metropolis(spin0,J,niter,displ);

% compare the states before and after burn-in
disp('** initial state (random) **')
ising_stats(spin0,displ);

disp('** final state (after burn-in) **')
ising_stats(spin,displ);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now a more realistic simulation with a large lattice (takes a few minutes)

displ = 50;     % display every 50th iteration
N = 500;        % number of spins per dimension
J = 0.8;        % dimensionless coupling strength
niter = 10000;  % number of iterations

spin0 = ising_initialisation(N,0.5);
[spin,Es,Ms] = ising_metropolis(spin0,J,niter,displ);

% watch how energy and magnetisation evolve in time
clf
subplot(211)
Emin = -4*J;
plot(1:niter,Es,".-",[1 niter],[Emin Emin],'--')
ylabel('Es')
xlabel("time steps")
title('Energy per spin')
text(niter*1.01,Emin,'E_{min}')
grid on

subplot(212)
plot(Ms,".-")
ylabel('Ms')
xlabel("time steps")
title('Magnetisation per spin')
grid on


