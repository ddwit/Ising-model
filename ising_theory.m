function [Es,Ms] = ising_theory(J)

% [Es,Ms] = ising_theory(J)
%
% ising_theory computes the theoretical values of the magnetisation Ms per spin
% and the energy per spin Es for a 2D Ising model in the limit of an infinitely 
% large lattice, assuming that the background magnetic B is zero. 
%
%	J 	dimensionless interaction strength 0 ≤ J (array), taking kT =1
%	Es	energy per site (vertical array)
%	Ms	magnetization per site (vertical array with two columns, one per branch)

% T. Dudok de Wit, 3/2019
% See for example https://notendur.hi.is/jeg1/Ising.pdf


if any(J<0)
	error('** J must be ≥0 **')
end

Jc = log(1+sqrt(2))/2;  % critical value below which Ms vanishes
J = J(:);

kappa =  2.*sinh(2*J)./(cosh(2*J)).^2;
ksi = 2*(tanh(2*J)).^2-1;
z = exp(-2*J);

% Calculating Energy & Magnetisation as a function of interaction strength

K1 = ellipke(kappa);
Es = -J.*(coth(2*J)).*(1+(2/pi)*ksi.*K1)*2;
z2 = z.*z;
Ms =  ((1+z2).^(1/4).*(1-6*z2+z2.*z2).^(1/8))./(1-z2).^(1/2);
Ms = [Ms -Ms];  % two branches of M

ind = find(J<Jc);
if ind
	Ms(ind,:) = 0; 
end



