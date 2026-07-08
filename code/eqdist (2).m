function [u,x] = eqdist(par)

% finds equilibrium distribution for the opinion dynamics model
%   u_t=f_x, where f=Du_x-g(x)*u

% u is the equilibrium solution values
% x is the list of spatial points

% global attractor is at x0 with strength a
% local attractors are at mu_j with width sigma_j and strength b_j

% D is the diffusion coefficient
% X is the value of x on the right boundary
% pts is the number of nodes in the spatial domain

% G is the integral of g
% eta is the PDE solution before normalization

% Glenn Ledder
% 2025/07/28

%% SAMPLE DATA

% X = 3;
% pts = 121;
% D = 0.1;
% x0 = 0;
% a = 0.1;
% b = [1,.5];
% mu = [1,-1];
% sigma = [1,.25];

%% DATA

X = par.X;
pts = par.pts;
x0 = par.x0;
D = par.D;
a = par.a;
b = par.b;
mu = par.mu;
sigma = par.sigma;

%% COMPUTATION

G = @(x) G_global(x,x0,a)+G_local(x,mu,sigma,b);
eta = @(x) exp(G(x)/D);

I = integral(eta,-Inf,Inf);

x = linspace(-X,X,pts);

u = eta(x)/I;

%% function G_global

function G=G_global(x,x0,a)

% Finds the integral function G for g = -2a(x-x0)

% The interval should have the form -X:dx:X where X is
% large enough that P->0.

G = 0.5*a*x.*(2*x0-x);

end

%% function G_local

function G=G_local(x,x0,s,b)

% Finds the integral function G for g = -2(b/s^2)(x-x0)exp(-(x-x0)^2/s^2)

% The interval should have the form -X:dx:X where X is
% large enough that P->0.

L = length(x0);
G = 0;
for j=1:L
    z = (x-x0(j))/s(j);
    z0 = -x0(j)/s(j);
    G = G+b(j)*(exp(-z.^2)-exp(-z0.^2));
end
G = 0.5*G;

end

end

