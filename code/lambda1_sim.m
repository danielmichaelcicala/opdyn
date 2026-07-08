% runs a simulation for the opinion dynamics model

% example 1 is the "rfk" scenario:
% the initially positive medical community local attractor splits off
%   a neutral government establishment local attractor

% the black dash-dot curve is the equilibrium solution for t -> infty

% Glenn Ledder

% uses eqdist.m, 2025/07/28
% uses pdesim.m, 2025/07/28

%% DATA

rvals = [0.08,0.09,0.10];

% location of right interval boundary
par.X = 3;
% number of spatial mesh points
par.pts = 241;

% total time
par.tf = 20;
% time increment
par.dt = .5;

% local attractor strengths at time 0
par.b = [1,.7];
% local attractor locations at time 0
par.mu = [1,-1];
% local attractor standard deviations at time 0
par.sigma = [1,.5];

% local attractor strengths at time t>0
b = [.7,.3,.7];
% local attractor locations at time t>0
mu = [1,0,-1];
% local attractor standard deviations at time t>0
sd = [1,1,.5];

% global attractor strength
par.a = 0.1;
% global attractor location
par.x0 = 0;
% diffusion coefficient
par.D = 0.1;

%% INITIALIZATION

% find the equilibrium solution for the t=0 parameters
% This will be the initial condition for the pde.

[u0,xx] = eqdist(par);
I = 1+par.tf/par.dt;
tt = linspace(0,par.tf,I);
Y = zeros(1,I);

%% COMPUTATION

% find the pde solution
% soln(i,j,1) is the value of u and time t_i and position x_j
% position x_j is determined from X and pts in eqdist.m
% time t_i is determined from tf and dt in pdesim.m

soln = pdesim(b,mu,sd,u0,xx,tt,par);

par.b = b;
par.mu = mu;
par.sigma = sd;
[ustar,xx] = eqdist(par);

%% PLOT SETUP

figure(3)
clf
for p=1:2
    subplot(1,2,p)
    hold on
    box on
    if p==1
        xlabel('$t$','interpreter','latex','FontSize',14)
        ylabel('$||Q||$','interpreter','latex','FontSize',14,'Rotation',0)
    else
        xlabel('$x$','interpreter','latex','FontSize',14)
        ylabel('$\frac{Q}{||Q||}$','interpreter','latex','FontSize',14,'Rotation',0)
        xlim([-par.X,par.X])
    end
end

%% EXPERIMENT

subplot(1,2,1)

for ir=1:length(rvals)
    r = rvals(ir);
    expfactor = exp(r*tt);
    for i=1:I
        deltau = soln(i,:,1)-ustar;
        Q = expfactor(i)*deltau;
        f = Q.^2;
        nnorm = sqrt(Simpson(f,2*par.X));
        Y(i) = nnorm;
    end
    plot(tt,Y,'linewidth',1)
    if ir==(1+length(rvals))/2
        phi1 = Q/nnorm;
    end
end

ylim([0,Y(1)])
legend('r=0.08','r=0.09','r=0.10','location','southwest')

subplot(1,2,2)
plot([-par.X,par.X],[0,0],'k:','linewidth',1)
plot(xx,phi1,'linewidth',1)


%% function Simpson

% computes an integral from function values f for equally spaced 
%   nodes with total interval width w.

function y=Simpson(f,w)
pts = length(f);
n = length(f)-1;
h = w/n;
c = 2*(1-mod(1:pts,2));
% convert f to a column vector if needed
if isscalar(f(1,:))
    ff = f;
else
    ff = f';
end
y = h*(2*sum(f)+c*ff-f(1)-f(end))/3;
end
