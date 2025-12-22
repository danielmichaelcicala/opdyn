% runs a simulation for the opinion dynamics model

% example 1 is the "rfk" scenario:
% the initially positive medical community local attractor splits off
%   a neutral government establishment local attractor

% the black dash-dot curve is the equilibrium solution for t -> infty

% Glenn Ledder


%% DATA

% R0 sign determines sign of eigenfunction slope at x=-X
R0 = -1;

% 1 to show eigenfunction phi2 or 0 not to
phi2 = 0;

first = 0.05;
last = 1;
pts = 200;
lambdavals = linspace(first,last,pts);

% location of right interval boundary
par.X = 3;

% local attractor strengths at time t>0
par.b = [.7,.3,.7];
% local attractor locations at time t>0
par.mu = [1,0,-1];
% local attractor standard deviations at time t>0
par.sigma = [1,1,.5];

% global attractor strength
par.a = 0.1;
% global attractor location
par.x0 = 0;
% diffusion coefficient
par.D = 0.1;

%% COMPUTATION

rhovals = zeros(1,pts);
for i=1:pts
    lambda = lambdavals(i);
    rho = rhofnc(lambda,R0,par);
    rhovals(i) = rho;
end

objfnc = @(lambda) rhofnc(lambda,R0,par);
lambda1 = fzero(objfnc,0.12)
[~,xx1,R1,~] = rhofnc(lambda1,R0,par);
lambda2 = fzero(objfnc,0.32)
[~,xx2,R2,~] = rhofnc(lambda2,R0,par);


%% OUTPUT

clf
for k=1:2
    subplot(1,2,k)
    hold on
    box on
end
C = colororder;
colour = C(2,:);

subplot(1,2,1)
plot(lambdavals,rhovals,'linewidth',1)
plot([0,last],[0,0],'k-.','linewidth',1)
xlabel('$\lambda$','interpreter','latex','FontSize',14)
ylabel('$S(X)$','interpreter','latex','FontSize',12,'Rotation',0)

subplot(1,2,2)
plot(xx1,R1,'linewidth',1,'color',colour)
if phi2==1
    plot(xx2,R2,'linewidth',1)
end
xlim([-par.X,par.X])
xlabel('$x$','interpreter','latex','FontSize',14)
ylabel('$\phi_1$','interpreter','latex','FontSize',14,'Rotation',0)

