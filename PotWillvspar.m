% This file plots the Potential Willing Population versus model parameters
% Pstar = int(p*u_equist,[-X,X])
% In this file, we take p = 1*(x>0)+0*(x<0)
% that is, P = int(p*u_equist,[0,X])
%
% j=1: professional, j=2: denier, and j=3: government
% x_j (mu_j): opinion position of influencer j 
% b_j: strength of x_j opinion
% sigma_j: standard deviation of x_j opinion
%
% Plot (I):
% Pstar vs x_g (mu_3)
% 
% Plot (II): contour plot
% Pstar vs b_p (b_1) and sigma_p (sigma_1)


% The integration is approximated by a Simpson's Rule. 

% Jacky Jiang
% 10/27/2025

close all

clc
% Default values for 2025 case 
X = 3; % interested interval [-X,X]
pts = 241;
dx = 2*X/pts; 
D = 0.1;
x0 = 0;
a = 0.1;
mu = [1,-1,0];  % x_j 
b = [.7, .7, .3]; 
sigma = [1,.5,1];

% input for equilibrium solu solver
par.X = X;
par.pts = pts;
par.D = D;
par.x0 = x0;
par.a = a;
par.mu = mu;
par.b = b;
par.sigma = sigma;


%%%%% Pstar vs b_p and sigma_p %%%%%

b_p = 0.4:0.01:1;
sigma_p = 0.5:0.01:1.5;
b1 = b;
sig1 = sigma;
Pstar =[];
for j = 1:length(sigma_p)
    sig1(1) = sigma_p(j);
    par.sigma = sig1;
    for i = 1:length(b_p)
        b1(1) = b_p(i);
        par.b = b1;
        [u,x] = eqdist(par);
        mid = ceil(par.pts/2);
        Pstar(i,j) = Simp(u(mid:end),X);
    end
end
figure(1)
% subplot(1,2,1)
[sY,bX]=meshgrid(sigma_p,b_p);
levels = 0.2:0.1:0.9;
contourf(bX,sY,Pstar,levels);
colorbar
hold on
[C,h] = contour(bX,sY,Pstar,levels, 'ShowText', 'on', 'LineColor', 'k');
clabel(C, h, 'FontSize', 20);
ax = gca;
ax.FontSize = 18;
xlabel('$b_p$','interpreter','latex','FontSize',30,'Rotation',0)
ylabel('$\sigma_p$','interpreter','latex','FontSize',30,'Rotation',0)
% set(gcf, 'Position',get(0,'Screensize'));
saveas(gcf,'Pstarvsbpsigp','epsc') % 'epsc'
min(Pstar,[],'all')
max(Pstar,[],'all')



%%%%% Pstar vs x_g %%%%%
par.b = b;% recover to default values
par.sigma = sigma;% recover to default values

x_g = -1:0.01:1;
mu1 = mu;
Pstar =[];
%P1 = [];
for i = 1:length(x_g)
    mu1(3) = x_g(i);
    par.mu = mu1;
    [u,x] = eqdist(par);
    mid = ceil(par.pts/2);
    Pstar(i) = Simp(u(mid:end),X);% using Simpson rule
    %P1(i) = sum(u(mid:end-1))*dx;% using left-end point
end
figure(2)
% subplot(1,2,2)
plot(x_g,Pstar,'LineWidth',1.5)
% hold on
% plot(x_g,P1,'r+','LineWidth',1.5)
% hold off
%title('$P^*_+$ vs Government Opinion','interpreter','latex','FontSize',20)
xlim([-1 1])
xlabel('$x_g$','interpreter','latex','FontSize',30,'Rotation',0)
ylabel('$P^*_+$','interpreter','latex','FontSize',30,'Rotation',0)
% set(gcf, 'Position',get(0,'Screensize'));
saveas(gcf,'Pstarvsxg','epsc')

% saveas(gcf,'paramstudies','epsc')


% Numerical Integration
function y=Simp(f,w)
pts = length(f);
n = length(f)-1;
h = w/n;
c = 2*(1-mod(1:pts,2));
% convert f to a column vector
if isscalar(f(1,:))
    ff = f;
else
    ff = f';
end
y = h*(2*sum(f)+c*ff-f(1)-f(end))/3;
end


