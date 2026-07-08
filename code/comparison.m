par.b = [1,.7];
par.mu = [1,-1];
par.sigma = [1,.5];
par.a = 0.1;
par.x0 = 0;
par.D = 0.1;
par.X = 3;
par.pts = 241;

figure(3)
clf
for k=1:2
    subplot(1,2,k)
    hold on
    box on
    xlabel('$x$','interpreter','latex','FontSize',16)
    xlim([-par.X par.X])
end

[u,x] = eqdist(par);

plot(x,u,'LineWidth',1.5)
ylabel('$u^*$','interpreter','latex','FontSize',16,'Rotation',0)

par.b = [.7,.7,.3];
par.mu = [1,-1,0];
par.sigma = [1,.5,1];
[u,x] = eqdist(par);
plot(x,u,'k-.','LineWidth',1.5)

subplot(1,2,1)

a = 0.1;
D = 0.1;
x0 = 0;
X = 3;
x = -X:.02:X;

ylabel('$\dot{x}$','interpreter','latex','FontSize',16,'Rotation',0)
for k=1:2
    if k==1
        b = [1,.7];
        sigma = [1,.5];
        mu = [1,-1];
    else
        b = [.7,.7,.3];
        sigma = [1,.5,1];
        mu = [1,-1,0];
    end

    g = @(x) g_global(x,x0,a)+g_local(x,mu,sigma,b);

    if k==2
        plot(x,g(x),'k-.','LineWidth',1.5)
    else
        plot(x,g(x),'LineWidth',1.5)
    end
    x0 = fzero(g,-1)
    plot(x0,0,'ko','MarkerSize',6)
    plot(x0,0,'b.','MarkerSize',15)
    x0 = fzero(g,1)
    plot(x0,0,'ko','MarkerSize',6)
    plot(x0,0,'b.','MarkerSize',15)
    x0 = fzero(g,-.3)
    for s=3:7
        plot(x0,0,'rs','MarkerSize',s)
    end

end


%% function g_global

function g=g_global(x,xj,a)

% Finds the integral function F for f = -2a(x-xj)

% The interval should have the form -X:dx:X where X is
% large enough that P->0.

g = -a*(x-xj);

end

%% function g_local

function g=g_local(x,x0,s,b)

% Finds the integral function F for f = -2(b/s^2)(x-x0)exp(-(x-x0)^2/s^2)

% The interval should have the form -X:dx:X where X is
% large enough that P->0.

L = length(x0);
g = 0;
for j=1:L
    z = (x-x0(j))/s(j);
    g = g-b(j)*z.*exp(-z.^2)/s(j);
end

end

