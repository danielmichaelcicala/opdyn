par.b = [1,.7];
par.mu = [1,-1];
par.sigma = [1,.25];
par.a = 0.1;
par.x0 = 0;
par.D = 0.1;
par.X = 3;
par.pts = 241;

figure(1)
clf
hold on
box on
xlabel('$x$','interpreter','latex','FontSize',14)
xlim([-par.X par.X])

[u,x] = eqdist(par);

plot(x,u,'LineWidth',1.5)
xlim([-3 3])
ylabel('$u^*$','interpreter','latex','FontSize',14,'Rotation',0)

par.b = [.6,.4,.6];
par.mu = [1,0,-1];
par.sigma = [1,1,.25];
[u,x] = eqdist(par);
plot(x,u,'k-.','LineWidth',1.5)



