par.b = [1,.7];
par.mu = [1,-1];
par.sigma = [1,.5];
par.a = 0.1;
par.x0 = 0;
par.D = 0.1;
par.X = 3;
par.pts = 241;

tt = [0 1 4 10];

[u0,xx] = eqdist(par);

b = [.7,.3,.7];
mu = [1,0,-1];
sd = [1,1,.5];

soln = pdesim(b,mu,sd,u0,xx,tt,par);

I = length(soln(:,1,1));

figure(1)
clf
hold on
box on
xlim([-par.X,par.X])
for i=1:I
    plot(xx,soln(i,:,1),'LineWidth',1.2)
end
par.b = b;
par.mu = mu;
par.sigma = sd;
[u1,xx] = eqdist(par);
plot(xx,u1,'k-.','LineWidth',1.2)
legend('$t = 0$','$t = 1$','$t = 4$','$t = 10$','$t \to \infty$','interpreter','latex','FontSize',6,'location','northwest')
xlabel('$x$','interpreter','latex','FontSize',14)
ylabel('$u$','interpreter','latex','FontSize',14,'Rotation',0)


