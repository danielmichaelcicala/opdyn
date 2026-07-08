%
% plot_equilibrium_vs_xg.m
%
% Sweep the global attractor location x_g from -1 to 1 and
% plot the equilibrium distribution u*(x) for each value.
%

%% Base parameters
params.D = 0.1;

% Global attractors (only the first is active here; others set to zero)
params.xg = [0, 0, 0];
params.a  = [0.3, 0, 0];

% Local attractors
params.xl    = [  1,   -1,   0];
params.b     = [0.7,  0.7, 0.3];
params.sigma = [  1, 0.5,   1];

%% Sweep xg over [-1, 1]
xl3_vals = linspace(-1, 1, 10);  

figure; hold on;
set(gcf, 'Color', 'w');
cmap = parula(length(xl3_vals));   % one color per curve

% Precompute x-grid once
xplot = linspace(-3, 3, 400);

for k = 1:length(xl3_vals)
    params.xl(3) = xl3_vals(k);  % move the govn't attractor
    [x, ueq] = function_equilibrium_distribution(params, xplot);
    plot(x, ueq, 'LineWidth', 1.6, 'Color', cmap(k,:));
end

ax=gca; 
ax.FontSize = 18;
xlabel('$x$', 'FontSize',30, 'Interpreter','latex', 'Rotation',0);
xlim([-2 2]);
ylabel('$u^*$', 'FontSize',30, 'Interpreter', 'latex','Rotation',0);
%ylim([0 1.5]);
set(gcf,'Position',[50 50 600 600])
% title('Equilibrium distributions for varying x_g \in [-1, 1]');
grid off; box on;

% Colorbar
colormap(cmap);
cb = colorbar;
cb.FontSize = 14;
clim([min(xl3_vals) max(xl3_vals)]);
ylabel(cb, '$x_g$','FontSize',30, 'Interpreter','latex', 'Rotation',0);

%% Export to EPS
print('-depsc', 'equilibrium_vs_xg.eps');
disp('EPS file saved as "equilibrium_vs_xg.eps" in current folder.');

%% ---------------- Helper functions ----------------

function [x,y] = function_equilibrium_distribution(params, xgrid)
    % Compute normalized equilibrium u*(x) ∝ exp(G(x)/D)
    if nargin < 2
        xgrid = linspace(-3, 3, 400);
    end
    D = params.D;

    % Antiderivative of g(x): G(x) = ∫ g(x) dx
    G = @(x) function_integrate_g(params, x);   % returns same shape as x

    % Unnormalized density (vectorized; preserves input shape)
    e = @(x) exp(G(x) / D);

    % Normalizer over plotting domain
    normalizer = integral(@(x) e(x), min(xgrid), max(xgrid), ...
                          'RelTol',1e-8,'AbsTol',1e-12);

    % Normalized equilibrium
    u = @(x) e(x) / normalizer;

    x = xgrid(:).';          % return row vector
    y = u(x);
end

function y = function_integrate_g(params,x)
    
    %% Parameters
    
    xg    = params.xg;    % global attractor locations
    a     = params.a;     % global attractor strengths
    xl    = params.xl;    % local attractor locations
    b     = params.b;     % local attractor strengths
    sigma = params.sigma; % local kernel sd's


    % Initialize output
    y = NaN(size(x)); 
    
    %% Calculate
    for i = 1:length(x)
    
        ex = exp(-( (x(i) - xl) ./ sigma).^2);
        ex0 = exp(-( xl ./ sigma ).^2);

        G_global = sum(0.5 * a .* x(i) .* (2 * xg - x(i)));
        G_local = sum(0.5 * b .* (ex - ex0));

        y(i) = G_global + G_local;
    end
end
