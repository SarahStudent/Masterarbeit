clc; close all;

%% Daten
r = [0.095; 0.105]; % Rendite
Q = [0.030 -0.012; -0.012 0.020]; % Kovarianz

alpha = 0.75; % Skalierungsfaktor der Unsicherheit

V = [ 0.010  -0.005; % Unsicherheitsrichtungen
      0.020   0.015;
     -0.006   0.012;
      0.000  -0.008 ];

M = 0.4;
D = [ 0   0; % -R^2_+ addieren
     -M   0;   
      0   M ]; 

% Rechteckgrenzen weil wir nicht den ganzen -R^2_+ brauchen
x_min = 0.089; x_max = 0.11;
y_min = 0;    y_max = 0.036;

% Startportfolio
x0 = [0;1];             
eps = 1e-6; % kleine Toleranz
n = length(x0);

%% Plot Figure 1

% Pareto-Front
lambda = linspace(0,1,200); % Diskretisierung
mu = lambda*r(1) + (1-lambda)*r(2); % erwartete Rendite
risk_plot = lambda.^2*Q(1,1) + 2*lambda.*(1-lambda)*Q(1,2) + (1-lambda).^2*Q(2,2); %Varianz

ax(1) = subplot(1,3,1);
plot(mu, risk_plot, 'k', 'LineWidth', 2); % Pareto-Front
hold on;
xlabel('Erwartete Rendite');
ylabel('Risiko (Varianz)');
grid on;

% Dominierte Flächen
mid_idx = 87;
idx_right = mid_idx:length(mu);
idx_left = 1:mid_idx;
X_fill = [mu(idx_left), fliplr(mu(idx_left))];
Y_fill = [risk_plot(idx_left), y_max*ones(size(idx_left))];
fill(X_fill, Y_fill, [0.85 0.85 0.85], 'EdgeColor','none', 'FaceAlpha',0.3);

[~, idx_min_risk] = min(risk_plot); % minimaler Varianzwert
mu_min_risk = mu(idx_min_risk);
risk_min_risk = risk_plot(idx_min_risk);

X_hor = [x_min, mu_min_risk, mu_min_risk, x_min];
Y_hor = [risk_min_risk, risk_min_risk, y_max, y_max];

fill(X_hor, Y_hor, [0.85 0.85 0.85], 'EdgeColor','none', 'FaceAlpha',0.3);

% Pareto_Front erneut zeichnen
plot(mu, risk_plot, 'k', 'LineWidth', 2);

% Ausgewählte Portfolios plotten (f(x))
lambdas_x = 0:0.2:1;
nSteps = length(lambdas_x);

hold on

for k = 1:length(lambdas_x)
    lam = lambdas_x(k);
    mu_x   = lam*r(1) + (1-lam)*r(2);
    risk_x = lam^2*Q(1,1) + 2*lam*(1-lam)*Q(1,2) + (1-lam)^2*Q(2,2);
    color_k = interp1(linspace(0,1,nSteps), parula(nSteps), lam);
    plot(mu_x, risk_x, 'o','MarkerFaceColor', color_k,'MarkerEdgeColor', color_k);
end

hold off

% Achsen und Labels
xlabel('Erwartete Rendite');
ylabel('Risiko (Varianz)');
set(gca,'YDir','reverse'); % invertierte Risikoachse

xlim([x_min + 0.001, x_max - 0.001]);
ylim([y_min + 0.001, y_max - 0.002]);

%% Plot figure 2
ax(2) = subplot(1,3,2);

plot(mu, risk_plot, 'k', 'LineWidth', 2); % Pareto-Front
hold on;
xlabel('Erwartete Rendite');
ylabel('Risiko (Varianz)');
grid on;
fill(X_hor, Y_hor, [0.85 0.85 0.85], 'EdgeColor','none', 'FaceAlpha',0.3); % dominierte PF
fill(X_fill, Y_fill, [0.85 0.85 0.85], 'EdgeColor','none', 'FaceAlpha',0.3);
plot(mu, risk_plot, 'k', 'LineWidth', 2); % Pareto-Front

% Ausgewählte Portfolios plotten (conv(f(x) + delta V)
lambdas_x = 0:0.2:1;
nSteps = length(lambdas_x);

hold on

for k = 1:length(lambdas_x)
    lam = lambdas_x(k);
    mu_x   = lam*r(1) + (1-lam)*r(2);
    risk_x = lam^2*Q(1,1) + 2*lam*(1-lam)*Q(1,2) + (1-lam)^2*Q(2,2);
    f_x = [mu_x; risk_x];

    % delta(x)
    delta = alpha * lam * (1-lam);

    Fx = f_x' + delta * V;
    Fx = unique(Fx,'rows');
    color_k = interp1(linspace(0,1,nSteps), parula(nSteps), lam);

    if size(Fx,1) >= 3
        % 2D: konvexe Hülle
        K = convhull(Fx(:,1), Fx(:,2));
        fill(Fx(K,1), Fx(K,2), color_k,'FaceAlpha', 0.3,'EdgeColor', color_k);
    elseif size(Fx,1) == 2
        % Liniensegment:
        plot(Fx(:,1), Fx(:,2), 'b-', 'LineWidth',2);
    else
        % Einzelpunkt:
        plot(Fx(1,1), Fx(1,2), 'yo', 'MarkerFaceColor',color_k, 'MarkerEdgeColor', color_k);
    end
end

hold off

% Achsen und Labels
xlabel('Erwartete Rendite');
ylabel('Risiko (Varianz)');
set(gca,'YDir','reverse'); % invertierte Risikoachse

xlim([x_min + 0.001, x_max - 0.001]);
ylim([y_min + 0.001, y_max - 0.002]);

%% Plot figure 3
ax(3) = subplot(1,3,3);

plot(mu, risk_plot, 'k', 'LineWidth', 2); % Pareto-Front
hold on;
xlabel('Erwartete Rendite');
ylabel('Risiko (Varianz)');
grid on;

% Ausgewählte Portfolios plotten (F(x))
lambdas_x = 0:0.2:1;
nSteps = length(lambdas_x);

hold on

for k = 1:nSteps
    lam = lambdas_x(k);

    mu_x   = lam*r(1) + (1-lam)*r(2);
    risk_x = lam^2*Q(1,1) + 2*lam*(1-lam)*Q(1,2) + (1-lam)^2*Q(2,2);

    x0 = [lam; (1-lam)];
    Fx = ComputeF(x0, r, Q, V, D, alpha, x_min, x_max, y_min, y_max)';
    
    color_k = interp1(linspace(0,1,nSteps), parula(nSteps), lam);

    if size(Fx,1) >= 3
        K = convhull(Fx(:,1), Fx(:,2));
        fill(Fx(K,1), Fx(K,2), color_k,'FaceAlpha', 0.3,'EdgeColor', color_k);
    end
end

hold off

% Achsen und Labels
xlabel('Erwartete Rendite');
ylabel('Risiko (Varianz)');
set(gca,'YDir','reverse'); % invertierte Risikoachse

xlim([x_min + 0.001, x_max - 0.001]);
ylim([y_min + 0.001, y_max - 0.002]);

% Colorbar einsetzen
nSteps = length(lambdas_x);
cmap = flipud(parula(nSteps));

colormap(cmap)
caxis(ax, [0 1])

for i = 1:3
    pos = ax(i).Position;
    pos(3) = pos(3) * 0.92;
    ax(i).Position = pos;
end

cb = colorbar('Position',[0.93 0.11 0.02 0.78]);
cb.Label.String = '$1-x_1$';
cb.Label.Interpreter = 'latex';
cb.Ticks = 0:0.2:1;

%% Funktionen

function Fx = ComputeF(x0, r, Q, V, D, alpha, x_min, x_max, y_min, y_max)
   lam = x0(1);
   mu_x0   = lam*r(1) + (1-lam)*r(2); % erwartete Rendite
   risk_x0 = lam^2*Q(1,1) + 2*lam*(1-lam)*Q(1,2) + (1-lam)^2*Q(2,2); % Varianz
   f_x0 = [mu_x0; risk_x0]; % f(x)
   delta_x0 = alpha * lam * (1-lam); % delta(x)

   Fx_base = f_x0' + delta_x0 * V;
   Fx = [];
   for i = 1:size(Fx_base,1)
       Fx = [Fx; Fx_base(i,:) + D];
   end
   Fx = unique(Fx,'rows');

   Fx(:,1) = min(max(Fx(:,1), x_min), x_max); % 'abschneiden' von F(x), weil wir in der Implementierung
   Fx(:,2) = min(max(Fx(:,2), y_min), y_max); % keinen unendlichen Kegel -R^2_+ benötigen, Andeutung reicht
   Fx = [Fx; [x_min,y_max]];
   Fx = Fx';
end
