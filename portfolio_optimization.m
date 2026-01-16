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

%% Plot 1 Ausgangsportfolio
% Pareto-Front
lambda = linspace(0,1,200); % Diskretisierung
mu = lambda*r(1) + (1-lambda)*r(2); % erwartete Rendite
risk_plot = lambda.^2*Q(1,1) + 2*lambda.*(1-lambda)*Q(1,2) + (1-lambda).^2*Q(2,2); % Varianz

subplot(1,2,1)
plot(mu, risk_plot, 'k', 'LineWidth', 2); % Pareto-Front
hold on;
xlabel('Erwartete Rendite');
ylabel('Risiko (Varianz)');
grid on;

% Ausgangsportfolio plotten
x_star = [0;1];
I = ComputeF(x_star, r, Q, V, D, alpha, x_min, x_max, y_min, y_max);
K = convhull(I(1,:), I(2,:));
fill(I(1,K), I(2,K), 'y', 'FaceAlpha',0.2, 'EdgeColor','y', 'LineWidth',1.5);
mu_star = x_star(1) * r(1) + x_star(2) * r(2); 
risk_star = x_star' * Q * x_star;
point_star = [mu_star; risk_star];
plot(point_star(1), point_star(2), 'yo','MarkerFaceColor','y');

% zugehörige Normalenvektoren
q = quiver(0.1, 0.02, 0, -0.0025, 0);
q.LineWidth = 2;
q.Color = 'k';
q.MaxHeadSize = 5;

q = quiver(0.105, 0.025, 0.0025, 0, 0);
q.LineWidth = 2;
q.Color = 'r';
q.MaxHeadSize = 5;

% Achsen und Labels
xlabel('Erwartete Rendite');
ylabel('Risiko (Varianz)');
set(gca,'YDir','reverse'); % invertierte Risikoachse

xlim([x_min + 0.001, x_max - 0.001]);
ylim([y_min + 0.001, y_max - 0.002]);

%% Plot 2 (Optimierer)

subplot(1,2,2)
plot(mu, risk_plot, 'k', 'LineWidth', 2); % Pareto-Front
hold on;
xlabel('Erwartete Rendite');
ylabel('Risiko (Varianz)');
grid on;

% Ausgangsportfolio
K = convhull(I(1,:), I(2,:));
fill(I(1,K), I(2,K), 'y', 'FaceAlpha',0.2, 'EdgeColor','y', 'LineWidth',1.5);
mu_star = x_star(1) * r(1) + x_star(2) * r(2); 
risk_star = x_star' * Q * x_star;
point_star = [mu_star; risk_star];
plot(point_star(1), point_star(2), 'yo','MarkerFaceColor','y');

% Bestimmung nächstes Fx (eins besser)
N = ComputeOuterNormals(I);
N = unique(N.', 'rows').';
w1 = [-1; 0]; % hier in diesem konkreten Beispiel nicht von Interesse
w2 = [ 0; 1]; % hier in diesem konkreten Beispiel nicht von Interesse
W = [w1 w2];

[w_star, y_star, x_star] = ChooseWorstNormalGeneral(N, W, I, r, Q, V, D, alpha, y_min, y_max, x_min, x_max);
W = [W w_star];
I = ComputeF(x_star, r, Q, V, D, alpha, x_min, x_max, y_min, y_max);
[N, midpoints] = ComputeOuterNormals(I);

% Außennormalen plotten, die uns interessieren
for i = 3:6
    q = quiver(midpoints(1, i), midpoints(2, i), N(1, i)*0.0025, N(2, i)*0.0025, 0, 'b');
    q.LineWidth = 2;
    q.Color = 'k';
    q.MaxHeadSize = 5;
end

% Plot des F(xbar)
K = convhull(I(1,:), I(2,:));
fill(I(1,K), I(2,K), 'b', 'FaceAlpha',0.2, 'EdgeColor','b', 'LineWidth',1.5);
mu_star = x_star(1) * r(1) + x_star(2) * r(2); 
risk_star = x_star' * Q * x_star;
point_star = [mu_star; risk_star];
plot(point_star(1), point_star(2), 'bo','MarkerFaceColor','b');

% Achsen und Labels
xlabel('Erwartete Rendite');
ylabel('Risiko (Varianz)');
set(gca,'YDir','reverse'); % invertierte Risikoachse

xlim([x_min + 0.001, x_max - 0.001]);
ylim([y_min + 0.001, y_max - 0.002]);

%% nochmaliges Durchlaufen zur Bestätigung des Optimierers

N = unique(N.', 'rows').';
[w_star, y_star, x_star] = ChooseWorstNormalGeneral(N, W, I, r, Q, V, D, alpha, y_min, y_max, x_min, x_max);

%% Funktionen

function Fx = ComputeF(x0, r, Q, V, D, alpha, x_min, x_max, y_min, y_max)
   lam = x0(1);
   mu_x0   = lam*r(1) + (1-lam)*r(2);
   risk_x0 = lam^2*Q(1,1) + 2*lam*(1-lam)*Q(1,2) + (1-lam)^2*Q(2,2);
   f_x0 = [mu_x0; risk_x0];
   delta_x0 = alpha * lam * (1-lam);

   Fx_base = f_x0' + delta_x0 * V;
   Fx = [];
   for i = 1:size(Fx_base,1)
       Fx = [Fx; Fx_base(i,:) + D];
   end
   Fx = unique(Fx,'rows');

   Fx(:,1) = min(max(Fx(:,1), x_min), x_max);
   Fx(:,2) = min(max(Fx(:,2), y_min), y_max);
   Fx = [Fx; [x_min,y_max]];
   Fx = Fx';
end

function [N, midpoints] = ComputeOuterNormals(I)
   kHull = convhull(I(1,:), I(2,:));
   kHull = kHull(1:end-1);
   nEdges = length(kHull);
   N = zeros(2, nEdges);
   midpoints = zeros(2, nEdges);
   centroid = mean(I,2);

   for i = 1:nEdges
       p1 = I(:,kHull(i));
       p2 = I(:,kHull(mod(i,nEdges)+1));
    
       e = p2 - p1;
       normal = [0 -1; 1 0] * e;
       normal = normal / norm(normal);
    
       if dot(normal, centroid - p1) > 0
           normal = -normal;
       end
       N(:,i) = normal;
       midpoints(:,i) = (p1+p2)/2;
   end
end

function [w_star, y_star, x_star] = ChooseWorstNormalGeneral(N, W, Fxbar, r, Q, V, D, alpha, y_min, y_max, x_min, x_max)
    best_value = -inf; w_star = []; y_star = []; x_star = [];
    nV = size(V,1);
    nD = size(D,1);
    
    for k = 1:size(N,2)
        w = N(:,k);

        if ~isempty(W)
            if any(vecnorm(W - w, 2, 1) < 1e-8)
                continue;
            end
        end

        fprintf('========================\nEvaluating normal w = [%g, %g]\n', w(1), w(2));
        for i = 1:nV
            for j = 1:nD
                % Funktion zur Maximierung entlang der Normalen
                fun = @(lam) w' * ([lam*r(1) + (1-lam)*r(2);
                                  -(lam^2*Q(1,1) + 2*lam*(1-lam)*Q(1,2) + (1-lam)^2*Q(2,2))] ...
                                  + alpha*lam*(1-lam)*V(i,:)' + D(j,:)');

                % i-Schleife: probiere jede Extremstörung V(i,:)
                % j-Schleife: probiere jede Basisverschiebung D(j,:)
                % fun(lam) = den Punkt entlang der λ-Linie, verschoben durch V(i,:) + D(j,:), projiziert auf die Richtung w
                
                % Finde lam, das maximal entlang w wirkt
                [lam_star, val_neg] = fminbnd(@(l) -fun(l), 0, 1); % positiv w
                val = -val_neg;

                % Berechne Fx-Punkte für diese lam_star
                mu_x = lam_star*r(1) + (1-lam_star)*r(2);
                risk_x = lam_star^2*Q(1,1) + 2*lam_star*(1-lam_star)*Q(1,2) + (1-lam_star)^2*Q(2,2);
                delta_x = alpha*lam_star*(1-lam_star);

                Fx = zeros(2, nV*nD);
                idx = 1;
                for vi = 1:nV
                    for dj = 1:nD
                        Fx_point = [mu_x; -risk_x] + delta_x*V(vi,:)' + D(dj,:)';
                        Fx_point(1) = min(max(Fx_point(1), x_min), x_max);
                        Fx_point(2) = min(max(Fx_point(2), y_min), y_max);
                        Fx(:,idx) = Fx_point;
                        idx = idx + 1;
                    end
                end
  
                % Prüfe, ob Fxbar ⊆ Fx (toolbox-frei, approximativ)
                Fx_min = min(Fx,[],2);
                Fx_max = max(Fx,[],2);
                is_in_hull = all(all(Fxbar >= Fx_min & Fxbar <= Fx_max));

                % Live-Tracking Ausgabe
                fprintf('i=%d, j=%d, lam_star=%.4f, val=%.4f, Fxbar_in_Fx=%d, current_best=%.4f\n', ...
                    i, j, lam_star, val, is_in_hull, best_value);

                % Speichere, falls Fxbar enthalten und Verbesserung besser
                if is_in_hull && val > best_value
                    best_value = val;
                    x_star = [lam_star; 1-lam_star];
                    y_star = Fx(:,1);
                    w_star = w;
                end
            end
        end
    end

    fprintf('========================\nAlgorithm finished.\n');
    if ~isempty(w_star)
        fprintf('Best normal w_star = [%g, %g]\n', w_star(1), w_star(2));
        fprintf('Best y_star = [%g, %g]\n', y_star(1), y_star(2));
        fprintf('Best x_star = [%g, %g]\n', x_star(1), x_star(2));
    else
        fprintf('No further strict improvement possible.');
    end
end