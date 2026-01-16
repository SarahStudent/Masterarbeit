% Hauptalgorithmus: ComputeEpsilonOptimizer
function [xbar, I] = ComputeEpsilonOptimizer(ybar, eps, n)
    xbar = ChooseInitialX(ybar, n); 
    I = ComputeInitialApproximation(xbar);
    [N, ~] = ComputeOuterNormals(I);
    W = [];

    while ~AllNormalsCovered(N, W)
        [w_star, y_star, x_star] = ChooseWorstNormalGeneral(N, I, W, n);
        W = [W w_star];

        if IsOutsideEpsilonNeighborhood(y_star, I, eps)
            xbar = x_star;
            I = ComputeUpdatedApproximation(I, y_star);
            [N, ~] = ComputeOuterNormals(I);
        end
    end
end

function xbar = ChooseInitialX(ybar, n)
    cvx_begin quiet
        variable x(n)
        minimize(0)
        subject to
            % zulässiges x
            0 <= x;
            x <= 2;
            % y in F(x)
            norm(ybar - [0;0], 2) <= 1 + 0.3*x;
    cvx_end
    xbar = x;
end

function I = ComputeInitialApproximation(xbar)
    I = [];
    dirs = [eye(2), -ones(2,1)];

    for k = 1:size(dirs,2)
        d = dirs(:,k);
        cvx_begin quiet
            variables t y(2)
            maximize(t)
            subject to
                y == t*d;
                t >= 0;
                % zulässiges x
                0 <= xbar;
                xbar <= 2;
                % y in F(xbar)
                norm(y - [0;0],2) <= 1 + 0.3*xbar;
        cvx_end
        I = append_point_to_hull(I, y);
    end
end

function [w_star, y_star, x_star] = ChooseWorstNormalGeneral(N, I, W, n)
    best_value = -inf; w_star = []; y_star = []; x_star = [];

    for k = 1:size(N,2)
        w = N(:,k);

        if ~isempty(W) % überspringe bereits getestete Richtungen
            if any(vecnorm(W - w, 2, 1) < 1e-8)
                continue;
            end
        end

        cvx_begin quiet
            variables x(n) y(2)
            maximize(w' * y)
            subject to
                % zulässiges x
                0 <= x;
                x <= 2;
                % y in F(x)
                norm(y - [0;0],2) <= 1 + 0.3*x;
                % I Teilmenge von F(x)
                for j = 1:size(I,2)
                    norm(I(:,j) - [0;0],2) <= 1 + 0.3*x;
                end
        cvx_end

        if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')
            if cvx_optval > best_value
                best_value = cvx_optval;
                w_star = w;
                y_star = y;
                x_star = x;
            end
        end
    end

    if isempty(w_star)
        error('Keine zulässige Normale gefunden.');
    end
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
        midpoints(:,i) = (p1 + p2)/2;
    end
end

function flag = IsOutsideEpsilonNeighborhood(y, I, eps)
    k = convhull(I(1,:), I(2,:));
    poly = I(:,k);
    d = point_to_polygon_distance(y, poly);
    flag = (d > eps);
end

function I = ComputeUpdatedApproximation(I, y)
    I = append_point_to_hull(I, y);
end

function flag = AllNormalsCovered(N, W)
    if isempty(W)
        flag = false;
        return;
    end
    flag = true;
    for k = 1:size(N,2)
        if ~any(vecnorm(W - N(:,k), 2, 1) < 1e-8)
            flag = false;
            return;
        end
    end
end

function d = point_to_polygon_distance(y, P)
    if norm(P(:,1)-P(:,end)) > 1e-12
        P = [P P(:,1)];
    end
    m = size(P,2);
    d = inf;
    for k = 1:m-1
        d = min(d, point_segment_distance(y, P(:,k), P(:,k+1)));
    end
end

function d = point_segment_distance(y, a, b)
    ab = b - a;
    t = dot(y-a,ab)/dot(ab,ab);
    t = max(0,min(1,t));
    d = norm(y - (a + t*ab));
end

function I = append_point_to_hull(I, y)
    if isempty(I)
        I = y;
        return;
    end
    Inew = [I y];
    if size(Inew,2) >= 3
        k = convhull(Inew(1,:), Inew(2,:));
        I = Inew(:,k(1:end-1));
    else
        I = Inew;
    end
end


%% Plot

clc; clear; close all;

n = 1;
ybar = [0; 0];
epsilons = [0.1, 0.01];
W = [];

xbar = ChooseInitialX(ybar,n);
I = ComputeInitialApproximation(xbar);
[N, midpoints] = ComputeOuterNormals(I); % initiales I
I_new = I;
Ystar = zeros(n, size(N,2));

[w_star, y_star, x_star] = ChooseWorstNormalGeneral(N, I, W, n);
W = [W, w_star];
I_new1 = ComputeUpdatedApproximation(I_new, y_star);
[N, midpoints] = ComputeOuterNormals(I_new1); 

[w_star, y_star1, x_star1] = ChooseWorstNormalGeneral(N, I, W, n);
W = [W, w_star];
I_new2 = ComputeUpdatedApproximation(I_new1, y_star1);
[N, midpoints] = ComputeOuterNormals(I_new2);

[w_star, y_star2, x_star2] = ChooseWorstNormalGeneral(N, I, W, n);
W = [W, w_star];
I_new3 = ComputeUpdatedApproximation(I_new2, y_star2);
[N, midpoints] = ComputeOuterNormals(I_new3);

[w_star, y_star3, x_star3] = ChooseWorstNormalGeneral(N, I, W, n);
W = [W, w_star];

% erster Plot: Ausgangslage und erste drei iterationen des Algorithmus
figure(1);

subplot(1,4,1); hold on; axis equal; grid on; axis([-2.5 2.5 -2.5 2.5]);

theta = linspace(0, 2 * pi, 200); % F(x0)
r0 = 1 + 0.3 * xbar;
cx0 = 0 + r0 * cos(theta);
cy0 = 0 + r0 * sin(theta);
plot(cx0, cy0, 'k-', 'LineWidth',1.5);

% Initiale Approximation I
I_plot = [I, I(:,1)];
fill(I_plot(1,:), I_plot(2,:), 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'b', 'LineWidth', 1.5);
plot(y_star(1), y_star(2), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); % erster Verbesserungsschritt

subplot(1,4,2); hold on; axis equal; grid on; axis([-2.5 2.5 -2.5 2.5]);
I_plot = [I_new1, I_new1(:,1)];
fill(I_plot(1,:), I_plot(2,:), 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'b', 'LineWidth', 1.5);
r0 = 1 + 0.3 * x_star; % F(xbar)
cx0 = 0 + r0 * cos(theta);
cy0 = 0 + r0 * sin(theta);
plot(cx0, cy0, 'k-', 'LineWidth',1.5);
plot(y_star1(1), y_star1(2), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); % zweiter

subplot(1,4,3); hold on; axis equal; grid on; axis([-2.5 2.5 -2.5 2.5]);
I_plot = [I_new2, I_new2(:,1)];
fill(I_plot(1,:), I_plot(2,:), 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'b', 'LineWidth', 1.5);
r0 = 1 + 0.3 * x_star1;
cx0 = 0 + r0 * cos(theta);
cy0 = 0 + r0 * sin(theta);
plot(cx0, cy0, 'k-', 'LineWidth',1.5);
plot(y_star2(1), y_star2(2), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); % dritter

subplot(1,4,4); hold on; axis equal; grid on; axis([-2.5 2.5 -2.5 2.5]);
I_plot = [I_new3, I_new3(:,1)];
fill(I_plot(1,:), I_plot(2,:), 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'b', 'LineWidth', 1.5);
r0 = 1 + 0.3 * x_star2;
cx0 = 0 + r0 * cos(theta);
cy0 = 0 + r0 * sin(theta);
plot(cx0, cy0, 'k-', 'LineWidth',1.5);
plot(y_star3(1), y_star3(2), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r'); % vierter

% zweiter Plot: Endapproximationen für zwei verschiedene epsilon
figure(2)
for i = 1:2
    epsilon = epsilons(i);
    [x, I_eps] = ComputeEpsilonOptimizer(ybar, epsilon,n);

    num_vertices = size(I_eps, 2); % Anzahl der Ecken
    fprintf("epsilon = %.3g → vertices = %d\n", epsilon, num_vertices);

    subplot(1,2,i); hold on; axis equal; grid on; axis([-2.5 2.5 -2.5 2.5]);
    r0 = 1 + 0.3 * x;
    cx0 = 0 + r0 * cos(theta);
    cy0 = 0 + r0 * sin(theta);
    plot(cx0, cy0, 'k-', 'LineWidth',1.5);
    I_plot = [I_eps, I_eps(:,1)];
    fill(I_plot(1,:), I_plot(2,:), 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'b', 'LineWidth', 1.5);
end