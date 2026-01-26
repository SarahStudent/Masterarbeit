Berechnung eines $\epsilon$-Optimierers kompakt-konvexer Mengenoptimierungsprobleme
# Masterarbeit
Ergänzendes und zusätzliches Material zur Masterarbeit "Ein Algorithmus für konvexe Mengenoptimierungsprobleme". Der Algorithmus basiert auf dem Algorithmus 1 aus Kapitel 3 des Artikels "A solution method for arbitrary polyhedral convex set optimization problems" von Andreas Löhne (https://arxiv.org/abs/2310.06602).

## Voraussetzungen
Damit die Codedateien ausgeführt werden können, muss MATLAB, CVX und GLPK installiert sein. Die Dateien sind mit MATLAB Version 25.2 (R2025b) erstellt worden.

## Hauptziel und Thema der Masterarbeit und diesem Repository
In der Masterarbeit wird eine **Lösungsmethode für kompakt-konvexe Mengenoptimierungsprobleme** vorgestellt. Das entwickelte Verfahren liefert einen $\epsilon$-Optimierer, der die Zielfunktion bzgl. der Mengeninklusion $\subseteq$ bis auf eine vorgegebene Toleranz maximiert. In Zuge dessen sind ein **Algorithmus** zur Konstruktion eines $\epsilon$-Optimierers und eine **beispielhafte Implementierung**, die für zwei ausgewählte Beispiele einen $\epsilon$-Optimierer berechnet, entwickelt worden. Die zentrale Funktionsweise des Algorithmus besteht darin iterativ die Facetten einer polyedrischen Innenapproximation möglichst weit nach außen zu erschieben und dadurch diejenige Menge $F(x)$ identifizieren, die die Funktion bezüglich der Mengeninklusion maximiert. Für eine ausführliche Einführung in das Thema wird auf die Masterarbeit verwiesen, die sich als PDF in dem Repository befindet.

## Codedateien (.m)
In allen Codedateien sind die Variablen wie folgt zu interpretieren:
- $x$ ist die Eingabevariable
- $F$ ist die konvexe, mengenwertige Funktion
- $y$ ist ein Funktionswert: $y\in F(x)$ 
- $eps$ ist die Toleranz $\epsilon$
- $n$ ist die Dimension von $x$
- $I$ ist eine Innenapproximation von $F(x)$
- $N$ ist die Menge der Außennormalen von $I$
- $W$ ist die Menge der bereits überprüften Außennormalen
- $w$ ist die Außennormale, die gerade getestet wird

### running_example.m
In diesem Beispiel wird eine kompakt-konvexe, mengenwertige Funktion maximiert, die in Abhängigkeit von $x$ eine Kreisscheibe um den Mittelpunkt erzeugt: $F(x)=$ \{ $y \in \mathbb{R}^2 \vert \Vert y \Vert \leq 1 + 0.3x$ \}  mit $x\in [0,2] \cap \mathbb{R}$. $\epsilon$-Optimierer ist $x=2$.


Funktionen:
- function [xbar, I] = ComputeEpsilonOptimizer(ybar, eps, n): Berechnet den $\epsilon$-Optimierer $xbar$ und die Endinnenapproximation, sodass $ybar \in F(xbar)$ und $I \subseteq  F(xbar) \subseteq I + \epsilon B$.
- function xbar = ChooseInitialX(ybar, n): Bestimmt ein beliebiges $xbar$, sodass $ybar \in F(xbar)$.
- function I = ComputeInitialApproximation(xbar): Berechnet eine initiale Innenapproximation $I$ von $F(xbar)$ (siehe Algorithmus 3.3 von Daniel Dörfler in https://www.db-thueringen.de/servlets/MCRFileNodeServlet/dbt_derivate_00060748/Dissertation_DanielDoerfler.pdf).
- function [N, midpoints] = ComputeOuterNormals(I): Berechnet die Außennormalen der polyedrischen Approximation $I$.
- function [w_star, y_star, x_star] = ChooseWorstNormalGeneral(N, I, W, n): Löst ein konvexes Optimierungsproblem zur Bestimmung einer Richtung, die potentiell ein Verbesserung erlaubt und noch nicht getestet wurde.
- function flag = IsOutsideEpsilonNeighborhood(y, I, eps): Prüft, ob die Verschiebung zum Punkt $y$ außerhalb von $I+\epsilon B$ liegt.
- function I = ComputeUpdatedApproximation(I, y): Aktualisiert die Innenapproximation $I$ (siehe Algorithmus 3.3 von Daniel Dörfler in https://www.db-thueringen.de/servlets/MCRFileNodeServlet/dbt_derivate_00060748/Dissertation_DanielDoerfler.pdf).
- function flag = AllNormalsCovered(N, W): Überprüft, ob alle Außennormalen getestet wurden.


Hilfsfunktionen:
- function d = point_to_polygon_distance(y, P): Bestimmt die Distanz von einem Punkt $y$ zu einem Polyeder $P$.
- function d = point_segment_distance(y, a, b): Bestimmt die Distanz von einem Punkt $y$ zu einer Strecke zwischen den Punkten $a$ und $b$.
- function I = append_point_to_hull(I, y): Fügt den Punkt $y$ der konvexen Hülle von $I$ hinzu.

### construct_Fx.m und portfolio_optimization.m
In diesem Beispiel wird die Portfoliooptimierung behandelt. Die herkömmliche multikriterielle Portfoliooptimierung (maximiere Rendite mit gleichzeitiger Minimierung des Riskikos) wird mengenwertig erweitert zu einem zweistufigen Modell, indem eine Unsicherheit hinzugefügt wird: $F(x) = conv((r^Tx,-x^TQx)^T + \delta(x)V)-\mathbb{R}^2_+$ mit $\sum_j x_j=1, x \geq 0$. Ausgangsportfolio ist $(0,1)^T$ und $\epsilon$-Optimierer ist $(\frac{1}{3}, \frac{2}{3})^T$.

Eingabegrößen
- $x0$ ist das Ausgangsportfolio
- $r$ ist die erwartete Rendite
- $Q$ ist die Kovarianzmatrix
- $V$ sind die Unsicherheitsrichtungen
- $\delta(x)$ ist die Skalierungsfunktion der Unsicherheit
- $alpha$ ist ein zusätzlicher Skalierungsfaktor in $\delta(x)$
- $D, x_{min}, x_{max}, y_{min}, y_{max}$ simulieren die Addition von $-\mathbb{R}^2_+$


Funktionen:
- function Fx = ComputeF(x0, r, Q, V, D, alpha, x_min, x_max, y_min, y_max): Berechnet $F(x)$.
- function [N, midpoints] = ComputeOuterNormals(I): Berechnet die Außennormalen der polyedrischen Menge $I$.
- function [w_star, y_star, x_star] = ChooseWorstNormalGeneral(N, W, Fxbar, r, Q, V, D, alpha, y_min, y_max, x_min, x_max): Bestimmt eine Richtung, die potentiell ein Verbesserung erlaubt und noch nicht getestet wurde.

## Bilddateien (.jpg)

Alle Bilddateien in diesem Repository werden durch die Dateien running_example.m, construct_Fx.m und portfolio_optimization.m erzeugt, wobei Erstere die Abbildungen 1 und 2 für das Beispiel 5.1 aus der Masterarbeit kreiert. Die anderen beiden Dateien erstellen jeweils Abbildung 3 und 4 des Portfoliooptimierungsbeispiels 5.2. Sie dienen der Veranschaulichung der Ergebnisse, Funktionsweise und Funktionalität des entwickelten Algorithmus. Alle Abbildungen sind in Kapitel 5 der Masterarbeit auffindbar.
