function main_fem_galerkin_1d()

% =========================================================================
% SCRIPT: main_fem_galerkin_1d
% DESCRIPCIÓN: 
%   Implementación computacional del Método de Elementos Finitos (FEM) 1D 
%   basado en la formulación variacional de Galerkin para un problema de 
%   contorno de tipo reacción-difusión.
%
% CARACTERÍSTICAS TÉCNICAS:
%   - Ensamblaje topológico con matrices dispersas (sparse).
%   - Espacios de interpolación: Lineal (P1) y Cuadrático (P2).
%   - Cuadratura de Gauss-Legendre exacta para evitar 'quadrature pollution'.
%   - Cálculo empírico del Orden de Convergencia (EOC) en norma L2.
% =========================================================================

clear; clc; close all;

% Secuencia de refinamiento de la partición espacial (h = 1/N).
% Rango acotado para evidenciar la convergencia asintótica sin 
% incurrir en saturación por precisión de máquina o exceso de coste computacional.
N = [4, 8, 16, 32, 64, 128];

% Almacenamiento de normas de error global
errores_P1 = zeros(length(N), 1); % Interpolación lineal P1
errores_P2 = zeros(length(N), 1); % Interpolación cuadrática P2
h_val   = zeros(length(N), 1);    % Parámetro de malla espacial

fprintf('============================================================\n');
fprintf('\n[INFO] Inicializando motor de resolución FEM 1D (Galerkin)...\n');
fprintf('[INFO] Problema de contorno: d2u/dx2 + u = f(x) con Dirichlet homogéneas.\n');
fprintf('[INFO] Ejecutando secuencia paramétrica de refinamiento espacial (h -> 0)...\n');
fprintf('============================================================\n');

% =========================================================================
%  BUCLE DE REFINAMIENTO DE LA MALLA
% =========================================================================
for i = 1:length(N)
    Ni = N(i);
    h_val(i) = 1/Ni;
    
    % --- Evaluación del subespacio de interpolación polinómica P1 --- 
    [u_h1, x_nodos1] = resolver_sistema_galerkin_P1(Ni); 
    errores_P1(i) = calcular_error_L2(u_h1, x_nodos1, 1);
    
    % --- Evaluación del subespacio de interpolación polinómica cuadrática P2 ---
    [u_h2, x_nodos2] = resolver_sistema_galerkin_P2(Ni);
    errores_P2(i) = calcular_error_L2(u_h2, x_nodos2, 2);
end

% Cálculo del Orden Empírico de Convergencia (EOC):
% Asumiendo Error ~ C*h^p, aplicando logaritmos: log(Error) ~ log(C) + p*log(h).
% El operador 'diff' extrae el orden asintótico p mediante diferencias progresivas.
% Se incluye un padding de ceros para preservar la dimensionalidad de los vectores.
ord_P1 = abs([0; diff(log(errores_P1)) ./ diff(log(h_val))]);
ord_P2 = abs([0; diff(log(errores_P2)) ./ diff(log(h_val))]);

% Visualización de resultados empíricos
fprintf('\n[INFO] Evaluando convergencia asintótica del error global en norma L2:\n');
fprintf('------------------------------------------------------------\n');
fprintf('|  N  | Error L2 (P1) | Orden P1 | Error L2 (P2) | Orden P2 |\n');
fprintf('------------------------------------------------------------\n');
for i = 1:length(N)
    fprintf('| %3d |   %1.4e  |   %1.2f   |   %1.4e  |   %1.2f   |\n', ...
        N(i), errores_P1(i), ord_P1(i), errores_P2(i), ord_P2(i));
end
fprintf('------------------------------------------------------------\n');

% =========================================================================
%  GRÁFICA DE CONVERGENCIA (Norma L2)
% =========================================================================
figure('Name', 'Convergencia L2 FEM', 'Color', 'w');

loglog(h_val, errores_P1, '-o', 'LineWidth', 2, 'DisplayName', 'Lineal (P1)');
hold on;
loglog(h_val, errores_P2, '-s', 'LineWidth', 2, 'DisplayName', 'Cuadrático (P2)');

% Curvas asintóticas de referencia para verificación teórica del error
loglog(h_val, h_val.^2 * (errores_P1(1)/h_val(1)^2), '--k', 'DisplayName', 'Pendiente 2 (O(h^2))');
loglog(h_val, h_val.^3 * (errores_P2(1)/h_val(1)^3), '-.k', 'DisplayName', 'Pendiente 3 (O(h^3))');

grid on;
xlabel('Tamaño de elemento h');
ylabel('Error L2 ||u - u_h||');
title('Convergencia Asintótica del Método de Elementos Finitos');
legend('Location', 'best');

% Inversión del eje de abscisas para representar el refinamiento 
% progresivo de la malla (h -> 0) de izquierda a derecha.
set(gca, 'XDir', 'reverse'); 
end

% =========================================================================
%   CÁLCULO DE SOLUCIONES Y ENSAMBLAJE
% =========================================================================

function [u, x] = resolver_sistema_galerkin_P1(Ni)
h = 1/Ni;
x = linspace(0, 1, Ni+1)';
n_nodos = Ni + 1;

% Ensamblaje mediante matrices dispersas (sparse) para optimizar 
% la complejidad espacial O(N) de los operadores globales.
K = spalloc(n_nodos, n_nodos, 3*n_nodos); 
M = spalloc(n_nodos, n_nodos, 3*n_nodos);
F = zeros(n_nodos, 1);

% Nodos y pesos de Gauss-Legendre en el dominio de referencia [-1,1].
% Para n=3, se integran exactamente polinomios de grado <= 2n-1 = 5.
gp = [-sqrt(3/5), 0, sqrt(3/5)]; 
gw = [5/9, 8/9, 5/9]; 

for i = 1:Ni
    idx = [i, i+1]; % Conectividad del elemento
    
    % Matrices elementales de Rigidez (difusión) y Masa (reacción)
    k_loc = (1/h) * [1 -1; -1 1]; 
    m_loc = (h/6) * [2 1; 1 2]; 
    
    K(idx, idx) = K(idx, idx) + k_loc; 
    M(idx, idx) = M(idx, idx) + m_loc; 
    
    f_loc = zeros(2, 1);
    x_a = x(i);
    x_b = x(i+1);
    
    for q = 1:length(gp)
        % Mapeo afín del dominio físico al dominio de referencia mediante el Jacobiano
        jac = (x_b - x_a)/2; 
        xi = gp(q); 
        x_r = jac * xi + (x_b + x_a)/2;
        
        % Evaluación de funciones de forma lineales nodales
        phi = [(1-xi)/2; (1+xi)/2];
        val_f = evaluar_termino_fuente(x_r); 
        
        % La cuadratura de Gauss-Legendre evalúa la integral de forma exacta para polinomios 
        % de grado <= 2n-1. Mapeamos el dominio físico [x_a, x_b] al dominio de referencia 
        % [-1, 1] mediante el Jacobiano (jac). 
        % Al usar 3 puntos de integración, garantizamos la evaluación exacta del sistema 
        % para elementos P1 sin contaminar el orden de convergencia espacial.
        f_loc = f_loc + val_f * phi * jac * gw(q);
    end
    F(idx) = F(idx) + f_loc; 
end

A = K + M;

% Imposición fuerte de condiciones de frontera de Dirichlet homogéneas
nodos_fijos = [1, n_nodos];
% Reducción al subespacio de grados de libertad (nodos interiores)
nodos_libres = setdiff(1:n_nodos, nodos_fijos);

u = zeros(n_nodos, 1);
% Resolución del sistema algebraico reducido
u(nodos_libres) = A(nodos_libres, nodos_libres) \ F(nodos_libres);
end

function [u, x] = resolver_sistema_galerkin_P2(Ni)
h = 1/Ni;
% Topología del elemento P2: Requiere 3 nodos por elemento 
% (extremos y punto medio), resultando en 2N+1 grados de libertad globales.
x = linspace(0, 1, 2*Ni+1)'; 
n_nodos = 2*Ni + 1;

% Ancho de banda ampliado debido al soporte de las funciones de forma cuadráticas
K = spalloc(n_nodos, n_nodos, 5*n_nodos);
M = spalloc(n_nodos, n_nodos, 5*n_nodos);
F = zeros(n_nodos, 1);

gp = [-sqrt(3/5), 0, sqrt(3/5)];
gw = [5/9, 8/9, 5/9];

for i = 1:Ni
    idx = [2*i-1, 2*i, 2*i+1];
    
    k_loc = (1/(3*h)) * [ 7, -8,  1; -8, 16, -8; 1, -8,  7]; 
    m_loc = (h/30) * [ 4,  2, -1; 2, 16,  2; -1,  2,  4];
    
    f_loc = zeros(3,1);
    x_a = x(idx(1));
    x_b = x(idx(3)); 
    
    for q = 1:length(gp)
        xi = gp(q);
        jac = (x_b - x_a)/2;
        x_r = (jac)*xi + (x_a + x_b)/2;
        
        % Funciones de forma cuadráticas en el elemento de referencia
        phi = [xi*(xi-1)/2;  (1-xi^2);  xi*(xi+1)/2];
        
        val_f = evaluar_termino_fuente(x_r);
        f_loc = f_loc + val_f * phi * jac * gw(q);
    end
    K(idx, idx) = K(idx, idx) + k_loc;
    M(idx, idx) = M(idx, idx) + m_loc;
    F(idx) = F(idx) + f_loc;
end

A = K + M; 

nodos_fijos = [1, n_nodos];
nodos_libres = setdiff(1:n_nodos, nodos_fijos);

u = zeros(n_nodos, 1);
u(nodos_libres) = A(nodos_libres, nodos_libres) \ F(nodos_libres);
end

% =========================================================================
%   FUNCIONES AUXILIARES
% =========================================================================

function y = evaluar_termino_fuente(x)
y = (1 + pi^2) * sin(pi*x);
end

function y = evaluar_solucion_analitica(x)
y = sin(pi*x);
end

function error_L2 = calcular_error_L2(u_num, x_nodes, type)
% Cuantificación del error en norma L2 mediante integración de alta precisión
suma_error = 0;

% Condición analítica necesaria para el ensamblaje de elementos de 
% interpolación cuadrática (P2), requiriendo un nodo central por elemento 
% (número total de nodos impar).
if type == 2 && mod(length(x_nodes), 2) == 0
    error('Error topológico: El método cuadrático (P2) exige un número impar de nodos globales.');
end

% Inyección directa de los nodos y pesos de Gauss-Legendre (n=5) para garantizar 
% un error de cuadratura casi nulo en la estimación del error.
gp = [-0.9061798459, -0.5384693101, 0, 0.5384693101, 0.9061798459];
gw = [0.2369268850, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268850];

if type == 1 % Elementos P1
    N_elem = length(x_nodes) - 1; 
    for i = 1:N_elem
        idx = [i, i+1];
        u_local = u_num(idx); 
        x_a = x_nodes(i); 
        x_b = x_nodes(i+1);
        jac = (x_b - x_a)/2; 
        
        for q = 1:length(gp)
            xi = gp(q);
            x_r = jac*xi + (x_a + x_b)/2;
            phi = [(1-xi)/2; (1+xi)/2];
            
            % Evaluación del interpolante numérico u_h mediante combinación lineal 
            % de las funciones de forma nodales (phi) y los grados de libertad.
            u_h_val = dot(u_local, phi);
            u_ex_val = evaluar_solucion_analitica(x_r);
            
            suma_error = suma_error + (u_ex_val - u_h_val)^2 * jac * gw(q);
        end
    end
    
elseif type == 2 % Elementos P2
    % Recorrido de la topología de la malla para elementos cuadráticos (salto de 2 nodos)
    N_elem = (length(x_nodes)-1)/2;
    for i = 1:N_elem
        idx = [2*i-1, 2*i, 2*i+1]; 
        
        u_local = u_num(idx);
        x_a = x_nodes(idx(1)); 
        x_b = x_nodes(idx(3));
        jac = (x_b - x_a)/2;
        
        for q = 1:length(gp)
            xi = gp(q);
            x_r = jac*xi + (x_a + x_b)/2;
           
            phi = [xi*(xi-1)/2; (1-xi^2); xi*(xi+1)/2];
            u_h_val = dot(u_local, phi);
            u_ex_val = evaluar_solucion_analitica(x_r);
            
            suma_error = suma_error + (u_ex_val - u_h_val)^2 * jac * gw(q);
        end
    end
end

error_L2 = sqrt(suma_error);
end