function main_fdm_parabolica_homogenea()

% =========================================================================
% SCRIPT: main_fdm_parabolica_homogenea
% DESCRIPCIÓN: 
%   Resolución numérica y perfilado de eficiencia computacional de la 
%   ecuación del calor 1D homogénea mediante esquemas de Diferencias Finitas.
%   Evalúa convergencia asintótica y estabilidad de Courant-Friedrichs-Lewy (CFL).
%
% ESQUEMAS IMPLEMENTADOS:
%   - Explícito (Condicionalmente estable)
%   - Implícito Puro (Incondicionalmente estable, O(h^2 + k))
%   - Crank-Nicolson (Incondicionalmente estable, O(h^2 + k^2))
% =========================================================================

clc; clear; close all; %Inicialización del entorno computacional y liberación de memoria

% Parámetros físicos y de malla
T_final = 0.5; % Límite temporal
L = 1.0; % Dominio espacial [0, L]
h_val = [0.1, 0.05, 0.025, 0.0125];

% Valores de k ajustados analíticamente para cumplir la restricción r <= 0.5
% y garantizar la estabilidad condicional del método explícito
k_val = [0.001, 0.0005, 0.00025];  

% Estructura de almacenamiento para perfilado de eficiencia
resultados = struct('h', {}, 'k', {}, 'method', {}, 'error', {}, 'time', {});
count = 1; % Inicialización del contador de iteraciones

% Cabecera de la tabla de resultados empíricos
fprintf('%-15s %-10s %-10s %-15s %-10s\n', 'Metodo', 'h', 'k', 'Error Max', 'Tiempo(s)');
fprintf('--------------------------------------------------------------\n');

% =========================================================================
% 1. BUCLE DE CÁLCULO Y ANÁLISIS DE EFICIENCIA
% =========================================================================
% Bucle paramétrico para evaluación de convergencia asintótica y 
% eficiencia computacional de los esquemas numéricos
for h = h_val
    for k = k_val
        
        % 1. Método EXPLÍCITO
        r = k/h^2; % Coeficiente de estabilidad (CFL)
        if r <= 0.5
            tic; % Inicialización de medición de tiempo CPU
            [u_ex, x_ex, ~] = esquema_explicito(L, T_final, h, k);
            t_cpu = toc;
            
            err = calcular_error(u_ex, x_ex, T_final);
            imprimir_datos('Explicito', h, k, err, t_cpu);
            resultados(count) = struct('h', h, 'k', k, 'method', 'Explicito', 'error', err, 'time', t_cpu);
            count = count + 1;
        else 
            % Restricción de estabilidad violada
            fprintf('%-15s %-10.4f %-10.5f %-15s\n', 'Explicito', h, k, 'INESTABLE (r>0.5)');
        end
        
        % 2. Método IMPLÍCITO PURO
        tic;
        [u_im, x_im, ~] = esquema_implicito(L, T_final, h, k);
        t_cpu = toc;
        err = calcular_error(u_im, x_im, T_final);
        imprimir_datos('Implicito', h, k, err, t_cpu);
        resultados(count) = struct('h', h, 'k', k, 'method', 'Implicito', 'error', err, 'time', t_cpu);
        count = count + 1;
        
        % 3. Método CRANK-NICOLSON
        tic;
        [u_cn, x_cn, ~] = esquema_crank_nicolson(L, T_final, h, k);
        t_cpu = toc;
        err = calcular_error(u_cn, x_cn, T_final);
        imprimir_datos('Crank-Nic', h, k, err, t_cpu);
        resultados(count) = struct('h', h, 'k', k, 'method', 'C-Nicolson', 'error', err, 'time', t_cpu);
        count = count + 1;
    end
    fprintf('--------------------------------------------------------------\n');
end

% --- GRÁFICAS Y ANÁLISIS ---
graficar_eficiencia_asintotica(resultados, h_val);

fprintf('\n[INFO] Ejecutando análisis empírico del Orden de Convergencia (EOC)...\n');
calcular_y_mostrar_oc(resultados);

fprintf('\n[INFO] Evaluando topología de esparcidad matricial...\n');
visualizar_matrices(L, 0.1, 0.001);

fprintf('\n[INFO] Renderizando la evolución espaciotemporal 3D de la dinámica térmica...\n');
graficar_evolucion_3d(L, T_final, 0.05, 0.005);

fprintf('\n[INFO] Forzando violación del criterio de estabilidad CFL (r > 0.5) para evidenciar divergencia numérica...\n');
demonstrar_inestabilidad(L, 0.2, 0.05);

fprintf('\n[INFO] Análisis computacional finalizado con éxito.\n');
end

% =========================================================================
%   ESQUEMAS NUMÉRICOS
% =========================================================================

function [u, x, t, U_hist] = esquema_explicito(L, T, h, k)
x = 0:h:L;
t = 0:k:T; % Discretización espaciotemporal

N = length(x);
M = length(t);
r = k / h^2; % Parámetro de red
u = sin(pi * x'); % Condición inicial analítica

if nargout > 3
    U_hist = zeros(N, M);
    U_hist(:,1) = u; 
end

u_nueva = zeros(N, 1);

% Integración temporal mediante aproximación de diferencias finitas centrales
for n = 1:M-1
    % Vectorización del operador para optimización de CPU
    u_nueva(2:N-1) = r*u(1:N-2) + (1 - 2*r)*u(2:N-1) + r*u(3:N);
    
    % Imposición fuerte de Condiciones de Frontera de Dirichlet Homogéneas
    u_nueva(1) = 0;
    u_nueva(N) = 0;
    
    u = u_nueva; 
    if nargout > 3, U_hist(:,n+1) = u; end 
end
end

function [u, x, t, A] = esquema_implicito(L, T, h, k)
x = 0:h:L;
N = length(x);
M = round(T/k);
r = k / h^2;

n_int = N - 2;

% Reducimos el sistema lineal al subespacio de los N-2 nodos interiores, 
% imponiendo las condiciones de frontera de Dirichlet de forma fuerte.
e = ones(n_int, 1); 

% Ensamblaje de la matriz de coeficientes (Topología Tridiagonal)
A = spdiags([-r*e, (1+2*r)*e, -r*e], -1:1, n_int, n_int);

u = sin(pi * x');
u_int = u(2:N-1); % Reducción al subespacio de nodos interiores

for n = 1:M  
    % Resolución del sistema algebraico
    u_int = A \ u_int;
end
u = [0; u_int; 0]; % Levantamiento de condiciones de frontera
t = 0:k:T;
end

function [u, x, t, U_history] = esquema_crank_nicolson(L, T, h, k)
x = 0:h:L;
t = 0:k:T;
N = length(x);
M = length(t);
r = k / h^2;
n_int = N - 2;
e = ones(n_int, 1);

% Ensamblaje de matrices del esquema incondicionalmente estable
A = spdiags([-r/2*e, (1+r)*e, -r/2*e], -1:1, n_int, n_int);
B = spdiags([ r/2*e, (1-r)*e,  r/2*e], -1:1, n_int, n_int);

u = sin(pi * x');
u_inter = u(2:N-1); % Reducción al subespacio de nodos interiores

if nargout > 3
    U_history = zeros(N, M);
    U_history(:,1) = u;
end

for n = 1:M-1 
    rhs = B * u_inter; 
    u_inter = A \ rhs; % Resolución del sistema
    if nargout > 3 
        U_history(:,n+1) = [0; u_inter; 0];
    end
end
u = [0; u_inter; 0];
end

% =========================================================================
%   FUNCIONES DE GRÁFICOS Y ANÁLISIS
% =========================================================================

function calcular_y_mostrar_oc(results)
h_ref = 0.025; 
metodos = {'Implicito', 'C-Nicolson'};
fprintf('   Analizando convergencia temporal para h fijo = %.4f\n', h_ref);

for m = metodos
    metodo = m{1};
    idx = strcmp({results.method}, metodo) & abs([results.h] - h_ref) < 1e-9;
    res_sub = results(idx);
    
    if length(res_sub) > 1
        errs = [res_sub.error];
        ks = [res_sub.k];
        
        [ks, sort_idx] = sort(ks, 'descend');
        errs = errs(sort_idx);
        
        % Computación empírica del orden de convergencia asintótico
        eoc = diff(log(errs)) ./ diff(log(ks)); 
        
        if strcmp(metodo, 'C-Nicolson')
            teorico = "~ 2.0";  
        else
            teorico = "~ 1.0";  
        end
        fprintf('   -> %-12s: EOC promedio = %.2f (Teórico: %s)\n', metodo, mean(eoc), teorico);
    end
end
end

function graficar_eficiencia_asintotica(results, h_val, k_val_all)
if nargin < 3
    k_val_all = unique([results.k]);
end

color = {'r', 'b', 'g'};
trazo = {'o', 's', '^', 'd', 'x', '*'};
metodos = {'Explicito', 'Implicito', 'C-Nicolson'};

% --- GRÁFICA 1: FIJANDO H  ---
figure('Name', 'Eficiencia - h Fijo (Varía k)', 'Color', 'w');
hold on; legend_str = {};
for i = 1:length(metodos)
    metodo_actual = metodos{i};
    col = color{i};
    for j = 1:length(h_val) 
        h_actual = h_val(j); 
        mark = trazo{min(j, length(trazo))};
        idx = strcmp({results.method}, metodo_actual) & abs([results.h] - h_actual) < 1e-9;
        res_m = results(idx);
        
        if ~isempty(res_m)
            [times, sort_idx] = sort([res_m.time]);
            errs = [res_m.error];
            errs = errs(sort_idx);
            loglog(times, errs, [col '-' mark], 'LineWidth', 1.5, ...
                'MarkerFaceColor', col, 'MarkerSize', 6);
            legend_str{end+1} = sprintf('%s (h=%.3f)', metodo_actual, h_actual);
        end
    end
end
grid on;
xlabel('Tiempo CPU (s)'); ylabel('Error Global');
title('Eficiencia con h fijo (Variando k)');
legend(legend_str, 'Location', 'bestoutside'); axis tight;

% --- GRÁFICA 2: FIJANDO K --- 
figure('Name', 'Eficiencia - k Fijo (Varía h)', 'Color', 'w');
hold on; legend_str2 = {};
k_val_all = sort(k_val_all, 'descend');
for i = 1:length(metodos)
    metodo_actual = metodos{i}; col = color{i};
    for j = 1:length(k_val_all)
        k_val = k_val_all(j);
        mark = trazo{min(j, length(trazo))};
        idx = strcmp({results.method}, metodo_actual) & abs([results.k] - k_val) < 1e-9;
        res_m = results(idx);
        if ~isempty(res_m)
            [times, sort_idx] = sort([res_m.time]);
            errs = [res_m.error]; errs = errs(sort_idx);
            loglog(times, errs, [col '--' mark], 'LineWidth', 1.5, ...
                'MarkerFaceColor', col, 'MarkerSize', 6);
            legend_str2{end+1} = sprintf('%s (k=%.4f)', metodo_actual, k_val);
        end
    end
end
grid on;
xlabel('Tiempo CPU (s)'); ylabel('Error Global');
title('Eficiencia con k fijo (Variando h)');
legend(legend_str2, 'Location', 'bestoutside'); axis tight;

end

function demonstrar_inestabilidad(L, T, h)
% Violación de la condición de estabilidad CFL (r > 0.5): Generación de oscilaciones 
% espurias de alta frecuencia por amplificación del error de truncamiento.
k_inestable = 0.60 * h^2; 
r = k_inestable / h^2;
fprintf('   -> Ejecutando Explícito con r=%.2f (Inestable > 0.5)\n', r);
[~, x, ~, U_inest] = esquema_explicito(L, T, h, k_inestable);

figure('Name', 'Inestabilidad Numérica', 'Color', 'w');
plot(x, U_inest(:, end), 'r-o', 'LineWidth', 1.5); hold on;
u_exacta = exp(-pi^2 * T) * sin(pi * x);
plot(x, u_exacta, 'b--', 'LineWidth', 2);
legend('Numérica (Inestable)', 'Exacta');
title(sprintf('Inestabilidad del Método Explícito (r = %.2f)', r));
grid on;
end

function visualizar_matrices(L, h, k)
% Inspección de la topología de la matriz del sistema (esparcidad)
[~, ~, ~, A] = esquema_implicito(L, 0.1, h, k);
figure('Name', 'Estructura Matriz', 'Color', 'w');
spy(A);
title('Estructura de la Matriz A (Implícito)');
end

function graficar_evolucion_3d(L, T, h, k)
[~, x, t, U_total] = esquema_crank_nicolson(L, T, h, k);
[T_mesh, X_mesh] = meshgrid(t, x);
figure('Name', 'Evolución 3D', 'Color', 'w');
surf(T_mesh, X_mesh, U_total);
shading interp; 
colormap jet; 
colorbar; 
xlabel('t'); ylabel('x'); zlabel('u'); title('Evolución de la Temperatura');
end

function imprimir_datos(metodo, h, k, err, t)
fprintf('%-15s %-10.4f %-10.5f %-15.2e %-10.4f\n', metodo, h, k, err, t);
end

function err = calcular_error(u_num, x, t_final)
u_exacta = exp(-pi^2 * t_final) * sin(pi * x');
err = max(abs(u_num - u_exacta));
end