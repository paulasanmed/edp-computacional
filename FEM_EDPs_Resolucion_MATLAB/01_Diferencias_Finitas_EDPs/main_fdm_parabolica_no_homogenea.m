function main_fdm_parabolica_no_homogenea()

% =========================================================================
% SCRIPT: main_fdm_parabolica_no_homogenea
% DESCRIPCIÓN: 
%   Resolución numérica de la ecuación del calor 1D no homogénea.
%   Incorpora un término fuente volumétrico f(x,t) y condiciones de frontera 
%   de Dirichlet dependientes del tiempo (borde activo).
%
% OBJETIVOS DEL ANÁLISIS:
%   - Validación del levantamiento de fronteras no homogéneas.
%   - Evaluación de la saturación del error global por dominancia espacial.
%   - Representación del flujo de calor asimétrico y disipación convectiva.
% =========================================================================

clc; clear; close all; 

% Parámetros físicos y de malla (Problema No Homogéneo)
L = 1.0;
T_final = 0.5;
h_val = [0.1, 0.05, 0.025];
k_val = [0.001, 0.0005, 0.00025];

% Estructura de almacenamiento para perfilado de eficiencia
resultados = struct('h', {}, 'k', {}, 'method', {}, 'error', {}, 'time', {});
count = 1;

fprintf('\n[INFO] Inicializando solver para EDP parabólica con término fuente y frontera dinámica...\n');
fprintf('%-15s %-10s %-10s %-15s %-10s\n', 'Metodo', 'h', 'k', 'Error Max', 'Tiempo(s)');

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
            [u_ex, x_ex, ~] = metodo_explicito_nh(L, T_final, h, k);
            t_cpu = toc;
            
            err = calcular_error_nh(u_ex, x_ex, T_final);
            imprimir_datos('Explicito', h, k, err, t_cpu);
            resultados(count) = struct('h', h, 'k', k, 'method', 'Explicito', 'error', err, 'time', t_cpu);
            count = count + 1;
        else
            fprintf('%-15s %-10.4f %-10.5f %-15s\n', 'Explicito', h, k, 'INESTABLE (r>0.5)');
        end
        
        % 2. Método IMPLÍCITO PURO
        tic;
        [u_im, x_im, ~] = metodo_implicito_nh(L, T_final, h, k);
        t_cpu = toc;
        err = calcular_error_nh(u_im, x_im, T_final);
        imprimir_datos('Implicito', h, k, err, t_cpu);
        resultados(count) = struct('h', h, 'k', k, 'method', 'Implicito', 'error', err, 'time', t_cpu);
        count = count + 1;
        
        % 3. Método CRANK-NICOLSON
        tic;
        [u_cn, x_cn, ~] = metodo_cn_nh(L, T_final, h, k);
        t_cpu = toc;
        err = calcular_error_nh(u_cn, x_cn, T_final);
        imprimir_datos('Crank-Nic', h, k, err, t_cpu);
        resultados(count) = struct('h', h, 'k', k, 'method', 'C-Nicolson', 'error', err, 'time', t_cpu);
        count = count + 1;
    end
    fprintf('--------------------------------------------------------------\n');
end

% =========================================================================
% 2. GENERACIÓN DE GRÁFICAS Y ANÁLISIS
% =========================================================================
% A) Gráficas de Eficiencia (Log-Log)
plot_resultados(resultados, h_val);

fprintf('\n[INFO] Evaluando Orden Empírico de Convergencia (EOC) y saturación del error espacial...\n');
calcular_y_mostrar_oc(resultados);

fprintf('\n[INFO] Renderizando topología térmica asimétrica (Evolución 3D)...\n');
fprintf('\n[INFO] Análisis computacional finalizado con éxito.\n');

% Validación Puntual del Perfil Espacial Espaciotemporal (T = T_final)
% La solución discreta nodal (marcadores rojos) converge puntualmente hacia 
% la solución analítica exacta (curva continua azul), validando la consistencia 
% del esquema numérico y el levantamiento de la condición de frontera.
figure('Name', 'Validación 2D (Final)', 'Color', 'w');
[u_final, x_final, ~] = metodo_cn_nh(L, T_final, 0.02, 0.001);
plot(x_final, u_final, 'ro', 'MarkerSize', 4, 'DisplayName', 'Numérica (C-N)');
hold on;

u_ex = (T_final/(T_final+1)) * cos(pi*x_final/2).^2; 
plot(x_final, u_ex, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Exacta');
title('Comparación en T_{final} (2D)'); legend('Location', 'best');
grid on; 

% Evolución Espaciotemporal de la Superficie Térmica (Dinámica Asimétrica)
% Dinámica térmica asimétrica. La topología refleja el equilibrio dinámico 
% entre la inyección de energía por la frontera de Dirichlet activa, la 
% generación volumétrica impuesta por la fuente y la disipación convectiva.
[~, x, t, U_hist] = metodo_cn_nh(L, T_final, 0.05, 0.005);
[T_mesh, X_mesh] = meshgrid(t, x);
figure('Name', 'Evolución 3D Completa', 'Color', 'w');
surf(T_mesh, X_mesh, U_hist);
shading interp;
colormap jet;
colorbar;
xlabel('Tiempo (t)'); ylabel('Espacio (x)'); zlabel('u(x,t)');
title('Evolución de la Solución No Homogénea');
view(-45, 30); 
end

% =========================================================================
%   DEFINICIÓN DE FUNCIONES DEL PROBLEMA
% =========================================================================

% --- TÉRMINO FUENTE f(x,t) ---
% Generación volumétrica de calor
function val = evaluar_termino_fuente_volumetrico(x, t)
term1 = (cos(pi*x/2).^2) ./ ((t+1)^2);
term2 = (pi^2 * t .* cos(pi*x)) ./ (2*(t+1));
val = term1 + term2;
end

% --- CONDICIÓN FRONTERA IZQUIERDA u(0,t) ---
% Imposición de condición de Dirichlet no homogénea variable en el tiempo
function val = evaluar_frontera_dirichlet_activa(t)
val = t / (t + 1);
end

% =========================================================================
%   ESQUEMAS NUMÉRICOS (Adaptados para f(x,t) y contornos variables)
% =========================================================================

function [u, x, t] = metodo_explicito_nh(L, T, h, k)
x = 0:h:L;
t = 0:k:T;
N = length(x);
M = length(t);
r = k / h^2;

u = zeros(N, 1); 
u_new = zeros(N, 1);

for n = 1:M-1
    t_act = t(n);
    
    % Evaluación del operador de difusión en el subespacio interior
    difusion = r*u(1:N-2) + (1-2*r)*u(2:N-1) + r*u(3:N);
    
    % Inyección del término fuente evaluado temporalmente
    fuente = k * evaluar_termino_fuente_volumetrico(x(2:N-1)', t_act);
    u_new(2:N-1) = difusion + fuente; 
    
    % Imposición fuerte de condiciones de frontera
    u_new(1) = evaluar_frontera_dirichlet_activa(t(n+1)); % Dirichlet no homogénea
    u_new(N) = 0;                 % Dirichlet homogénea
    u = u_new;
end
end

function [u, x, t] = metodo_implicito_nh(L, T, h, k)
x = 0:h:L;
t = 0:k:T;
N = length(x);
M = length(t);
r = k / h^2;
n_int = N - 2;

e = ones(n_int, 1);
A = spdiags([-r*e, (1+2*r)*e, -r*e], -1:1, n_int, n_int);

u = zeros(N, 1);
u_int = u(2:N-1); % Reducción al subespacio de nodos interiores

for n = 1:M-1
    t_fut = t(n+1); 
    
    % Ensamblaje del vector de cargas: estado previo + fuente nodal
    b = u_int + k * evaluar_termino_fuente_volumetrico(x(2:N-1)', t_fut);
    
    % Levantamiento de frontera (boundary lifting)
    b(1) = b(1) + r * evaluar_frontera_dirichlet_activa(t_fut);
    
    u_int = A \ b; % Resolución del sistema
end
u = [evaluar_frontera_dirichlet_activa(T); u_int; 0]; 
end

function [u, x, t, U_hist] = metodo_cn_nh(L, T, h, k)
x = 0:h:L;
t = 0:k:T;
N = length(x);
M = length(t);
r = k / h^2;
n_int = N - 2;

e = ones(n_int, 1);
A = spdiags([-r/2*e, (1+r)*e, -r/2*e], -1:1, n_int, n_int);
B = spdiags([ r/2*e, (1-r)*e,  r/2*e], -1:1, n_int, n_int);

u = zeros(N, 1);
u_int = u(2:N-1); 

if nargout > 3
    U_hist = zeros(N, M);
    U_hist(:, 1) = u;
end

for n = 1:M-1
    t_act = t(n);
    t_futu = t(n+1);
    
    % Integración de la fuente mediante promedio temporal (Crank-Nicolson)
    f_promedio = 0.5 * (evaluar_termino_fuente_volumetrico(x(2:N-1)', t_act) + evaluar_termino_fuente_volumetrico(x(2:N-1)', t_futu));
    b = B * u_int + k * f_promedio;
    
    % Levantamiento de frontera promediado
    front_prom = 0.5 * (evaluar_frontera_dirichlet_activa(t_act) + evaluar_frontera_dirichlet_activa(t_futu));
    b(1) = b(1) + r * front_prom;
    
    u_int = A \ b; 
    
    if nargout > 3
        U_hist(:, n+1) = [evaluar_frontera_dirichlet_activa(t_futu); u_int; 0];
    end
end
u = [evaluar_frontera_dirichlet_activa(T); u_int; 0]; 
if ~(nargout > 3)
    U_hist = [];
end 
end

% =========================================================================
%   FUNCIONES AUXILIARES DE CÁLCULO Y PLOT
% =========================================================================

function err = calcular_error_nh(u_num, x, t_final)
u_exacta = (t_final / (t_final + 1)) * (cos(pi * x' / 2).^2);
err = max(abs(u_num - u_exacta)); 
end

function imprimir_datos(metodo, h, k, err, t)
fprintf('%-15s %-10.4f %-10.5f %-15.2e %-10.4f\n', metodo, h, k, err, t);
end

function plot_resultados(results, h_val, k_val_all)
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
        
        eoc = diff(log(errs)) ./ diff(log(ks)); 
        
        if strcmp(metodo, 'C-Nicolson')
            teorico = "~ 2.0";  
        else
            teorico = "~ 1.0";  

        end
        fprintf('   -> %-12s: EOC promedio = %.2f (Teórico: %s)\n', metodo, mean(eoc), teorico);
    end

    % =========================================================================
% [NOTA ANALÍTICA DE CONVERGENCIA]
% El orden empírico resultante (~0.32) refleja la saturación del error 
% global por la componente de discretización espacial O(h^2). Al mantener 
% 'h' constante, el principio de superposición del error impide visualizar 
% la convergencia asintótica puramente temporal O(k^2) característica del 
% esquema incondicionalmente estable de Crank-Nicolson.
% =========================================================================

end
end