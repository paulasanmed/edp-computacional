# Resolución Computacional de EDPs: Diferencias y Elementos Finitos

Repositorio técnico que contiene la implementación en MATLAB y la fundamentación analítica para la resolución de Ecuaciones en Derivadas Parciales (EDPs) parabólicas y elípticas 1D.

## 📌 Arquitectura del Proyecto

### 📁 `01_Diferencias_Finitas_EDPs`
Implementación optimizada de discretización clásica y análisis de estabilidad:
* **Scripts:** Esquemas Explícito, Implícito y Crank-Nicolson para ecuaciones parabólicas homogéneas y no homogéneas.
* **Validación Analítica (`Analisis_Estabilidad_FDM.pdf`):** Perfilado de eficiencia computacional, flujo asimétrico e inestabilidades bajo el criterio CFL.

### 📁 `02_Elementos_Finitos_1D`
Formulación variacional pesada y resolución algebraica:
* **Script:** Solvers de Galerkin para elementos de Lagrange $P_1$ y $P_2$. Ensamblaje topológico mediante matrices dispersas (*sparse*) y cuadratura exacta de Gauss-Legendre.
* **Validación Analítica (`Validacion_Convergencia_FEM.pdf`):** Cálculo empírico del Orden de Convergencia (EOC) en norma $L^2$.

### 📄 Documentación Maestra (Directorio Raíz)
* **`Manual_Teorico_EDPs_FDM_FEM.pdf`:** Monografía analítica integral. Abarca desde la estabilidad espectral de Von Neumann hasta formulaciones variacionales fundamentadas en el Teorema de Lax-Milgram. 

## ⚙️ Requisitos de Ejecución
* MATLAB R2021a o superior.
* No requiere *toolboxes* adicionales (implementación y ensamblaje matricial construidos desde cero).

## 🚀 Uso
Ejecutar los scripts principales directamente desde la consola de MATLAB:
```matlab
>> main_fdm_parabolica_homogenea
>> main_fem_galerkin_1d
