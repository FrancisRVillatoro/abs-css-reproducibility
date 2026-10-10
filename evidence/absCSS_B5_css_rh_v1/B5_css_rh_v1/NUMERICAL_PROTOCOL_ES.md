# Protocolo numérico reducido para la ecuación absoluta CSS

Versión revisada tras auditoría independiente: 2026-09-04.

## 1. Objetivo científico y límites de interpretación

La campaña no pretende repetir un catálogo general de colisiones. Debe responder, en este orden, a cinco preguntas falsables:

1. ¿Un compactón y el dipolo exacto se propagan de forma convergente y conservan las cantidades correctas?
2. ¿Las formulaciones conservativas en `u` y en la variable de presión `v` convergen a la misma evolución?
3. ¿Dos regularizaciones hamiltonianas suaves distintas seleccionan el mismo límite cuando `epsilon -> 0`?
4. ¿Ese límite coincide con el obtenido mediante hiperviscosidad disipativa cuando `eta -> 0`?
5. ¿En qué difieren dinámicamente `|K*|(p+2,p)` y una ecuación absoluta Rosenau-Hyman escalada para poseer exactamente el mismo compactón inicial?

Ningún resultado numérico debe presentarse como prueba de buena posición, unicidad, existencia de frontera libre o estabilidad matemática. Las leyes de escala del cruce transversal son predicciones **de tiempo inicial** para un perfil impuesto; no deben extrapolarse sin prueba a tiempos `O(1)`.

## 2. Modelos que deben implementarse

### M1. Forma conservativa CSS en `u`

```text
u_t + Dx(|u|^(l-2) u/(l-1))
    - alpha p Dx(|u|^(p-2) u u_x^2)
    + 2 alpha Dx^2(|u|^p u_x) = 0.
```

Ésta es la formulación primaria. No debe evaluarse la forma expandida en los ceros.

Diagnósticos suaves esperados:

```text
C_CSS = integral u dx,
P_CSS = (1/2) integral u^2 dx,
H_CSS = integral [alpha |u|^p u_x^2-|u|^(p+2)/((p+1)(p+2))] dx.
```

### M2. Forma de presión sobre `l=p+2`

```text
v = sign(u) |u|^(p/2),
u = sign(v) |v|^(2/p),
```

```text
v_t + v^2 v_x
    + 2 alpha [v^2 v_xxx
      + (1+6/p) v v_x v_xx
      + (4/p^2) v_x^3] = 0.
```

La ecuación es polinómica, pero sigue siendo degenerada: el coeficiente de `v_xxx` es `2 alpha v^2`.

### M3 y M4. Dos regularizaciones hamiltonianas

Para `r=1,2`, usar

```text
A_{epsilon,r}(u) = (u^(2r)+epsilon^(2r))^(p/(2r)),
G_{epsilon,r}(u) = [(u^(2r)+epsilon^(2r))^((p+2)/(2r))-epsilon^(p+2)]/[(p+1)(p+2)],
H_{epsilon,r} = integral [alpha A_{epsilon,r}(u) u_x^2-G_{epsilon,r}(u)] dx,
u_t = Dx(delta H_{epsilon,r}/delta u).
```

Los casos `r=1` y `r=2` suavizan el mismo límite mediante kernels locales distintos. `epsilon` no debe confundirse con la disipación artificial.

### M5. Ruta disipativa por hiperviscosidad

Usar la forma conservativa M1 con

```text
u_t + N_CSS[u] + eta u_xxxx = 0,   eta>0,
```

donde `N_CSS[u]` denota todos los términos no disipativos del miembro izquierdo de M1. Con esta convención, el término `+ eta u_xxxx` produce amortiguamiento de Fourier `exp(-eta k^4 t)`.

Esta ruta no es hamiltoniana. Debe extrapolarse `eta -> 0` de forma independiente y compararse con los límites `epsilon -> 0` de M3 y M4. Cuando se use la maquinaria de 2023-2024, conservar exactamente su convención de signo y documentar cualquier filtro adicional.

### M6. Comparador Rosenau-Hyman de igual perfil

Con `q=p+1`, resolver

```text
u_t + a_RH Dx(|u|^(q-1)u) + b_RH Dx^3(|u|^(q-1)u) = 0,
a_RH = 4/(p+2)^2,
b_RH = 4 alpha/[(p+1)(p+2)].
```

M1 y M6 tienen el mismo compactón exacto, con la misma amplitud, anchura, velocidad y signo. La comparación aísla la dinámica, no la geometría del dato inicial.

Para RH no debe tratarse la variación de `(1/2) integral u^2` como error numérico. Los diagnósticos correctos son

```text
H_RH = - integral |u|^(q+1)/(q+1) dx,
C0_RH = integral u dx,
Cc_RH = integral u cos(kappa_RH x) dx,
Cs_RH = integral u sin(kappa_RH x) dx,
kappa_RH = sqrt(a_RH/b_RH),
```

si el dominio y las condiciones de contorno permiten definir los Casimires armónicos. En un dominio periódico, elegir una longitud compatible con `2 pi/kappa_RH`. Además debe comprobarse la identidad

```text
d/dt [(1/2) integral u^2 dx]
 = (b_RH q(q-1)/2) integral sign(u)|u|^(q-2) u_x^3 dx.
```

## 3. Normalización y conjunto de parámetros

Campaña principal:

```text
c = 1,
alpha = p^2/[(p+1)(p+2)],
l = p+2,
p in {1, 1.25, 1.5, 1.75, 2}.
```

Con esta elección, todos los compactones tienen semisupporte `L=pi`:

```text
A = [(p+1)(p+2)/2]^(1/p),
k = 1/2,
u_c(x,t) = sigma A cos((x-t)/2)^(2/p),  |x-t|<=pi.
```

Los CSV exactos del paquete son referencias independientes del código de evolución.

## 4. Prerrequisito T0: residuales y consistencia algebraica

Antes de integrar en tiempo:

- evaluar el residual discreto de M1 sobre el compactón en puntos interiores;
- comprobar el residual de M2 sobre `v_c = sigma A^(p/2) cos((x-t)/2)`;
- comprobar que M6 acepta exactamente el mismo perfil;
- verificar el orden espacial esperado lejos de los bordes;
- excluir inicialmente una banda de varios stencils alrededor de cada borde y después estudiar por separado cómo decrece esa banda al refinar;
- verificar en M6 la identidad instantánea para `dP/dt` anterior.

T0 es una prueba de implementación, no uno de los cinco experimentos científicos.

## 5. Los cinco benchmarks obligatorios

### B1. Propagación de un compactón

Dato: un compactón positivo aislado.

Medir:

```text
error_shape(t) = min_delta ||u_h(.,t)-u_c(.-delta,t)||_2 / ||u_c||_2,
error_phase(t) = |delta_best(t)|,
drift_C, drift_P, drift_H,
negative_mass = integral max(-u_h,0) dx.
```

Criterios:

- convergencia monótona del error de forma al refinar;
- orden observado compatible con la regularidad del borde, no necesariamente con el orden nominal del stencil;
- deriva de invariantes CSS que disminuya con `Delta x` y `Delta t`;
- ausencia de un ripple cuya amplitud tienda a un valor no nulo al refinar sólo la malla;
- repetición en M1, M2 y, al menos, una de M3-M5.

### B2. Dipolo compacto exacto

Usar dos lóbulos adyacentes, positivo y negativo, con la misma velocidad. Debe cumplirse

```text
C = 0,
P = 2 P_c,
H = 2 H_c.
```

Es el benchmark más importante para cambios de signo porque no contiene una colisión.

Criterios:

- `|C(t)|` converge a cero;
- las derivas relativas de `P` y `H` convergen a cero;
- `||u_M1-sign(v_M2)|v_M2|^(2/p)||_2` converge a cero;
- el cero interno se desplaza con velocidad `c`;
- no aparece una capa de amplitud o anchura fija dependiente de malla.

### B3. Cruce transversal regularizado: prueba de tiempo inicial y relajación observada

#### B3a. Kernel exacto en `t=0`

Imponer localmente `u=a x`, inicialmente con `a=1`, y evaluar `D_x B_{p,epsilon,r}` sin integrar la PDE. Para `r=1,2`, la referencia exacta es

```text
D_x B_{p,epsilon,r}
 = p a^3 (a x)^(2r-2) [(a x)^(2r)+epsilon^(2r)]^(p/(2r)-2)
   * [(2r-1)epsilon^(2r)+(p-1)(a x)^(2r)].
```

Predicciones de tiempo inicial:

```text
peak_abs ~ epsilon^(p-2),
transition/core width ~ epsilon/|a|,
integral over |x|<=Z epsilon/|a| ~ epsilon^(p-1),
```

con `Z` fijo. Para `p=1`, la masa total del kernel tiende a `2 a |a|`.

Usar inicialmente

```text
epsilon in {1e-1, 3e-2, 1e-2, 3e-3, 1e-3}.
```

Valores menores sólo se aceptan cuando la capa está resuelta. Objetivo: al menos 24 puntos en la escala `epsilon/|a|`; mínimo absoluto: 16. Por tanto,

```text
Delta x <= epsilon/(24 |a|)   (objetivo),
Delta x <= epsilon/(16 |a|)   (mínimo).
```

Para mallas uniformes globales, no reducir `epsilon` si esto exige un coste desproporcionado: B3a puede ejecutarse en coordenada local reescalada `z=a x/epsilon`.

Ajustar exponentes mediante regresión `log-log`; error absoluto objetivo menor que `0.05` y estabilidad menor que `0.03` al retirar el mayor `epsilon`.

El script `generate_zero_crossing_benchmarks.py` genera los datos y figuras de referencia para B3a.

#### B3b. Evolución del cero preparado

Construir un dato global suave con

```text
u_0(x)=a(x-x0)+O((x-x0)^2)
```

cerca de `x0`, evitando que el cutoff contamine la ventana local. Registrar:

- `u_t(x,0)` y/o el gradiente del flujo en `t=0`;
- posición del cero `x0(t)`;
- pendiente local;
- exponente ajustado `r_loc(t)` en `|u| ~ C |x-x0(t)|^r`;
- la misma cantidad en presión, donde un borde compatible tendría un cero aproximadamente lineal.

No asumir de antemano una escala de relajación. Determinarla empíricamente mediante tiempos de muestreo densos al comienzo y comprobar su convergencia con `Delta x`, `Delta t`, `epsilon` o `eta`.

### B4. Test de selección entre formulaciones y regularizaciones

Sólo después de B1-B3:

1. comparar M1 con M2;
2. comparar M3 (`r=1`) con M4 (`r=2`);
3. comparar el límite hamiltoniano con M5 (`eta -> 0`);
4. repetir para compactón, dipolo y cero preparado.

No variar todos los parámetros a la vez. Secuencia:

1. `Delta t -> 0` con `Delta x`, `epsilon` o `eta` fijos;
2. `Delta x -> 0` con el parámetro de regularización fijo;
3. `eta -> 0` para M5 a resolución convergida;
4. `epsilon -> 0` separadamente para M3 y M4;
5. comparar los tres límites y la formulación M2.

Magnitudes:

```text
D_12(t) = ||u_M1-sign(v_M2)|v_M2|^(2/p)||_2,
D_H(t;epsilon) = ||u_{r=1}-u_{r=2}||_2,
D_HD(t) = ||u_Hamiltonian-limit-u_hyperviscous-limit||_2.
```

Selección universal: las diferencias extrapolan a cero dentro de incertidumbre. Dependencia de regularización: alguna extrapola a un valor no nulo estable bajo refinamiento y cambios razonables de stencil.

### B5. Comparación CSS-RH con el mismo perfil

Primero propagar un lóbulo y el dipolo en M1 y M6. Después realizar una interacción de alcance con la misma separación y velocidades positivas.

Diagnósticos CSS:

```text
C_CSS, P_CSS, H_CSS.
```

Diagnósticos RH:

```text
C0_RH, Cc_RH, Cs_RH, H_RH,
residual_Pdot = numerical_dPdt - exact_rhs_Pdot.
```

Comparar:

- fase y velocidad;
- pérdida de coherencia de cada lóbulo;
- `L2` del residuo fuera de ventanas ligadas a los pulsos;
- signo, velocidad y número de estructuras emitidas;
- evolución física de `P` en RH frente a conservación de `P` en CSS;
- dependencia respecto de `eta`, `Delta x`, `Delta t` y formulación conservativa.

Una diferencia es física sólo si sobrevive a todas las extrapolaciones pertinentes.

## 6. Colisiones opcionales posteriores

Sólo si B1-B5 pasan:

- compactón-compactón de alcance, velocidades `1` y `0.5`;
- compactón-anticompactón de alcance, signos opuestos y velocidades `1` y `0.5`;
- opcionalmente relación de velocidades `1:0.25`.

No empezar con colisiones frontales: para `alpha>0` y `c>0`, las ramas relacionadas por `u -> -u` se propagan en la misma dirección.

## 7. Refinamiento y resolución

Para propagación de compactones sin capa `epsilon`, una parrilla inicial puede ser

```text
N in {1024, 2048, 4096, 8192},
Delta t = C_t Delta x^3,
```

con al menos dos valores de `C_t`.

Para regularizaciones, la resolución se decide por la escala local, no por una lista fija de `N`. Registrar siempre

```text
points_in_layer = epsilon/(|a| Delta x).
```

No incluir en un ajuste puntos con menos de 16 celdas por capa. Para M5, documentar la longitud disipativa asociada a `eta` y verificar que está resuelta.

## 8. Métricas e incertidumbres

Para cada observable `Q`, usar extrapolaciones escalonadas, no un ajuste simultáneo ciego:

```text
Q_h = Q_* + C_x Delta x^r + C_t Delta t^s
                + C_eta eta^gamma + C_eps epsilon^beta.
```

Incertidumbre mínima:

- mitad de la diferencia entre las dos extrapolaciones de mayor resolución;
- sensibilidad al subconjunto de puntos del ajuste;
- sensibilidad al alineamiento de fase;
- sensibilidad a la ventana usada para medir el cero o el ripple.

Para `r_loc(t)`, informar la ventana de ajuste y repetirla con al menos tres tamaños. Para un ripple, publicar la ventana espacial exacta, normas `L1`, `L2`, `Linf` y espectro de Fourier filtrado y sin filtrar.

## 9. Reglas de parada

Detener la campaña y no pasar a colisiones si ocurre cualquiera de estas situaciones:

- el compactón o dipolo no converge bajo refinamiento;
- M1 y M2 no se aproximan en B1-B2;
- los exponentes de B3a dependen de puntos no resueltos;
- una diferencia entre regularizaciones desaparece al cambiar sólo `Delta x`;
- una supuesta estructura física cambia de signo, velocidad o número al variar el filtro;
- en RH, se interpreta erróneamente la variación física de `P` como deriva numérica.

## 10. Productos de datos

Cada corrida debe guardar:

```text
config.json
invariants.csv
observables.csv
snapshots.npz o HDF5
convergence_summary.json
stdout.log
```

B3a añade:

```text
initial_time_scaling.csv
fitted_exponents.json
zero_crossing_peak_scaling.png
zero_crossing_width_scaling.png
zero_crossing_normalized_kernels.png
```

El artículo debe citar el commit, el entorno, las tolerancias, los scripts de postproceso y los datos de extrapolación.
