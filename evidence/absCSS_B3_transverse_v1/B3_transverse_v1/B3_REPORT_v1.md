# B3 — Cruce transversal regularizado

## Veredicto

**B3a: PASS.**  **B3b: PASS como benchmark diagnóstico de límite singular; no se observa relajación espontánea a un borde de compactón en la ventana local limpia.**

## B3a

Para

\[A_{\varepsilon,r}(u)=(u^{2r}+\varepsilon^{2r})^{p/(2r)},\qquad r=1,2,\]

la evaluación exacta sobre `u=a x` reproduce para todos `p=1,5/4,3/2,7/4,2`

\[\max|D_xB|\sim\varepsilon^{p-2},\qquad w\sim\varepsilon,\qquad I_Z\sim\varepsilon^{p-1}.\]

En `p=1` la masa normalizada total es `2.000000000006` (r=1) y `2.000000000000` (r=2), compatible con `2 a |a| delta_0` para `a=1`.

## B3b: dato global y control de contaminación

Se usa `u0=(a/k) sin(k x)`, periódico, con `a=1`. Es transversal y suave, sin cutoff. En `x=0`, `u=x+O(x^3)`. Para `L=2`, `k=pi`; el término de curvatura es conocido exactamente. Los ensayos `L=2,4,8` muestran acuerdo hasta aproximadamente `t=0.01`; la fuerte inversión del cero que aparece en `L=2` después de `t~0.015` cambia al alejar el segundo cero y no se interpreta como dinámica local.

## Velocidad instantánea

Para el seno global, la velocidad inicial exacta del cero central es

\[s'_{\varepsilon,1}(0)=\frac{\varepsilon^p}{p+1}+\alpha p a^2\varepsilon^{p-2}-2\alpha k^2\varepsilon^p,\]

\[s'_{\varepsilon,2}(0)=-2\alpha k^2\varepsilon^p.\]

En el límite local exactamente lineal (`k=0`), la segunda es cero. Por tanto, para `1<=p<2`, `r=1` tiene velocidad puntual divergente `~eps^(p-2)` mientras `r=2` tiene velocidad puntual cero. Esto no contradice que ambos kernels tengan la misma clase de escala y, en `p=1`, la misma masa delta: el límite es no uniforme en `t=0` y en el punto `u=0`.

La tabla `B3_INITIAL_SPEED.csv` verifica las fórmulas globales con errores relativos de `1e-8` o menores.

## Convergencia numérica

En `p=3/2, eps=0.06`, 16/24/32 puntos por capa dejan `x_zero(T=0.04)` y la pendiente prácticamente invariantes. Halvar `dt/dx` de `0.1` a `0.05` cambia la posición final en `1.37e-5` relativo (r=1) y `1.12e-4` (r=2). Las derivas de H disminuyen al refinar.

## Evolución limpia en L=8

Con `eps=0.1`, 16 puntos por capa y `T=0.1`, los cinco exponentes y ambas regularizaciones mantienen dos cruces simples. Para `p<=3/2`, el exponente de presión ajustado cerca del cero sigue próximo al valor transversal `p/2` y lejos de `1`; no hay evidencia de conversión a un borde pressure-regular. En `p=3/2` se extendió hasta `T=0.15` y no aparece el giro visto en `L=2`.

## Límite r=1 frente a r=2 a tiempo fijo

En el dominio limpio `L=8`, a `t=0.01`:

- `p=1`: `D12 = 9.149e-4, 4.550e-4, 2.107e-4` para `eps=0.1,0.06,0.03`, ajuste `D12~eps^1.213`.
- `p=3/2`: `D12 = 1.169e-3, 7.484e-4, 3.286e-4`, ajuste `D12~eps^1.062`.

Así, la diferencia de campo no permanece O(1) a tiempo fijo, pese a la incompatibilidad de las velocidades instantáneas. Este resultado es evidencia de una **capa inicial singular**; B3 no pretende sustituir el test de selección B4.

## Lectura para M2

Un cruce transversal `u~a(x-s)` corresponde en presión a `v~sgn(x-s)|x-s|^(p/2)`. Para `p<2` no es pressure-regular. Los ensayos L=8 no muestran que el cruce genérico se transforme automáticamente en un cero lineal de `v`. Por eso la variable de presión no puede considerarse una regularización inocua de cruces transversales. Esto es consistente con la sensibilidad encontrada en B2, aunque no prueba que sea su única causa.

## Decisión

B3 queda cerrado. El siguiente bloque correcto es B4: comparar, en orden separado, `dx->0`, después `eps->0` para M3/M4, y `eta->0` para M5, incluyendo la formulación M2. Las colisiones siguen bloqueadas.
