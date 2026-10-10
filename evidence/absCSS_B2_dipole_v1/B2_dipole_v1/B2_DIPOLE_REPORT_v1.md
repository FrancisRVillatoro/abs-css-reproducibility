# B2 — Propagación del dipolo compacto exacto

## Objetivo

Benchmark B2 de la campaña numérica para |K*|(p+2,p), con

- c = 1,
- alpha = p^2/((p+1)(p+2)),
- p in {1,5/4,3/2,7/4,2},
- semisoporte de cada compactón = pi.

El dipolo exacto es un lóbulo positivo adyacente a uno negativo, ambos con velocidad c. Su soporte es [-2pi,2pi] respecto al cero interno y sus invariantes exactos son C=0, P=2 P_c, H=2 H_c.

La malla principal usa h=2pi/m. Al desplazar el cero interno media celda, las tres interfaces (borde izquierdo, cero interno y borde derecho) mantienen la misma fase subcelda en toda la secuencia. Esta elección evita el defecto semisuave que aparece cuando p=1 y u=0 cae exactamente sobre un nodo.

## 1. Ruta M5: PASS para todo 1 <= p <= 2

Se usa la misma retirada hiperviscosa validada en B1,

    eta = 0.05 h^3 -> 0.

Para m=63,126,188 (h ~ 0.10,0.05,0.033), T=1, todos los casos completan. Los órdenes ajustados del error de forma son:

- p=1: q_shape=2.575, q_radiation=2.168, q_internal_max=2.050.
- p=1.25: q_shape=2.149, q_radiation=1.798, q_internal_max=1.609.
- p=1.5: q_shape=1.795, q_radiation=1.640, q_internal_max=1.201.
- p=1.75: q_shape=1.619, q_radiation=1.759, q_internal_max=0.805.
- p=2: q_shape=1.480, q_radiation=1.274, q_internal_max=1.661.

La amplitud del error local alrededor del cero interno y la radiación fuera del soporte disminuyen con h. La anchura definida por un error mayor que 10^-3 A se hace nula en las mallas refinadas en todos los casos. C permanece a precisión de redondeo por la simetría del dipolo; las derivas de P y H decrecen con h.

## 2. Velocidad del cero interno

En marco de laboratorio (c0=0), M5 con m=126 produce los siguientes ajustes lineales x_zero(t)=x0+c_num t:

- p=1: c_num=0.993662738, c_num-1=-6.337e-03.
- p=1.25: c_num=0.996880075, c_num-1=-3.120e-03.
- p=1.5: c_num=0.998507394, c_num-1=-1.493e-03.
- p=1.75: c_num=0.999422917, c_num-1=-5.771e-04.
- p=2: c_num=1.000056171, c_num-1=+5.617e-05.

Para p=1 y 5/4, donde el cero es de orden mayor y su localización puntual es peor condicionada, el refinamiento m=188 mejora los ajustes a c_num=0.995678423 y 0.997968381, respectivamente. La cinemática c=1 queda confirmada.

## 3. M1 sin regularización

La fase del cero es crítica en p=1: con un cero exactamente nodal el Newton tropieza con la no diferenciabilidad de sign(u); con las tres interfaces a media celda, M1 completa y converge. En particular:

- p=1: E = 1.9867e-4 (m=63) -> 8.4434e-5 (m=126) -> 2.5591e-5 (m=188).
- p=5/4: E = 5.0196e-3 -> 1.5207e-3 -> 5.3678e-4.
- p=3/2: E = 5.2912e-2 -> 1.4851e-2 -> 5.9408e-3 -> 2.3680e-3 (m=251).
- p=7/4: las mallas gruesas entran en una transición tardía; m=188 y 251 completan con E=4.6130e-2 y 2.2971e-2.
- p=2: T=1 sigue preasintótico y M1 falla antes de T=1 hasta m=251. Sin embargo, a T=0.5 converge limpiamente: E=6.3957e-2, 1.7582e-2, 8.1435e-3, 4.6307e-3 para m=63,126,188,251, con q=1.897 y radiación q=1.811.

Por tanto, M1 es consistente pero extremadamente ineficiente cerca de p=2 para T=1; no debe confundirse este umbral temporal con ausencia de convergencia a tiempo corto.

## 4. Formulación de presión M2

Con eta_v=0.1 h^2, M2 converge para p >= 3/2. Los órdenes de forma son:

- p=1.5: q_shape=0.912, q_radiation=1.121, q_C=0.947.
- p=1.75: q_shape=0.978, q_radiation=0.890, q_C=0.800.
- p=2: q_shape=0.985, q_radiation=0.752, q_C=2.216.

P y H también presentan derivas decrecientes. En las mallas m=188 se usa tolerancia Newton absoluta 1e-8 porque el residual bruto de las ecuaciones auxiliares contiene D4~h^-4 y se estanca en unos pocos 10^-9; esta tolerancia sigue muy por debajo del error espacial y no cambia los observables a la precisión reportada.

### Régimen bajo p

Para p=5/4, eta_v=0.1h^2 funciona en m=63 pero falla al refinar. La trayectoria estabilizante eta_v ~= 0.01h (mu=0.1,0.2,0.3 para m=63,126,188) produce errores 2.954e-3, 2.752e-3, 6.007e-3: el punto fino empeora y C_abs crece hasta 1.86e-2.

Para p=1 la situación es más fuerte. Regularizaciones crecientes permiten completar corridas individuales (por ejemplo mu=0.5 en m=63, mu=1.5 en m=126, mu=2.5 en m=188), pero la intensidad necesaria crece con el refinamiento y los errores estabilizados no forman una secuencia convergente al dipolo exacto. En m=188, mu=2.5 completa con E=2.372e-2 y |C|=9.845e-2.

Esto es dependencia de selección/formulación, no un fallo del dipolo exacto: M1/M5 sí lo propagan en ese mismo régimen.

## 5. Comparación directa M1–M2 en p=3/2

La norma requerida por el protocolo, sin realineamiento adicional, es:

- m=63, h=0.099733: ||u_M1-u_M2||_2/||u_M1||_2 = 5.292008e-02; fase relativa=1.369e-03.
- m=126, h=0.049867: ||u_M1-u_M2||_2/||u_M1||_2 = 1.491796e-02; fase relativa=3.299e-04.
- m=188, h=0.033421: ||u_M1-u_M2||_2/||u_M1||_2 = 6.015329e-03; fase relativa=1.942e-04.

Ajuste: direct_relL2 ~ h^1.971. La misma potencia se obtiene tras realinear fase. Este criterio B2 pasa limpiamente en p=3/2.


## 6. Refinamiento temporal

En p=3/2, m=126, reducir Delta t/h de 0.1 a 0.05 cambia el error de forma sólo en:

- M1: cambio relativo=1.056e-04.
- M5: cambio relativo=1.146e-05.
- M2: cambio relativo=2.369e-05.

Por tanto el benchmark está dominado por error espacial y por la selección de regularización, no por el integrador R_{1,2}.

## 7. Veredicto B2

- M5: PASS para los cinco p.
- Cero interno: PASS cinemático, c_num -> 1.
- Ausencia de capa fija: PASS en M5 y en M2 para p>=3/2.
- M1–M2: PASS demostrado en p=3/2; M2 converge además para p=7/4 y 2, pero M1 sigue fuertemente preasintótico allí a T=1.
- p=1 y 5/4: FAIL de equivalencia global M1–M2 bajo las rutas de presión ensayadas; aparece dependencia de selección al retirar la regularización.

Por las reglas de parada del protocolo, B2 no puede etiquetarse PASS global en [1,2]. El estado correcto es

    B2 = FAIL_GLOBAL_FORMULATION_AGREEMENT
         + PASS_EXACT_DIPOLE_M5_ALL_P
         + PASS_REGULAR_SUBREGIME_p>=3/2

Esto NO invalida la proposición analítica de que el dipolo es una solución débil exacta. Al contrario: el desacuerdo numérico cerca de p=1 indica que la evolución débil alrededor de cambios de signo puede depender de la realización/regularización, precisamente la cuestión científica central del proyecto.

No deben iniciarse colisiones. El siguiente benchmark lógico sigue siendo B3 (cruce transversal regularizado), porque proporciona el diagnóstico local necesario antes de B4 (selección dinámica completa).
