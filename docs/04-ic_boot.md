# Intervalos de confianza bootstrap {#icboot}





Consideremos el problema de construcción, mediante bootstrap, de un
intervalo de confianza bilateral, con nivel de confianza $1-\alpha$,
para un parámetro $\theta$ de la distribución $F$. Una vez elegido el
método bootstrap adecuado a la información disponible en el contexto del
que se trate, un aspecto importante es el de la posible corrección de
los intervalos de confianza bootstrap, aproximados por el método de
Monte Carlo, al objeto de que la probabilidad de cobertura sea lo más
parecida posible al nivel nominal $1-\alpha$. Comenzaremos analizando
el error de cobertura de los intervalos de confianza clásicos, los
basados en la distribución normal asintótica.

## Intervalos basados en la distribución normal asintótica {#icboot-norm}

Consideremos primeramente el caso más sencillo (y poco realista) de
construcción de un intervalo de confianza para la media, $\mu \,$, con
desviación típica, $\sigma$, conocida. El estadístico usado para
construir el intervalo de confianza es
$$R=\sqrt{n}\frac{\bar{X}-\mu }{\sigma }$$que cuando
$n\rightarrow \infty$ tiende en distribución a una $N\left(
0,1 \right)$. El intervalo de confianza basado en dicha aproximación
normal es
$\hat{I}=\left( \bar{X}-\frac{\sigma }{\sqrt{n}}z_{\alpha /2},
\bar{X}+\frac{\sigma }{\sqrt{n}}z_{\alpha /2} \right)$. Mediante el
desarrollo de Edgeworth dado por el Teorema de Cramer es fácil obtener
una cota para el error de cobertura de dicho intervalo:

$$\begin{aligned}
P\left( \mu \in \hat{I} \right) -\left( 1-\alpha \right) =&\  P\left(
R<z_{\alpha /2} \right) -P\left( R\leq -z_{\alpha /2} \right) \\
&-\left( \Phi \left( z_{\alpha /2} \right) -\Phi \left( -z_{\alpha
/2} \right) \right) \\
=&\ \ n^{-\frac{1}{2}}p_1\left( z_{\alpha /2} \right) \phi \left( z_{\alpha
/2} \right) +O\left( n^{-1} \right) \\
&-\left( n^{-\frac{1}{2}}p_1\left( -z_{\alpha /2} \right) \phi \left(
-z_{\alpha /2} \right) +O\left( n^{-1} \right) \right) \\
=&\  O\left( n^{-1} \right),\end{aligned}$$

ya que por la simetría de las funciones $p_1\left( u \right)$ y $\phi
\left( u \right)$ se tiene
$$p_1\left( -z_{\alpha /2} \right) \phi \left( -z_{\alpha /2} \right)
=p_1\left( z_{\alpha /2} \right) \phi \left( z_{\alpha /2} \right).$$
De esta forma, el orden del error de cobertura del intervalo de confianza
bilateral con desviación típica poblacional conocida es $O\left(
n^{-1} \right)$. Puede obtenerse fácilmente el orden del error de
cobertura de los intervalos unilaterales que resulta ser
$O\left( n^{-\frac{1}{2}} \right)$.

En el caso más realista en que la desviación típica, $\sigma$, sea
desconocida, el intervalo de confianza resulta 
$$\hat{I}_{0}=\left( 
\bar{X}-\frac{S_n}{\sqrt{n}}z_{\alpha /2},\bar{X}+\frac{S_n}{
\sqrt{n}}z_{\alpha /2} \right).$$
Ahora, el estadístico en el que se basa la inferencia resulta:
$$R_1=\sqrt{n}\frac{\bar{X}-\mu }{S_n}.$$
Un desarrollo de Edgeworth del tipo del obtenido en el Teorema de Bhattacharya-Ghosh
permite acotar el error de cobertura de este intervalo:
$$\begin{aligned}
P\left( \mu \in \hat{I}_{0} \right) -\left( 1-\alpha \right) =&\ P\left(
R_1<z_{\alpha /2} \right) -P\left( R_1\leq -z_{\alpha /2} \right) \\
&-\left( \Phi \left( z_{\alpha /2} \right) -\Phi \left( -z_{\alpha
/2} \right) \right) \\
=&\ n^{-\frac{1}{2}}q_1\left( z_{\alpha /2} \right) \phi \left( z_{\alpha
/2} \right) +O\left( n^{-1} \right) \\
&-\left( n^{-\frac{1}{2}}q_1\left( -z_{\alpha /2} \right) \phi \left(
-z_{\alpha /2} \right) +O\left( n^{-1} \right) \right) \\
=&\ O\left( n^{-1} \right).
\end{aligned}$$

De esta forma, el orden del error de cobertura del intervalo de
confianza bilateral con desviación típica desconocida es $O\left(
n^{-1} \right)$. El orden del error de cobertura para el intervalo de
confianza unilateral resulta $O\left( n^{-\frac{1}{2}} \right)$.

Si el parámetro de interés fuese otro arbitrario: $\theta =\theta
\left( F \right)$, no necesariamente la media, puede obtenerse,
análogamente un intervalo de confianza basado en la normal asintótica:
$$\hat{I}_{0}=\left( \hat{\theta}-\frac{\hat{\sigma}_{\theta }}{\sqrt{n}}
z_{\alpha /2},\hat{\theta}+\frac{\hat{\sigma}_{\theta }}{\sqrt{n}}
z_{\alpha/2} \right),$$
que está basado en el estadístico
$$R_1=\sqrt{n}\frac{\hat{\theta}-\theta }{\hat{\sigma}_{\theta }}.$$
De forma análoga a lo ya razonado para la media muestral, puede
deducirse que el error de cobertura del intervalo de confianza bilateral
tiene un orden de $O\left( n^{-1} \right)$, mientras que para intervalos
unilaterales el orden es $O\left( n^{-\frac{1}{2}} \right)$.

<!-- 
describir `type="norm"` en la llamada a la función `boot.ci()` 
para el caso general
$\frac{\hat{\sigma}_{\theta}}{\sqrt{n}}$ err. std.
-->


## Método percentil (básico) {#icboot-basic}

Este método se basa en la construcción del intervalo de confianza,
mediante bootstrap, a partir del estadístico no estandarizado
$$R_2=\sqrt{n}\left( \hat{\theta}-\theta \right).$$
Una vez realizado el correspondiente remuestreo (uniforme, suavizado,
simetrizado, ...), a partir de cierto estimador, $\hat{F}\,$, de la
distribución poblacional, $F$, la distribución en el muestreo de
$R_2$ se aproxima mediante la distribución bootstrap de
$$R_2^{\ast}=\sqrt{n}\left( \hat{\theta}^{\ast}-\theta \left( \hat{F}
 \right) \right).$$
Así se obtienen valores $x_{\alpha /2}$ y
$x_{1-\alpha /2}$, siendo $x_{\beta }$, tal que
$P^{\ast}\left( R_2^{\ast }\leq x_{\beta } \right) =\beta$, 
y a partir de ellos sabemos que
$$\begin{aligned}
1-\alpha &= 1-\frac{\alpha }{2}-\frac{\alpha }{2}=P^{\ast}\left(
R_2^{\ast}<x_{1-\alpha /2} \right) -P^{\ast}\left( R_2^{\ast}\leq
x_{\alpha /2} \right) \\
&= P^{\ast}\left( x_{\alpha /2}<R_2^{\ast}<x_{1-\alpha /2} \right),
\end{aligned}$$
con lo cual decimos que también ha de ser aproximadamente igual a
$1-\alpha$ la siguiente probabilidad$$\begin{aligned}
P\left( x_{\alpha /2}<R_2<x_{1-\alpha /2} \right) &= P\left( x_{\alpha /2}<
\sqrt{n}\left( \hat{\theta}-\theta \right) <x_{1-\alpha /2} \right) \\
&= P\left( \hat{\theta}-\frac{x_{1-\alpha /2}}{\sqrt{n}}<\theta <\hat{\theta}
-\frac{x_{\alpha /2}}{\sqrt{n}} \right).\end{aligned}$$
Ello da pie a definir el intervalo de confianza bootstrap calculado 
por el método percentil como
$$\hat{I}_1=\left( \hat{\theta}-\frac{x_{1-\alpha /2}}{\sqrt{n}},\hat{\theta}
-\frac{x_{\alpha /2}}{\sqrt{n}} \right).$$

Para estudiar el error de cobertura de este intervalo de confianza
conviene ver antes qué grado de aproximación existe entre la
distribución en el muestreo de $R_2$ y la distribución bootstrap de
$R_2^{\ast}$. A partir del Teorema de Bhattacharya-Ghosh se
tiene$$\begin{aligned}
P\left( R_2\leq v \right) &= P\left( \sqrt{n}\left( \hat{\theta}-\theta
 \right) \leq v \right) =P\left( \sqrt{n}\frac{\hat{\theta}-\theta }{\sigma
_{\theta }}\leq \frac{v}{\sigma _{\theta }} \right) \\
&= \Phi \left( \frac{v}{\sigma _{\theta }} \right) +O\left( n^{-\frac{1}{2}
} \right) \\
P^{\ast}\left( R_2^{\ast}\leq v \right) &= P^{\ast}\left( \sqrt{n}\left( 
\hat{\theta}^{\ast}-\theta \left( \hat{F} \right) \right) \leq v \right)
\\
&=  P^{\ast}\left( \sqrt{n}\frac{\hat{\theta}^{\ast}-\theta \left( \hat{F}
 \right)}{\sigma _{\hat{\theta}}}\leq \frac{v}{\sigma _{\hat{\theta}}} \right)
=\Phi \left( \frac{v}{\sigma _{\hat{\theta}}} \right) +O_{P}\left( n^{-
\frac{1}{2}} \right).
\end{aligned}$$

Como consecuencia
$$P^{\ast}\left( R_2^{\ast}\leq v \right) -P\left( R_2\leq v \right) =\Phi
\left( \frac{v}{\sigma _{\hat{\theta}}} \right) -\Phi \left( \frac{v}{\sigma
_{\theta }} \right) +O_{P}\left( n^{-\frac{1}{2}} \right) =O_{P}\left( n^{-
\frac{1}{2}} \right),$$
ya que, típicamente, $\sigma _{\hat{\theta}}-\sigma _{\theta}
= O_{P}\left( n^{-\frac{1}{2}} \right)$. En resumen, la distribución en
el muestreo de $R_2$ y la distribución bootstrap de $R_2^{\ast}$ se
aproximan, una a la otra, a la velocidad $O_{P}\left( n^{-\frac{1}{2}
} \right)$, cuando $n\rightarrow \infty$.

El error de cobertura del intervalo de confianza bilateral calculado
mediante bootstrap por el método percentil es
$$\begin{aligned}
P\left( \theta \in \hat{I}_1 \right) -\left( 1-\alpha \right) =&\ P\left(
R_2<x_{1-\alpha /2} \right) -P\left( R_2\leq x_{\alpha /2} \right) \\
&-\left[ P^{\ast}\left( R_2^{\ast}<x_{1-\alpha /2} \right) -P^{\ast
}\left( R_2^{\ast}\leq x_{\alpha /2} \right) \right] \\
=&\ P\left( R_2<x_{1-\alpha /2} \right) -P^{\ast}\left( R_2^{\ast
}<x_{1-\alpha /2} \right) \\
&-\left[ P\left( R_2\leq x_{\alpha /2} \right) -P^{\ast}\left(
R_2^{\ast}\leq x_{\alpha /2} \right) \right] \\
=&\ O_{P}\left( n^{-\frac{1}{2}} \right)
\end{aligned}$$

De esta forma el error de cobertura para los intervalos de confianza
bilaterales bootstrap obtenidos mediante el método percentil es $O\left(
n^{-\frac{1}{2}} \right)$. Puede deducirse que ese es también el orden
para los intervalos unilaterales obtenidos por este método. Así pues el
orden del error de cobertura para el método percentil cuando se
construyen intervalos de confianza unilaterales coincide con el de los
construidos usando la normal asintótica pero el orden del error de
cobertura de los intervalos bilaterales bootstrap constuidos por el
método percentil es peor que el de los basados en la normal asintótica,
que es del orden $O\left( n^{-1} \right)$.


\BeginKnitrBlock{example}\iffalse{-91-73-110-102-101-114-101-110-99-105-97-32-115-111-98-114-101-32-108-97-32-109-101-100-105-97-32-99-111-110-32-118-97-114-105-97-110-122-97-32-100-101-115-99-111-110-111-99-105-100-97-44-32-99-111-110-116-105-110-117-97-99-105-243-110-93-}\fi{}
<span class="example" id="exm:media-dt-desconocida-perc"><strong>(\#exm:media-dt-desconocida-perc)  \iffalse (Inferencia sobre la media con varianza desconocida, continuación) \fi{} </strong></span><br> \vspace{0.5cm}

Continuando con el ejemplo de los tiempos de vida de microorganismos
(sin asumir varianza conocida),
supongamos que queremos obtener una estimación por intervalo de confianza 
de su vida media empleando este método. 
El código necesario sería muy similar al del Ejemplo \@ref(exm:media-dt-desconocida):
\EndKnitrBlock{example}


```r
muestra <- c(0.143, 0.182, 0.256, 0.26, 0.27, 0.437, 0.509, 
             0.611, 0.712, 1.04, 1.09, 1.15, 1.46, 1.88, 2.08)
n <- length(muestra)
alfa <- 0.05
x_barra <- mean(muestra)

# Remuestreo
set.seed(1)
B <- 1000
remuestra <- numeric(n)
x_barra_boot <- numeric(B) 
estadistico_boot <- numeric(B)
for (k in 1:B) {
  remuestra <- sample(muestra, n, replace = TRUE)
  x_barra_boot[k] <- mean(remuestra)
  estadistico_boot[k] <- sqrt(n) * (x_barra_boot[k] - x_barra)
}

# Aproximación bootstrap de los ptos críticos
pto_crit <- quantile(estadistico_boot, c(alfa/2, 1 - alfa/2))

# Construcción del IC
ic_inf_boot <- x_barra - pto_crit[2]/sqrt(n)
ic_sup_boot <- x_barra - pto_crit[1]/sqrt(n)
IC_boot <- c(ic_inf_boot, ic_sup_boot)
names(IC_boot) <- paste0(100*c(alfa/2, 1-alfa/2), "%")
IC_boot
```

```
##      2.5%     97.5% 
## 0.4837233 1.1025650
```

Aunque en este caso también podemos obtener el intervalo a
partir de las réplicas bootstrap del estimador:

```r
pto_crit <- quantile(x_barra_boot, c(alfa/2, 1 - alfa/2))
ic_inf_boot <- 2*x_barra - pto_crit[2]
ic_sup_boot <- 2*x_barra - pto_crit[1]
IC_boot <- c(ic_inf_boot, ic_sup_boot)
names(IC_boot) <- paste0(100*c(alfa/2, 1-alfa/2), "%")
IC_boot
```

```
##      2.5%     97.5% 
## 0.4837233 1.1025650
```

Esta forma de proceder es la que emplea el paquete `boot` para obtener
el que denomina intervalo de confianza *bootstrap básico*
(estableciendo `type="basic"` en la llamada a la función `boot.ci()`):

```r
library(boot)
statistic <- function(data, i) mean(data[i])

set.seed(1)
res.boot <- boot(muestra, statistic, R = 1000)
res <- boot.ci(res.boot, type = "basic")
res
```

```
## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
## Based on 1000 bootstrap replicates
## 
## CALL : 
## boot.ci(boot.out = res.boot, type = "basic")
## 
## Intervals : 
## Level      Basic         
## 95%   ( 0.4825,  1.0980 )  
## Calculations and Intervals on Original Scale
```

```r
IC_boot <- res$basic[4:5]
IC_boot
```

```
## [1] 0.4824717 1.0980120
```

Además del paquete `boot`, otros autores también denominan a este método
*bootstrap básico* (*bootstrap percentil básico* o incluso *bootstrap natural*), 
y utilizan la terminología *bootstrap percentil* cuando se emplea 
directamente el estimador como estadístico ($R = \hat \theta$) para 
realizar inferencia. Con el paquete `boot` habrá que establecer `type="perc"` 
en la llamada a la función `boot.ci()` para obtener el intervalo
correspondiente:

```r
boot.ci(res.boot, type = "perc")
```

```
## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
## Based on 1000 bootstrap replicates
## 
## CALL : 
## boot.ci(boot.out = res.boot, type = "perc")
## 
## Intervals : 
## Level     Percentile     
## 95%   ( 0.5127,  1.1282 )  
## Calculations and Intervals on Original Scale
```

En este método se emplean directamente los cuantiles de las
réplicas bootstrap del estadístico:

```r
# IC_boot <-  quantile(res.boot$t, c(alfa/2, 1 - alfa/2)) # type = 7
IC_boot <-  quantile(res.boot$t, c(alfa/2, 1 - alfa/2), type = 6)
IC_boot 
```

```
##      2.5%     97.5% 
## 0.5126517 1.1281950
```

Asintóticamente ambos métodos son equivalentes, aunque en general
es preferible (evita sesgos) el bootstrap percentil básico.

<!-- 
ejercicio 
Obtener los ICs mediante los método básico y percentil
de la varianza y cuasivarianza muestrales.

Recorrido intercuartílico?
-->

## Método percentil-*t* {#icboot-perc-t}

Este método bootstrap, construye un intervalo de confianza bootstrap a
partir del estadístico studentizado:
$$R_1=\sqrt{n}\frac{\hat{\theta}-\theta }{\hat{\sigma}_{\theta }}.$$
Su distribución en el muestreo se aproxima mediante la distribución
bootstrap de
$$R_1^{\ast}=\sqrt{n}\frac{\hat{\theta}^{\ast}-\theta \left( \hat{F}
 \right)}{\hat{\sigma}_{\theta }^{\ast}}.$$
En este caso, los valores $x_{\alpha /2}$ y
$x_{1-\alpha /2}$, se obtienen a partir de esta última distribución
bootstrap, es decir, $x_{\beta }$ se define a partir de $P^{\ast}\left(
R_1^{\ast}\leq x_{\beta } \right) =\beta$.
Como
$$1-\alpha =P^{\ast}\left( x_{\alpha /2}<R_1^{\ast}<x_{1-\alpha /2} \right),$$
razonamos que también ha de ser aproximadamente igual a $1-\alpha$ la
siguiente probabilidad
$$\begin{aligned}
P\left( x_{\alpha /2}<R_1<x_{1-\alpha /2} \right) &= P\left( x_{\alpha /2}<
\sqrt{n}\frac{\hat{\theta}-\theta }{\hat{\sigma}_{\theta }}<x_{1-\alpha
/2} \right) \\
&= P\left( \hat{\theta}-\frac{\hat{\sigma}_{\theta }}{\sqrt{n}}x_{1-\alpha
/2}<\theta <\hat{\theta}-\frac{\hat{\sigma}_{\theta }}{\sqrt{n}}x_{\alpha
/2} \right).
\end{aligned}$$
Con lo cual, el intervalo de confianza bootstrap calculado 
por el método percentil-$t$ se define como
$$\hat{I}_2=\left( \hat{\theta}-\frac{\hat{\sigma}_{\theta }}{\sqrt{n}}
x_{1-\alpha /2},\hat{\theta}-\frac{\hat{\sigma}_{\theta }}{\sqrt{n}}
x_{\alpha /2} \right).$$

Utilizando el Teorema de Bhattacharya-Ghosh puede acotarse el error de
aproximación entre la distribución en el muestreo de $R_1$ y la
distribución bootstrap de $R_1^{\ast}$:
$$\begin{aligned}
P^{\ast}\left( R_1^{\ast}\leq u \right) -P\left( R_1\leq u \right)
=&\ P^{\ast}\left( \sqrt{n}\frac{\hat{\theta}^{\ast}-\theta \left( \hat{F}
 \right)}{\hat{\sigma}_{\theta }^{\ast}}\leq u \right)-P\left( \sqrt{n}
\frac{\hat{\theta}-\theta }{\hat{\sigma}_{\theta }}\leq u \right)  \\
=&\  \Phi \left( u \right) +n^{-\frac{1}{2}}\hat{q}_1\left( u \right) \phi
\left( u \right) +O_{P}\left( n^{-1} \right) \\
&-\left[ \Phi \left( u \right) +n^{-\frac{1}{2}}q_1\left( u \right) \phi
\left( u \right) +O\left( n^{-1} \right) \right] \\
=&\  n^{-\frac{1}{2}}\left[ \hat{q}_1\left( u \right) -q_1\left( u \right) 
\right] \phi \left( u \right) +O_{P}\left( n^{-1} \right) \\
=&\ O_{P}\left( n^{-1} \right),
\end{aligned}$$
ya que los coeficientes que aparecen en el polinomio $\hat{q}_1\left(
u \right)$ son estimadores $\sqrt{n}$-consistentes de los coeficientes
del polinomio $q_1\left( u \right)$. Los de éste último dependen de
los momentos poblacionales y los del primero son sus correspondientes
versiones empíricas.

Así pues, el error de cobertura del intervalo de confianza bootstrap
bilateral calculado por el método percentil-$t$ es
$$\begin{aligned}
P\left( \theta \in \hat{I}_2 \right) -\left( 1-\alpha \right) =&\  P\left(
R_1<x_{1-\alpha /2} \right) -P\left( R_1\leq x_{\alpha /2} \right) \\
& -\left[ P^{\ast}\left( R_1^{\ast}<x_{1-\alpha /2} \right) -P^{\ast
}\left( R_1^{\ast}\leq x_{\alpha /2} \right) \right] \\
=&\  P\left( R_1<x_{1-\alpha /2} \right) -P^{\ast}\left( R_1^{\ast
}<x_{1-\alpha /2} \right) \\
&-\left[ P\left( R_1\leq x_{\alpha /2} \right) -P^{\ast}\left(
R_1^{\ast}\leq x_{\alpha /2} \right) \right] \\
=&\  O_{P}\left( n^{-1} \right)
\end{aligned}$$

Se tiene entonces que el error de cobertura para los intervalos de
confianza bilaterales bootstrap obtenidos mediante el método
percentil-$t$ es
$O\left( n^{-1} \right)$. Puede deducirse que ese es también el orden
para los intervalos unilaterales obtenidos por este método. Así pues el
orden del error de cobertura para el método percentil-$t$ cuando se
construyen intervalos de confianza unilaterales mejora al de los
intervalos unilaterales basados en la normal asintótica, con un error de
cobertura de orden $O\left( n^{-\frac{1}{2}} \right)$. En el caso de los
intervalos de confianza bilaterales, el orden del error de cobertura
usando la normal asintótica o bien el bootstrap por el método
percentil-$t$ es el mismo, $O\left( n^{-1} \right)$ en ambos casos.

En el Ejemplo \@ref(exm:media-dt-desconocida), se implementó
este método para obtener una estimación por intervalo de confianza
del tiempo de vida medio de microorganismos.
En el Ejemplo \@ref(exm:media-dt-desconocida-boot) se mostró como
calcular este intervalo empleando el paquete `boot` (haciendo que 
la función `statistic` devuelva también la varianza del estadístico 
y estableciendo `type="stud"` en la llamada a la función `boot.ci()`).


## Método percentil-*t* simetrizado {#icboot-perc-t-sim}

Es un método análogo al percentil-$t$. Sólo difiere de él en la forma de
seleccionar los cuantiles de la distribución bootstrap. En lugar de
tomar cuantiles que dejen colas iguales ($\frac{\alpha }{2}$ a la
izquierda y a la derecha, respectivamente), se eligen los cuantiles de
forma que sean simétricos. Así, dado el estadístico $R_1$ y su versión
bootstrap $R_1^{\ast}$, se considera el valor $y_{1-\alpha }$ que
cumple $P^{\ast}\left( \left\vert R_1^{\ast}\right\vert
\leq y_{1-\alpha } \right) =1-\alpha$. 
Así se tiene que
$$1-\alpha =P^{\ast}\left( -y_{1-\alpha }\leq R_1^{\ast}
\leq y_{1-\alpha} \right).$$
De esa forma se razona que también ha de ser aproximadamente
igual a $1-\alpha$ la siguiente probabilidad
$$\begin{aligned}
P\left( -y_{1-\alpha }<R_1<y_{1-\alpha } \right) &= P\left( -y_{1-\alpha }
< \sqrt{n}\frac{\hat{\theta}-\theta }{\hat{\sigma}_{\theta }}
< y_{1-\alpha} \right) \\
&= P\left( \hat{\theta}-\frac{\hat{\sigma}_{\theta }}{\sqrt{n}}y_{1-\alpha
}<\theta <\hat{\theta}+\frac{\hat{\sigma}_{\theta }}{\sqrt{n}}y_{1-\alpha
} \right).
\end{aligned}$$

Con lo cual, el intervalo de confianza bootstrap calculado por el método
percentil-$t$ simetrizado se define como
$$\hat{I}_3=\left( \hat{\theta}-\frac{\hat{\sigma}_{\theta }}{\sqrt{n}}
y_{1-\alpha },\hat{\theta}+\frac{\hat{\sigma}_{\theta }}{\sqrt{n}}
y_{1-\alpha } \right).$$

Utilizando el Teorema de Bhattacharya-Ghosh con desarrollos hasta el
orden $n^{-1}$ se tiene:
$$\begin{aligned}
P^{\ast}\left( R_1^{\ast}\leq u \right) -P\left( R_1\leq u \right) 
=&\ P^{\ast}\left( \sqrt{n}\frac{\hat{\theta}^{\ast}-\theta \left( \hat{F}
 \right)}{\hat{\sigma}_{\theta }^{\ast}}\leq u \right) -P\left( \sqrt{n}
\frac{\hat{\theta}-\theta }{\hat{\sigma}_{\theta }}\leq u \right) \\
=&\ \Phi \left( u \right) +n^{-\frac{1}{2}}\hat{q}_1\left( u \right) \phi
\left( u \right) +n^{-1}\hat{q}_2\left( u \right) \phi \left( u \right)
+O_{P}\left( n^{-\frac{3}{2}} \right) \\
& -\left[ \Phi \left( u \right) +n^{-\frac{1}{2}}q_1\left( u \right) \phi
\left( u \right) +n^{-1}q_2\left( u \right) \phi \left( u \right) +O\left(
n^{-\frac{3}{2}} \right) \right] \\
=&\ n^{-\frac{1}{2}}\left[ \hat{q}_1\left( u \right) -q_1\left( u \right) 
\right] \phi \left( u \right) 
+n^{-1}\left[ \hat{q}_2\left( u \right) -q_2\left( u \right) \right] \phi
\left( u \right) +O_{P}\left( n^{-\frac{3}{2}} \right).
\end{aligned}$$

Como consecuencia, el error de cobertura del intervalo de confianza
bootstrap bilateral calculado por el método percentil-$t$ simetrizado
es
$$\begin{aligned}
P\left( \theta \in \hat{I}_3 \right) -\left( 1-\alpha \right) 
=&\ P\left(-y_{1-\alpha }<R_1<y_{1-\alpha } \right) -P^{\ast}\left( -y_{1-\alpha
}<R_1^{\ast}<y_{1-\alpha } \right) \\
=&\ P\left( R_1<y_{1-\alpha } \right) -P^{\ast}\left( R_1^{\ast
}<y_{1-\alpha } \right) - \\
&-\left[ P\left( R_1\leq -y_{1-\alpha } \right) -P^{\ast}\left( R_1^{\ast
}\leq -y_{1-\alpha } \right) \right] \\
=&\ n^{-\frac{1}{2}}\left[ q_1\left( y_{1-\alpha } \right) -\hat{q}_1\left(
y_{1-\alpha } \right) \right] \phi \left( y_{1-\alpha } \right) \\
&+n^{-1}\left[ q_2\left( y_{1-\alpha } \right) -\hat{q}_2\left(
y_{1-\alpha } \right) \right] \phi \left( y_{1-\alpha } \right) \\
&-n^{-\frac{1}{2}}\left[ q_1\left( -y_{1-\alpha } \right) -\hat{q}_1\left(
-y_{1-\alpha } \right) \right] \phi \left( -y_{1-\alpha } \right) \\
&-n^{-1}\left[ q_2\left( -y_{1-\alpha } \right) -\hat{q}_2\left(
-y_{1-\alpha } \right) \right] \phi \left( -y_{1-\alpha } \right) +O_{P}\left(
n^{-\frac{3}{2}} \right) \\
=&\ 2n^{-1}\left[ q_2\left( y_{1-\alpha } \right) -\hat{q}_2\left(
y_{1-\alpha } \right) \right] \phi \left( y_{1-\alpha } \right) +O_{P}\left(
n^{-\frac{3}{2}} \right) \\
=&\ O_{P}\left( n^{-\frac{3}{2}} \right)
\end{aligned}$$

ya que los polinomios $q_1\left( u \right)$ y
$\hat{q}_1\left( u \right)$ son simétricos, $q_2\left( u \right)$ y
$\hat{q}_2\left( u \right)$ son antisimétricos, la función
$\phi \left( u \right)$ es simétrica y los coeficientes que aparecen en
el polinomio $\hat{q}_2\left(
u \right)$ son estimadores $\sqrt{n}$-consistentes de los coeficientes
del polinomio $q_2\left( u \right)$. Como consecuencia, el error de
cobertura para los intervalos de confianza bilaterales bootstrap
obtenidos mediante el método percentil-$t$ simetrizado es
$O\left( n^{-\frac{3}{2}} \right)$. Puede deducirse que el orden para
los intervalos unilaterales obtenidos por este método es
$O\left( n^{-1} \right)$. Así pues el orden del error de cobertura para
el método percentil-$t$ simetrizado cuando se construyen intervalos de
confianza unilaterales mejora al de los intervalos unilaterales basados
en la normal asintótica, con un error de cobertura de orden
$O\left( n^{-\frac{1}{2}} \right)$, e iguala al orden del error de
cobertura de los obtenidos mediante el percentil-$t$. En el caso de los
intervalos de confianza bilaterales, el orden del error de cobertura
usando el método percentil-$t$ simetrizado es
$O\left( n^{-\frac{3}{2}} \right)$, el cual mejora el orden
$O\left( n^{-1} \right)$, que es el que presentan los intervalos basados
en la normal asintótica o bien en el método percentil-$t$.

## Tabla resumen de los errores de cobertura {#cap-err-cober}

  Tipo de I.C.                Unilateral                              Bilateral 
  ------------                --------------------------              -------------------------- 
  Percentil                   $O_{P}\left( n^{-\frac{1}{2}}\right)$   $O_{P}\left(n^{-\frac{1}{2}} \right)$     
  Percentil-$t$               $O_{P}\left( n^{-1} \right)$            $O_{P}\left( n^{-1} \right)$
  Percentil-$t$ simetrizado   $O_{P}\left( n^{-1} \right)$            $O_{P}\left( n^{-\frac{3}{2}} \right)$
          
                

## Ejemplos {#icboot-ejem}

### IC bootstrap para la media mediante el método percentil-*t* simetrizado {#media-dt-desconocida-persim}

Continuando con el ejemplo de los tiempos de vida de microorganismos
(sin asumir varianza conocida),
para obtener una estimación por intervalo de confianza 
de su vida media empleando el método bootstrap percentil-*t* simetrizado 
se podría utilizar (ver Ejemplo \@ref(exm:media-dt-desconocida)):


```r
muestra <- c(0.143, 0.182, 0.256, 0.26, 0.27, 0.437, 0.509, 
             0.611, 0.712, 1.04, 1.09, 1.15, 1.46, 1.88, 2.08)
n <- length(muestra)
alfa <- 0.05
x_barra <- mean(muestra)
cuasi_dt <- sd(muestra)

# Remuestreo
set.seed(1)
B <- 1000
remuestra <- numeric(n)
estadistico_boot <- numeric(B)
for (k in 1:B) {
  remuestra <- sample(muestra, n, replace = TRUE)
  x_barra_boot <- mean(remuestra)
  cuasi_dt_boot <- sd(remuestra)
  estadistico_boot[k] <- sqrt(n) * abs(x_barra_boot - x_barra)/cuasi_dt_boot
}

# Aproximación bootstrap del pto crítico
pto_crit <- quantile(estadistico_boot, 1 - alfa)

# Construcción del IC
ic_inf_boot <- x_barra - pto_crit * cuasi_dt/sqrt(n)
ic_sup_boot <- x_barra + pto_crit * cuasi_dt/sqrt(n)
IC_boot <- c(ic_inf_boot, ic_sup_boot)
names(IC_boot) <- paste0(100*c(alfa/2, 1-alfa/2), "%")
IC_boot
```

```
##      2.5%     97.5% 
## 0.4334742 1.1771924
```

<!-- 
Tabla resumen de los ics del tiempo de vida medio de microorganismos 
-->

### Estudio de simulación {#estudio-sim-exp}

El siguiente código permite realizar estudios de
simulación comparando las probabilidades de cobertura y las longitudes
de los intervalos de confianza clásicos (basados en normalidad),
bootstrap percentil, bootstrap percentil-$t$ y bootstrap percentil-$t$
simetrizado para la media, en el caso de muestras procedentes de una
distribución $\exp \left( \lambda \right)$. 
En este caso se obtienen las estimaciones Monte Carlo a partir de 500
simulaciones con $\lambda = 0.01$, tamaño muestral $n=100$ y $B=1000$ réplicas
bootstrap para un nivel de confianza nominal del 90\% ($\alpha =0.10$).



```r
t.ini <- proc.time()
rate <- 0.01
mu <- 1/rate
n <- 100

alfa <- 0.1
namesI <- paste0(100*c(alfa/2, 1-alfa/2), "%")

B <- 1000
percentil <- numeric(B)
percentilt <- numeric(B)
percentilts <- numeric(B)

nsim <- 500
resultados <- array(dim = c(nsim, 2, 4))
dimnames(resultados) <- list(NULL, c("Cobertura", "Longitud"),
        c("Normal", "Percentil", "Percentil-t", "Percentil-t simetrizado"))
# Bucle simulación
set.seed(1)
for (isim in 1:nsim) {
    # Aproximación clásica
    muestra <- rexp(n, rate = 0.01)
    media <- mean(muestra)
    desv <- sd(muestra)
    z <- qnorm(1 - alfa/2)
    ic_inf <- media - z*desv/sqrt(n)
    ic_sup <- media + z*desv/sqrt(n)
    I0 <- c(ic_inf, ic_sup)
    # names(I0) <- namesI
    resultados[isim, 1, 1] <- (I0[1] < mu) && (mu < I0[2])
    resultados[isim, 2, 1] <- I0[2] - I0[1]
    
    # Remuestreo bootstrap
    for (k in 1:B) {
        remuestra <- sample(muestra, n, replace = TRUE)
        percentil[k] <- sqrt(n) * (mean(remuestra) - media)
        percentilt[k] <- percentil[k]/sd(remuestra)
        percentilts[k] <- abs(percentilt[k])
    }
    
    # Aproximación bootstrap percentil
    pto_crit <- quantile(percentil, c(alfa/2, 1 - alfa/2))
    # Construcción del IC
    ic_inf_boot <- media - pto_crit[2]/sqrt(n)
    ic_sup_boot <- media - pto_crit[1]/sqrt(n)
    I1 <- c(ic_inf_boot, ic_sup_boot)
    # names(I1) <- namesI
    resultados[isim, 1, 2] <- (I1[1] < mu) && (mu < I1[2])
    resultados[isim, 2, 2] <- I1[2] - I1[1]
    
    # Aproximación bootstrap percentil-t
    pto_crit <- quantile(percentilt, c(alfa/2, 1 - alfa/2))
    # Construcción del IC
    ic_inf_boot <- media - pto_crit[2] * desv/sqrt(n)
    ic_sup_boot <- media - pto_crit[1] * desv/sqrt(n)
    I2 <- c(ic_inf_boot, ic_sup_boot)
    # names(I2) <- namesI
    resultados[isim, 1, 3] <- (I2[1] < mu) && (mu < I2[2])
    resultados[isim, 2, 3] <- I2[2] - I2[1]
    
    # Aproximación bootstrap percentil-t simetrizado
    pto_crit <- quantile(percentilts, 1 - alfa)
    # Construcción del IC
    ic_inf_boot <- media - pto_crit * desv/sqrt(n)
    ic_sup_boot <- media + pto_crit * desv/sqrt(n)
    I3 <- c(ic_inf_boot, ic_sup_boot)
    # names(I3) <- namesI
    resultados[isim, 1, 4] <- (I3[1] < mu) && (mu < I3[2])
    resultados[isim, 2, 4] <- I3[2] - I3[1]
}

t.fin <- proc.time() - t.ini
t.fin
```

```
##    user  system elapsed 
##   12.47    0.03   12.50
```

```r
apply(resultados, c(2, 3), mean)
```

```
##            Normal Percentil Percentil-t Percentil-t simetrizado
## Cobertura  0.8800   0.86600      0.8880                 0.88800
## Longitud  32.5022  32.13928     33.5653                33.49959
```

```r
# knitr::kable(t(apply(resultados, c(2, 3), mean)), digits = 3)
```

  Aproximación               Cobertura   Longitud
  ------------------------  ----------  ---------
  Normal                         0.892     32.243
  Percentil                      0.886     32.051
  Percentil-t                    0.912     33.395
  Percentil-t simetrizado        0.904     33.342

La siguiente tabla recoge las probabilidades de cobertura, estimadas por
Monte Carlo, en una ejecución con tamaño muestral $n=100$, $N=10000$
trials y $B=1000$ réplicas bootstrap para un nivel de confianza nominal
del 90\% ($\alpha =0.10$).

  Aproximación                        Cobertura IC  
  -------------------------------     ---------------------
  Normal                              88.60\%  
  Boot. percentil                     88.60\%
  Boot. percentil-$t$                 89.76\%
  Boot. percentil-$t$ simetrizado     89.46\%


<!-- 
n = 30, n = 100, t.fin
Example 5.12 (Exponential mean) 
-->

En la Sección \@ref(estudio-sim-boot) del Apéndice \@ref(intro-hpc) se incluye un estudio similar, empleando computación en paralelo para comparar las probabilidades de cobertura y las longitudes de los intervalos de confianza implementados en la función `boot.ci()`.


### IC bootstrap para el coeficiente de correlación {#icboot-trans}

Supongamos que queremos estudiar la correlación entre dos variables $X$ e $Y$ a partir del coeficiente de correlación lineal de Pearson:
$$\rho =\frac{ Cov \left( X, Y \right) }
{ \sigma \left( X \right) \sigma \left( Y \right) },$$
cuyo estimador natural es el coeficiente de correlación muestral:
$$r=\frac{\sum_{i=1}^{n}(x_i-\overline{x})(y_i-\overline{y})}
{\sqrt{ \sum_{i=1}^{n}(x_i-\overline{x})^{2}} 
\sqrt{\sum_{i=1}^{n}(y_i-\overline{y})^{2}}},$$
que podemos calcular en `R` empleando la función `cor()`.

Para realizar inferencias sobre el coeficiente de correlación, como aproximación más simple, se puede considerar que la distribución muestral de $r$ es aproximadamente normal y emplear el estadístico:

\begin{equation} 
\frac{r -\rho}{\sqrt{\frac{1 - r^2}{n - 2}}} \underset{aprox}{\sim } t_{n-2}
(\#eq:cor-t)
\end{equation} 

<!-- \mathcal{t}_{n-2} error en LaTeX-->

Pero esta aproximación solo sería válida en el caso de muestras grandes (o si la distribución bivariante de $(X, Y)$ es aproximadamente normal) cuando la correlación entre las variables es débil o moderada. 
En caso contrario la distribución muestral de $r$ puede ser muy asimétrica y los resultados obtenidos con el estadístico anterior no ser muy adecuados (esto concuerda con lo observado en la Sección \@ref(boot-unif-multi), al emplear bootstrap uniforme multidimensional para hacer inferencia sobre $R = r -\rho$).
Para evitar este problema se suelen obtener intervalos de confianza para $\rho$ empleando la transformación $Z$ de Fisher (1915):
$$Z = \frac{1}{2}\ln \left( \frac{1+r}{1-r} \right) = \operatorname{arctanh}(r),$$
que es una transformación (aprox.) normalizadora y estabilizadora de la varianza.
Suponiendo que $(X, Y)$ es normal bivariante y que hay independencia entre las observaciones:
$$Z \sim \mathcal{N}\left( \frac{1}{2}\ln \left( \frac{1+\rho}{1-\rho} \right), \frac{1}{n-3} \right).$$
El intervalo de confianza asintótico se obtiene empleando la aproximación normal tradicional en la escala $Z$ y aplicando posteriormente la transformación inversa:
$$r = \frac{\exp(2Z)-1}{\exp(2Z)+1} = \operatorname{tanh}(Z).$$

Esta aproximación está implementada en la función `cor.test()` del paquete base `stat` de R^[Se puede obtener el código tecleando en la consola `stats:::cor.test.default`.], además de que también realiza el contraste $H_0: \rho = 0$ empleando el estadístico \@ref(eq:cor-t).

Continuando con el ejemplo de la Sección \@ref(boot-unif-multi), para obtener un intervalo de confianza para el coeficiente de correlación lineal entre las variables `income` y `prestige` del conjunto de datos `Prestige`, podríamos emplear el siguiente código:


```r
data(Prestige, package="carData")
# with(Prestige, cor.test(income, prestige))
cor.test(Prestige$income, Prestige$prestige)
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  Prestige$income and Prestige$prestige
## t = 10.224, df = 100, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.6044711 0.7983807
## sample estimates:
##       cor 
## 0.7149057
```

La función `boot.ci()` del paquete `boot` también permite obtener intervalos de confianza calculados en una escala transfomada del estadístico,
mediante los parámetros:

- `h`: función vectorial que define la transformación. 
  Los intervalos se calculan en la escala de $h(t)$ y se aplica la función inversa (si se especifica) para transformarlos a la escala original.

- `hinv`: (opcional) función inversa de la transformación 
  (si no se especifica solo se calculan los intervalos en la escala transformada). 

- `hdot`: (opcional) función derivada de la transformación 
  (empleada por algunos métodos para aproximar la varianza en la escala transformada mediante el método delta).

Por ejemplo, para considerar la transformación $Z$ de Fisher en este caso, se podría emplear el siguiente código:

```r
library(boot)

statistic <- function(data, i){
  remuestra <- data[i, ]
  cor(remuestra$income, remuestra$prestige)
}

set.seed(1)
res.boot <- boot(Prestige, statistic, R = 1000)

h <- function(t) atanh(t)
hdot <- function(t) 1/(1 - t^2)
hinv <- function(t) tanh(t)

# boot.ci(res.boot, type = "norm", h = h)
# boot.ci(res.boot, type = "norm", h = h, hinv = hinv)
boot.ci(res.boot, type = "norm", h = h, hdot = hdot, hinv = hinv)
```

```
## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
## Based on 1000 bootstrap replicates
## 
## CALL : 
## boot.ci(boot.out = res.boot, type = "norm", h = h, hdot = hdot, 
##     hinv = hinv)
## 
## Intervals : 
## Level      Normal        
## 95%   ( 0.6016,  0.7858 )  
## Calculations on Transformed Scale;  Intervals on Original Scale
```

Esto sería en principo preferible a trabajar en la escala original, ya que la distribución bootstrap en la escala transformada se aproximaría más a la normalidad:


```r
ht <- h(res.boot$t)
hist(ht, freq = FALSE, breaks = "FD",
     main = "Distribución bootstrap en la escala transformada")
curve(dnorm(x, mean=mean(ht), sd=sd(ht)), lty = 2, add = TRUE)
```



\begin{center}\includegraphics[width=0.7\linewidth]{04-ic_boot_files/figure-latex/unnamed-chunk-10-1} \end{center}






