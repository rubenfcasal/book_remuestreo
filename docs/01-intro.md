# Motivación del principio Bootstrap  {#intro}




Etimología: bootstrap = cinta de la bota (oreja lateral para calzarse
las botas). Modismo anglosajón: to pull oneself up by one’s bootstraps.


## Introducción

El bootstrap es un procedimiento estadístico que sirve para aproximar la
distribución en el muestreo (normalmente) de un estadístico. Para ello
procede mediante remuestreo, es decir, obteniendo muestras mediante
algún procedimiento aleatorio que utilice la muestra original.

Su ventaja principal es que no requiere hipótesis sobre el mecanismo
generador de los datos. Sí las requiere, aunque suelen ser más
relajadas, para obtener propiedades asintóticas del mismo. Por otra
parte, su implementación en ordenador suele ser sencilla, en comparación
con otros métodos. Su principal inconveniente es la necesidad de
computación intensiva, debido a la fuerza bruta del método de Monte
Carlo. Con la capacidad computacional actual, esta mayor carga
computacional del bootstrap no suele ser un problema hoy en día. En
raras ocasiones el bootstrap no necesita del uso de técnicas de Monte
Carlo.

### Breve nota histórica

Precursores teóricos remotos:

-   Laplace (1810). Teoría límite de primer orden.

-   Chebychev (final siglo XIX). Teoría límite de segundo orden.

Primeras contribuciones:

-   Hubback (1878-1968). Esquemas de muestreo espacial para ensayos
    agrícolas.

-   Mahalanobis (años 1930 y segunda guerra mundial). Precursor del
    bootstrap por bloques.

Otras contribuciones:

-   Gurney, McCarthy, Hartigan (años 1960, 1970). Métodos de
    half-sampling para estimación de varianzas (U.S. Bureau of the
    Census).

-   Maritz, Jarret, Simon (años 1970, 1980). Métodos de permutaciones
    relacionados con el bootstrap.

En la actualidad:

-   Bradley Efron (Stanford University, 1979). Creador oficial del
    método. Acuñó su nombre. Fusionó la potencia de Monte Carlo con la
    resolución de problemas planteados de forma muy general.

-   Peter Hall (1951-2016). Fue uno de los estadísticos contemporáneos
    más prolíficos. Dedicó al bootstrap gran parte de su producción a
    partir de los años 1980.

### Paradigma inferencial y análogo bootstrap

***Paradigma inferencial***

Suponemos que $\mathbf{X}=\left( X_1,\ldots ,X_n \right)$ es una m.a.s. 
de una población con distribución $F$ y que 
estamos interesados en hacer inferencia sobre $\theta =\theta \left(F \right)$.
Para ello nos gustaría conocer la distribución en el muestreo de
$R\left( \mathbf{X},F \right)$, cierto estadístico función
de la muestra y de la distribución poblacional.
Por ejemplo: 
$$R=R\left( \mathbf{X},F \right) =\theta \left( F_n \right) 
-\theta \left( F \right) = \hat \theta - \theta,$$
siendo $F_n$ la función de distribución empírica.

A veces podemos calcular directamente la distribución de $R\left( \mathbf{X},F \right)$, 
aunque suele depender de cantidades poblacionales, 
no conocidas en la práctica.
Por ejemplo, bajo normalidad $X_i \overset{i.i.d.}{\sim} 
\mathcal{N}\left( \mu ,\sigma^2 \right)$, si estamos interesados en 
$$\theta \left( F \right) =\mu =\int x~dF\left( x \right) =\int xf\left( x \right) ~dx$$
como $\theta \left( F_n \right) = \int x~dF_n\left( x \right) = 
\sum \frac{1}{n}X_i = \bar{X}$, podríamos considerar el estadístico:
$$R=R\left( \mathbf{X},F \right) = \bar{X} - \mu \sim \mathcal{N}\left( 0 ,\frac{\sigma^2}{n} \right).$$
Aunque en la práctica la varianza no es normalmente conocida y habría que aproximarla
(sería preferible considerar como estadístico la media estudentizada).


Otras veces sólo podemos llegar a aproximar la distribución de 
$R\left( \mathbf{X},F \right)$ cuando $n \rightarrow \infty$.
Por ejemplo, cuando estamos interesados en la media pero desconocemos la 
distribución de los datos.

***Análogo bootstrap***

El primer paso es reemplazar la distribución poblacional (desconocida) $F$ por una
estimación, $\hat{F}$, de la misma. Por ejemplo, podríamos considerar la
distribución empírica $\hat{F}=F_n$ (bootstrap uniforme; Sección \@ref(intro-unif)), 
o una aproximación paramétrica $\hat{F}=F_{\hat \theta}$ (bootstrap paramétrico; Sección \@ref(modunif-boot-par)).

Como ejemplo ilustrativo consideramos los datos simulados [Figura \@ref(fig:muestra-sim)]:

```r
set.seed(1)
muestra <- rnorm(100)
hist(muestra, freq = FALSE, xlim = c(-3, 3),
     main = '', xlab = 'x', ylab = 'densidad')
curve(dnorm, lty = 2, add = TRUE)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{01-intro_files/figure-latex/muestra-sim-1} 

}

\caption{Distribución de la muestra simulada.}(\#fig:muestra-sim)
\end{figure}
Como aproximación de la distribución poblacional, desconocida en la práctica,
siempre podemos considerar la distribución empírica 
(o una versión suavizada: bootstrap suavizado; Sección \@ref(modunif-boot-suav)). 
Alternativamente podríamos asumir un modelo paramétrico y estimar los parámetros a partir de la muestra [Figura \@ref(fig:muestra-sim-aprox)].

```r
# Distribución bootstrap uniforme
curve(ecdf(muestra)(x), xlim = c(-3, 3), ylab = "F(x)", type = "s")
# Distribución bootstrap paramétrico (asumiendo normalidad)
curve(pnorm(x, mean(muestra), sd(muestra)), lty = 2, add = TRUE)
# Distribución teórica
curve(pnorm, lty = 3, add = TRUE)
legend("bottomright", legend = c("Empírica", "Aprox. paramétrica", "Teórica"), lty = 1:3)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{01-intro_files/figure-latex/muestra-sim-aprox-1} 

}

\caption{Distribución teórica de la muestra simulada y distintas aproximaciones.}(\#fig:muestra-sim-aprox)
\end{figure}


A partir de la aproximación $\hat{F}$ podríamos generar, condicionalmente a la muestra observada, 
remuestras 
$$\mathbf{X}^{\ast}=\left( X_1^{\ast},\ldots ,X_n^{\ast} \right)$$ 
con distribución $X_i^{\ast} \sim \hat{F}$, que demoninaremos remuestras bootstrap.
Por lo que podemos hablar de la distribución en el remuestreo de 
$$R^{\ast}=R\left( \mathbf{X}^{\ast},\hat{F} \right),$$ 
llamada distribución bootstrap.

La idea original (Efron, 1979) es que la distribución
de $\hat{\theta}_{b}^{\ast }$ en torno a $\hat{\theta}$ aproxima la
distribución de $\hat{\theta}$ en torno a $\theta$. 
Por tanto se pretende aproximar la distribución en el muestreo de $R$ por la
distribución bootstrap de $R^{\ast}$. 

En raras ocasiones la distribución bootstrap de $R^{\ast}$ es
calculable directamente, pero siempre suele poder aproximarse por
Monte Carlo.

<!-- 
Selección del estadístico
Interesaría un estadístico pivotal
Bootstrap percentil
Bootstrap básico o natural
Boostrap estudentizado
-->

### Implementación en la práctica {#intro-implementacion}

En el caso i.i.d., si empleamos como aproximación la distribución empírica $\hat{F}=F_n$,
la generación de las muestras bootstrap puede hacerse mediante remuestreo 
(manteniendo el tamaño muestral). Habría que simular una muestra de tamaño $n$ 
de una variable aleatoria discreta que toma los valores
$X_1,\ldots ,X_n$ todos ellos con probabilidad $\frac{1}{n}$:

* Para cada $i=1,\ldots, n$,
  $P^{\ast}\left( X_i^{\ast}=X_j \right) =
  \frac{1}{n}$, $j=1,\ldots ,n$.

Existen multitud de algoritmos para simular variables discretas, pero en
este caso de equiprobabilidad hay un procedimiento muy eficiente (método
de la transformación cuantil con búsqueda directa) que se reduce a
simular un número aleatorio $U$, con distribución $\mathcal{U}\left( 0,1 \right)$,
y hacer $X^{\ast}=X_{\left\lfloor nU\right\rfloor +1}$, donde $\left\lfloor x\right\rfloor$
representa la parte entera de $x$, es decir, el mayor número entero que
sea menor o igual que $x$. 
Empleando ese método, el procedimiento para generar la muestra bootstrap sería:

* Para cada $i=1,\ldots ,n$ 
  generar $U_i\sim \mathcal{U}\left( 0,1 \right)$ y
  hacer $X_i^{\ast}=X_{\left\lfloor nU_i\right\rfloor +1}$.
  

```r
set.seed(1)
n <- length(muestra)
u <- runif(n)
muestra_boot <- muestra[floor(n*u) + 1]
head(muestra_boot)
```

```
## [1] -0.1557955 -0.0593134 -1.0441346 -0.5425200  0.9189774  0.2670988
```


En `R` es recomendable^[De esta forma se evitan posibles problemas numéricos 
al emplear el método de la transformación cuantil cuando $n$ es extremadamente grande
(e.g. <https://stat.ethz.ch/pipermail/r-devel/2018-September/076817.html>).] 
emplear la función `sample` para generar muestras aleatorias con reemplazamiento 
del conjunto de datos original:

```r
muestra_boot <- sample(muestra, replace = TRUE)
head(muestra_boot)
```

```
## [1] -1.4707524  0.7685329  0.3876716 -0.6887557  0.9189774  1.3586796
```

En el caso multidimensional, cuando trabajamos con un conjunto de datos
con múltiples variables, 
podríamos emplear un procedimiento análogo, a partir de remuestras del
vector de índices. Por ejemplo:

```r
data(iris)
str(iris)
```

```
## 'data.frame':	150 obs. of  5 variables:
##  $ Sepal.Length: num  5.1 4.9 4.7 4.6 5 5.4 4.6 5 4.4 4.9 ...
##  $ Sepal.Width : num  3.5 3 3.2 3.1 3.6 3.9 3.4 3.4 2.9 3.1 ...
##  $ Petal.Length: num  1.4 1.4 1.3 1.5 1.4 1.7 1.4 1.5 1.4 1.5 ...
##  $ Petal.Width : num  0.2 0.2 0.2 0.2 0.2 0.4 0.3 0.2 0.2 0.1 ...
##  $ Species     : Factor w/ 3 levels "setosa","versicolor",..: 1 1 1 1 1 1 1 1 1 1 ...
```

```r
n <- nrow(iris)
# i_boot <- floor(n*runif(n)) + 1
# i_boot <- sample.int(n, replace = TRUE)
i_boot <- sample(n, replace = TRUE)
data_boot <- iris[i_boot, ]
str(data_boot)
```

```
## 'data.frame':	150 obs. of  5 variables:
##  $ Sepal.Length: num  6.9 7.2 4.3 6.4 5.7 5.8 7.2 5.6 5.8 5.4 ...
##  $ Sepal.Width : num  3.1 3.2 3 3.2 4.4 4 3 2.9 2.7 3.9 ...
##  $ Petal.Length: num  5.4 6 1.1 5.3 1.5 1.2 5.8 3.6 5.1 1.3 ...
##  $ Petal.Width : num  2.1 1.8 0.1 2.3 0.4 0.2 1.6 1.3 1.9 0.4 ...
##  $ Species     : Factor w/ 3 levels "setosa","versicolor",..: 3 3 1 3 1 1 3 2 3 1 ...
```

Esta forma de proceder es la que emplea por defecto el paquete `boot` que 
describiremos más adelante (Sección \@ref(intro-pkgboot)).

\BeginKnitrBlock{example}\iffalse{-91-73-110-102-101-114-101-110-99-105-97-32-115-111-98-114-101-32-108-97-32-109-101-100-105-97-32-99-111-110-32-118-97-114-105-97-110-122-97-32-99-111-110-111-99-105-100-97-93-}\fi{}<div class="example"><span class="example" id="exm:media-dt-conocida"><strong>(\#exm:media-dt-conocida)  \iffalse (Inferencia sobre la media con varianza conocida) \fi{} </strong></span>
<br> \vspace{0.5cm}

Hemos observado 15 tiempos de vida de microorganismos: 
0.143, 0.182, 0.256, 0.260, 0.270, 0.437, 0.509, 
0.611, 0.712, 1.04, 1.09, 1.15, 1.46, 1.88, 2.08.
A partir de los cuales queremos 
obtener una estimación por intervalo de confianza de su vida media,
suponiendo que la desviación típica es conocida e igual a 0.6
(en el Capítulo \@ref(icboot) se tratará con más detalle la construcción de intervalos de confianza).</div>\EndKnitrBlock{example}

```r
muestra <- c(0.143, 0.182, 0.256, 0.26, 0.27, 0.437, 0.509, 
    0.611, 0.712, 1.04, 1.09, 1.15, 1.46, 1.88, 2.08)
sigma <- 0.6
summary(muestra)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.1430  0.2650  0.6110  0.8053  1.1200  2.0800
```

```r
sd(muestra)
```

```
## [1] 0.6237042
```

```r
hist(muestra)
rug(muestra)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{01-intro_files/figure-latex/microorganismos-1} 

}

\caption{Distribución del tiempo de vida de microorganismos.}(\#fig:microorganismos)
\end{figure}
[Figura \@ref(fig:microorganismos)]

***Contexto clásico***

Suponemos que los datos $\mathbf{X}=\left( X_1,\ldots ,X_n \right)$ son una m.a.s. 
de una población con distribución $F$, con $\mu$ desconocida y $\sigma$ conocida,
y que estamos interesados en hacer inferencia sobre:
$$\theta \left( F \right) =\mu =\int x~dF\left( x \right)$$
Para ello, un estadístico adecuado para este caso es:
$$R=R\left( \mathbf{X},F \right) =\sqrt{n}\frac{\bar{X}-\mu }{\sigma},$$
con $\theta \left( F_n \right) =\int x~dF_n\left( x \right) = \bar{X}$.

Bajo normalidad $\left( X\sim \mathcal{N}\left( \mu ,\sigma^2 \right) \right)$,
$R\sim N\left( 0,1 \right)$. Si $F$ no es normal, tan sólo sabemos que,
bajo ciertas condiciones,
$R\overset{d}{\rightarrow }\mathcal{N}\left( 0, 1 \right)$.

A partir de esta última aproximación, se obtiene el intervalo de
confianza asintótico (de nivel $1-\alpha$) para la media $\mu$:
$$\hat{IC}_{1-\alpha}\left(  \mu\right)  = 
\left(  \overline{X}-z_{1-\alpha/2}\dfrac{\sigma}{\sqrt{n}},\ \overline{X} 
+ z_{1-\alpha/2}\dfrac{\sigma}{\sqrt{n}} \right).$$


```r
alfa <- 0.05
x_barra <- mean(muestra)
z <- qnorm(1 - alfa/2)
n <- length(muestra)
ic_inf <- x_barra - z*sigma/sqrt(n)
ic_sup <- x_barra + z*sigma/sqrt(n)
IC <- c(ic_inf, ic_sup)
IC
```

```
## [1] 0.501697 1.108970
```

***Contexto bootstrap*** 

Consideramos la función de distribución empírica $\hat{F}=F_n$
como aproximación de la distribución poblacional (bootstrap uniforme).
Para aproximar la distribución bootstrap del estadístico por Monte Carlo,
se generan $B=1000$ muestras bootstrap 
$\mathbf{X}^{\ast (b)}=\left( X_1^{\ast (b)},\ldots ,X_n^{\ast (b)} \right)$
de forma que 
$P^{\ast (b)}\left( X_i^{\ast}=X_j \right) =
\frac{1}{n}$, $j=1,\ldots ,n$, para $i=1,\ldots, n$ y $b=1,\ldots, B$.
A partir de las cuales se obtienen las $B$ réplicas bootstrap del estadístico:
$$R^{\ast (b)}=R\left( \mathbf{X}^{\ast (b)},\hat{F} \right) =\sqrt{n}\frac{
\bar{X}^{\ast  (b)}-\bar{X}}{\sigma }, \ b=1,\ldots, B, $$
con $\bar{X}^{\ast (b)} = \frac{1}{n}\sum X_i^{\ast (b)}$.


```r
set.seed(1)
B <- 1000
estadistico_boot <- numeric(B)
for (k in 1:B) {
    remuestra <- sample(muestra, n, replace = TRUE)
    x_barra_boot <- mean(remuestra)
    estadistico_boot[k] <- sqrt(n) * (x_barra_boot - x_barra)/sigma
}
```

Las características de interés de la distribución en el muestreo de $R$ 
se aproximan por las correspondientes de la distribución bootstrap de $R^{\ast}$.
En este caso nos interesa aproximar los puntos críticos $x_{\alpha /2}$ y
$x_{1-\alpha /2}$, tales que:
$$P\left( x_{\alpha /2} < R < x_{1-\alpha /2} \right) = 1-\alpha.$$
Para lo que podemos emplear los cuantiles muestrales^[
Se podrían considerar distintos estimadores del cuantil $x_{\alpha}$ 
(ver p.e. la ayuda de la función `quantile()`).
Si empleamos directamente la distribución empírica, el cuantil se 
correspondería con la observación ordenada en la posición $B \alpha$ 
(se suele hacer una interpolación lineal si este valor no es entero), 
lo que equivale a emplear la función `quantile()` de `R` con el parámetro 
`type = 1`. Esta función considera por defecto la posición 
$1 + (B - 1) \alpha$ (`type = 7`).
En el libro de Davison y Hinkley (1997), y en el paquete `boot`, se emplea $(B + 1) \alpha$ (equivalente a `type = 6`; lo que justifica que
consideren habitualmente 99, 199 ó 999 réplicas bootstrap).]:



```r
# Empleando la distribución empírica del estadístico bootstrap: 
estadistico_boot_ordenado <- sort(estadistico_boot)
indice_inf <- floor(B * alfa/2)
indice_sup <- floor(B * (1 - alfa/2))
pto_crit <- estadistico_boot_ordenado[c(indice_inf, indice_sup)]
# Empleando la función `quantile`:
# pto_crit <- quantile(estadistico_boot, c(alfa/2, 1 - alfa/2), type = 1)
pto_crit <- quantile(estadistico_boot, c(alfa/2, 1 - alfa/2))
pto_crit
```

```
##      2.5%     97.5% 
## -1.918622  2.075984
```
A partir de los cuales obtenemos la correspondiente estimación por IC
boostrap:
$$\hat{IC}^{boot}_{1-\alpha}\left(  \mu\right)  = 
\left(  \overline{X}-x_{1-\alpha/2}\dfrac{\sigma}{\sqrt{n}},\ \overline{X} 
- x_{\alpha/2}\dfrac{\sigma}{\sqrt{n}} \right).$$


```r
# Construcción del IC
ic_inf_boot <- x_barra - pto_crit[2] * sigma/sqrt(n)
ic_sup_boot <- x_barra - pto_crit[1] * sigma/sqrt(n)
IC_boot <- c(ic_inf_boot, ic_sup_boot)
names(IC_boot) <- paste0(100*c(alfa/2, 1-alfa/2), "%") # rev(names(IC_boot))
IC_boot
```

```
##      2.5%     97.5% 
## 0.4837233 1.1025650
```
Nótese que este intervalo de confianza no está centrado en la media,
al contrario que el obtenido con la aproximación tradicional.
Aunque en este caso no se observan grandes diferencias ya que 
la distribución bootstrap obtenida es muy similar a la aproximación normal
(ver Figura \@ref(fig:estad-boot)).


```r
hist(estadistico_boot, freq = FALSE)
lines(density(estadistico_boot))
abline(v = pto_crit)
curve(dnorm, lty = 2, add = TRUE)
abline(v = c(-z, z), lty = 2)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{01-intro_files/figure-latex/estad-boot-1} 

}

\caption{Distribución del estadístico boostrap y aproximaciones de los cuantiles. Con línea discontinua se muestra la distribución normal asintótica.}(\#fig:estad-boot)
\end{figure}


## El Bootstrap uniforme {#intro-unif}

Como ya se comentó anteriormente el bootstrap uniforme es aquel en el que
se reemplaza la distribución poblacional (desconocida) por la distribución
empírica:
$$F_n\left( x \right) =\frac{1}{n}\sum_{i=1}^{n}\mathbf{1}\left\{ X_i\leq x\right\}.$$

Es decir $\hat{F}=F_n$ y, por lo tanto,
$R^{\ast}=R\left( \mathbf{X}^{\ast},F_n \right)$. 

Conviene recordar algunas propiedades de la distribución empírica:
$$\begin{aligned}
nF_n\left( x \right) &= \sum_{i=1}^{n}\mathbf{1}\left\{ X_i\leq x\right\}
\sim \mathcal{B}\left( n,F\left( x \right) \right), \\
E\left( nF_n\left( x \right) \right) &= nF\left( x \right) \implies E\left(
F_n\left( x \right) \right) =F\left( x \right), \\
Var\left( nF_n\left( x \right) \right) &=  nF\left( x \right) \left(
1-F\left( x \right) \right) \\
&\implies  Var\left( F_n\left( x \right) \right) =\frac{F\left( x \right) \left( 1-F\left( x \right) \right)}{n}
\end{aligned}$$

Así pues, en este caso el algoritmo bootstrap uniforme (también llamado
bootstrap naïve) es el siguiente:

1. Para cada $i=1,\ldots ,n$ arrojar $X_i^{\ast}$ a partir de
$F_n$, es decir
$P^{\ast}\left( X_i^{\ast}=X_j \right) =\frac{1}{n}$, $j=1,\ldots
,n$

2. Obtener $\mathbf{X}^{\ast}=\left( X_1^{\ast},\ldots
,X_n^{\ast} \right)$

3. Calcular $R^{\ast}=R\left( \mathbf{X}^{\ast},F_n \right)$

Como veremos más adelante, a veces (muy poco frecuentemente) es posible
calcular exactamente la distribución bootstrap de $R^{\ast}$. Cuando
eso no es posible, esa distribución es fácilmente aproximable por Monte
Carlo, arrojando una gran cantidad, $B$, de réplicas de $R^{\ast}$. En
ese caso, el algoritmo se convierte en:

1. Para cada $i=1,\ldots ,n$ arrojar $X_i^{\ast}$ a partir de $F_n$

2. Obtener $\mathbf{X}^{\ast}=\left( X_1^{\ast},\ldots
,X_n^{\ast} \right)$

3. Calcular $R^{\ast}=R\left( \mathbf{X}^{\ast},F_n \right)$

4. Repetir $B$ veces los pasos 1-3 para obtener las réplicas bootstrap
$R^{\ast (1)}, \ldots, R^{\ast (B)}$

5. Utilizar esas réplicas bootstrap para aproximar la distribución en el
muestreo de $R$


### Ejemplos

\BeginKnitrBlock{example}\iffalse{-91-73-110-102-101-114-101-110-99-105-97-32-115-111-98-114-101-32-108-97-32-109-101-100-105-97-32-99-111-110-32-118-97-114-105-97-110-122-97-32-99-111-110-111-99-105-100-97-44-32-99-111-110-116-105-110-117-97-99-105-243-110-93-}\fi{}<div class="example"><span class="example" id="exm:media-dt-conocida-perturbando"><strong>(\#exm:media-dt-conocida-perturbando)  \iffalse (Inferencia sobre la media con varianza conocida, continuación) \fi{} </strong></span></div>\EndKnitrBlock{example}

En el Ejemplo \@ref(exm:media-dt-conocida) anteriormente visto de inferencia para la media con
varianza conocida, el algoritmo bootstrap (basado en Monte Carlo) para
aproximar la distribución en el muestreo de $R$ empleado fue:

1. Para cada $i=1,\ldots ,n$ arrojar $U_i\sim \mathcal{U}\left( 0,1 \right)$ y
hacer $X_i^{\ast}=X_{\left\lfloor nU_i\right\rfloor +1}$

2. Obtener $\bar{X}^{\ast}=\frac{1}{n}\sum X_i^{\ast}$

3. Calcular
$R^{\ast}=\sqrt{n}\frac{\bar{X}^{\ast}-\bar{X}}{
\sigma }$

4. Repetir $B$ veces los pasos 1-3 para obtener las réplicas bootstrap
$R^{\ast (1)}, \ldots, R^{\ast (B)}$

5. Aproximar la distribución en el muestreo de $R$ mediante la empírica
de $R^{\ast (1)}, \ldots, R^{\ast (B)}$

Como curiosidad podemos calcular la esperanza y la varianza de $R$ y la
esperanza y varianza bootstrap de $R^{\ast}$. Para $R$ tenemos:
$$\begin{aligned}
E\left( R \right) &=\sqrt{n}\frac{E\left( \bar{X} \right) -\mu }{\sigma }
=0, \\
Var\left( R \right) &=n\frac{Var\left( \bar{X} \right)}{\sigma^2}=n
\frac{\frac{1}{n}\sigma^2}{\sigma^2}=1.
\end{aligned}$$

Para calcular esos mismos momentos de $R^{\ast}$, resultará útil
obtener previamente la esperanza y varianza bootstrap de
$\bar{X}^{\ast}$:
$$\begin{aligned}
E^{\ast}\left( \bar{X}^{\ast} \right) &= \frac{1}{n}
\sum_{i=1}^{n}E^{\ast}\left( X_i^{\ast} \right) =\frac{1}{n}
\sum_{i=1}^{n}E^{\ast}\left( X_1^{\ast} \right) =E^{\ast}\left(
X_1^{\ast} \right) =\bar{X}, \\
Var^{\ast}\left( \bar{X}^{\ast} \right) &= \frac{1}{n^2}
\sum_{i=1}^{n}Var^{\ast}\left( X_i^{\ast} \right) =\frac{1}{n^2}
\sum_{i=1}^{n}Var^{\ast}\left( X_1^{\ast} \right) =\frac{1}{n}Var^{\ast
}\left( X_1^{\ast} \right) =\frac{S_n^2}{n},
\end{aligned}$$
ya que
$$\begin{aligned}
E^{\ast}\left( X_1^{\ast} \right) &= \sum_{j=1}^{n}X_jP^{\ast}\left(
X_1^{\ast}=X_j \right) =\sum_{j=1}^{n}\frac{1}{n}X_j=\bar{X}, \\
Var^{\ast}\left( X_1^{\ast} \right) &= E^{\ast}\left( X_1^{\ast
2} \right) -\left[ E^{\ast}\left( X_1^{\ast} \right) \right]
^2=\sum_{j=1}^{n}X_j^2P^{\ast}\left( X_1^{\ast}=X_j \right) -\bar{X}
^2 \\
&= \frac{1}{n}\sum_{j=1}^{n}X_j^2-\bar{X}^2=\frac{1}{n}
\sum_{j=1}^{n}\left( X_j-\bar{X} \right)^2=S_n^2
\end{aligned}$$

Así pues, la esperanza y la varianza bootstrap de $R^{\ast}$
resultan:
$$\begin{aligned}
E^{\ast}\left( R^{\ast} \right) &= \sqrt{n}\frac{E^{\ast}\left( \bar{X}^{\ast} \right) -\bar{X}}{\sigma }=0, \\
Var^{\ast}\left( R^{\ast} \right) &= n\frac{Var^{\ast}\left( \bar{X}^{\ast} \right)}{\sigma^2}=n\frac{\frac{1}{n}S_n^2}{\sigma^2}=
\frac{S_n^2}{\sigma^2}.
\end{aligned}$$

Es curioso observar que la esperanza de $R$ y la esperanza bootstrap de
$R^{\ast}$ coinciden (son ambas cero), pero no ocurre lo mismo con sus
varianzas: la de $R$ es $1$ y la varianza bootstrap de $R^{\ast}$ es
$S_n^2/\sigma^2$, que, aunque tiende a $1$ (en probabilidad o de
forma casi segura, bajo las condiciones adecuadas) cuando $n\rightarrow
\infty$, no es igual a $1$. Eso nos lleva a intuir que el método de
remuestreo bootstrap propuesto quizá podría modificarse ligeramente para
que imitase exactamente al caso no bootstrap también en la varianza.
Puede comprobarse que eso se consigue remuestreando $X^{\ast}$ de la
distribución empírica de la muestra modificada: 
$\left( \tilde{X}_1,\ldots ,\tilde{X}_n \right)$, siendo
$$\tilde{X}_i=\bar{X}+\frac{\sigma }{S_n}\left( X_i-\bar{X}
 \right) \text{, }i=1,\ldots ,n.$$

Efectivamente, bajo ese nuevo remuestreo, se tiene 
$$\begin{aligned}
E^{\ast}\left( \bar{X}^{\ast} \right) &= E^{\ast}\left( X_1^{\ast
} \right) =\overline{\tilde{X}}=\frac{1}{n}\sum_{i=1}^{n}\left[ \bar{X}+
\frac{\sigma }{S_n}\left( X_i-\bar{X} \right) \right] \\
&= \bar{X}+\frac{1}{n}\frac{\sigma }{S_n}\sum_{i=1}^{n}\left( X_i-
\bar{X} \right) =\bar{X}, \\
Var^{\ast}\left( \bar{X}^{\ast} \right) &= \frac{1}{n}Var^{\ast
}\left( X_1^{\ast} \right) =\frac{\sigma^2}{n},
\end{aligned}$$
ya que
$$\begin{aligned}
Var^{\ast}\left( X_1^{\ast} \right) &= E^{\ast}\left( X_1^{\ast
2} \right) -\left[ E^{\ast}\left( X_1^{\ast} \right) \right]
^2=\sum_{j=1}^{n}\tilde{X}_j^2P^{\ast}\left( X_1^{\ast}=\tilde{X}
_j \right) -\overline{\tilde{X}}^2 \\
&= \frac{1}{n}\sum_{j=1}^{n}\tilde{X}_j^2-\overline{\tilde{X}}^2=\frac{
1}{n}\sum_{j=1}^{n}\left( \tilde{X}_j-\overline{\tilde{X}} \right)^2=
\frac{1}{n}\sum_{j=1}^{n}\left[ \frac{\sigma }{S_n}\left( X_j-\bar{X} \right) \right]^2 \\
&= \frac{\sigma^2}{S_n^2}\frac{1}{n}\sum_{j=1}^{n}\left( X_i-
\bar{X} \right)^2=\frac{\sigma^2}{S_n^2}S_n^2=\sigma^2.
\end{aligned}$$
Como consecuencia 
$$\begin{aligned}
E^{\ast}\left( R^{\ast} \right) &= \sqrt{n}\frac{E^{\ast}\left( 
\bar{X}^{\ast} \right) -\bar{X}}{\sigma }=0, \\
Var^{\ast}\left( R^{\ast} \right) &= n\frac{Var^{\ast}\left( 
\bar{X}^{\ast} \right)}{\sigma^2}=n\frac{\frac{\sigma^2}{n}}{\sigma^2}
=1.
\end{aligned}$$

Esto es muy coherente con lo que nos diría la intuición pues, si la
varianza poblacional, $\sigma^2$, es conocida (ese es el motivo de
que podamos usarla directamente en la definición del estadístico $R$),
el plan de remuestreo bootstrap también ha de conocer $\sigma^2$, es
decir ha de diseñarse de modo que la distribución bootstrap de
$X^{\ast}$ tenga también varianza bootstrap $\sigma^2$. Eso ocurre
con el remuestreo uniforme de la muestra transformada
$\left( \tilde{X}_1,\ldots ,\tilde{X}_n \right)$, pero no ocurre con
el remuestreo naïve (a partir de la distribución empírica de la muestra
original). Esto da pie a una de las consideraciones más importantes a la
hora de diseñar un buen método de remuestreo bootstrap: ha de procurarse
que **el bootstrap imite todas las condiciones que cumple la población
original**.

El código para realizar remuestreo bootstrap uniforme sobre la empírica de la
muestra perturbando es análogo:

```r
# Remuestreo
B <- 1000
estadistico_boot <- numeric(B)
coeficiente <- sigma/sd(muestra)
muestra_perturbada <- x_barra + coeficiente * (muestra - x_barra)
for (k in 1:B) {
  remuestra <- sample(muestra_perturbada, n, replace = TRUE)
  x_barra_boot <- mean(remuestra)
  estadistico_boot[k] <- sqrt(n) * (x_barra_boot - x_barra)/sigma
}

# Aproximación bootstrap de los ptos críticos
pto_crit <- quantile(estadistico_boot, c(alfa/2, 1 - alfa/2))
# Construcción del IC
ic_inf_boot <- x_barra - pto_crit[2] * sigma/sqrt(n)
ic_sup_boot <- x_barra - pto_crit[1] * sigma/sqrt(n)
IC_boot <- c(ic_inf_boot, ic_sup_boot)
names(IC_boot) <- paste0(100*c(alfa/2, 1-alfa/2), "%")
IC_boot
```

```
##      2.5%     97.5% 
## 0.5024398 1.0827052
```



\BeginKnitrBlock{example}\iffalse{-91-73-110-102-101-114-101-110-99-105-97-32-115-111-98-114-101-32-108-97-32-109-101-100-105-97-110-97-93-}\fi{}<div class="example"><span class="example" id="exm:mediana"><strong>(\#exm:mediana)  \iffalse (Inferencia sobre la mediana) \fi{} </strong></span><br> \vspace{0.5cm}

Continuando con el ejemplo de los tiempos de vida de microorganismos,
supongamos que queremos obtener una estimación por intervalo de confianza 
de su vida mediana a partir de los 15 valores observados.</div>\EndKnitrBlock{example}

Consideramos la mediana poblacional como parámetro de interés:
$$\theta = \theta \left( F \right) = F^{-1}\left( \frac{1}{2} \right) 
= \inf \left\{ x\in \mathbb{R} : F\left( x \right) \geq \frac{1}{2}\right\}.$$
Dada una muestra $\mathbf{X}=\left( X_1,\ldots ,X_n \right) \sim
F$, $\theta$ puede estimarse mediante la mediana muestral
$$\begin{aligned}
\hat{\theta} &= \theta \left( F_n \right) =F_n^{-1}\left( \frac{1}{2} \right) 
=\inf \left\{ x\in \mathbb{R} : F_n\left( x \right) \geq \frac{1}{2}
\right\} \\
&= \left\{ 
\begin{array}{ll}
X_{(m)} & \text{si } n=2m-1 \text{ es impar} \\ 
\frac{X_{(m)}+X_{\left( m+1 \right)}}{2} & \text{si } n=2m \text{ es par}
\end{array}
\right.
\end{aligned}$$
siendo $X_{(1)},\ldots ,X_{(n)}$ los estadísticos ordenados.

El estadístico interesante para realizar inferencia en este contexto es
$R=\sqrt{n}\left( \hat{\theta}-\theta \right)$. Si la población de
partida es continua, puede demostrarse que su distribución asintótica
(i.e., cuando $n \rightarrow \infty$) viene dada por
$$R=\sqrt{n}\left( \hat{\theta}-\theta \right) \overset{d}{\rightarrow }
\mathcal{N}\left( 0,\frac{1}{f\left( \theta \right)^2} \right),$$donde $f$ es
la función de densidad de la población. Como consecuencia, la
utilización de esta distribución límite, 
$\mathcal{N}\left( 0, 1/f\left( \theta \right)^2 \right)$, para realizar
inferencia sobre la mediana, además de comportar una aproximación de la
distribución en el muestreo real, no puede utilizarse directamente
porque la densidad (desconocida) aparece en la expresión de la varianza
asintótica. Para ser utilizable en la práctica deberíamos estimar $f$,
lo cual es un problema añadido.

Esta es pues una situación muy natural en la que usar un método
bootstrap para aproximar la distribución de $R$. Consideremos como
estimador de $F$ la distribución empírica, $F_n$, y procedamos según
un bootstrap uniforme (supongamos $n=2m-1$, impar, por simplicidad):

1. Para cada $i=1,\ldots ,n$ arrojar $U_i\sim \mathcal{U}\left( 0,1 \right)$ y
hacer $X_i^{\ast}=X_{\left\lfloor nU_i\right\rfloor +1}$

2. Obtener $X_{(1)}^{\ast},\ldots ,X_{(n)}^{\ast}$ 
los estadísticos ordenados de la remuestra bootstrap y quedarse con
el que ocupa lugar central:
$\hat{\theta}^{\ast}=\theta \left( F_n^{\ast} \right) =X_{(m)}^{\ast}$

3. Calcular
$R^{\ast}=\sqrt{n}\left( X_{(m)}^{\ast}-X_{\left(m \right)} \right)$

4. Repetir $B$ veces los pasos 1-3 para obtener las réplicas bootstrap
$R^{\ast (1)}, \ldots, R^{\ast (B)}$

5. Aproximar la distribución en el muestreo de $R$ mediante la empírica
de $R^{\ast (1)}, \ldots, R^{\ast (B)}$

El código implementando este algoritmo sería muy similar al de los casos anteriores:

```r
x_mediana<- median(muestra)

# Remuestreo
B <- 1000
estadistico_boot <- numeric(B)
coeficiente <- sigma/sd(muestra)
for (k in 1:B) {
  remuestra <- sample(muestra, n, replace = TRUE)
  x_mediana_boot <- median(remuestra)
  estadistico_boot[k] <- sqrt(n) * (x_mediana_boot - x_mediana)
}

# Aproximación bootstrap de los ptos críticos
pto_crit <- quantile(estadistico_boot, c(alfa/2, 1 - alfa/2))
# Construcción del IC
ic_inf_boot <- x_mediana - pto_crit[2]/sqrt(n)
ic_sup_boot <- x_mediana - pto_crit[1]/sqrt(n)
IC_boot <- c(ic_inf_boot, ic_sup_boot)
names(IC_boot) <- paste0(100*c(alfa/2, 1-alfa/2), "%")
IC_boot
```

```
##  2.5% 97.5% 
## 0.132 0.962
```

Sin embargo, como veremos más adelante, este caso de inferencia de la
mediana es uno de los pocos casos en los que la distribución bootstrap
se puede calcular de forma exacta, siendo dicha expresión utilizable en
la práctica.

## Cálculo de la distribución Bootstrap: exacta y aproximada 

### Distribución bootstrap exacta

En principio siempre es posible calcular la distribución en el
remuestreo del estadístico bootstrap de forma exacta. Al menos para el
bootstrap uniforme, que es el más habitual. El motivo es que la
distribución de probabilidad de la que se remuestrea en el universo
bootstrap es discreta y con un número finito de valores: $X_1,\ldots
,X_n$. Así pues, cada observación bootstrap, $X_i^{\ast}$, ha de
tomar necesariamente alguno de esos $n$ valores y, por tanto, el número
de posibles remuestras, $\mathbf{X}^{\ast}=\left( X_1^{\ast
},\ldots ,X_n^{\ast} \right)$, obtenibles mediante el bootstrap
uniforme es finito, concretamente $n^{n}$. Aún siendo finito, este
número es gigantescamente grande incluso para tamaños muestrales
pequeños (salvo casos extremos del tipo $n=2,\ldots ,9$). Por ejemplo,
para $n=10$, tenemos $10^{10}$ (diez mil millones de) posibles
remuestras bootstrap y para $n=20$, tendríamos
$20^{20}\simeq 10.4857\cdot 10^{25}$ (algo más de cien
cuatrillones). Incluso para estos tamaños muestrales el problema de
cálculo de la distribución bootstrap exacta de
$\mathbf{X}^{\ast}$ es inabordable.

### Vectores de remuestreo

Una forma alternativa de representar las posibles remuestras bootstrap
es mediante los llamados vectores de remuestreo. Son utilizables en el
caso de que el estadístico de interés sea funcional, es decir, cuando
$R$ depende de la muestra sólo a través de la distribución empírica o,
lo que es lo mismo, el valor de $R$ no cambia cuando realizamos una
permutación arbitraria sobre los elementos de la muestra (los cambiamos
de orden). Consideremos la remuestra bootstrap
$\mathbf{X}^{\ast}
=\left( X_1^{\ast},\ldots ,X_n^{\ast} \right)$ y denotemos por
$$N_j=\#\left\{ i\in \left\{ 1,\ldots ,n\right\} : 
X_i^{\ast}=X_j\right\}.$$
Obviamente, si el orden en el que se han obtenido los
elementos de la muestra no es importante, entonces el vector 
$\mathbf{N}=\left( N_1,\ldots ,N_n \right)$ contiene la misma
información que la remuestra bootstrap $\mathbf{X}^{\ast}$.
Esencialmente lo que hace el vector $\mathbf{N}$ es contabilizar
cuantas veces se repite cada elemento de la muestra original en la
remuestra bootstrap. Con esta notación, el vector de remuestreo
bootstrap, 
$\mathbf{P}^{\ast}=\left( P_1^{\ast},\ldots ,P_n^{\ast} \right)$, 
se define como $P_i^{\ast}=\frac{N_i}{n}$, $i=1,\ldots ,n$.

La distribución en el remuestreo de $\mathbf{N}$, bajo el
bootstrap uniforme, es multinomial: $\mathbf{N}\sim \mathcal{M}_n\left(
n,\left( \frac{1}{n},\ldots ,\frac{1}{n} \right) \right)$. Así que su
masa de probabilidad, y por tanto la de $\mathbf{P}^{\ast}$, es
fácilmente calculable:
$$\begin{aligned}
P\left( N_1=m_1,\ldots ,N_n=m_n \right) &= \frac{n!}{m_1!\cdots
m_n!}\left( \frac{1}{n} \right)^{m_1}\cdots \left( \frac{1}{n} \right)
^{m_n} \\
&= \frac{n!}{m_1!\cdots m_n!n^{n}}\text{, } \\
\end{aligned}$$
para $m_1,\ldots ,m_n$ enteros con $\sum_{i=1}^{n}m_i = n$,
donde el número de átomos de probabilidad de $\mathbf{N}$ es ahora
$\binom{n+n-1}{n}=\binom{2n-1}{n}$. 

En general $\binom{2n-1}{n}<n^{n}$,
pues el hecho de que no importe el orden de las componentes de las
remuestras bootstrap provoca un menor número de átomos de probabilidad.
Aún así dicho cardinal es prohibitivamente grande incluso para tamaños
muestrales pequeños: para $n=10$, resultaría abordable pues
$\binom{19}{10}= 92\,378$, pero para $n=20$ tendríamos
$\binom{39}{20}= 68\, 923\,264\,410$. De toda esa
enorme cantidad de átomos, el de más grande probabilidad resulta tener
una probabilidad de $\frac{n!}{n^{n}}$, que es insignificantemente
pequeña para tamaños pequeños como $n=20$, con
$\frac{20!}{20^{20}}\sim 2.\, 320\,2\times 10^{-8}$. De todas
formas, existen raras ocasiones en las que el número de átomos de
probabilidad de $R^{\ast}$ resulta ser mucho menor que el de
$\mathbf{P}^{\ast}$.

Para tamaños muestrales realmente pequeños es posible encontrar todos
los átomos de probabilidad de la distribución bootstrap. Un ejemplo es
la media muestral con, por ejemplo, $n=3$.


\BeginKnitrBlock{example}\iffalse{-91-77-101-100-105-97-32-109-117-101-115-116-114-97-108-32-112-97-114-97-32-117-110-97-32-109-117-101-115-116-114-97-32-100-101-32-116-97-109-97-241-111-32-51-93-}\fi{}<div class="example"><span class="example" id="exm:media3"><strong>(\#exm:media3)  \iffalse (Media muestral para una muestra de tamaño 3) \fi{} </strong></span><br> \vspace{0.5cm}

Consideremos una muestra aleatoria simple de tamaño $n=3$ de una
población con distribución $F$ y tomemos como parámetro de interés la
media poblacional
$\theta \left( F \right) =\mu =\int xdF\left( x \right)$. 
Tomemos como estadístico de interés 
$R=R\left( \mathbf{X},F \right) =\bar{X}$. 
El análogo bootstrap de esta estadístico es
$R^{\ast}=R\left( \mathbf{X}^{\ast},F_n \right) =\bar{X}^{\ast}$, 
cuya distribución en el remuestreo se puede calcular de forma exacta 
debido al reducido número de átomos de probabilidad que tiene. 
Esta es una distribución discreta con 10 posibles valores, cuyo
valor más probable es precisamente $\bar{X}$ que tiene una
probabilidad bootstrap de $\frac{2}{9}$, como puede verse en la
siguiente tabla.

| $\mathbf{X}^{\ast}$ (salvo permutaciones)  |    $\mathbf{N}$ = $\left( m_1,m_2,m_3 \right)$  |    $\mathbf{P}^{\ast}$ = $\left(p_1,p_2,p_3 \right)$   |    $\frac{3!}{m_1!m_2!m_3!3^{3}}$   |   $\bar{X}^{\ast}$ | 
| -----------------------------  | -------------------------- | -------------------------- | --------------------------- | ----------------------------- |
| $\left( X_1,X_1,X_1 \right)$   |    $\left( 3,0,0 \right)$   |    $\left(1,0,0 \right)$   |    $\frac{1}{27}$   |    $X_1$    | 
| $\left( X_2,X_2,X_2 \right)$   |    $\left( 0,3,0 \right)$   |    $\left(0,1,0 \right)$   |    $\frac{1}{27}$   |    $X_2$    | 
| $\left( X_3,X_3,X_3 \right)$   |    $\left( 0,0,3 \right)$   |    $\left(0,0,1 \right)$   |    $\frac{1}{27}$   |    $X_3$    | 
| $\left( X_1,X_1,X_2 \right)$   |    $\left( 2,1,0 \right)$   |    $\left( \frac{2}{3},\frac{1}{3},0 \right)$   |    $\frac{1}{9}$   |    $\frac{2X_1+X_2}{3}$    | 
| $\left( X_1,X_1,X_3 \right)$   |    $\left( 2,0,1 \right)$   |    $\left( \frac{2}{3},0,\frac{1}{3} \right)$   |    $\frac{1}{9}$   |    $\frac{2X_1+X_3}{3}$    | 
| $\left( X_1,X_2,X_2 \right)$   |    $\left( 1,2,0 \right)$   |    $\left( \frac{1}{3},\frac{2}{3},0 \right)$   |    $\frac{1}{9}$   |    $\frac{X_1+2X_2}{3}$    | 
| $\left( X_2,X_2,X_3 \right)$   |    $\left( 0,2,1 \right)$   |    $\left( 0,\frac{2}{3},\frac{1}{3} \right)$   |    $\frac{1}{9}$   |    $\frac{2X_2+X_3}{3}$    | 
| $\left( X_1,X_3,X_3 \right)$   |    $\left( 1,0,2 \right)$   |    $\left( \frac{1}{3},0,\frac{2}{3} \right)$   |    $\frac{1}{9}$   |    $\frac{X_1+2X_3}{3}$    | 
| $\left( X_2,X_3,X_3 \right)$   |    $\left( 0,1,2 \right)$   |    $\left( 0,\frac{1}{3},\frac{2}{3} \right)$   |    $\frac{1}{9}$   |    $\frac{X_2+2X_3}{3}$    | 
| $\left( X_1,X_2,X_3 \right)$   |    $\left( 1,1,1 \right)$   |    $\left( \frac{1}{3},\frac{1}{3},\frac{1}{3} \right)$   |    $\frac{2}{9}$   |    $\frac{X_1+X_2+X_3}{3}$    | 
</div>\EndKnitrBlock{example}

En algunas ocasiones es factible encontrar expresiones cerradas
para la distribución de $R^{\ast}$, más allá de las obvias que
consisten en enumerar el ingente número de átomos de probabilidad de
$\mathbf{P}^{\ast}$:
$$P^{\ast}\left( R^{\ast}=R\left( \left( m_1,\ldots ,m_n \right)
,F_n \right) \right)$$ 
Veamos un ejemplo.

### Inferencia sobre la mediana

En el caso de la mediana, consideremos, por simplicidad el caso de
tamaño muestral impar, $n=2m-1$. Supongamos también que no hay empates
en los valores de la muestra (si los hubiese las expresiones serían más
farragosas pero también calculables). La versión bootstrap del
estadístico sobre el cual pivota la inferencia es $R^{\ast}=X_{\left(
m \right)}^{\ast}-X_{(m)}$. Su distribución bootstrap
podría calcularse si se obtuviese la de $X_{(m)}^{\ast}$.
Pero ésta es factible de calcular por los pocos posibles valores que
puede tomar el estadístico $X_{(m)}^{\ast}$ (tan sólo los
valores de la muestra original) y por la sencillez del bootstrap
uniforme.
Veámoslo:
$$P^{\ast}\left( X_{(m)}^{\ast}>X_{(j)} \right)
=P^{\ast}\left( \#\left\{ X_i^{\ast}\leq X_{(j)}\right\}
\leq m-1 \right),$$
pero 
$$\#\left\{ X_i^{\ast}\leq X_{(j)}\right\} \sim \mathcal{B}\left(
n,\frac{j}{n} \right),$$
con lo cual
$$P^{\ast}\left( X_{(m)}^{\ast}>X_{(j)} \right)
=\sum_{k=0}^{m-1}\binom{n}{k}\left( \frac{j}{n} \right)^{k}
\left( \frac{n-j}{n} \right)^{n-k}$$
y, por lo tanto, si $j\geq 2$, 
$$\begin{aligned}
P^{\ast}\left( X_{(m)}^{\ast}=X_{(j)} \right)
=&\ P^{\ast}\left( X_{(m)}^{\ast}>X_{\left( j-1 \right)} \right)
-P^{\ast}\left( X_{(m)}^{\ast}>X_{(j)} \right) \\
=&\ \sum_{k=0}^{m-1}\binom{n}{k}\left( \frac{j-1}{n} \right)^{k}\left( \frac{
n-j+1}{n} \right)^{n-k} \\
&-\sum_{k=0}^{m-1}\binom{n}{k}\left( \frac{j}{n} \right)^{k}\left( \frac{n-j}{
n} \right)^{n-k} \\
=&\ \sum_{k=0}^{m-1}\binom{n}{k}\left[ \left( \frac{j-1}{n} \right)^{k}\left( 
\frac{n-j+1}{n} \right)^{n-k}-\left( \frac{j}{n} \right)^{k}\left( \frac{n-j
}{n} \right)^{n-k}\right] .
\end{aligned}$$

Cuando $j=1$, entonces 
$$\begin{aligned}
P^{\ast}\left( X_{(m)}^{\ast} = X_{(1)} \right)
&= 1-P^{\ast}\left( X_{(m)}^{\ast}>X_{(1)} \right) \\
&=  1-\sum_{k=0}^{m-1}\binom{n}{k}\left( \frac{1}{n} \right)^{k}
\left( \frac{n-1}{n} \right)^{n-k}.
\end{aligned}$$

<!-- 
Se podría comparar la distribución exacta 
con la aproximación por Monte carlo 
-->

### Distribución Bootstrap aproximada por Monte Carlo

Como ya se comentó anteriormente, 
al conocer el mecanismo que genera los datos en el bootstrap,
siempre se podrá simular dicho mecanismo mediante el método de Monte
Carlo. Por lo que el algoritmo general para la aproximación de Monte Carlo del
bootstrap uniforme es:

1. Para cada $i=1,\ldots ,n$ arrojar $X_i^{\ast}$ a partir de $F_n$

2. Obtener $\mathbf{X}^{\ast}=\left( X_1^{\ast},\ldots
,X_n^{\ast} \right)$

3. Calcular $R^{\ast}=R\left( \mathbf{X}^{\ast},F_n \right)$

4. Repetir $B$ veces los pasos 1-3 para obtener las réplicas bootstrap
$R^{\ast (1)}$, $\ldots$, $R^{\ast (B)}$

5. Utilizar esas réplicas bootstrap para aproximar la distribución en el
muestreo de $R$

Como se mostró en la Sección \@ref(intro-implementacion), el paso 1 se puede llevar a cabo simulando una distribución uniforme discreta
mediante el método de la transformación cuantil:

1. Para cada $i=1,\ldots ,n$ arrojar $U_i\sim \mathcal{U}\left( 0,1 \right)$ y
hacer $X_i^{\ast}=X_{\left\lfloor nU_i\right\rfloor +1}$

Aunque en `R` se recomienda emplear la función `sample`.

\BeginKnitrBlock{example}\iffalse{-91-73-110-102-101-114-101-110-99-105-97-32-115-111-98-114-101-32-108-97-32-109-101-100-105-97-32-99-111-110-32-118-97-114-105-97-110-122-97-32-100-101-115-99-111-110-111-99-105-100-97-93-}\fi{}<div class="example"><span class="example" id="exm:media-dt-desconocida"><strong>(\#exm:media-dt-desconocida)  \iffalse (Inferencia sobre la media con varianza desconocida) \fi{} </strong></span><br> \vspace{0.5cm}

Continuando con el ejemplo de los tiempos de vida de microorganismos,
supongamos que queremos obtener una estimación por intervalo de confianza 
de su vida media a partir de los 15 valores observados pero en la
situación mucho más realista de que la varianza sea desconocida.</div>\EndKnitrBlock{example}

Tenemos pues
$\mathbf{X}=\left( X_1,\ldots ,X_n \right) \sim F\,$, con
$\mu$ y $\sigma$ desconocidas

El parámetro de interés es
$$\theta \left( F \right) =\mu =\int x~dF\left( x \right)$$
que se estima mediante
$$\theta \left( F_n \right) =\int x~dF_n\left( x \right) =\bar{X}.$$
Así pues, el estadístico en el que basar la inferencia es
$$R=R\left( \mathbf{X},F \right) =\sqrt{n}\frac{\bar{X}-\mu }{S_{n-1}},$$
donde $S_{n-1}^2$ es la cuasivarianza muestral:
$$S_{n-1}^2=\frac{1}{n-1}\sum_{j=1}^{n}\left( X_j-\bar{X} \right)^2.$$

Bajo normalidad $\left( X\sim \mathcal{N}\left( \mu ,\sigma^2 \right) \right)$,
se sabe que $R\sim t_{n-1}$ y, en particular,
$R\overset{d}{\rightarrow } \mathcal{N}\left( 0,1 \right)$ cuando $n\rightarrow \infty$. 
Si $F$ no es normal entonces la distribución de $R$ ya no es una $t_{n-1}$, 
pero también es cierto que, bajo ciertas condiciones,
$R\overset{d}{\rightarrow}\mathcal{N}\left(0,1 \right)$.

En el contexto bootstrap elegimos $\hat{F}=F_n\,$, con lo cual se
trata de un bootstrap naïve o uniforme. El análogo bootstrap del
estadístico $R$ será
$$R^{\ast}=R\left( \mathbf{X}^{\ast},F_n \right) =\sqrt{n}\frac{
\bar{X}^{\ast}-\bar{X}}{S_{n-1}^{\ast}},$$
siendo
$$\begin{aligned}
\bar{X}^{\ast} &= \frac{1}{n}\sum_{i=1}^{n}X_i^{\ast}, \\
S_{n-1}^{\ast 2} &= \frac{1}{n-1}\sum_{i=1}^{n}\left( X_i^{\ast}-
\bar{X}^{\ast} \right)^2.
\end{aligned}$$

El algoritmo bootstrap (aproximado por Monte Carlo) procedería así:

1. Para cada $i=1,\ldots ,n$ arrojar $U_i\sim \mathcal{U}\left( 0,1 \right)$ y
hacer $X_i^{\ast}=X_{\left\lfloor nU_i\right\rfloor +1}$

2. Obtener $\bar{X}^{\ast}$ y $S_{n-1}^{\ast 2}$

3. Calcular
$R^{\ast}=\sqrt{n}\frac{\bar{X}^{\ast}-\bar{X}}{
S_{n-1}^{\ast}}$

4. Repetir $B$ veces los pasos 1-3 para obtener las réplicas bootstrap
$R^{\ast (1)}, \ldots, R^{\ast (B)}$

5. Aproximar la distribución en el muestreo de $R$ mediante la
distribución empírica de $R^{\ast (1)}, \ldots, R^{\ast (B)}$

El código para implementar este método es similar al del caso de varianza conocida
del Ejemplo \@ref(exm:media-dt-conocida):

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
  estadistico_boot[k] <- sqrt(n) * (x_barra_boot - x_barra)/cuasi_dt_boot
}

# Aproximación bootstrap de los ptos críticos
pto_crit <- quantile(estadistico_boot, c(alfa/2, 1 - alfa/2))

# Construcción del IC
ic_inf_boot <- x_barra - pto_crit[2] * cuasi_dt/sqrt(n)
ic_sup_boot <- x_barra - pto_crit[1] * cuasi_dt/sqrt(n)
IC_boot <- c(ic_inf_boot, ic_sup_boot)
names(IC_boot) <- paste0(100*c(alfa/2, 1-alfa/2), "%")
IC_boot
```

```
##      2.5%     97.5% 
## 0.5030131 1.2888063
```

Este procedimiento para la construcción de intervalos de confianza
se denomina *método percentil-t* y se tratará en la Sección \@ref(icboot-perc-t).

Como ejemplo adicional podemos comparar la aproximación de la distribución bootstrap del estadístico con la aproximación $t_{n-1}$ basada en normalidad.


```r
hist(estadistico_boot, freq = FALSE, ylim = c(0, 0.4))
abline(v = pto_crit)
curve(dt(x, n-1), add=TRUE, lty = 2)
pto_crit_t <- qt(1 - alfa/2, n-1)
abline(v = c(-pto_crit_t, pto_crit_t), lty = 2)
```



\begin{center}\includegraphics[width=0.7\linewidth]{01-intro_files/figure-latex/unnamed-chunk-11-1} \end{center}

En este caso la distribución bootstrap del estadístico es más asimétrica, lo que se traducirá en diferencias entre las estimaciones por intervalos de confianza.
Por ejemplo, podemos obtener la estimación basada en normalidad mediante la función `t.test()`:


```r
t.test(muestra)$conf.int
```

```
## [1] 0.4599374 1.1507292
## attr(,"conf.level")
## [1] 0.95
```


### Elección del número de réplicas Monte Carlo

Normalmente el valor de $B$ se toma del orden de varias centenas o
incluso millares. En los casos en los que el bootstrap se utiliza para
estimar el sesgo o la varianza de un estimador, bastará tomar un número,
$B$, de réplicas bootstrap del orden de $B = 100, 200, 500$. Sin embargo,
cuando se trata de utilizar el bootstrap para realizar contrastes de
hipótesis o construir intervalos de confianza son necesarios valores
mayores, del tipo $B = 500, 1000, 2000, 5000$.

<!-- 
considerar distintos estadísticos: 
media, error estándar y cuantiles
-->

Evidentemente, la función de distribución del estadístico de interés,
$\psi \left( u \right) =P\left( R\leq u \right)$, se estimaría mediante
la distribución empírica de las $B$ realizaciones de la aproximación de
Monte Carlo, 
$$\hat{\psi}_{B}\left( u \right) =
\frac{1}{B}\sum_{i=1}^{B}\mathbf{1}\left\{ R^{\ast (i)}\leq u\right\},$$
de la verdadera distribución bootstrap exacta: 
$\hat{\psi}\left(u \right) =P^{\ast}\left( R^{\ast}\leq u \right)$. 
El error de MonteCarlo de $\hat{\psi}_{B}\left( u \right)$ con respecto 
a $\hat{\psi}\left( u \right)$ viene dado por su varianza Monte Carlo, 
pues su sesgo Monte Carlo es cero:

$$\begin{aligned}
E^{MC}\left( \hat{\psi}_{B}\left( u \right) \right) &= \frac{1}{B}
\sum_{i=1}^{B}E^{MC}\left( \mathbf{1}\left\{ R^{\ast (i)}\leq
u\right\} \right) =\frac{1}{B}\sum_{i=1}^{B}P^{\ast}\left( R^{\ast \left(
i \right)}\leq u \right) \\
&= \frac{1}{B}\sum_{i=1}^{B}\hat{\psi}\left( u \right) =\hat{\psi}\left(
u \right), \\
Var^{MC}\left( \hat{\psi}_{B}\left( u \right) \right) &= \frac{1}{B^2}
\sum_{i=1}^{B}Var^{MC}\left( \mathbf{1}\left\{ R^{\ast (i)}\leq
u\right\} \right) \\
&= \frac{1}{B^2}\sum_{i=1}^{B}P^{\ast}\left( R^{\ast (i)
}\leq u \right) \left[ 1-P^{\ast}\left( R^{\ast (i)}\leq
u \right) \right] = \\
&= \frac{1}{B^2}\sum_{i=1}^{B}\hat{\psi}\left( u \right) \left( 1-\hat{\psi}
\left( u \right) \right) =\frac{1}{B}\hat{\psi}\left( u \right) \left( 1-\hat{
\psi}\left( u \right) \right) \leq \frac{1}{4B}
\end{aligned}$$

Así, el error de la aproximación de Monte Carlo al bootstrap exacto
(raíz cuadrada de la varianza del Monte Carlo), puede acotarse por
$\frac{1}{2\sqrt{B}}$
(para más detalles sobre la convergencia de una aproximación Monte Carlo ver p.e. el [Capítulo 4](https://rubenfcasal.github.io/simbook/cap4.html) de Fernández-Casal y Cao, 2020).

## Herramientas disponibles en R sobre bootstrap  {#intro-paquetes}

En `R` hay una gran cantidad de paquetes que implementan métodos bootstrap.
Por ejemplo, al ejecutar el comando `??bootstrap` (o `help.search('bootstrap')`)
se mostrarán las funciones de los paquetes instalados que incluyen este término
en su documentación (se puede realizar la búsqueda en todos los paquetes disponibles
de `R` a través de <https://www.rdocumentation.org>).

De entre todos estas herramientas destacan dos librerías 
como las más empleadas:

* `bootstrap`: contiene las rutinas (bootstrap, cross-validation,
  jackknife) y los datos del libro "An Introduction to the Bootstrap" de B.
  Efron y R. Tibshirani, 1993, Chapman and Hall. La librería fue
  desarrollada originalmente en `S` por Rob Tibshirani y exportada a `R` por
  Friedrich Leisch. Es útil para desarrollar los ejemplos que se citan en
  ese libro.

* `boot`: incluye las funciones y conjuntos de datos utilizados en el libro 
  "Bootstrap Methods and Their Applications" de A. C. Davison y D. V. Hinkley, 1997,
  Cambridge University Press. Esta librería fue desarrollada originalmente 
  en `S` por Angelo J. Canty y posteriormente exportada a `R` (ver [Canty, 2002](http://cran.fhcrc.org/doc/Rnews/Rnews_2002-3.pdf)).
  Este paquete es mucho más completo que el paquete `bootstrap`, forma parte de la distribución estándar de `R` y es el que emplearemos como referencia en este libro (ver Sección \@ref(intro-pkgboot)).

Por otra parte existen numerosas rutinas (scripts) realizadas en `R` por
diversos autores, que están disponibles en Internet 
(por ejemplo, puede ser interesante realizar una búsqueda en 
<https://rseek.org>). 

El bootstrap uniforme se puede implementar fácilmente. Por ejemplo,
una rutina general para el caso univariante sería la siguiente:


```r
#' @param x vector que contiene la muestra.
#' @param B número de réplicas bootstrap.
#' @param statistic función que calcula el estadístico.
boot.strap0 <- function(x, B=1000, statistic=mean){
  ndat <- length(x)
  x.boot <- sample(x, ndat*B, replace=TRUE)
  x.boot <- matrix(x.boot, ncol=B, nrow=ndat)
  stat.boot <- apply(x.boot, 2, statistic)
}
```

Podríamos aplicar esta función a la muestra de tiempos de vida de
microorganismos con el siguiente código:

```r
fstatistic <- function(x){
  #  mean(x)
  #  mean(x, trim=0.2)
  median(x)
  #  max(x)
}

B <- 1000
set.seed(1)
stat.dat <- fstatistic(muestra)
stat.boot <- boot.strap0(muestra, B, fstatistic)

res.boot <- c(stat.dat, mean(stat.boot)-stat.dat, sd(stat.boot))
names(res.boot) <- c("Estadístico", "Sesgo", "Error Std.")
res.boot
```

```
## Estadístico       Sesgo  Error Std. 
##   0.6110000   0.0609260   0.2481362
```

La función `boot.strap0()` anterior no es adecuada para el caso multivariante
(por ejemplo cuando estamos interesados en regresión).
Como se mostró en la Sección \@ref(intro-implementacion)
sería preferible emplear remuestras del vector de índices. Por ejemplo:


```r
#' @param datos vector, matriz o data.frame que contiene los datos.
#' @param B número de réplicas bootstrap.
#' @param statistic función con al menos dos parámetros, 
#' los datos y el vector de índices de remuestreo, 
#' y que devuelve el vector de estadísticos.
#' @param ... parámetros adicionales de la función statistic.
boot.strap <- function(datos, B=1000, statistic, ...) {
  ndat <- NROW(datos)
  i.boot <- sample(ndat, ndat*B, replace=TRUE)
  i.boot <- matrix(i.boot, ncol=B, nrow=ndat)
  stat.boot <- drop(apply(i.boot, 2, function(i) statistic(datos, i, ...)))
}
```

El paquete `boot`, descrito a continuación, emplea una implementación similar.

### El paquete `boot` {#intro-pkgboot}

La función principal de este paquete es la función `boot()` que implementa distintos métodos de remuestreo para datos i.i.d..
En su forma más simple permite realizar bootstrap uniforme (que en la práctica también se denomina habitualmente *bootstrap noparamétrico*):

```r
boot(data, statistic, R)
```
donde `data` es un vector, matriz o `data.frame` que contiene los datos, 
`R` es el número de réplicas bootstrap, y `statistic` es una función 
con al menos dos parámetros (con las opciones por defecto), 
los datos y el vector de índices de remuestreo, 
y que devuelve el vector de estadísticos.

Por ejemplo, para hacer inferencia sobre la mediana del tiempo de microorganismos,
podríamos emplear el siguiente código:

```r
library(boot)
muestra <- c(0.143, 0.182, 0.256, 0.26, 0.27, 0.437, 0.509, 
             0.611, 0.712, 1.04, 1.09, 1.15, 1.46, 1.88, 2.08)

statistic <- function(data, i){
  # remuestra <- data[i]; median(remuestra)
  median(data[i])
}

set.seed(1)
res.boot <- boot(muestra, statistic, R = 1000)
```

El resultado que devuelve esta función es un objeto de clase `boot`, una lista con los siguientes componentes:

```r
names(res.boot)
```

```
##  [1] "t0"        "t"         "R"         "data"      "seed"      "statistic"
##  [7] "sim"       "call"      "stype"     "strata"    "weights"
```
Además de los parámetros de entrada (incluyendo los valores por defecto), contiene tres componentes adicionales:

* `tO`: el valor observado del estadístico 
  (su evaluación en los datos originales).
  
* `t`: la matriz de réplicas bootstrap del estadístico
  (cada fila se corresponde con una remuestra).
  
* `seed`: el valor inicial de la semilla (`.Random.seed`)
  empleada para la generación de las réplicas.

Este tipo de objetos dispone de dos métodos principales:
el método `print()` que muestra un resumen de los resultados
(incluyendo  aproximaciones bootstrap del sesgo y del error
estándar de los estadísticos; ver Capítulo \@ref(prec-sesgo)):

```r
res.boot
```

```
## 
## ORDINARY NONPARAMETRIC BOOTSTRAP
## 
## 
## Call:
## boot(data = muestra, statistic = statistic, R = 1000)
## 
## 
## Bootstrap Statistics :
##     original   bias    std. error
## t1*    0.611 0.058523   0.2526519
```
y el método `plot()` que genera gráficas básicas de diagnosis
de los resultados (correspondientes al estadístico determinado por el parámetro `index`, por defecto `= 1`): [Figura \@ref(fig:plot-res-boot)]


```r
plot(res.boot)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{01-intro_files/figure-latex/plot-res-boot-1} 

}

\caption{Gráficos de diagnóstico de los resultados bootstrap de la mediana de los tiempos de vida de microorganismos.}(\#fig:plot-res-boot)
\end{figure}

Es recomendable examinar la distribución bootstrap del estimador (o estadístico) para detectar posibles problemas.
Como en este caso puede ocurrir que el estadístico bootstrap tome pocos valores distintos, lo que indicaría que el número de réplicas bootstrap es insuficiente o que hay algún problema con método de remuestreo empleado (en este caso la distribución objetivo es continua). 
Se darán más detalles sobre los posibles problemas del bootstrap uniforme en la Sección \@ref(deficien-unif).

Además de estos métodos, las principales funciones de interés serían:

* `jack.after.boot()`: genera un gráfico para diagnósticar la inluencia 
  de las observaciones individuales en los resultados bootstrap 
  (se representan los cuantiles frente a las diferencias en el estadístico 
  al eliminar una observación; este gráfico también se puede obtener estableciendo
  `jack = TRUE` en `plot.boot()`).
  
* `boot.array()`: genera la matriz de índices a partir de la que se obtuvieron las remuestras (permite reconstruir las remuestras bootstrap).
  
* `boot.ci()`: construye distintos tipos de intervalos de confianza 
  (se tratarán en el Capítulo \@ref(icboot)) dependiendo del parámetro `type`:
  
    - `"norm"`: utiliza la distribución asintótica normal considerando las
        aproximaciones bootstrap del sesgo y de la varianza.
        
    - `"basic"`: emplea el estadístico $R = \hat \theta - \theta$ para la
      construcción del intervalo de confianza.
      
    - `"stud"`: calcula el intervalo a partir del estadístico estudentizado 
        $R = \left( \hat \theta - \theta \right) / \sqrt{Var(\hat \theta)}$.
    
    - `"perc"`: utiliza directamente la distribución bootstrap del estadístico
      ($R = \hat \theta$).
    
    - `"bca"`: emplea el método $BCa$ ("bias-corrected and accelerated") 
        propuesto por Efron (1987) (ver Sección 5.3.2 de Davison y Hinkley, 1997).
    
    - `"all"`: calcula los cinco tipos de intervalos anteriores.


Como ya se comentó, la función `boot()` admite estadísticos multivariantes 
(haciendo que la función `statistic` devuelva un vector en lugar de un escalar),
pero por defecto las funciones anteriores consideran el primer componente
como el estadístico principal. 
Para obtener resultados de otros componentes del vector de estadísticos
habrá que establecer el parámetro `index` igual al índice deseado.
Además, en algunos casos (por ejemplo para la obtención de intevalos de confianza
estudentizados con la función `boot.ci()`) se supone, por defecto, que el segundo
componente del vector de estadísticos contiene estimaciones de la varianza del
estadístico para cada réplica boostrap.

\BeginKnitrBlock{example}\iffalse{-91-73-110-102-101-114-101-110-99-105-97-32-115-111-98-114-101-32-108-97-32-109-101-100-105-97-32-99-111-110-32-118-97-114-105-97-110-122-97-32-100-101-115-99-111-110-111-99-105-100-97-44-32-99-111-110-116-105-110-117-97-99-105-243-110-93-}\fi{}<div class="example"><span class="example" id="exm:media-dt-desconocida-boot"><strong>(\#exm:media-dt-desconocida-boot)  \iffalse (Inferencia sobre la media con varianza desconocida, continuación) \fi{} </strong></span><br> \vspace{0.5cm}

Continuando con el Ejemplo \@ref(exm:media-dt-desconocida) de
inferencia sobre la media con varianza desconocida. 
Para obtener la estimación por intervalo de confianza del tiempo de vida medio 
de los microorganismos con el paquete `boot`, podríamos emplear
el siguiente código:</div>\EndKnitrBlock{example}

```r
library(boot)
muestra <- c(0.143, 0.182, 0.256, 0.26, 0.27, 0.437, 0.509, 
             0.611, 0.712, 1.04, 1.09, 1.15, 1.46, 1.88, 2.08)

statistic <- function(data, i){
  remuestra <- data[i]
  c(mean(remuestra), var(remuestra)/length(remuestra))
}

set.seed(1)
res.boot <- boot(muestra, statistic, R = 1000)
res.boot
```

```
## 
## ORDINARY NONPARAMETRIC BOOTSTRAP
## 
## 
## Call:
## boot(data = muestra, statistic = statistic, R = 1000)
## 
## 
## Bootstrap Statistics :
##      original       bias    std. error
## t1* 0.8053333  0.003173267 0.158330646
## t2* 0.0259338 -0.002155755 0.007594682
```

```r
boot.ci(res.boot)
```

```
## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
## Based on 1000 bootstrap replicates
## 
## CALL : 
## boot.ci(boot.out = res.boot)
## 
## Intervals : 
## Level      Normal              Basic             Studentized     
## 95%   ( 0.4918,  1.1125 )   ( 0.4825,  1.0980 )   ( 0.4715,  1.2320 )  
## 
## Level     Percentile            BCa          
## 95%   ( 0.5127,  1.1282 )   ( 0.5384,  1.1543 )  
## Calculations and Intervals on Original Scale
```

El intervalo marcado como `Studentized` se obtuvo empleando el mismo estadístico
del Ejemplo \@ref(exm:media-dt-desconocida).


***Modificaciones del bootstrap uniforme*** 

Establecenciendo parámetros adicionales de la función `boot` se pueden llevar 
a cabo modificaciones del bootstrap uniforme. 
Algunos de estos parámetros son los siguientes:

* `strata`: permite realizar remuestreo estratificado estableciendo este parámetro
  como un vector numérico o factor que defina los grupos.

* `sim = c("ordinary" , "parametric", "balanced", "permutation", "antithetic")`:
  permite establecer distintos tipos de remuestreo. 
  Por defecto es igual a `"ordinary"` que se corresponde con el bootstrap uniforme,
  descrito anteriormente. Entre el resto de opciones destacaríamos 
  `sim = "permutation"`, que permite realizar contrastes de
  permutaciones (remuestreo sin reemplazamiento), y `sim = "parametric"`,
  que permite realizar bootstrap paramétrico (Sección \@ref(modunif-boot-par)). 
  En este último caso también habrá que establecer los parámetros `ran.gen` y
  `mle`, y la función `statistics` no empleará el segundo parámetro de índices.

* `ran.gen`: función que genera los datos. El primer argumento será el conjunto de datos
  original y el segundo un vector de parámetros adicionales 
  (normalmente los valores de los parámetros de la distribución).

* `mle`: parámetros de la distribución (típicamente estimados por máxima verosimilitud)
  o parámetros adicionales para `ran.gen` ó `statistics`.

Además hay otros parámetros para el procesamiento en paralelo: `parallel = c("no", "multicore", "snow")`, `ncpus`, `cl`. 
En el Apéndice \@ref(intro-hpc) se incluye una pequeña introducción al procesamiento en paralelo y se muestran algunos ejemplos sobre el uso de estos parámetros.
También se puede consultar la ayuda de la función `boot()` (`?boot`).

El paquete `boot` también incluye otras funciones que implementan métodos
boostrap para otros tipos de datos, como la función `censboot()` para datos 
censurados (Capítulo \@ref(bootcen)) o la función `tsboot()` para series de tiempo (Capítulo \@ref(bootdep)).

Finalmente destacar que hay numerosas extensiones implementadas en otros paquetes utilizando el paquete `boot` (ver *Reverse dependencies* en la [web de CRAN](https://cran.r-project.org/package=boot)).
Por ejemplo en la Sección \@ref(boot-reg) se ilustrará el uso de la función `Boot()` del paquete `car` para hacer inferencia sobre modelos de regresión.


### Ejemplo: Bootstrap uniforme multidimensional {#boot-unif-multi}

Como ya se mostró en las Secciones \@ref(intro-implementacion) y \@ref(intro-paquetes) podemos implementar el bootstrap uniforme en el caso multidimensional (denominado también *remuestreo de casos* o *bootstrap de las observaciones*) de modo análogo al unidimensional.

Consideraremos como ejemplo el conjunto de datos `Prestige` del paquete `carData`, y supongamos que queremos realizar inferencias sobre el coeficiente de correlación entre `prestige` (puntuación de ocupaciones obtenidas a partir de una encuesta) e`income` (media de ingresos en la ocupación).


```r
data(Prestige, package = "carData")
str(Prestige)
```

```
## 'data.frame':	102 obs. of  6 variables:
##  $ education: num  13.1 12.3 12.8 11.4 14.6 ...
##  $ income   : int  12351 25879 9271 8865 8403 11030 8258 14163 11377 11023 ...
##  $ women    : num  11.16 4.02 15.7 9.11 11.68 ...
##  $ prestige : num  68.8 69.1 63.4 56.8 73.5 77.6 72.6 78.1 73.1 68.8 ...
##  $ census   : int  1113 1130 1171 1175 2111 2113 2133 2141 2143 2153 ...
##  $ type     : Factor w/ 3 levels "bc","prof","wc": 2 2 2 2 2 2 2 2 2 2 ...
```

```r
# with(Prestige, cor(income, prestige))
cor(Prestige$income, Prestige$prestige)
```

```
## [1] 0.7149057
```

En el siguiente código se emplea el paquete `boot` para realizar bootstrap uniforme multidimensional sobre este estadístico:


```r
library(boot)

statistic <- function(data, i){
  remuestra <- data[i, ]
  cor(remuestra$income, remuestra$prestige)
}

set.seed(1)
B <- 1000
res.boot <- boot(Prestige, statistic, R = B)
res.boot
```

```
## 
## ORDINARY NONPARAMETRIC BOOTSTRAP
## 
## 
## Call:
## boot(data = Prestige, statistic = statistic, R = B)
## 
## 
## Bootstrap Statistics :
##      original      bias    std. error
## t1* 0.7149057 0.006306905  0.04406473
```

```r
plot(res.boot)
```



\begin{center}\includegraphics[width=0.7\linewidth]{01-intro_files/figure-latex/unnamed-chunk-22-1} \end{center}

En este caso podemos observar que la distribución bootstrap del estimador es asimétrica, por lo que asumir que su distribución es normal podría no ser adecuado (por ejemplo para la construcción de intervalos de confianza, que se tratarán en la Sección \@ref(icboot-trans)).

Como comentario final, nótese que en principio el paquete boot está diseñado para obtener réplicas bootstrap de un estimador, por lo que si lo que nos interesa es emplear otro estadístico habría que construirlo a partir de ellas (como hacen otras funciones secundarias como `boot.ci()`).
Por ejemplo, si queremos emplear el estadístico $R = \hat \theta - \theta$
(bootstrap percentil básico o natural), podemos obtener la correspondiente distribución bootstrap (aproximada por Monte Carlo) con el siguiente código:


```r
estadistico_boot <- res.boot$t - res.boot$t0 
hist(estadistico_boot)
```



\begin{center}\includegraphics[width=0.7\linewidth]{01-intro_files/figure-latex/unnamed-chunk-23-1} \end{center}

A partir de la distribución empírica del estadístico bootstrap $R^{\ast} = \hat \theta^{\ast} - \hat \theta$ aproximaríamos la característica de interés de la distribución en el muestreo de $R = \hat \theta - \theta$.
Por ejemplo, para aproximar $\psi \left( u \right) =P\left( R\leq u \right)$ emplearíamos la frecuencia relativa: 
$$\hat{\psi}_{B}\left( u \right) =
\frac{1}{B}\sum_{i=1}^{B}\mathbf{1}\left\{ R^{\ast (i)}\leq u\right\}.$$


```r
u <- 0
sum(estadistico_boot <= u)/B
```

```
## [1] 0.427
```

```r
# Equivalentemente:
mean(estadistico_boot <= u)
```

```
## [1] 0.427
```

