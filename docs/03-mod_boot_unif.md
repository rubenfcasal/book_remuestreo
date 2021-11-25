# Modificaciones del Bootstrap uniforme {#modunif}





El bootstrap uniforme (o naïve) es aquel en el que remuestreamos a
partir de la función de distribución empírica. Eso es muy razonable
cuando no tenemos ninguna información adicional sobre la función de
distribución poblacional, ya que la distribución empírica es el
estimador máximo verosímil no paramétrico de la función de distribución
poblacional. Sin embargo, cuando en el contexto en el que nos
encontremos conozcamos alguna propiedad adicional de dicha distribución
poblacional, entonces debemos incorporarla en el método de remuestreo,
dando lugar a otro método bootstrap que ya no debemos llamar uniforme o
naïve. Veremos algunos de ellos.

## Bootstrap paramétrico {#modunif-boot-par}

Supongamos que sabemos que la función de distribución poblacional
pertenece a cierta familia paramétrica. Es decir $F=F_{\theta }$ para
algún vector $d$-dimensional $\theta \in \Theta$. En ese caso parece
lógico estimar $\theta$ a partir de la muestra (denotemos
$\hat{\theta}$ un estimador de $\theta$, por ejemplo el de máxima
verosimilitud) y obtener remuestras de $F_{\hat{\theta}}$ no de $F_n$.
Entonces, el bootstrap uniforme se modifica de la siguiente forma, dando
lugar al llamado bootstrap paramétrico:

1.  Dada la muestra
    $\mathbf{X}=\left( X_1,\ldots ,X_n \right)$, calcular
    $\hat{\theta}$

2.  Para cada $i=1,\ldots ,n$ arrojar $X_i^{\ast}$ a partir de
    $F_{\hat{\theta}}$

3.  Obtener $\mathbf{X}^{\ast}=\left( X_1^{\ast},\ldots
    ,X_n^{\ast} \right)$

4.  Calcular $R^{\ast}=R\left( \mathbf{X}^{\ast},F_{\hat{\theta}} \right)$

Así utilizaremos las distribución en el remuestreo de $R^{\ast}$ para
aproximar la distribución en el muestreo de $R$. Lógicamente, cuando no
sea posible obtener una expresión explícita para la distribución
bootstrap de $R^{\ast}$ utilizaremos una aproximación de Monte Carlo de
la misma:

1. Dada la muestra
$\mathbf{X}=\left( X_1,\ldots ,X_n \right)$, calcular
$\hat{\theta}$

2. Para cada $i=1,\ldots ,n$ arrojar $X_i^{\ast}$ a partir de
$F_{\hat{\theta}}$

3. Obtener $\mathbf{X}^{\ast}=\left( X_1^{\ast},\ldots
,X_n^{\ast} \right)$

4. Calcular
$R^{\ast}=R\left( \mathbf{X}^{\ast},F_{\hat{\theta}
} \right)$

5. Repetir $B$ veces los pasos 2-4 para obtener las réplicas bootstrap
$R^{\ast (1)}$, $\ldots$, $R^{\ast (B)}$

6. Utilizar esas réplicas bootstrap para aproximar la distribución en el
muestreo de $R$

En general, para llevar a cabo el paso 2, debemos poder simular valores
de la distribución $F_{\hat{\theta}}$ (en el caso del bootstrap uniforme
se trataba de simular valores de la distribución empírica, lo cual es
muy sencillo y rápido). Para ello podemos utilizar el método de
inversión, que consiste en simular un valor $U$ procedente de una
distribución $\mathcal{U}\left( 0,1 \right)$ (es decir, $U$ es un número aleatorio
uniforme) y devolver $X^{\ast}=F_{\hat{\theta}}^{-1}\left(
U \right)$. Así, podríamos escribir el paso 2 de una forma más
detallada:

2. Para cada $i=1,\ldots ,n$ arrojar $U_i\sim \mathcal{U}\left( 0,1 \right)$ y
hacer $X_i^{\ast}=F_{\hat{\theta}}^{-1}\left( U_i \right)$

No en todos los modelos paramétricos es fácil de calcular la inversa
$F_{\hat{\theta}}^{-1}$. En algunos modelos paramétricos (como el caso
de la distribución normal) ni siquiera tenemos una fórmula explícita
para $F_{\theta }\left( x \right)$, con lo cual difícilmente podremos
calcular explícitamente su inversa. En casos como esos es frecuente
recurrir a otros métodos para simular la distribución en cuestión.
Normalmente existen rutinas incorporadas a la mayoría de los lenguajes
de programación y software estadístico (como `R`) que permiten simular
directamente la mayoría de las distribuciones paramétricas habituales.


\BeginKnitrBlock{example}\iffalse{-91-73-110-102-101-114-101-110-99-105-97-32-115-111-98-114-101-32-108-97-32-109-101-100-105-97-32-99-111-110-32-118-97-114-105-97-110-122-97-32-99-111-110-111-99-105-100-97-44-32-99-111-110-116-105-110-117-97-99-105-243-110-93-}\fi{}
<span class="example" id="exm:media-dt-conocida-par"><strong>(\#exm:media-dt-conocida-par)  \iffalse (Inferencia sobre la media con varianza conocida, continuación) \fi{} </strong></span>
\EndKnitrBlock{example}
Continuando con el ejemplo de tiempo de vida de microorganismos,
podemos modificar fácilmente el código mostrado en el Ejemplo \@ref(exm:media-dt-conocida), de forma que se emplee bootstrap
paramétrico (normal), con desviación típica conocida, para
calcular un intervalo de confianza para la media poblacional.


```r
muestra <- c(0.143, 0.182, 0.256, 0.26, 0.27, 0.437, 0.509, 
             0.611, 0.712, 1.04, 1.09, 1.15, 1.46, 1.88, 2.08)
n <- length(muestra)
sigma <- 0.6

alfa <- 0.05
x_barra <- mean(muestra)

# Remuestreo
set.seed(1)
B <- 1000
estadistico_boot <- numeric(B)
for (k in 1:B) {
    # u <- rnorm(n)
    # remuestra <- u * sigma + x_barra
    remuestra <- rnorm(n, x_barra, sigma)
    x_barra_boot <- mean(remuestra)
    estadistico_boot[k] <- sqrt(n) * (x_barra_boot - x_barra)/sigma
}

# Aproximación Monte Carlo de los ptos críticos
# Empleando la distribución empírica del estadístico bootstrap: 
    # estadistico_boot_ordenado <- sort(estadistico_boot)
    # indice_inf <- floor(B * alfa/2)
    # indice_sup <- floor(B * (1 - alfa/2))
    # pto_crit <- estadistico_boot_ordenado[c(indice_inf, indice_sup)]
# Empleando la función `quantile`:
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
## 0.5236922 1.1217871
```

En este caso concreto la distribución bootstrap del estadístico sería conocida (normal estándar) y realmente no sería necesario emplear la aproximación Monte Carlo:


```r
hist(estadistico_boot, freq = FALSE)
abline(v = pto_crit)
curve(dnorm, add=TRUE, lty = 2)
pto_crit_teor <- qnorm(1 - alfa/2)
abline(v = c(-pto_crit_teor, pto_crit_teor), lty = 2)
```



\begin{center}\includegraphics[width=0.7\linewidth]{03-mod_boot_unif_files/figure-latex/unnamed-chunk-2-1} \end{center}

```r
ic_inf_boot_teor <- x_barra - pto_crit_teor * sigma/sqrt(n)
ic_sup_boot_teor <- x_barra + pto_crit_teor * sigma/sqrt(n)
IC_boot_teor <- c(ic_inf_boot_teor, ic_sup_boot_teor)
names(IC_boot_teor) <- paste0(100*c(alfa/2, 1-alfa/2), "%")
IC_boot_teor
```

```
##     2.5%    97.5% 
## 0.501697 1.108970
```


Para emplear el paquete `boot`, como se comentó en la Sección 
\@ref(intro-pkgboot), habría que establecer en la llamada a la 
función `boot()` los argumentos: `sim = "parametric"`, 
`mle` igual a los parámetros necesarios para la simulación y 
`ran.gen = function(data, mle)`, una función de los datos originales 
y de los parámetros que devuelve los datos generados. 
En este caso además, la función `statistic` no necesita el vector 
de índices como segundo parámetro. 
Por ejemplo, para calcular el intervalo de confianza para la media del 
tiempo de vida de los microorganismos, podríamos utilizar el siguiente código:


```r
library(boot)
ran.gen.norm <- function(data, mle) {
    # Función para generar muestras aleatorias normales
    # con desviación típica sigma = 0.6,
    # mle contendrá la media de los datos originales
    out <- rnorm(length(data), mle, sigma)
    out
}

statistic <- function(data){
    c(mean(data), sigma^2/length(data))
}

set.seed(1)
res.boot <- boot(muestra, statistic, R = B, sim = "parametric",
                 ran.gen = ran.gen.norm, mle = mean(muestra))

boot.ci(res.boot, type = "stud")
```

```
## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
## Based on 1000 bootstrap replicates
## 
## CALL : 
## boot.ci(boot.out = res.boot, type = "stud")
## 
## Intervals : 
## Level    Studentized     
## 95%   ( 0.5208,  1.1232 )  
## Calculations and Intervals on Original Scale
```

Aunque los resultados dependerán en gran medida de que el modelo paramétrico
sea adecuado para describir la variabilidad de los datos 
(en este caso no es muy razonable que el modelo admita tiempos de vida negativos).
Si, por ejemplo, consideramos que un modelo exponencial es más adecuado: [Figura \@ref(fig:boot-par-aprox)]


```r
# Distribución bootstrap uniforme
curve(ecdf(muestra)(x), xlim = c(-.5, 3), ylab = "F(x)", type = "s")
# Distribución bootstrap paramétrico normal
curve(pnorm(x, mean(muestra), 0.6), lty = 2, add = TRUE)
# Distribución bootstrap paramétrico exponencial
curve(pexp(x, 1/mean(muestra)), lty = 3, add = TRUE)
legend("bottomright", legend = c("Empírica", "Aprox. normal", "Aprox. exponencial"), lty = 1:3)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{03-mod_boot_unif_files/figure-latex/boot-par-aprox-1} 

}

\caption{Distribución empírica de la muestra de tiempos de vida de microorganismos y aproximaciones paramétricas.}(\#fig:boot-par-aprox)
\end{figure}

Solo tendríamos que cambiar la función que genera los datos:


```r
ran.gen.exp <- function(data, mle) {
    # Función para generar muestras aleatorias exponenciales
    # mle contendrá la media de los datos originales
    out <- rexp(length(data), 1/mle)
    out
}
```

Una de las principales aplicaciones del bootstrap paramétrico 
es el contraste de hipótesis que se tratará en la Sección \@ref(contrastes-parametricos).

<!-- 
https://rubenfcasal.github.io/simbook/aplicaciones-de-la-simulacion-en-inferencia-estadistica.html#contrastes-de-hipotesis
-->


## Bootstrap simetrizado

Supongamos que conocemos que la función de distribución poblacional es
simétrica entorno a cierto valor. Eso significa que existe un valor $c$
tal que $F\left( c-h \right) =1-F\left( c+h \right)$ para todo $h>0$.
Equivalentemente, una variable aleatoria es simétrica entorno a $c$ si
su función de distribución verifica 
$$F\left( x \right) = 1 - F\left( 2c - x \right)$$ 
para todo $x\in \mathbb{R}$. Puede demostrarse que dicho
centro de simetría, $c$, ha de ser la media de la distribución, 
$\mu$, en caso de que exista. Esa información (la simetría) sobre la
distribución poblacional también se debe incorporarse en el bootstrap.
Así, para estimar la función de distribución poblacional, $F$, supuesto
que es simétrica entorno a $\mu$, es razonable utilizar una versión
simetrizada de la distribución empírica, $F_n^{sim}$. Ese estimador
empírico simetrizado de la función de distribución es el que otorga
igual masa de probabilidad a una muestra artificialmente construida
simetrizando, alrededor de la media muestral, la muestra original:

$$Y_i=\left\{ 
\begin{array}{ll}
X_i & \text{si } i=1,\ldots ,n \\ 
2\bar{X}-X_{i-n} &\text{si } i=n+1,\ldots ,2n
\end{array}
\right.$$

con lo cual
$$F_n^{sim}\left( x \right) =\frac{1}{2n}\sum_{i=1}^{2n}\mathbf{1}\left( Y_i\leq x \right).$$
Puede demostrarse fácilmente que
$$F_n^{sim}\left( x \right) =\frac{1}{2}\left( F_n\left( x \right)
+1-F_n\left( 2\bar{X}-x \right) \right).$$

Al diseñar el plan de remuestreo debemos utilizar $F_n^{sim}$
(bootstrap simetrizado), en lugar de $F_n$ (bootstrap uniforme).

1. Para cada $i=1,\ldots ,n$ arrojar $X_i^{\ast}$ a partir de
$F_n^{sim}$, es decir
$P^{\ast}\left( X_i^{\ast}=Y_j \right) =\frac{1
}{2n}$, $j=1,\ldots ,2n$

2. Obtener $\mathbf{X}^{\ast}=\left( X_1^{\ast},\ldots
,X_n^{\ast} \right)$

3. Calcular
$R^{\ast}=R\left( \mathbf{X}^{\ast},F_n^{sim} \right)$

Como veremos más adelante, a veces (muy poco frecuentemente) es posible
calcular exactamente la distribución bootstrap de $R^{\ast}$. Cuando
eso no es posible, esa distribución es fácilmente aproximable por Monte
Carlo, arrojando una gran cantidad, $B$, de réplicas de $R^{\ast}$. En
ese caso, el algoritmo se convierte en:

1. Para cada $i=1,\ldots ,n$ arrojar $X_i^{\ast}$ a partir de
$F_n^{sim}$

2. Obtener $\mathbf{X}^{\ast}=\left( X_1^{\ast},\ldots
,X_n^{\ast} \right)$

3. Calcular
$R^{\ast}=R\left( \mathbf{X}^{\ast},F_n^{sim} \right)$

4. Repetir $B$ veces los pasos 1-3 para obtener las réplicas bootstrap
$R^{\ast (1)}$, $\ldots$, $R^{\ast (B)}$

5. Utilizar esas réplicas bootstrap para aproximar la distribución en el
muestreo de $R$

Para llevar a cabo el paso 1 podemos proceder de dos formas
equivalentes. La primera consiste en definir explícitamente la muestra
simetrizada en torno a la media, $\mathbf{Y}$, y luego obtener
uno de los valores de dicha muestra con equiprobabilidad. El paso 1
quedaría de la siguiente forma:

1. Para cada $i=1,\ldots ,n$ arrojar $U_i\sim \mathcal{U}\left( 0,1 \right)$ y
hacer $X_i^{\ast}=Y_{\left\lfloor 2nU_i\right\rfloor +1}$

Alternativamente podemos proceder con el paso 1 utilizando el hecho de que la función de distribución $F_n^{sim}\left( x \right)$ resultar ser la distribución de una variable aleatoria obtenida en dos etapas: 
en la primera etapa se genera un valor según la empírica, $F_n\left( x \right)$, y en la segunda se decide (con equiprobabilidad) si el valor obtenido no se altera o bien si se refleja alrededor de la media muestral, $\bar{X}$ 
(equivalentemente, la distribución simetrizada es una mixtura de la distribución empírica $F_n\left( x \right)$ y de su versión "reflejada" $1-F_n\left( 2\bar{X}-x \right)$ y se puede simular mediante el método de composición; ver p.e. Fernández-Casal y Cao, 2020, [Sección 5.4](https://rubenfcasal.github.io/simbook/m%C3%A9todo-de-composici%C3%B3n.html)). Así el paso 1 resulta:

1. Para cada $i=1,\ldots ,n$ arrojar
$U_i,V_i\sim \mathcal{U}\left( 0,1 \right)$. Si $V_i\leq \frac{1}{2}$
entonces hacer $X_i^{\ast}=X_{\left\lfloor
nU_i\right\rfloor +1}$ y en caso contrario hacer
$X_i^{\ast}=2\overline{X
}-X_{\left\lfloor nU_i\right\rfloor +1}$

La utilización de $F_n^{sim}\left( x \right)$ en lugar de $F_n\left(
x \right)$ altera las propiedades conocidas de la distribución
(empírica) de la que se remuestrea en el bootstrap uniforme. Así, en
primer lugar, $F_n^{sim}\left( x \right)$ es simétrica (como se desea)
con lo cual todos los momentos impares de esta distribución con respecto
a $\bar{X}$ son cero. 
En particular la media de $F_n^{sim}\left(x \right)$ es
$$\begin{aligned}
\int x~dF_n^{sim}\left( x \right) &= \frac{1}{2n}\sum_{i=1}^{2n}Y_i=\frac{
1}{2n}\left[ \sum_{i=1}^{n}X_i+\sum_{i=1}^{n}\left( 2\bar{X}
-X_i \right) \right] \\
&= \frac{1}{2n}\left( n\bar{X}+2n\bar{X}-n\bar{X} \right) =
\bar{X}.\end{aligned}$$
También se conservan los momentos centrales de orden par:

$$\begin{aligned}
\int \left( x-\bar{X} \right)^{2k}~dF_n^{sim}\left( x \right) &= \frac{
1}{2n}\sum_{i=1}^{2n}\left( Y_i-\bar{X} \right)^{2k} \\
&= \frac{1}{2n}\left[ \sum_{i=1}^{n}\left( X_i-\bar{X} \right)
^{2k}+\sum_{i=1}^{n}\left[ \left( 2\bar{X}-X_i \right) -\bar{X}
\right]^{2k}\right] \\
&= \frac{1}{2n}\left[ \sum_{i=1}^{n}\left( X_i-\bar{X} \right)
^{2k}+\sum_{i=1}^{n}\left( \bar{X}-X_i \right)^{2k}\right] \\
&= \frac{1}{n}\sum_{i=1}^{n}\left( X_i-\bar{X} \right)^{2k}.
\end{aligned}$$

En particular, la varianza de $F_n^{sim}\left( x \right)$ coincide con
la de $F_n\left( x \right)$, que es $S_n^2$.

En general, cuando la distribución de partida es simétrica, es más
adecuado utilizar el bootstrap simetrizado que el bootstrap uniforme.
Aún así, cuando se realiza inferencia sobre algún estadístico (como
$\sqrt{n}(\bar{X}-\mu)/\sigma$) cuya distribución
asintótica ya es simétrica (como la normal), la aproximación bootstrap
uniforme para distribuciones de partida simétricas, ya es especialmente
buena y, por tanto, la ganancia del bootstrap simétrizado aporta una
mejora difícil de detectar en la práctica. Ese no es el caso de otros
estadísticos (como los asociados a inferencia sobre la varianza) con
distribución más alejada de la simetría.


\BeginKnitrBlock{exercise}\iffalse{-91-73-110-102-101-114-101-110-99-105-97-32-115-111-98-114-101-32-108-97-32-109-101-100-105-97-32-99-111-110-32-118-97-114-105-97-110-122-97-32-99-111-110-111-99-105-100-97-32-101-109-112-108-101-97-110-100-111-32-98-111-111-116-115-116-114-97-112-32-115-105-109-101-116-114-105-122-97-100-111-93-}\fi{}
<span class="exercise" id="exr:media-dt-conocida-sim"><strong>(\#exr:media-dt-conocida-sim)  \iffalse (Inferencia sobre la media con varianza conocida empleando bootstrap simetrizado) \fi{} </strong></span>
Modificar adecuadamente el código del Ejemplo \@ref(exm:media-dt-conocida), para
implementar un método bootstrap simetrizado, con el objeto de calcular
un intervalo de confianza para la media con desviación típica conocida.
Qué diferencias se observan entre los intervalos obtenidos por el
bootstrap uniforme y por el simetrizado?

\EndKnitrBlock{exercise}


## Bootstrap suavizado {#modunif-boot-suav}

Cuando la distribución poblacional, $F$, es continua es lógico
incorporar dicha información al bootstrap. Eso significa que la función
de distribución tiene una función de densidad asociada, relacionadas
mediante la expresión: $f\left( x \right) =F^{\prime}\left(
x \right)$. Para ello, debemos utilizar un método bootstrap que
remuestree de un universo bootstrap continuo. En otras palabras debemos
utilizar un estimador de la función de densidad y remuestrear de él.

Pasamos a considerar brevemente el problema de estimar, no
paramétricamente, la función de densidad, $f$, de una población, a
partir de una muestra, $\left( X_1,X_2,\ldots ,X_n \right)$,
procedente de la misma. En ese contexto es bien conocido el método
histograma (basado en el cual sería posible idear un método bootstrap)
aunque es más recomendable utilizar el estimador tipo núcleo propuesto
por Parzen (1962) y Rosenblatt (1956), que viene dado por

$$\hat{f}_{h}\left( x \right) =\frac{1}{nh}\sum_{i=1}^{n}K\left( \frac{x-X_i}{
h} \right) =\frac{1}{n}\sum_{i=1}^{n}K_{h}\left( x-X_i \right),$$

donde 
$$K_{h}\left( u \right) =\frac{1}{h}K\left( \frac{u}{h} \right),$$
$K$ es una función núcleo (normalmente una densidad simétrica en torno
al cero) y $h>0$ es una parámetro de suavizado, llamado ventana, que
regula el tamaño del entorno que se usa para llevar a cabo la
estimación. Este estimador generaliza el bien conocido histograma y, más
concretamente, su versión histograma móvil. Así, eligiendo como función
$K$ la densidad de una $\mathcal{U}\left( -1,1 \right)$, el estimador de
Parzen-Rosenblatt resulta:
$$\begin{aligned}
\frac{1}{nh}\sum_{i=1}^{n}\frac{1}{2}\mathbf{1}\left\{ \frac{x-X_i}{h}\in
\left( -1,1 \right) \right\} &= \frac{1}{2nh}\sum_{i=1}^{n}\mathbf{1}\left\{
X_i\in \left( x-h,x+h \right) \right\} \\
&= \frac{\#\left\{ X_i\in \left( x-h,x+h \right) \right\} }{2nh},
\end{aligned}$$
que no es más que la frecuencia relativa de datos $X_i$ en el
intervalo $\left( x-h,x+h \right)$ dividida entre la longitud del
intervalo en cuestión ($2h$).

Es habitual exigir que la función núcleo $K$ sea no negativa y su
integral sea uno:
$$K\left( u \right) \geq 0,~\forall u,~\int_{-\infty }^{\infty }K\left(u \right) du=1.$$
Además también es frecuente exigir que $K$ sea una
función simétrica ($K\left( -u \right) =K\left( u \right)$).

Aunque la elección de la función $K$ no tiene gran impacto en las
propiedades del estimador (salvo sus condiciones de regularidad:
continuidad, diferenciabilidad, etc.) la elección del parámetro de
suavizado sí es muy importante para una correcta estimación. En otras
palabras, el tamaño del entorno usado para la estimación no paramétrica
debe ser adecuado (ni demasiado grande ni demasiado pequeño). 

En `R` podemos emplear la función `density()` del paquete base para obtener
una estimación tipo núcleo de la densidad (con la ventana determinada
por el parámetro `bw`), aunque podríamos emplear implementaciones de otros
paquetes (en la Sección \@ref(npden-r) se incluyen más detalles).
Por ejemplo, considerando el conjunto de datos `precip` (que contiene el promedio de precipitación, en pulgadas de lluvia, de 70 ciudades de Estados Unidos), podríamos utilizar el siguiente código [Figura \@ref(fig:density)]:


```r
x <- precip
npden <- density(x)
# npden <- density(x, bw = "SJ")

# plot(npden)
bandwidth <- npden$bw
hist(x, freq = FALSE, main = "Kernel density estimation", 
     xlab = paste("Bandwidth =", formatC(bandwidth)), lty = 2, 
     border = "darkgray", xlim = c(0, 80), ylim = c(0, 0.08))
lines(npden, lwd = 2)
rug(x, col = "darkgray")
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{03-mod_boot_unif_files/figure-latex/density-1} 

}

\caption{Estimación tipo núcleo de la densidad de `precip`. }(\#fig:density)
\end{figure}

La sensibilidad del estimador tipo núcleo al parámetro de suavizado puede
observarse ejecutando el siguiente código (ver Figura \@ref(fig:bandwidth-movie), [bandwidth-movie.gif](./bandwidth-movie.gif)):

```r
bws <- 2^seq(log2(bandwidth * 0.01), log2(bandwidth * 20), len = 50)
bws <- c(bws, rev(bws))
for (bw in bws)
  plot(density(x, bw = bw) , main = "Kernel density estimation", 
         xlab = paste("Bandwidth =", formatC(bw)), 
         xlim = c(0, 80), ylim = c(0, 0.08))
```


\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{03-mod_boot_unif_files/figure-latex/bandwidth-movie-1} 

}

\caption{Efecto de cambio en la ventana en la estimación tipo núcleo de la densidad.}(\#fig:bandwidth-movie)
\end{figure}


La función de distribución asociada al estimador tipo núcleo de la
función de densidad viene dada por
$$\begin{aligned}
\hat{F}_{h}\left( x \right) &= \int_{-\infty }^{x}\hat{f}_{h}\left( y \right) dy
=\int_{-\infty }^{x}\frac{1}{n}\sum_{i=1}^{n}\frac{1}{h}
K\left( \frac{y-X_i}{h} \right) dy \\
&= \frac{1}{nh}\sum_{i=1}^{n}\int_{-\infty }^{x}
K\left( \frac{y-X_i}{h} \right) dy \\
&= \frac{1}{n}\sum_{i=1}^{n}\int_{-\infty }^{\frac{x-X_i}{h}}K\left( u \right) du
=\frac{1}{n}\sum_{i=1}^{n}\mathbb{K}\left( \frac{x-X_i}{h} \right)
\end{aligned}$$
donde $\mathbb{K}$ es la función de distribución
asociada al núcleo $K$, es decir
$$\mathbb{K}\left( t \right) =\int_{-\infty }^{t}K\left(
u \right) du.$$

Por ejemplo, en el caso de del conjunto de datos de precipitaciones, el siguiente código compara la estimación tipo núcleo de la distribución con la empírica [Figura \@ref(fig:pnp)]:


```r
Fn <- ecdf(precip)
curve(Fn, xlim = c(0, 75), ylab = "F(x)", type = "s")
Fnp <- function(x) sapply(x, function(y) mean(pnorm(y, precip, bandwidth)))
curve(Fnp, lty = 2, add = TRUE) 
legend("bottomright", legend = c("Empírica", "Tipo núcleo"), lty = 1:2)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{03-mod_boot_unif_files/figure-latex/pnp-1} 

}

\caption{Estimación empírica y tipo núcleo de la función de distribución de `precip`. }(\#fig:pnp)
\end{figure}


El método bootstrap suavizado procede de la siguiente forma:

1. A partir de la muestra $\left( X_1,X_2,\ldots ,X_n \right)$ y
utilizando un valor $h>0$ como parámetro de suavizado, se calcula el
estimador de Parzen-Rosenblatt $\hat{f}_{h}$

2. Se arrojan remuestras bootstrap $\mathbf{X}^{\ast}=\left(
X_1^{\ast},X_2^{\ast},\ldots ,X_n^{\ast} \right)$ a partir de
la densidad $\hat{f}_{h}$

3. Calcular
$R^{\ast}=R\left( \mathbf{X}^{\ast},\hat{F}_{h} \right)$

4. Repetir $B$ veces los pasos 2-3 para obtener las réplicas bootstrap
$R^{\ast (1)}$, $\ldots$, $R^{\ast (B)}$

Para llevar a cabo un bootstrap que remuestree a partir del estimador
$\hat{f}_{h}\left( x \right)$ es útil pensar en dicho estimador como una
combinación lineal convexa de funciones de densidad, $K_{h}\left(
x-X_i \right)$, cada una con coeficiente $\frac{1}{n}$ en dicha
combinación lineal. Gracias a esa representación podemos simular
valores,
$X^{\ast}$, procedentes de $\hat{f}_{h}\left( x \right)$ en dos pasos (empleando el denominado método de composición; ver p.e. Fernández-Casal y Cao, 2020, [Sección 5.4](https://rubenfcasal.github.io/simbook/m%C3%A9todo-de-composici%C3%B3n.html)).

En un primer paso elegiremos (aleatoriamente y con equiprobabilidad)
cuál de los índices $i\in \left\{ 1,\ldots ,n\right\}$ vamos a
considerar y en un segundo paso simularemos $X^{\ast}$ a partir de la
densidad $K_{h}\left(
\cdot -X_i \right)$. Esta última fase puede relacionarse fácilmente
con la simulación de un valor, $V$, con densidad $K$, sin más que hacer
$X_i+hV$. Así, el paso 2 del algoritmo previo puede llevarse a cabo
mediante el siguiente procedimiento:

2. Para cada $i=1,\ldots ,n$ arrojar $U_i\sim \mathcal{U}\left( 0,1 \right)$ y
$V_i$ con densidad $K$ y hacer $X_i^{\ast}=X_{\left\lfloor
nU_i\right\rfloor +1}+hV_i$

La equivalencia de ambas presentaciones del paso 2 viene dada por el
siguiente razonamiento. Denotando $U\sim \mathcal{U}\left( 0,1 \right)$, 
$I=\left\lfloor nU\right\rfloor +1$ y $V\sim K$, independiente de $U$, se
tiene:

$$\begin{aligned}
P^{\ast}\left( X^{\ast}\leq x \right) &= \sum_{i=1}^{n}P^{\ast}\left(
\left. X^{\ast}\leq x\right\vert _{I=i} \right) P^{\ast}\left( I=i \right) \\
&= \sum_{i=1}^{n}P^{\ast}\left( \left. X_i+hV\leq x\right\vert
_{I=i} \right) P^{\ast}\left( I=i \right) \\
&= \sum_{i=1}^{n}P^{\ast}\left( \left. V\leq \frac{x-X_i}{h}
\right\vert _{X_i} \right) \frac{1}{n}=\frac{1}{n}\sum_{i=1}^{n}\mathbb{K}
\left( \frac{x-X_i}{h} \right),
\end{aligned}$$

cuya función de densidad es, como ya sabemos, $\hat{f}_{h}\left(
x \right)$. Esto justifica la presentación alternativa del paso 2, de
forma que el bootstrap suavizado puede pensarse, a partir del bootstrap
uniforme ($X_i^{\ast}=X_{\left\lfloor nU_i\right\rfloor +1}$)
añadiendo al mismo una perturbación ($hV_i$) cuya magnitud viene dada
por el parámetro de suavizado ($h$) y cuya forma imita a la de una
variable aleatoria ($V_i$) con densidad $K$.

Por ejemplo, la función `density()` emplea por defecto un núcleo
gaussiano, y como se muestra en la ayuda de esta función,
podemos emplear un código como el siguiente para obtener
`nsim` simulaciones (ver Figura \@ref(fig:density-sim)):

```r
## simulation from a density() fit:
# a kernel density fit is an equally-weighted mixture.
nsim <- 1e6
set.seed(1)
# x_boot <- sample(x, nsim, replace = TRUE)
# x_boot <- x_boot + bandwidth * rnorm(nsim)
x_boot <- rnorm(nsim, sample(x, nsim, replace = TRUE), bandwidth)

plot(npden, main = "")
lines(density(x_boot), col = "blue", lwd = 2, lty = 2)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{03-mod_boot_unif_files/figure-latex/density-sim-1} 

}

\caption{Estimaciónes tipo núcleo de las densidades de `precip` y de una simulación.}(\#fig:density-sim)
\end{figure}

Es fácil percatarse de que los posibles valores que puede tomar una
observación $X_i^{\ast}$ de cada remuestra bootstrap son infinitos,
pues la variable $V$ puede tomar infinitos posibles valores, según la
densidad de probabilidad $K$. Esto significa que la distribución en el
remuestreo de la remuestra bootstrap, $\mathbf{X}^{\ast}$, es
mucho más complicada que para el bootstrap uniforme. En particular no es
discreta y por tanto no puede caracterizarse a partir de vectores de
remuestreo sobre la muestra original. Un problema importante es la
elección del parámetro de suavizado, $h$, en este procedimiento de
remuestreo. En la práctica es razonable elegir $h$ como un valor
bastante pequeño, en relación con la desviación típica de la muestra. Es
fácil observar que en el caso extremo $h=0$ este método de remuestreo se
reduce al bootstrap uniforme.


\BeginKnitrBlock{example}\iffalse{-91-73-110-102-101-114-101-110-99-105-97-32-115-111-98-114-101-32-108-97-32-109-101-100-105-97-32-99-111-110-32-118-97-114-105-97-110-122-97-32-99-111-110-111-99-105-100-97-44-32-99-111-110-116-105-110-117-97-99-105-243-110-93-}\fi{}
<span class="example" id="exm:media-dt-conocida-suav"><strong>(\#exm:media-dt-conocida-suav)  \iffalse (Inferencia sobre la media con varianza conocida, continuación) \fi{} </strong></span>
Continuando con el ejemplo de tiempo de vida de microorganismos,
podemos modificar fácilmente el código mostrado en el Ejemplo \@ref(exm:media-dt-conocida), para implementar bootstrap suavizado 
con función núcleo gaussiana, para calcular un intervalo de confianza 
para la media poblacional con desviación típica conocida:

\EndKnitrBlock{example}



```r
muestra <- c(0.143, 0.182, 0.256, 0.26, 0.27, 0.437, 0.509, 
             0.611, 0.712, 1.04, 1.09, 1.15, 1.46, 1.88, 2.08)
n <- length(muestra)
sigma <- 0.6

alfa <- 0.05
x_barra <- mean(muestra)

# Remuestreo
set.seed(1)
B <- 1000
# h <- 1e-08
h <- bw.SJ(muestra)/2
estadistico_boot <- numeric(B)
for (k in 1:B) {
    # remuestra <- sample(muestra, n, replace = TRUE)
    # remuestrasu <- remuestra + h * rnorm(n, 0, 1)
    remuestrasu <- rnorm(n, sample(muestra, n, replace = TRUE), h)
    x_barra_boot <- mean(remuestrasu)
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
## 0.4668897 1.0798549
```

Con el paquete `boot`, la recomendación es implementarlo como
un bootstrap paramétrico:

```r
library(boot)
ran.gen.smooth <- function(data, mle) {
    # Función para generar muestras aleatorias mediante 
    # bootstrap suavizado con función núcleo gaussiana,
    # mle contendrá la ventana.
    n <- length(data)
    h <- mle
    out <- rnorm(n, sample(data, n, replace = TRUE), h)
    out
}

statistic <- function(data){
    c(mean(data), sigma^2/length(data))
}

set.seed(1)
res.boot <- boot(muestra, statistic, R = B, sim = "parametric",
                 ran.gen = ran.gen.smooth, mle = h)

boot.ci(res.boot, type = "stud")
```

```
## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
## Based on 1000 bootstrap replicates
## 
## CALL : 
## boot.ci(boot.out = res.boot, type = "stud")
## 
## Intervals : 
## Level    Studentized     
## 95%   ( 0.4664,  1.0830 )  
## Calculations and Intervals on Original Scale
```

## Bootstrap ponderado y bootstrap sesgado

Mediante el nombre bootstrap ponderado se incluyen todos aquellos
métodos de remuestreo bootstrap en los que la distribución de la que se
remuestrea es discreta y asigna probabilidades sólo a los datos de la
muestra:
$$\hat{F}\left( X_i \right) -\hat{F}\left( X_i^{-} \right) = p_i
\text{, para }i=1,\ldots, n$$
siendo $p_i\geq 0$ y $\sum_{i=1}^{n}p_i=1$. En el caso particular $p_i=
\frac{1}{n}$ para todo $i=1,\ldots ,n$, se tiene el bootstrap uniforme.
Veremos más adelante casos particulares de métodos bootstrap ponderados
en el contexto de datos censurados y también para datos dependientes.

El bootstrap ponderado da lugar al bootstrap sesgado cuando los pesos,
$p_i$, se eligen de forma que el vector $\mathbf{p}$ minimice
la distancia al vector de pesos del bootstrap uniforme
$\left( \frac{1}{n},\ldots ,\frac{1}{n} \right)$, 
sujeto a una serie de restricciones inherentes al problema en estudio. 
Este método fue propuesto por Hall (1998).


## Deficiencias del bootstrap uniforme {#deficien-unif}

Supongamos un contexto paramétrico en el que la distribución
poblacional, $F$, es la $\mathcal{U}\left( 0,\theta \right)$. Nuestro interés
será hacer inferencia acerca del parámetro $\theta$, para lo cual, dada
una muestra observada, $\mathbf{X}=\left( X_1,X_2,\ldots
,X_n \right)$, consideraremos el estimador máximo verosímil en este
contexto: $\hat{\theta}=X_{(n)}$. Para realizar dicha
inferencia estaremos interesados en aproximar la distribución de
$R\left( \mathbf{X},F \right) =\hat{\theta}-\theta$

La función de distribución en el muestreo, $G\left( x \right)$, de
$\hat{\theta}$ puede calcularse de forma sencilla:$$\begin{aligned}
G\left( x \right) &= P\left( \hat{\theta}\leq x \right) =P\left( X_{\left(
n \right)}\leq x \right) =P\left( X_i\leq x\,,\forall i\in \left\{ 1,\ldots
n\right\} \right) \\
&= \prod_{i=1}^{n}P\left( X_i\leq x \right) =F\left( x \right)^{n}=\left( 
\frac{x}{\theta } \right)^{n},\text{ si }x\in \left[ 0,\theta \right]\end{aligned}$$

con lo cual su función de densidad viene dada por
$$g\left( x \right) =\frac{n}{\theta }\left( \frac{x}{\theta } \right)^{n-1},
\text{ si }x\in \left[ 0,\theta \right] .$$
Tomando, por ejemplo, $\theta =1$ y $n=50$, esta función de densidad resulta
[Figura \@ref(fig:den-max)]:

```r
theta <- 1
n <- 50
curve(n/theta * (x/theta)^(n - 1), 0, theta, ylab = "Density")
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{03-mod_boot_unif_files/figure-latex/den-max-1} 

}

\caption{Función de densidad del máximo de una muestra procedente de una uniforme.}(\#fig:den-max)
\end{figure}

Como consecuencia podemos hallar fácilmente el sesgo del estimador
$\hat{\theta}$, ya que
$$E\left( \hat{\theta} \right) =\int_{0}^{\theta }x\frac{n}{\theta }\left( 
\frac{x}{\theta } \right)^{n-1}dx=\left[ \frac{n}{n+1}\frac{x^{n+1}}{\theta
^{n}}\right] _{x=0}^{x=\theta }=\frac{n}{n+1}\theta ,$$
con lo cual
$$Sesgo\left( \hat{\theta} \right) =E\left( \hat{\theta} \right)
-\theta = -\frac{\theta }{n+1}.$$
Se ve claramente que $\hat{\theta}$
es un estimador sesgado de $\theta$, puesto que se tiene que
$\hat{\theta}\leq
\theta$ con probabilidad 1.

Si deseamos aproximar mediante bootstrap la distribución en el muestreo
de $\hat{\theta}$ (o la de $R$) y utilizamos un bootstrap uniforme
(naïve), la versión bootstrap del estimador resulta ser
$\hat{\theta}^{\ast
}=X_{(n)}^{\ast}$, siendo
$\mathbf{X}^{\ast}=\left(
X_1^{\ast}\text{, }X_2^{\ast}\text{, }\ldots \text{, }X_n^{\ast
} \right)$ una remuestra bootstrap obtenida a partir de la distribución
empírica $F_n$. La distribución en el remuestreo de $\hat{\theta}
^{\ast}\,$ resulta un poco más complicada pues es discreta y sólo puede
tomar cualquiera de los valores de la muestra.

Suponiendo que no hay empates en las observaciones de la muestra, es
fácil darse cuenta de que
$$P^{\ast}\left( \hat{\theta}^{\ast}\leq X_{(j)} \right)
=P^{\ast}\left( X_{(n)}^{\ast}\leq X_{(j)
} \right) =P^{\ast}\left( X_i^{\ast}\leq X_{(j)}\,,
 1 \leq i \leq n \right) =\left( \frac{j}{n} \right)^{n}$$
y, por tanto, su masa de probabilidad viene dada por
$$P^{\ast}\left( \hat{\theta}^{\ast}=X_{(j)} \right) =\left( 
\frac{j}{n} \right)^{n}-\left( \frac{j-1}{n} \right)^{n}\text{, }j=1,\ldots,n.$$

En particular,
$$P^{\ast}\left( \hat{\theta}^{\ast}=X_{(n)} \right) =1-\left( 1-
\frac{1}{n} \right)^{n}\rightarrow 1-\frac{1}{e}\simeq 0.6321,$$
con lo cual la distribución en remuestreo de $R^{\ast}=R\left( 
\mathbf{X}^{\ast},F_n \right) =\hat{\theta}^{\ast}-X_{\left(
n \right)}$ tiene un átomo de probabilidad en el valor $0$ cuya
probabilidad tiende a $1-\frac{1}{e}$ cuando el tamaño muestral tiende a
infinito, es decir

$$\lim_{n\rightarrow \infty }P^{\ast}\left( R^{\ast}=0 \right) =1-\frac{1}{e},$$

cosa que no ocurre con la distribución en el muestreo de $R$, que es continua
con densidad:
$$g_R\left( x \right) =\frac{n}{\theta }\left( \frac{x + \theta}{\theta } \right)^{n-1},
\text{ si }x\in \left[ -\theta, 0\right].$$ 
De esta forma vemos que el bootstrap uniforme (no paramétrico) es inconsistente.

\BeginKnitrBlock{example}\iffalse{-91-73-110-102-101-114-101-110-99-105-97-32-115-111-98-114-101-32-101-108-32-109-225-120-105-109-111-32-100-101-32-117-110-97-32-100-105-115-116-114-105-98-117-99-105-243-110-32-117-110-105-102-111-114-109-101-93-}\fi{}
<span class="example" id="exm:boot-maximo-uniforme"><strong>(\#exm:boot-maximo-uniforme)  \iffalse (Inferencia sobre el máximo de una distribución uniforme) \fi{} </strong></span>
\EndKnitrBlock{example}

El siguiente código implementa el método
bootstrap uniforme (también llamado naïve) para aproximar la
distribución del estadístico $R=\hat{\theta}-\theta$, para una muestra
de tamaño $n=50$, proveniente de una población con distribución
$\mathcal{U}\left( 0,1\right)$ [Figura \@ref(fig:boot-uniforme-maximo)]:


```r
theta <- 1
n <- 50
set.seed(1)
muestra <- runif(50) * theta
theta_est <- max(muestra)
# Remuestreo
B <- 2000
maximo <- numeric(B)
estadistico <- numeric(B)
for (k in 1:B) {
    remuestra <- sample(muestra, n, replace = TRUE)
    maximo[k] <- max(remuestra)
    estadistico[k] <- maximo[k] - theta_est
}
# Distribución estadístico
xlim <- c(-theta/2, 0) # c(-theta, 0)
hist(estadistico, freq = FALSE, main = "", lty = 2, 
     border = "darkgray", xlim = xlim)
lines(density(estadistico))
rug(estadistico, col = "darkgray")
curve(n/theta * ((x + theta)/theta)^(n - 1), col = "blue", lty = 2, lwd = 2, add = TRUE)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{03-mod_boot_unif_files/figure-latex/boot-uniforme-maximo-1} 

}

\caption{Distribución de las réplicas bootstrap (uniforme) del estadístico y distribución poblacional.}(\#fig:boot-uniforme-maximo)
\end{figure}

### Ejemplo (método alternativo)

En este contexto, al conocer la familia paramétrica ($\mathcal{U}\left( 0,\theta 
\right)$) a la cual pertenece la distribución de la población de
partida, lo natural sería utilizar un bootstrap paramétrico, consistente
en obtener las remuestras bootstrap a partir de una distribución
uniforme con parámetro estimado:
$$\mathbf{X}^{\ast}=\left( X_1^{\ast}\text{, }X_2^{\ast}\text{, 
}\ldots \text{, }X_n^{\ast} \right), \text{ con } X_i^{\ast} \sim \mathcal{U}\left( 0,\hat{\theta}\right).$$
En estas circunstancias es muy sencillo obtener la distribución en el
remuestreo de $\hat{\theta}^{\ast}$, ya que su deducción es totalmente
paralela a la de la distribución en el muestreo de $\hat{\theta}$. Así,
la función de densidad de $\hat{\theta}^{\ast}$ es
$$\hat{g}\left( x \right) =\frac{n}{\hat{\theta}}\left( \frac{x}{\hat{\theta}}
 \right)^{n-1},\text{ si }x\in \left[ 0,\hat{\theta}\right] .$$

Con lo cual, al utilizar un bootstrap paramétrico, la distribución en el
remuestreo de $R^{\ast}=R\left( \mathbf{X}^{\ast},F_{\hat{
\theta}} \right) =\hat{\theta}^{\ast}-\hat{\theta}$ imita a la
distribución en muestreo de
$R=R\left( \mathbf{X},F \right) =\hat{\theta}-\theta$.


\BeginKnitrBlock{example}\iffalse{-91-73-110-102-101-114-101-110-99-105-97-32-115-111-98-114-101-32-101-108-32-109-225-120-105-109-111-32-100-101-32-117-110-97-32-100-105-115-116-114-105-98-117-99-105-243-110-32-117-110-105-102-111-114-109-101-44-32-99-111-110-116-105-110-117-97-99-105-243-110-93-}\fi{}
<span class="example" id="exm:boot-maximo-parametrico"><strong>(\#exm:boot-maximo-parametrico)  \iffalse (Inferencia sobre el máximo de una distribución uniforme, continuación) \fi{} </strong></span>
\EndKnitrBlock{example}

Para emplear el bootstrap paramétrico (que remuestrea de una distribución
uniforme con parámetro estimado) podríamos emplear un código muy similar al 
del Ejemplo \@ref(exm:boot-maximo-uniforme) [Figura \@ref(fig:boot-parametrico-maximo)]:


```r
# Remuestreo
B <- 2000
maximo <- numeric(B)
estadistico <- numeric(B)
for (k in 1:B) {
    remuestra <- runif(n) * theta_est
    maximo[k] <- max(remuestra)
    estadistico[k] <- maximo[k] - theta_est
}
# Distribución estadístico
xlim <- c(-theta/2, 0) # c(-theta, 0)
hist(estadistico, freq = FALSE, main = "", lty = 2, 
     border = "darkgray", xlim = xlim)
lines(density(estadistico))
rug(estadistico, col = "darkgray")
curve(n/theta * ((x + theta)/theta)^(n - 1), col = "blue", lty = 2, lwd = 2, add = TRUE)
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{03-mod_boot_unif_files/figure-latex/boot-parametrico-maximo-1} 

}

\caption{Distribución bootstrap paramétrica y distribución poblacional.}(\#fig:boot-parametrico-maximo)
\end{figure}

<!-- 
Ejercicio con el máximo de otra distribución?
bootstrap suavizado?
-->

## Validez de la aproximación Bootstrap

Trataremos ahora de dar una justificación teórica del buen
funcionamiento del bootstrap uniforme. Para ello, por simplicidad, nos
centraremos en el problema de aproximar la distribución en el muestreo
del estadístico
$$R=R\left( \mathbf{X},F \right) =\sqrt{n}\frac{\bar{X}-\mu }{\sigma },$$
donde $\mathbf{X}=\left( X_1,X_2,\ldots ,X_n \right)$ es una
m.a.s. procedente de una distribución $F$, con media $\mu$ y desviación
típica $\sigma$. Sabemos que, bajo ciertas condiciones, el teorema
central del límite permite obtener la distribución asintótica de $R$,
que es una $\mathcal{N}\left( 0,1 \right)$, es decir
$$\lim_{n\rightarrow \infty }P\left( R\leq u \right) =\Phi \left( u \right),
\quad\forall u\in \mathbb{R},$$
siendo $\Phi$ la función de distribución de una normal estándar, cuya
función de densidad denotaremos por $\phi$.

Para ver cómo de buena es la aproximación por normal del estadístico
$R$, debemos razonar cómo de rápida es la convergencia en el límite
anteriormente expuesto. La respuesta a esa pregunta viene dada por el
Teorema de Cramer que usa los llamados desarrollos de Edgeworth de un
estadístico para aproximarlo por una suma de términos, el primero es la
función de distribución normal estándar y los siguientes irán tendiendo
a cero sucesivamente más rápido cuando el tamaño muestral tiende a
infinito. Enunciemos ese resultado.

\BeginKnitrBlock{theorem}\iffalse{-91-67-114-97-109-101-114-93-}\fi{}
<span class="theorem" id="thm:aprox-cramer"><strong>(\#thm:aprox-cramer)  \iffalse (Cramer) \fi{} </strong></span><br> \vspace{0.5cm}

Consideremos variables aleatorias
$X_1,X_2,\ldots ,X_n,\ldots$ independientes e idénticamente
distribuidas procedentes de una distribución $F$, con media $\mu$ y
desviación típica $\sigma$. Supongamos que existe cierto $j$, natural,
para el cual $E\left( \left\vert X\right\vert^{j+2} \right) <\infty
\,$, y que $\lim_{\left\vert t\right\vert \rightarrow \infty }\left\vert
\alpha \left( t \right) \right\vert <1$, siendo $\alpha \left( t \right)
=E\left( e^{itX} \right)$ la función característica de la población.
Entonces:
$$\begin{aligned}
P\left( R\leq u \right) &=P\left( \sqrt{n}\frac{\bar{X}-\mu }{\sigma }
\leq u \right) \\
&= \Phi \left( u \right) +n^{-\frac{1}{2}}p_1\left( u \right) \phi \left(
u \right) +\cdots +n^{-\frac{j-1}{2}}p_{j-1}\left( u \right) \phi \left(
u \right) +O\left( n^{-\frac{j}{2}} \right),
\end{aligned}$$
siendo los $p_i\left( u \right)$ polinomios de grado $3i-1$ cuyos coeficientes
dependen de los momentos de $X$ de orden menor o igual que $i+2$. 
En particular
$$\begin{aligned}
p_1\left( u \right) &= -\frac{1}{6}\frac{k_3}{\sigma^{3}}\left(
u^2-1 \right), \\
p_2\left( u \right) &= -u\left[ \frac{1}{24}\frac{k_4}{\sigma^{4}}\left(
u^2-3 \right) +\frac{1}{72}\left( \frac{k_3}{\sigma^{3}} \right)
^2\left( u^{4}-10u^2+15 \right) \right] ,
\end{aligned}$$

siendo $k_j$ el $j$-ésimo cumulante, es decir el términos que acompaña
a $\frac{\left( it \right)^{j}}{j!}$ en el desarrollo en serie del
logaritmo de la función característica:
$$\log \alpha \left( t \right) =\sum_{j=1}^{\infty }k_j\frac{\left( it \right)
^{j}}{j!}.$$
  
Además dichos polinomios tienen paridad alternada, es
decir, $p_1$ es simétrico, $p_2$ es antisimétrico, $p_3$ es
simétrico, y así sucesivamente:
$$p_1\left( -u \right) = p_1\left( u \right),\quad p_2\left( -u \right)
= -p_2\left( u \right),\quad p_3\left( -u \right) = p_3\left( u \right)
,\cdots$$
\EndKnitrBlock{theorem}

Existen ecuaciones que relacionan todos los cumulantes hasta cierto
orden con todos los momentos poblacionales hasta ese mismo orden. Dichas
ecuaciones permiten expresar los cumulantes en función de los momentos y
viceversa.

Como consecuencia de este resultado teórico, el grado de aproximación
entre la distribución de $R$ y la normal estándar límite es
$O (n^{-\frac{1}{2}})$. Sin embargo, puede razonarse
fácilmente que este orden de aproximación mejorará cuando utilizamos el
bootstrap uniforme, en lugar de la normal estándar, para aproximar la
distribución de $R$. Un desarrollo de Edgeworth para la distribución en
el remuestreo de $R^{\ast}$ permite obtener la siguiente expresión:
$$\begin{aligned}
P^{\ast}\left( R^{\ast}\leq u \right) &= \Phi \left( u \right) +n^{-\frac{1}{2}
}\hat{p}_1\left( u \right) \phi \left( u \right) +\cdots +n^{-\frac{j-1}{2}}
\hat{p}_{j-1}\left( u \right) \phi \left( u \right) \\
&+ O_{P}\left( n^{-\frac{j}{2}} \right),
\end{aligned}$$

donde los polinomios $\hat{p}_i\left( u \right)$ tienen la misma
estructura que los $p_i\left( u \right)$ pero reemplazando los
cumulantes teóricos por los empíricos y la desviación típica teórica por
la empírica. Así pues el grado de aproximación entre cada polinomio
$\hat{p}_i( u )$ y su análogo teórico $p_i( u )$ es
$\hat{p}_i( u ) -p_i( u ) = O_{P}( n^{-\frac{1}{2}} )$. 
Como consecuencia, puede obtenerse el orden de aproximación entre la distribución 
en el muestreo de $R$ y la distribución en el remuestreo de $R^{\ast}$:
$$\begin{aligned}
P\left( R\leq u \right) -P^{\ast}\left( R^{\ast}\leq u \right) &=  n^{-\frac{1}{
2}}\left[ p_1\left( u \right) -\hat{p}_1\left( u \right) \right] \phi
\left( u \right) +O_{P}\left( n^{-1} \right) \\
&=  O_{P}\left( n^{-1} \right),\end{aligned}$$que es mejor que el orden
de aproximación de la normal estándar límite. Dichos órdenes pueden
resumirse en la siguiente tabla.

  Aproximación      Orden                       
  --------------    ----------------------------------
  Normal límite     $O\left( n^{-\frac{1}{2}} \right)$   
  Boot. uniforme    $O_{P}\left( n^{-1} \right)$

Usando razonamiento similares pueden encontrarse los órdenes de
aproximación, tanto de la normal límite, como del bootstrap uniforme y
del bootstrap simetrizado, cuando la distribucional de partida es
simétrica. En ese caso, $p_1\left( u \right) =0$, ya que $k_3$ es
cero debido a la simetría de la distribución poblacional. Sin embargo
$\hat{p}_1\left( u \right)$ no es cero cuando se usa el bootstrap
uniforme, aunque sí lo es en el caso del bootstrap simetrizado. La
siguiente tabla recoge los órdenes de las distintas aproximaciones.

  Aproximación       Orden                       
  --------------     ----------------------------------
  Normal límite      $O\left( n^{-1} \right)$   
  Boot. uniforme     $O_{P}\left( n^{-1} \right)$
  Boot. simetrizado  $O_{P}\left( n^{-\frac{3}{2}} \right)$

El siguiente resultado permite generalizar los desarrollos de Edgeworth
(Teorema \@ref(thm:aprox-cramer))
a otros estadísticos (estandarizados o studentizados) obtenidos para
otros estimadores arbitrarios, $\hat{\theta}$, no necesariamente iguales
a la media muestral.

\BeginKnitrBlock{theorem}\iffalse{-91-66-104-97-116-116-97-99-104-97-114-121-97-45-71-104-111-115-104-93-}\fi{}
<span class="theorem" id="thm:aprox-bhat-gho"><strong>(\#thm:aprox-bhat-gho)  \iffalse (Bhattacharya-Ghosh) \fi{} </strong></span><br> \vspace{0.5cm}

Consideremos variables aleatorias
$X_1,X_2,\ldots ,X_n,\ldots$ independientes e idénticamente
distribuidas procedentes de una distribución $F$. Sea $\theta =\theta
\left( F \right)$ un parámetro de dicha distribución y $\hat{\theta}$ un
estimador de dicho parámetro. Supongamos además
que$$\sqrt{n}\left( \hat{\theta}-\theta \right) \rightarrow \mathcal{N}\left( 0,\sigma
_{\theta }^2 \right),$$
en distribución. Entonces, bajo ciertas condiciones de regularidad 
(pueden verse en Bhattacharya y Ghosh, 1978) se tiene:
$$\begin{aligned}
P\left( \sqrt{n}\frac{\hat{\theta}-\theta }{\sigma _{\theta }}\leq u \right)
= &\ \Phi \left( u \right) +n^{-\frac{1}{2}}p_1\left( u \right) \phi \left(
u \right) +\cdots \\
& +n^{-\frac{j-1}{2}}p_{j-1}\left( u \right) \phi \left( u \right) +O\left(
n^{-\frac{j}{2}} \right), \\
P\left( \sqrt{n}\frac{\hat{\theta}-\theta }{\hat{\sigma}_{\theta }}\leq
u \right) = &\ \Phi \left( u \right) +n^{-\frac{1}{2}}q_1\left( u \right) \phi
\left( u \right) +\cdots \\
& +n^{-\frac{j-1}{2}}q_{j-1}\left( u \right) \phi \left( u \right) +O\left(
n^{-\frac{j}{2}} \right),\end{aligned}$$

siendo los $p_i\left( u \right)$ y $q_i\left( u \right)$ polinomios
de grado $3i-1$ con paridad alternada, es decir, $p_1$ y $q_1$ son
simétricos, $p_2$ y $q_2$ son antisimétricos, $p_3$ y $q_3$ son
simétricos y así sucesivamente.
\EndKnitrBlock{theorem}


## Bootstrap semiparamétrico y bootstrap residual {#boot-reg}

En ocasiones nos pueden interesar modelos semiparamétricos, en los que se asume una componente paramétrica pero no se especifica por completo la distribución de los datos. 
Una de las situaciones más habituales es en regresión, donde se puede considerar un modelo para la tendencia pero sin asumir una forma concreta para la distribución del error.

Nos centraremos en el caso de regresión y consideraremos como base el siguiente modelo general: 
\begin{equation} 
  Y = m(\mathbf{X}) + \varepsilon,
  (\#eq:modelogeneral)
\end{equation}
donde $Y$ es la respuesta, $\mathbf{X}=(X_1, X_2, \ldots, X_p)$ es el vector de variables explicativas, $m(\mathbf{x}) = E\left( \left. Y\right\vert_{\mathbf{X}=\mathbf{x}} \right)$ es la media condicional, denominada función de regresión (o tendencia), y $\varepsilon$ es un error aleatorio de media cero y varianza $\sigma^2$, independiente de $\mathbf{X}$ (errores homocedásticos independientes).

Supondremos que el objetivo es, a partir de una muestra:
$$\left\{ \left( X_{1i}, \ldots, X_{pi}, Y_{i} \right)  : i = 1, \ldots, n \right\},$$
realizar inferencias sobre la distribución condicional 
$\left.Y \right\vert_{\mathbf{X}=\mathbf{x}}$.

El modelo \@ref(eq:modelogeneral) se corresponde con el denominado *diseño aleatorio*, mas general.
Alternativamente se podría asumir que los valores de las variables explicativas no son aleatorios (por ejemplo han sido fijados por el experimentador), hablaríamos entonces de *diseño fijo*.
Para realizar inferencias sobre modelos de regresión con errores homocedásticos se podrían emplear dos algoritmos bootstrap (e.g. [Canty, 2002](http://cran.fhcrc.org/doc/Rnews/Rnews_2002-3.pdf), y subsecciones siguientes).
El primero consistiría en utilizar directamente bootstrap uniforme, remuestreando las observaciones, y sería adecuado para el caso de diseño aleatorio.
La otra alternativa, que podría ser más adecuada para el caso de diseño fijo, sería lo que se conoce como *remuestreo residual*, *remuestreo basado en modelos* o *bootstrap semiparamétrico*.
En esta aproximación se mantienen fijos los valores de las variables explicativas y se remuestrean los residuos.
Una de las aplicaciones del bootstrap semiparamétrico es el contraste de hipótesis en regresión, que se tratará en la Sección \@ref(contrastes-semiparametricos). 

Se puede generalizar el modelo \@ref(eq:modelogeneral) de diversas formas, por ejemplo asumiendo que la distribución del error depende de $X$ únicamente a través de la varianza (error heterocedástico independiente).
En este caso se suele reescribir como:
$$Y = m(\mathbf{X}) + \sigma(\mathbf{X}) \varepsilon,$$
siendo $\sigma^2(\mathbf{x}) = Var\left( \left. Y\right\vert_{\mathbf{X}=\mathbf{x}} \right)$ la varianza condicional y suponiendo adicionalmente que $\varepsilon$ tiene varianza uno.
Se podría modificar el bootstrap residual para este caso pero habría que modelizar y estimar la varianza condicional.
Alternativamente se podría emplear el denominado  *Wild Bootstrap* que se describirá en la Sección \@ref(wild-bootstrap) para el caso de modelos de regresión no paramétricos.

En esta sección nos centraremos en el caso de regresión lineal:
$$m_{\boldsymbol{\beta}}(\mathbf{x}) =  \beta_{0} + \beta_{1}X_{1} + \beta_{2}X_{2} + \cdots + \beta_{p}X_{p},$$ 
siendo $\boldsymbol{\beta} = \left(  \beta_{0}, \beta_{1}, \ldots, \beta_{p} \right)^{T}$ el vector de parámetros (desconocidos).
Su estimador mínimo cuadrático es:
$$\boldsymbol{\hat{\beta}} = \left( X^{T}X\right)^{-1}X^{T}\mathbf{Y},$$
siendo $\mathbf{Y} = \left( Y_{1}, \ldots, Y_{n} \right)^{T}$ el vector de observaciones de la variable $Y$ y $X$ la denominada *matriz del diseño* de las variables regresoras, cuyas filas son los valores observados de las variables explicativas.


En regresión lineal múltiple, bajo las hipótesis estructurales del modelo de normalidad y homocedásticidad, se dispone de resultados teóricos que permiten realizar inferencias sobre características de la distribución condicional. Si alguna de estas hipótesis no es cierta se podrían emplear aproximaciones basadas en resultados asintóticos, pero podrían ser poco adecuadas para tamaños muestrales no muy grandes. Alternativamente se podría emplear bootstrap.
Con otros métodos de regresión, como los modelos no paramétricos descritos en el Capítulo \@ref(npreg), es habitual emplear bootstrap para realizar inferencias sobre la distribución condicional.

En esta sección se empleará el conjunto de datos `Prestige` del paquete `carData`, considerando como variable respuesta `prestige` (puntuación de ocupaciones obtenidas a partir de una encuesta) y como variables explicativas: `income` (media de ingresos en la ocupación) y `education` (media de los años de educación).
Para ajustar el correspondiente modelo de regresión lineal podemos emplear el siguiente código:


```r
data(Prestige, package = "carData")
# ?Prestige
modelo <- lm(prestige ~ income + education, data = Prestige)
summary(modelo)
```

```
## 
## Call:
## lm(formula = prestige ~ income + education, data = Prestige)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -19.4040  -5.3308   0.0154   4.9803  17.6889 
## 
## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -6.8477787  3.2189771  -2.127   0.0359 *  
## income       0.0013612  0.0002242   6.071 2.36e-08 ***
## education    4.1374444  0.3489120  11.858  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 7.81 on 99 degrees of freedom
## Multiple R-squared:  0.798,	Adjusted R-squared:  0.7939 
## F-statistic: 195.6 on 2 and 99 DF,  p-value: < 2.2e-16
```

Como ejemplo, consideraremos que el objetivo es realizar inferencias sobre el coeficiente de determinación ajustado:


```r
res <- summary(modelo)
names(res)
```

```
##  [1] "call"          "terms"         "residuals"     "coefficients" 
##  [5] "aliased"       "sigma"         "df"            "r.squared"    
##  [9] "adj.r.squared" "fstatistic"    "cov.unscaled"
```

```r
res$adj.r.squared
```

```
## [1] 0.7939201
```

### Remuestreo de las observaciones {#boot-unif-reg}

Como ya se comentó, en regresión podríamos emplear bootstrap uniforme multidimensional para el caso de diseño aleatorio, aunque hay que tener en cuenta que con este método la distribución en el remuestreo de $\left. Y^{\ast}\right\vert _{X^{\ast}=X_i}$ es degenerada.

En este caso, podríamos realizar inferencias sobre el coeficiente de determinación ajustado empleando el siguiente código:


```r
library(boot)

case.stat <- function(data, i) {
  fit <- lm(prestige ~ income + education, data = data[i, ])
  summary(fit)$adj.r.squared
}

set.seed(1)
boot.case <- boot(Prestige, case.stat, R = 1000)
boot.case
```

```
## 
## ORDINARY NONPARAMETRIC BOOTSTRAP
## 
## 
## Call:
## boot(data = Prestige, statistic = case.stat, R = 1000)
## 
## 
## Bootstrap Statistics :
##      original      bias    std. error
## t1* 0.7939201 0.002495631   0.0315275
```

```r
# plot(boot.case)
boot.ci(boot.case, type = c("basic", "perc", "bca"))
```

```
## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
## Based on 1000 bootstrap replicates
## 
## CALL : 
## boot.ci(boot.out = boot.case, type = c("basic", "perc", "bca"))
## 
## Intervals : 
## Level      Basic              Percentile            BCa          
## 95%   ( 0.7331,  0.8570 )   ( 0.7308,  0.8547 )   ( 0.7203,  0.8497 )  
## Calculations and Intervals on Original Scale
```


### Bootstrap residual {#boot-residual}

Como ya se comentó, en el caso de diseño fijo podemos realizar un remuestreo de los residuos:
$$\mathbf{r} = \mathbf{Y} - X\hat{\mathbf{\beta}} = \mathbf{Y} - \hat{\mathbf{Y}}$$
obteniéndose las réplicas bootstrap:
$$\mathbf{Y}^{\ast} = \hat{\mathbf{Y}} + \mathbf{r}^{\ast}.$$
Por ejemplo, adaptando el código en Canty (2002) para este conjunto de datos, podríamos emplear:


```r
pres.dat <- Prestige
pres.dat$fit <- fitted(modelo)
pres.dat$res <- residuals(modelo)

mod.stat <- function(data, i) {
    data$prestige <- data$fit + data$res[i]
    fit <- lm(prestige ~ income + education, data = data)
    summary(fit)$adj.r.squared
}

set.seed(1)
boot.mod <- boot(pres.dat, mod.stat, R = 1000)
boot.mod
```

```
## 
## ORDINARY NONPARAMETRIC BOOTSTRAP
## 
## 
## Call:
## boot(data = pres.dat, statistic = mod.stat, R = 1000)
## 
## 
## Bootstrap Statistics :
##      original      bias    std. error
## t1* 0.7939201 0.004401997  0.02671996
```

```r
# plot(boot.mod)
boot.ci(boot.mod, type = c("basic", "perc", "bca"))
```

```
## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
## Based on 1000 bootstrap replicates
## 
## CALL : 
## boot.ci(boot.out = boot.mod, type = c("basic", "perc", "bca"))
## 
## Intervals : 
## Level      Basic              Percentile            BCa          
## 95%   ( 0.7407,  0.8464 )   ( 0.7415,  0.8471 )   ( 0.7244,  0.8331 )  
## Calculations and Intervals on Original Scale
## Some BCa intervals may be unstable
```

Sin embargo, la variabilidad de los residuos no reproduce la de los verdaderos errores, por lo que podría ser preferible (especialmente si el tamaño muestral es pequeño) emplear la modificación descrita en Davison y Hinkley (1997, Alg. 6.3, p. 271).
Teniendo en cuenta que:
$$\mathbf{r} = \left( I - H \right)\mathbf{Y},$$
siendo $H = X\left( X^{T}X\right)^{-1}X^{T}$ la matriz de proyección.
La idea es remuestrear los residuos reescalados (de forma que su varianza sea constante) y centrados $e_i - \bar{e}$, siendo:
$$e_i = \frac{r_i}{\sqrt{1 - h_{ii}}},$$
donde $h_{ii}$ es el valor de influencia o leverage, el elemento $i$-ésimo de la diagonal de $H$.

En `R` podríamos obtener estos residuos mediante los comandos:


```r
pres.dat$sres <- residuals(modelo)/sqrt(1 - hatvalues(modelo))
pres.dat$sres <- pres.dat$sres - mean(pres.dat$sres)
```

Sin embargo puede ser más cómodo emplear la función `Boot()` del paquete `car` (que internamente llama a la función `boot()`), 
como se describe en el apéndice "Bootstrapping Regression Models in R" del libro "An R Companion to Applied Regression" de Fox y Weisberg (2018), disponible [aquí](https://socialsciences.mcmaster.ca/jfox/Books/Companion/appendices/Appendix-Bootstrapping.pdf).

Esta función es de la forma:

```r
Boot(object, f = coef, labels = names(f(object)), R = 999, 
     method = c("case", "residual"))
```
donde:

- `object`: es un objeto que contiene el ajuste de un modelo de regresión.

- `f`: es la función de estadísticos (utilizando el ajuste como argumento).

- `method`: especifíca el tipo de remuestreo: remuestreo de observaciones (`"case"`)
  o de residuos (`"residual"`), empleando la modificación descrita anteriormente.


\BeginKnitrBlock{exercise}
<span class="exercise" id="exr:boot-car"><strong>(\#exr:boot-car) </strong></span>
\EndKnitrBlock{exercise}
Emplear la función `Boot()` del paquete `car` para hacer inferencia sobre 
el coeficiente de determinación ajustado del modelo de regresión lineal 
que explica `prestige` a partir de `income` y `education` 
(obtener una estimación del sesgo y de la predicción,
y una estimación por intervalo de confianza de este estadístico).


```r
library(car)

# set.seed(DNI)
# ...
```

