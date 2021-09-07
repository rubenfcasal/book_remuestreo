# Estimación de la precisión y el sesgo de un estimador {#prec-sesgo}




Uno de los problemas más interesantes que pueden ser abordados desde la
perspectiva de los métodos de remuestreo es el de la estimación del
sesgo y la precisión de un estimador. En dicho contexto surgió el método
Jackknife (bastante antes que el bootstrap), que, en ese sentido, puede
considerarse el método de remuestreo más antiguo como tal.

## Estimación bootstrap de la precisión y el sesgo de un estimador {#prec-sesgo-boot}

Consideremos $\mathbf{X}=\left( X_1,\ldots ,X_n \right)$ una
m.a.s. de una población con distribución $F$ y supongamos que tenemos
interés en realizar inferencia sobre un parámetro de la población
$\theta =\theta \left( F \right)$. Consideremos un estimador, 
$\hat{\theta}=T\left( \mathbf{X} \right)$, de dicho parámetro y
definamos el estadístico 
$$R=R\left( \mathbf{X}, F \right) = T\left( \mathbf{X} \right) 
- \theta \left( F \right) = \hat{\theta} - \theta.$$ 
El sesgo del estimador no es más que la esperanza del
estadístico $R$ y la varianza de $\hat{\theta}$ es también la varianza
de $R$ (pues $\theta$ no es aleatorio). Además, el error cuadrático
medio del estimador también se puede escribir como el momento de orden 2
de $R$:
$$\begin{aligned}
Sesgo\left( \hat{\theta} \right) &= E\left( \hat{\theta}-\theta \right)
=E\left( R \right), \\
Var\left( \hat{\theta} \right) &= Var\left( \hat{\theta}-\theta \right)
=Var\left( R \right), \\
MSE\left( \hat{\theta} \right) &= E\left[ \left( \hat{\theta}-\theta \right)
^2\right] =E\left( R^2 \right).
\end{aligned}$$

Dado que el principio bootstrap es útil para aproximar la distribución
en muestreo del estadístico $R$, entonces también permitirá aproximar
sus momentos (su esperanza, su varianza, la esperanza de su cuadrado) y
así proceder según sigue:

1. Estimar la función de distribución de probabilidad mediante $\hat{F}$

2. Para cada $i=1,\ldots ,n$ arrojar $X_i^{\ast}$ a partir de
$\hat{F}$ y obtener
$\mathbf{X}^{\ast}=\left( X_1^{\ast}, \ldots ,X_n^{\ast} \right)$

3. Calcular $R^{\ast}=R\left( \mathbf{X}^{\ast},\hat{F} \right)
=T\left( \mathbf{X}^{\ast} \right) -\theta \left( \hat{F} \right) =
\hat{\theta}^{\ast}- \hat \theta$

4. Repetir $B$ veces los pasos 2-3 para obtener las réplicas bootstrap
$R^{\ast (1)}, \ldots, R^{\ast (B)}$

5. Calcular el estimador bootstrap del sesgo:
$$Sesgo^{\ast}\left( \hat{\theta}^{\ast} \right) =\bar{R}^{\ast}=\frac{1
}{B}\sum_{b=1}^{B}R^{\ast (b)}.$$

El algoritmo anterior es útil para aproximar por bootstrap el sesgo. Si
se desea aproximar la varianza puede sustituirse al paso 5 por:

5. Calcular el estimador bootstrap de la varianza:
$$Var^{\ast}\left( \hat{\theta}^{\ast} \right) =\frac{1}{B}
\sum_{b=1}^{B}\left( R^{\ast (b)}-\bar{R}^{\ast} \right)^2$$

Si se trata de aproximar por bootstrap el error cuadrático medio, el
paso 5 pasaría a ser:

5. Calcular el estimador bootstrap del error cuadrático medio: 
$MSE^{\ast}\left( \hat{\theta}^{\ast} \right) =\frac{1}{B}\sum_{b=1}^{B}R^{\ast (b) 2}$

En el caso de la varianza podría ahorrarse algunos cálculos definiendo
directamente $R=T\left( \mathbf{X} \right) =\hat{\theta}$ y,
consecuentemente,
$R^{\ast}=T\left( \mathbf{X}^{\ast} \right) =
\hat{\theta}^{\ast}$. Así otro algoritmo de cálculo algo menos
intensivo sería:

1. Estimar la función de distribución de probabilidad mediante $\hat{
F}$

2. Para cada $i=1,\ldots ,n$ arrojar $X_i^{\ast}$ a partir de
$\hat{F}$ y obtener
$\mathbf{X}^{\ast}=\left( X_1^{\ast}, \ldots, X_n^{\ast} \right)$

3. Calcular
$T\left( \mathbf{X}^{\ast} \right) = \hat{\theta}^{\ast}$

4. Repetir $B$ veces los pasos 2-3 para obtener las réplicas bootstrap
$\hat{\theta}^{\ast (1)}, \ldots, \hat{\theta}^{\ast(B)}$

5. Calcular
$\overline{\hat{\theta}^{\ast}}=\frac{1}{B}\sum_{b=1}^{B}\hat{
\theta}^{\ast (b)}$ y, con ello,
$Var^{\ast}\left( \hat{\theta}^{\ast} \right) =\frac{1}{B}\sum_{b=1}^{B}\left( \hat{\theta}^{\ast \left(b \right)}-\overline{\hat{\theta}^{\ast}} \right)^2$

Es interesante mencionar que esto permite aproximar por bootstrap
(mediante Monte Carlo) la varianza de un estimador sin conocer una
expresión explícita para dicha varianza teórica.
En general, el estimador, $\hat{F}$, de $F$ a utilizar en los pasos 1-2
de estos algoritmos se elije según proceda al caso. En el caso del
bootstrap uniforme sería $\hat{F}=F_n$ y se puede proceder
como se  mostró en la Sección \@ref(intro-implementacion).

### Ejemplo: la media muestral

Consideremos como parámetro de interés la media de la población,
$\theta =\theta \left( F \right) =\mu =\int xdF\left( x \right)$, y
tomemos como estimador la media muestral:
$\hat{\theta}=\hat{\mu}=T\left( 
\mathbf{X} \right) =\bar{X}=\frac{1}{n}\sum_{j=1}^{n}X_j$.
Supongamos que deseamos estudiar la varianza de este estimador:
$Var\left( 
\hat{\theta} \right) =Var\left( \bar{X} \right) =\frac{\sigma^2}{n}$.

A la hora de aproximar por bootstrap $Var\left( \hat{\theta} \right)$,
si no disponemos de ninguna otra información adicional (como que la
distribución sea de cierta familia paramétrica o que sea continua),
parece razonable elegir como método de remuestreo el bootstrap uniforme.
En tal caso el algoritmo bootstrap de Monte Carlo procedería de esta
forma:

1. Para cada $i=1,\ldots ,n$ arrojar $U_i\sim \mathcal{U}\left( 0,1 \right)$ y
hacer $X_i^{\ast}=X_{\left\lfloor nU_i\right\rfloor +1}$

2. Calcular
$T\left( \mathbf{X}^{\ast} \right) =\bar{X}^{\ast}=
\frac{1}{n}\sum_{i=1}^{n}X_i^{\ast}$

3. Repetir $B$ veces los pasos 1-2 para obtener las réplicas bootstrap
$\bar{X}^{\ast (1)}, \ldots, \bar{X}^{\ast
(B)}$

4. Calcular $\overline{\bar{X}^{\ast}}=\frac{1}{B}\sum_{b=1}^{B}$
$\bar{X}^{\ast (b)}$ y
$Var^{\ast}\left( \bar{X}^{\ast} \right) =\frac{1}{B} \sum_{b=1}^{B}\left( \bar{X}^{\ast (b)}-\overline{\bar{X}^{\ast}} \right)^2$

De todas formas, en este caso puede verse fácilmente que no es necesario
realizar Monte Carlo. En efecto,
$$\begin{aligned}
Var^{\ast}\left( \bar{X}^{\ast} \right) & = \frac{1}{n^2}
\sum_{i=1}^{n}Var^{\ast}\left( X_i^{\ast} \right) =\frac{1}{n}Var^{\ast}\left( X_1^{\ast} \right) \\
& = \frac{1}{n}\left\{ E^{\ast}\left( X_1^{\ast 2} \right) -
\left[ E^{\ast}\left( X_1^{\ast} \right) \right]^2\right\} =\frac{1}{n}\left[ 
\frac{1}{n}\sum_{j=1}^{n}X_j^2-\left( \frac{1}{n}\sum_{j=1}^{n}X_j \right)^2 \right] \\
& = \frac{1}{n^2}\sum_{j=1}^{n}\left( X_j-\bar{X} \right)^2=\frac{S_n^2}{n},
\end{aligned}$$
que es precisamente el estimador
plug-in de la varianza de la media muestral.


\BeginKnitrBlock{example}\iffalse{-91-65-112-114-111-120-105-109-97-99-105-243-110-32-98-111-111-116-115-116-114-97-112-32-100-101-32-108-97-32-112-114-101-99-105-115-105-243-110-32-100-101-32-101-115-116-105-109-97-99-105-111-110-101-115-32-100-101-108-32-116-105-101-109-112-111-32-100-101-32-118-105-100-97-32-109-101-100-105-111-32-100-101-32-109-105-99-114-111-111-114-103-97-110-105-115-109-111-115-93-}\fi{}
<span class="example" id="exm:estimacion-boot-precision"><strong>(\#exm:estimacion-boot-precision)  \iffalse (Aproximación bootstrap de la precisión de estimaciones del tiempo de vida medio de microorganismos) \fi{} </strong></span><br> \vspace{0.5cm}

Continuando con el ejemplo de los tiempos de vida de microorganismos,
supongamos que queremos estimar
la precisión de dos estimadores de su vida media: media muestral y
mediana muestral, a partir de los datos observados: 0.143, 0.182, 0.256, 0.260, 0.270,
0.437, 0.509, 0.611, 0.712, 1.04, 1.09, 1.15, 1.46, 1.88, 2.08.
\EndKnitrBlock{example}

La estimación media muestral resulta $\bar{X}=0.8053333$. Por su
parte la estimación mediana muestral es $x_{\left( 8 \right)}=0.611$.

La varianza del estimador media muestral, $\bar{X}$ es $Var\left( 
\bar{X} \right) =\frac{\sigma^2}{n}$, desconocida en este caso.
Su estimación bootstrap (idéntica a la plug-in) mediante un remuestreo
uniforme es calculable sin necesidad de realizar Monte Carlo y
resulta:
$$Var^{\ast}\left( \bar{X}^{\ast} \right) =\frac{1}{n^2}
\sum_{j=1}^{n}\left( x_j-\bar{X} \right)^2=0.024204877.$$
Con lo cual $\sqrt{Var^{\ast}\left( \bar{X}^{\ast} \right)}=\sqrt{
0.024204877}= 0.155\,58$

Si consideramos ahora la mediana muestral (como estimador de la media),
también sabemos que su distribución bootstrap puede calcularse de forma
explícita, sin necesidad de realizar Monte Carlo. Su masa de
probabilidad viene dada por:

$$\begin{aligned}
P^{\ast}\left( X_{\left( 8 \right)}^{\ast}=x_{(1)} \right)
&= 1-\sum_{k=0}^{m-1}\binom{n}{k}\left( \frac{1}{n} \right)^{k}\left( \frac{
n-1}{n} \right)^{n-k} \\
P^{\ast}\left( X_{\left( 8 \right)}^{\ast}=x_{(j)} \right)
&= \sum_{k=0}^{m-1}\binom{n}{k}\left[ \left( \frac{j-1}{n} \right)^{k}\left( 
\frac{n-j+1}{n} \right)^{n-k} - \left( \frac{j}{n} \right)^{k}\left( \frac{n-j}{n} \right)^{n-k}
\right] \\
\text{para }j &= 2,\ldots ,n
\end{aligned}$$

Con los datos concretos del ejemplo resulta: 
$$\begin{array}{ll}
P^{\ast}\left( X_{\left( 8 \right)}^{\ast} = 0.143 \right) = 1.639\times 10^{-6}\text{, } 
& P^{\ast}\left( X_{\left( 8 \right)}^{\ast} = 0.182 \right) = 2.655\times 10^{-4},\\
P^{\ast}\left( X_{\left( 8 \right)}^{\ast} = 0.256 \right) = 3.973\times 10^{-3}\text{, } 
& P^{\ast}\left( X_{\left( 8 \right)}^{\ast} = 0.260 \right) = 2.121\times 10^{-2}, \\
P^{\ast}\left( X_{\left( 8 \right)}^{\ast} = 0.270 \right) = 6.278\times 10^{-2}\text{, }
& P^{\ast}\left( X_{\left( 8 \right)}^{\ast} = 0.437 \right) = 0.1249, \\
P^{\ast}\left( X_{\left( 8 \right)}^{\ast} = 0.509 \right) = 0.1832\text{, }
& P^{\ast}\left( X_{\left( 8 \right)}^{\ast} = 0.611 \right) = 0.2073,\\
P^{\ast}\left( X_{\left( 8 \right)}^{\ast} = 0.712 \right) = 0.1832\text{, }
& P^{\ast}\left( X_{\left( 8 \right)}^{\ast} = 1.04 \right) = 0.1249,\\
P^{\ast}\left( X_{\left( 8 \right)}^{\ast} = 1.09 \right) = 6.278\times 10^{-2}\text{, }
& P^{\ast}\left( X_{\left( 8 \right)}^{\ast} = 1.15 \right) = 2.121\times 10^{-2}, \\
P^{\ast}\left( X_{\left( 8 \right)}^{\ast} = 1.46 \right) = 3.973\times 10^{-3}\text{, }
& P^{\ast}\left( X_{\left( 8 \right)}^{\ast} = 1.88 \right) = 2.655\times 10^{-4}, \\
P^{\ast}\left( X_{\left( 8 \right)}^{\ast} = 2.08 \right) = 1.639\times 10^{-6}. 
\end{array}$$

Como consecuencia, 
$$\begin{aligned}
E^{\ast}\left( X_{\left( 8 \right)}^{\ast} \right)
&= \sum_{j=1}^{15}x_{(j)}P^{\ast}\left( X_{\left( 8 \right)
}^{\ast}=x_{(j)} \right) =0.65749924 \\
E^{\ast}\left( X_{\left( 8 \right)}^{\ast 2} \right)
&= \sum_{j=1}^{15}x_{(j)}^2P^{\ast}\left( X_{\left( 8 \right)
}^{\ast}=x_{(j)} \right) =0.49500381 \\
Var^{\ast}\left( X_{\left( 8 \right)}^{\ast} \right)
&= 0.49500381-0.65749924^2=6.\, 2699\times 10^{-2} \\
\sqrt{Var^{\ast}\left( X_{\left( 8 \right)}^{\ast} \right)} &= \sqrt{
6.\, 2699\times 10^{-2}}=0.250\,40
\end{aligned}$$

Las estimaciones bootstrap de los errores cuadráticos medios de ambos
estimadores (como estimadores de la media poblacional) son:
$$\begin{aligned}
MSE^{\ast}\left( \bar{X}^{\ast} \right) &= \left( E^{\ast}\left( 
\bar{X}^{\ast} \right) -\bar{X} \right)^2+Var^{\ast}\left( 
\bar{X}^{\ast} \right) \\
&= Var^{\ast}\left( \bar{X}^{\ast} \right) =0.024204877, \\
MSE^{\ast}\left( X_{\left( 8 \right)}^{\ast} \right) &= \left( E^{\ast
}\left( X_{\left( 8 \right)}^{\ast} \right) -\bar{X} \right)
^2+Var^{\ast}\left( X_{\left( 8 \right)}^{\ast} \right) = \\
&= \left( 0.65749924-0.8053333 \right)^2+6.2699 \times 10^{-2}\\
&=  0.084554.
\end{aligned}$$

Una aproximación de Monte Carlo de estas varianzas bootstrap se puede
llevar a cabo mediante el siguiente código:


```r
# Para la muestra de TIEMPOS DE VIDA, estima (plug-in) la precisión 
# de la media muestral (desvmedia) y también aproxima por Monte Carlo
# la estimación bootstrap de dicha precisión (desvmediaboot) y
# también la precisión de la mediana muestral, de la cual no se
# conoce su expresión, (desvmedianaboot).
# También estima el sesgo bootstrap de esos dos estimadores 
# (sesgomediaboot y sesgomedianaboot, respectivamente).

muestra <- c(0.143, 0.182, 0.256, 0.26, 0.27, 0.437, 0.509,
    0.611, 0.712, 1.04, 1.09, 1.15, 1.46, 1.88, 2.08)
n <- length(muestra)
varmedia <- (1/(n^2)) * sum((muestra - mean(muestra))^2)
# Alternativamente: varmedia <- var(muestra)/n
desvmedia <- sqrt(varmedia)

# Remuestreo
B <- 1e+04
media <- numeric(B)
mediana <- numeric(B)
for (k in 1:B) {
    remuestra <- sample(muestra, n, replace = TRUE)
    media[k] <- mean(remuestra)
    # remordenada <- sort(remuestra)
    # mediana[k] <- remordenada[8]
    mediana[k] <- median(remuestra)
}

# Aproximaciones precisión
varmediaboot <- (1/B) * sum((media - mean(media))^2)
desvmediaboot <- sqrt(varmediaboot)
varmedianaboot <- (1/B) * sum((mediana - mean(mediana))^2)
desvmedianaboot <- sqrt(varmedianaboot)
desvmedia
```

```
## [1] 0.1555792
```

```r
desvmediaboot
```

```
## [1] 0.1544907
```

```r
desvmedianaboot
```

```
## [1] 0.2475729
```

```r
# Aproximaciones sesgo
sesgomediaboot <- mean(media) - mean(muestra)
sesgomedianaboot <- mean(mediana) - muestra[8]
sesgomediaboot
```

```
## [1] 0.0009302
```

```r
sesgomedianaboot
```

```
## [1] 0.045781
```

Empleando el paquete `boot` el código sería más simple:


```r
library(boot)
statistic <- function(data, i){
  remuestra <- data[i]
  c(mean(remuestra), median(remuestra))
}

set.seed(1)
res.boot <- boot(muestra, statistic, R = B)
res.boot
```

```
## 
## ORDINARY NONPARAMETRIC BOOTSTRAP
## 
## 
## Call:
## boot(data = muestra, statistic = statistic, R = B)
## 
## 
## Bootstrap Statistics :
##      original       bias    std. error
## t1* 0.8053333 7.115333e-05   0.1572396
## t2* 0.6110000 4.529410e-02   0.2511022
```

Lamentablemente la función `print.boot()` calcula las aproximaciones
bootstrap del sesgo y de la precisión pero no las almacena. 
En el caso más simple podríamos obtenerlas con el siguiente código:


```r
op <- with(res.boot, cbind(
  t0, apply(t, 2, mean, na.rm = TRUE) -  t0,
  apply(t, 2, sd, na.rm = TRUE)
  ))
rownames(op) <- paste0("t", 1:ncol(res.boot$t), "*")
colnames(op) <- c("original", "bias  ", " std. error")
op
```

```
##      original       bias    std. error
## t1* 0.8053333 7.115333e-05   0.1572396
## t2* 0.6110000 4.529410e-02   0.2511022
```

## Motivación del método Jackknife {#jackknife}

El jackknife es probablemente el método de remuestreo, propiamente
dicho, más antiguo. Fue propuesto por Quenouille (1949) para estimar el
sesgo de un estimador. Tukey (1958) bautiza el método y lo utiliza para
estimar la varianza de un estimador. En realidad el jackknife no suele
utilizarse para aproximar la distribución de $R\left( \mathbf{X},F \right)$, 
sino más bien para estimar características de dicha
variable aleatoria, como su esperanza o su varianza.

La diferencia entre el bootstrap y el jackknife es muy fácil de expresar
en términos de los vectores de remuestreo. Así, el bootstrap uniforme
utiliza vectores de remuestreo de la forma 
$\mathbf{P}^{\ast}=\left( \frac{m_1}{n},\ldots ,\frac{m_n}{n} \right)$, con
$m_i\in \mathbb{Z}^{+}$, $i=1,\ldots ,n$, mientras que el jackknife considera
vectores de remuestreo de la forma 
$$\mathbf{P}_{(i)}^{\ast}=\left( \frac{1}{n-1},\ldots ,\underset{(i)}{0}
,\ldots ,\frac{1}{n-1} \right).$$
En otras palabras todas las remuestras
jackknife posibles son tantas como el tamaño muestral y cada una
consiste en eliminar una observación de la muestra, quedándose con una
remuestra de tamaño $n-1$ en la que las demás observaciones aparecen
exactamente con frecuencia $1$.

Evidentemente, el número de posibles remuestras jackknife, $n$, es
muchísimo más pequeño que el número de remuestras bootstrap, 
$\binom{2n-1}{n}$, lo cual permite calcular con rapidez las realizaciones
del estadístico de interés en todas las posibles remuestras jackknife.

## Estimación Jackknife de la precisión y el sesgo de un estimador

Cuando estamos interesados en el sesgo o la varianza de un estimador
$\hat{\theta}=\theta \left( \mathbf{X} \right)$ de un parámetro 
$\theta =\theta \left( F \right)$, el estadístico de interés suele
definirse como
$R=R\left( \mathbf{X},F \right) =\hat{\theta}-\theta$. 
En este caso
$$\begin{aligned}
Sesgo\left( \hat{\theta} \right) &= E\left( \hat{\theta} \right) -\theta
=E\left( R \right), \\
Var\left( \hat{\theta} \right) &= Var\left( \hat{\theta}-\theta \right)
=Var\left( R \right).
\end{aligned}$$
Así pues trataremos de usar el
jackknife para aproximar la esperanza y varianza de $R$, o,
equivalentemente, el sesgo y la varianza de $\hat{\theta}$.

El conjunto de remuestras jackknife es
$$\mathcal{X}_{jackk}=\left\{ \mathbf{X}^{\ast}=
\mathbf{X}_{(i)}=\left( X_1,\ldots ,X_{i-1},X_{i+1},\ldots
,X_n \right) : i=1,\ldots ,n\right\}$$
y todas ellas se consideran con
equiprobabilidad en el universo jackknife. Como primera tentativa
estimaríamos el sesgo y la varianza jackknife mediante:
$$\begin{aligned}
E^{\ast}\left( R^{\ast} \right) &= E_{jackk}^{\ast}\left( \hat{\theta}
^{\ast} \right) -\theta \left( \mathbf{X} \right) =\frac{1}{n}
\sum_{i=1}^{n}\theta \left( \mathbf{X}_{(i)} \right) -
\hat{\theta}=\overline{\theta \left( \mathbf{X}_{(\cdot)} \right)}-\hat{\theta}, \\
Var^{\ast}\left( R^{\ast} \right) &= Var_{jackk}^{\ast}\left( \hat{\theta}
^{\ast} \right) =\frac{1}{n}\sum_{i=1}^{n}\left[ \theta \left( 
\mathbf{X}_{(i)} \right) -\overline{\theta \left( 
\mathbf{X}_{(\cdot)} \right)}\right]^2,
\end{aligned}$$
con $\overline{\theta \left( \mathbf{X}_{(\cdot)} \right)}
= \frac{1}{n}\sum_{j=1}^{n}\theta \left( \mathbf{X}_{(j)} \right)$.

Sin embargo es evidente que las réplicas jackknife son mucho más
parecidas a la muestra original de lo que lo son las remuestras
bootstrap, en general. De hecho se puede demostrar que el valor absoluto
de ese estimador jackknife del sesgo es siempre menor que el valor
absoluto del sesgo bootstrap y que la estimación jackknife de la
varianza que se acaba de proponer también es menor que la varianza
bootstrap. En resumen, el método jackknife necesita de un **factor de
elevación** para que las estimaciones que proporciona sean consistentes.
La idea es elegir dicho factor de elevación como aquel que provoca que,
al multiplicar los estadísticos anteriores por él, y considerando como
parámetro a estimar la media o la varianza poblacional, el estimador
jackknife finalmente resultante sea insesgado. Así, el factor de
elevación resulta ser $n-1$ y las estimaciones jackknife finales son

$$\begin{aligned}
Sesgo_{jackk}^{\ast}\left( \hat{\theta}^{\ast} \right) &= \left( n-1 \right)
\left( \overline{\theta \left( \mathbf{X}_{(\cdot)}
 \right)}-\hat{\theta} \right) =\frac{n-1}{n}\sum_{i=1}^{n}\left( \theta
\left( \mathbf{X}_{(i)} \right) -\hat{\theta} \right), \\
Var_{jackk}^{\ast}\left( \hat{\theta}^{\ast} \right) &= \frac{n-1}{n}
\sum_{i=1}^{n}\left[ \theta \left( \mathbf{X}_{(i)}
 \right) -\overline{\theta \left( \mathbf{X}_{(\cdot)}
 \right)}\right]^2.
\end{aligned}$$

Tomando como parámetro de interés la media, $\theta =\mu$, tenemos que
$$\begin{aligned}
\overline{\theta \left( \mathbf{X}_{(\cdot)} \right)}-
\hat{\theta} &= \frac{1}{n}\sum_{i=1}^{n}\theta \left( 
\mathbf{X}_{(i)} \right) -\hat{\theta}=\frac{1}{n}\sum_{i=1}^{n}
\overline{X_{(i)}}-\bar{X}= \\
&= \frac{1}{n}\sum_{i=1}^{n}\frac{1}{n-1}\sum_{j=1,j\neq i}^{n}X_j-
\bar{X}=\frac{1}{n\left( n-1 \right)}\sum_{i,j=1,i\neq j}^{n}X_j-
\bar{X} \\
&= \frac{1}{n\left( n-1 \right)}\sum_{j=1}^{n}\left( n-1 \right) X_j-
\bar{X}=\frac{1}{n}\sum_{j=1}^{n}X_j-\bar{X}=0,\end{aligned}$$

así que $\overline{\theta \left( \mathbf{X}_{(\cdot)}
\right)}-\hat{\theta}$ es un estimador insesgado del sesgo de
$\bar{X}$ (que es cero). De esta forma, utilizando un factor de
elevación arbitrario, $c$, se tiene igualmente que
$$c\left( \overline{\theta \left( \mathbf{X}_{(\cdot)} \right)}
- \hat{\theta} \right) = 0,$$
así que es también un estimador insesgado de 
$Sesgo\left( \bar{X} \right) =0.$$
Determinaremos el valor de $c$
imponiendo que el estimador jackknife de la varianza de dicho estimador
($\hat{\theta}=\bar{X}$) es una estimador insesgado de la varianza de
dicho estimador. Por una parte, es bien conocido que la varianza de
$\hat{\theta}$ es $Var\left( \bar{X} \right) = \sigma^2 /n$. 
Por otra parte, la estimación jackknife de la varianza de $\hat{\theta}$,
con factor de elevación $c$ es

$$\begin{aligned}
Var_{jackk}^{\ast}\left( \bar{X} \right) &= \frac{c}{n}\sum_{i=1}^{n}
\left[ \overline{X_{(i)}}-\overline{\overline{X_{\left( \cdot
 \right)}}}\right]^2 \\
&= \frac{c}{n}\sum_{i=1}^{n}\left[ \frac{1}{n-1}\sum_{j=1,j\neq i}^{n}X_j-
\frac{1}{n}\sum_{k=1}^{n}\frac{1}{n-1}\sum_{j=1,j\neq k}^{n}X_j\right]^2
\\
&= \frac{c}{n}\sum_{i=1}^{n}\left[ \frac{1}{n-1}\sum_{j=1,j\neq i}^{n}X_j-
\frac{1}{n\left( n-1 \right)}\sum_{k,j=1,j\neq k}^{n}X_j\right]^2 \\
&= \frac{c}{n}\sum_{i=1}^{n}\left[ \frac{1}{n\left( n-1 \right)}
\sum_{j=1,j\neq i}^{n}X_j-\frac{1}{n}X_i\right]^2.\end{aligned}$$

La esperanza de esta cantidad resulta:
$$\begin{aligned}
&E\left[ \frac{c}{n}\sum_{i=1}^{n}\left( \frac{1}{n\left( n-1 \right)}
\sum_{j=1,j\neq i}^{n}X_j-\frac{1}{n}X_i \right)^2\right] \\
&= \frac{c}{n}\sum_{i=1}^{n}E\left[ \left( \frac{1}{n\left( n-1 \right)}
\sum_{j=1,j\neq i}^{n}X_j-\frac{1}{n}X_i \right)^2\right] \\
&= \frac{c}{n}\sum_{i=1}^{n}E\left[ \left( \frac{1}{n\left( n-1 \right)}
\sum_{j=1,j\neq i}^{n}\left( X_j-\mu \right) 
-\frac{1}{n}\left( X_i-\mu \right) \right)^2\right] \\
&= \frac{c}{n}\sum_{i=1}^{n}Var\left[ \frac{1}{n\left( n-1 \right)}
\sum_{j=1,j\neq i}^{n}\left( X_j-\mu \right) 
-\frac{1}{n}\left( X_i-\mu \right) \right]\\
&= \frac{c}{n}\sum_{i=1}^{n}\left( \frac{1}{n^2\left( n-1 \right)^2}
\sum_{j=1,j\neq i}^{n}\sigma^2+\frac{1}{n^2}\sigma^2 \right) \\
&= \frac{c}{n}\sum_{i=1}^{n}\left( \frac{1}{n^2\left( n-1 \right)}\sigma^2+\frac{1}{n^2}\sigma^2 \right) 
=\frac{c\sigma^2}{n\left(n-1 \right)}.
\end{aligned}$$

Así pues, el sesgo del estimador jackknife de la varianza de la media
muestral es
$$E\left[ Var_{jackk}^{\ast}\left( \bar{X} \right) \right] -\frac{\sigma
^2}{n}=\frac{c\sigma^2}{n\left( n-1 \right)}-\frac{\sigma^2}{n}=
\frac{\sigma^2}{n}\left( \frac{c}{n-1}-1 \right),$$
que vale cero si y solamente si $c=n-1$. 
Dicho en otras palabras, tomando como factor de
elevación $c=n-1$, entonces, tanto el estimador jackknife del sesgo de
$\bar{X}$ como el estimador jackknife de la varianza de 
$\bar{X}$ son estimadores insesgados, respectivamente, del sesgo y
la varianza de $\bar{X}$.

Dichos estimadores resultan$$\begin{aligned}
Sesgo_{jackk}^{\ast}\left( \bar{X} \right) &= \left( n-1 \right) \left( 
\overline{\overline{X_{(\cdot)}}}-\bar{X} \right), \\
Var_{jackk}^{\ast}\left( \bar{X} \right) &= \frac{n-1}{n}\sum_{i=1}^{n}
\left[ \overline{X_{(i)}}-\overline{\overline{X_{\left( \cdot
 \right)}}}\right]^2.\end{aligned}$$

Podemos realizar un razonamiento análogo cuando el parámetro de interés
es la varianza poblacional, $\theta =\sigma^2$. En ese caso,
considerando el estimador varianza muestral: $\hat{\theta}=S_n^2$,
se tiene que su esperanza viene dada por

$$\begin{aligned}
E\left( S_n^2 \right) &= E\left[ \frac{1}{n}\sum_{i=1}^{n}\left( X_i-
\bar{X} \right)^2\right] =\frac{1}{n}\sum_{i=1}^{n}E\left[ \left(
X_i-\bar{X} \right)^2\right] \\
&= E\left[ \left( X_1-\frac{1}{n}\sum_{j=1}^{n}X_j \right)^2\right] =E
\left[ \left( \left( X_1-\mu \right) -\frac{1}{n}\sum_{j=1}^{n}\left(
X_j-\mu \right) \right)^2\right] \\
&= Var\left[ \left( X_1-\mu \right) -\frac{1}{n}\sum_{j=1}^{n}\left(
X_j-\mu \right) \right] \\
&= Var\left[ \frac{n-1}{n}\left( X_1-\mu \right) -\frac{1}{n}
\sum_{j=1,j\neq 1}^{n}\left( X_j-\mu \right) \right] \\
&= \left( \frac{n-1}{n} \right)^2\sigma^2+\frac{1}{n^2}
\sum_{j=2}^{n}\sigma^2 \\
&= \frac{\left( n-1 \right)}{n^2}^2\sigma^2+\frac{n-1}{n^2}\sigma
^2=\frac{n\left( n-1 \right)}{n^2}\sigma^2=\frac{n-1}{n}\sigma^2,\end{aligned}$$

así que su sesgo es
$$Sesgo\left( S_n^2 \right) =E\left( S_n^2 \right) -\sigma^2=-\frac{1
}{n}\sigma^2.$$Para un factor de elevación, $c$, el estimador
jackknife del sesgo de este estimador
es$$Sesgo_{jackk}^{\ast}\left( S_n^2 \right) =c\left( \overline{\theta
\left( \mathbf{X}_{(\cdot)} \right)}-\hat{\theta}
 \right) =c\left( \overline{S_{n,(\cdot)}^2}-S_n^2 \right)
.$$Con lo cual la esperanza de este estimador resulta
$$E\left( c\left( \overline{S_{n,(\cdot)}^2}-S_n^2 \right)
 \right) =c\left[ E\left( \overline{S_{n,(\cdot)}^2} \right)
-E\left( S_n^2 \right) \right]$$

Estudiemos por separado cada término:
$$\begin{aligned}
\overline{S_{n,(\cdot)}^2} &= \frac{1}{n}\sum_{i=1}^{n}S_{n,
(i)}^2=\frac{1}{n}\sum_{i=1}^{n}\frac{1}{n-1}\sum_{j=1,j\neq
i}^{n}\left( X_j-\overline{X_{(i)}} \right)^2 \\
&= \frac{1}{n}\sum_{i=1}^{n}\frac{1}{n-1}\sum_{j=1,j\neq i}^{n}\left( X_j-
\frac{1}{n-1}\sum_{k=1,k\neq i}^{n}X_{k} \right)^2 \\
&= \frac{1}{n\left( n-1 \right)}\sum_{i=1}^{n}\sum_{j=1,j\neq i}^{n}\left( 
\frac{n-2}{n-1}X_j-\frac{1}{n-1}\sum_{k=1,k\neq i,k\neq j}^{n}X_{k} \right)^2,
\end{aligned}$$
con lo cual

$$\begin{aligned}
E\left( \overline{S_{n,(\cdot)}^2} \right) 
&=\frac{1}{n\left(
n-1 \right)}\sum_{i=1}^{n}\sum_{j=1,j\neq i}^{n}E\left[ \left( \frac{n-2}{n-1
}X_j-\frac{1}{n-1}\sum_{k=1,k\neq i,k\neq j}^{n}X_{k} \right)^2\right] \\
&=\frac{1}{n\left( n-1 \right)}\sum_{i=1}^{n}\sum_{j=1,j\neq i}^{n}E\left[
\left( \frac{n-2}{n-1}\left( X_j-\mu \right) -\frac{1}{n-1}\sum_{k=1,k\neq
i,k\neq j}^{n}\left( X_{k}-\mu \right) \right)^2\right] \\
&=\frac{1}{n\left( n-1 \right)}\sum_{i=1}^{n}\sum_{j=1,j\neq i}^{n}Var\left[ 
\frac{n-2}{n-1}\left( X_j-\mu \right) -\frac{1}{n-1}\sum_{k=1,k\neq
i,k\neq j}^{n}\left( X_{k}-\mu \right) \right] \\
&=\frac{1}{n\left( n-1 \right)}\sum_{i=1}^{n}\sum_{j=1,j\neq i}^{n}\left[ 
\frac{\left( n-2 \right)^2}{\left( n-1 \right)^2}\sigma^2+\frac{1}{
\left( n-1 \right)^2}\left( n-2 \right) \sigma^2\right] \\
&=\frac{1}{n\left( n-1 \right)}\sum_{i=1}^{n}\sum_{j=1,j\neq i}^{n}\left[ 
\frac{\left( n-2 \right) \left( n-1 \right)}{\left( n-1 \right)^2}\sigma
^2\right] =\frac{n-2}{n-1}\sigma^2
\end{aligned}$$

Además, ya hemos visto anteriormente que $E\left( S_n^2 \right) =
(n-1)\sigma^2/n$, con o cual la esperanza del estimador
jackknife del sesgo es
$$\begin{aligned}
c\left( \frac{n-2}{n-1}\sigma^2-\frac{n-1}{n}\sigma^2 \right)
&= c\sigma^2\left( \frac{n-2}{n-1}-\frac{n-1}{n} \right) \\
&= c\sigma^2\frac{n\left( n-2 \right) -\left( n-1 \right)^2}{n\left( n-1 \right)} \\
&= c\sigma^2\frac{-1}{n\left( n-1 \right)}
=-\frac{c\sigma^2}{n\left(n-1 \right)}.
\end{aligned}$$

De esta forma, el sesgo del estimador jackknife del sesgo de $S_n^2$
resulta ser
$$\begin{aligned}
E\left[ Sesgo_{jackk}^{\ast}\left( S_n^2 \right) \right] 
-Sesgo\left(S_n^2 \right) &= -\frac{c\sigma^2}{n\left( n-1 \right)}
-\left( -\frac{1}{n}\sigma^2 \right) \\
&= -\frac{\sigma^2}{n}\left( \frac{c}{n-1}-1 \right),
\end{aligned}$$

con lo cual este sesgo será cero si y sólo si $c=n-1$.

Esto da pie al estimador jackknife del sesgo de la varianza muestral:
$$Sesgo_{jackk}^{\ast}\left( S_n^2 \right) =\left( n-1 \right) \left( 
\overline{S_{n,(\cdot)}^2}-S_n^2 \right).$$

Así pues, queda justificado, en el caso de estimación de los parámetros
media y varianza, la razón de la elección del factor de elevación
$c=n-1$.

## Relación Bootstrap/Jackknife en dicha estimación

Consideremos un parámetro de interés $\theta \left( F \right)$ y su
correspondiente estimador que supondremos funcional, $\theta \left(
F_n \right)$. En realidad, cuando calculamos cantidades como $\theta
\left( \mathbf{X}_{(i)} \right)$, lo que estamos
haciendo es evaluar el funcional $\theta$ en otra función de
distribución$$F_{n,(i)}\left( x \right) =\frac{1}{n-1}\sum_{j=1,j\neq i}^{n}
\mathbf{1}\left( X_j\leq x \right).$$Dicho en terminología de vectores
de remuestreo, el estimador habitual consiste en evaluar $\theta$ en el
vector $\mathbf{p}=\left( \frac{1}{n},\ldots ,\frac{1}{n} \right)$,
$\theta \left( \mathbf{p} \right)$, mientras que el estimador
construido con toda la muestra excepto el dato $i$-ésimo es la evaluación 
$\theta \left( \mathbf{p}_{(i)} \right)$, siendo
$$\mathbf{p}_{(i)}=\left( \frac{1}{n-1},\ldots ,\frac{1}{n-1},
\underset{(i)}{0},\frac{1}{n-1},\ldots ,\frac{1}{n-1} \right).$$

En lo que sigue, consideraremos funcionales $\theta$ que definen
estimadores lineales o cuadráticos en los vectores de remuestreo:

-   Estimadores lineales: 
    $$\theta \left( \mathbf{p} \right) = a+\mathbf{b}^{T} \mathbf{p}$$

-   Estimadores cuadráticos:
    $$\theta \left( \mathbf{p} \right) = a+\mathbf{b}^{T}
    \mathbf{p}+\mathbf{p}^{T}C\mathbf{p}$$

Una forma alternativa de definir estos estimadores es

-   Estimadores lineales: 
    $$\theta \left( \mathbf{p} \right) = \mathbf{b}^{T}
    \left( \mathbf{p}-\mathbf{p}_{0} \right) $$

-   Estimadores cuadráticos:
    $$\theta \left( \mathbf{p} \right) = \left( \mathbf{p}-\mathbf{p}_{0} \right)^{T}
    C\left( \mathbf{p}-\mathbf{p}_{0} \right)$$

Por ejemplo, puede demostrarse fácilmente que la media, $\bar{X}$,
es un estimador lineal en el vector de remuestreo y que la varianza
muestral, $S_n^2$, es un estimador cuadrático en el vector de
remuestreo.

Existen dos resultados que relacionan el sesgo bootstrap y el sesgo
jackknife de cualquier estimador cuadrático y la varianza bootstrap y la
varianza jackknife de cualquier estimador lineal.


\BeginKnitrBlock{theorem}
<span class="theorem" id="thm:jack-boot-sesgo"><strong>(\#thm:jack-boot-sesgo) </strong></span><br> \vspace{0.5cm}

Si $\hat{\theta}$ es un estimador cuadrático, entonces
$$Sesgo_{jackk}\left( \hat{\theta} \right) =\frac{n}{n-1}Sesgo_{boot}\left( 
\hat{\theta} \right)$$

\EndKnitrBlock{theorem}

\BeginKnitrBlock{theorem}
<span class="theorem" id="thm:jack-boot-precision"><strong>(\#thm:jack-boot-precision) </strong></span><br> \vspace{0.5cm}

Si $\hat{\theta}$ es un estimador lineal, entonces
$$Var_{jackk}\left( \hat{\theta} \right) =\frac{n}{n-1}Var_{boot}\left( \hat{
\theta} \right)$$
\EndKnitrBlock{theorem}


Dicho en otras palabras, para cualquier estimador cuadrático, el sesgo
jackknife es mayor que el sesgo bootstrap. Si el estimador es lineal, la
varianza jackknife es mayor que la varianza bootstrap. En ambos casos,
el factor multiplicador es $n/(n-1)$.


\BeginKnitrBlock{example}\iffalse{-91-65-112-114-111-120-105-109-97-99-105-243-110-32-106-97-99-107-107-110-105-102-101-32-100-101-32-108-97-32-112-114-101-99-105-115-105-243-110-32-100-101-32-101-115-116-105-109-97-99-105-111-110-101-115-32-100-101-108-32-116-105-101-109-112-111-32-100-101-32-118-105-100-97-32-109-101-100-105-111-32-100-101-32-109-105-99-114-111-111-114-103-97-110-105-115-109-111-115-93-}\fi{}
<span class="example" id="exm:estimacion-jack-precision"><strong>(\#exm:estimacion-jack-precision)  \iffalse (Aproximación jackknife de la precisión de estimaciones del tiempo de vida medio de microorganismos) \fi{} </strong></span><br> \vspace{0.5cm}

Consideremos la muestra de tiempos de vida de microorganismos ya
tratada. El siguiente código permite calcular los
estimadores jackknife del sesgo y de la precisión tanto de la media como
de la mediana muestral.
\EndKnitrBlock{example}


```r
# Para la muestra de TIEMPOS DE VIDA, estima (plug-in) la precisión 
# de la media muestral (desvmedia) y también estima mediante el jackknife 
# dicha precisión (desvmediajackk) y también la precisión de la mediana 
# muestral, de la cual no se conoce su expresión, (desvmedianajackk). 
# También estima el sesgo jackknife de esos dos estimadores
# (sesgomediajackk y sesgomedianajackk, respectivamente).

muestra <- c(0.143, 0.182, 0.256, 0.26, 0.27, 0.437, 0.509, 
    0.611, 0.712, 1.04, 1.09, 1.15, 1.46, 1.88, 2.08)
n <- length(muestra)
varmedia <- (1/(n^2)) * sum((muestra - mean(muestra))^2)
desvmedia <- sqrt(varmedia)

# Jackknife
media <- numeric(n)
mediana <- numeric(n)
for (i in 1:n) {
    imuestra <- muestra[-i]
    media[i] <- mean(imuestra)
    # remordenada <- sort(imuestra)
    # mediana[i] <- (remordenada[7] + remordenada[8])/2
    mediana[i] <- median(imuestra)
}

# Aproximaciones precisión
varmediajackk <- ((n - 1)/n) * sum((media - mean(media))^2)
desvmediajackk <- sqrt(varmediajackk)
varmedianajackk <- ((n - 1)/n) * sum((mediana - mean(mediana))^2)
desvmedianajackk <- sqrt(varmedianajackk)
desvmedia
```

```
## [1] 0.1555792
```

```r
desvmediajackk
```

```
## [1] 0.1610397
```

```r
desvmedianajackk
```

```
## [1] 0.1834505
```

```r
# Aproximaciones sesgo
sesgomediajackk <- (n - 1) * (mean(media) - mean(muestra))
# sesgomedianajackk <- (n - 1) * (mean(mediana) - muestra[8])
sesgomedianajackk <- (n - 1) * (mean(mediana) - median(muestra))
sesgomediajackk
```

```
## [1] 0
```

```r
sesgomedianajackk
```

```
## [1] -0.003733333
```

También podríamos emplear la la función `jackknife` del paquete `bootstrap`:


```r
library(bootstrap)
resmedia <- jackknife(muestra, mean)
resmedia
```

```
## $jack.se
## [1] 0.1610397
## 
## $jack.bias
## [1] 0
## 
## $jack.values
##  [1] 0.8526429 0.8498571 0.8445714 0.8442857 0.8435714 0.8316429 0.8265000
##  [8] 0.8192143 0.8120000 0.7885714 0.7850000 0.7807143 0.7585714 0.7285714
## [15] 0.7142857
## 
## $call
## jackknife(x = muestra, theta = mean)
```

```r
resmediana <- jackknife(muestra, median)
resmediana
```

```
## $jack.se
## [1] 0.1834505
## 
## $jack.bias
## [1] -0.003733333
## 
## $jack.values
##  [1] 0.6615 0.6615 0.6615 0.6615 0.6615 0.6615 0.6615 0.6105 0.5600 0.5600
## [11] 0.5600 0.5600 0.5600 0.5600 0.5600
## 
## $call
## jackknife(x = muestra, theta = median)
```

Estos resultados pueden compararse con los obtenidos en el Ejemplo \@ref(exm:estimacion-boot-precision) empleando bootstrap.
En general las aproximaciones jackknife son adecuadas para el caso de estadísticos
"suaves", como la media, pero pueden ser inconsistentes cuando no lo son,
como es el caso de la mediana.


<!-- 
| Estimación | $\bar{X}=0.804$   |    $x_{\left( 8 \right)}=0.611$    | 
| ------------- | ------------- | ------------- |
| Precisión (plugin)   |    $0.155579$   |    ?    | 
| Precisión (jackknife)   |    $0.161040$   |    $0.183451$    | 
| Sesgo (jackknife)   |    $0$   |    $-0.003733$ |

Recordemos que, para el bootstrap teníamos (aproximado por Monte Carlo):

| Estimación | $\bar{X}=0.804$   |    $x_{\left( 8 \right)}=0.611$    | 
| ------------- | ------------- | ------------- |
| Precisión (plugin)   |    $0.155579$   |    ?    | 
| Precisión (bootstrap)   |    $0.155391$   |    $0.2501718$    | 
| Sesgo (bootstrap)   |    $0.000613$   |    $0.0473460$ |
 -->
 
