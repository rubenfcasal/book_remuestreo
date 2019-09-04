# Estimación de la precisión y el sesgo de un estimador {#cap2}





Uno de los problemas más interesantes que pueden ser abordados desde la
perspectiva de los métodos de remuestreo es el de la estimación del
sesgo y la precisión de un estimador. En dicho contexto surgió el método
Jackknife (bastante antes que el bootstrap), que, en ese sentido, puede
considerarse el método de remuestreo más antiguo como tal.

## Estimación bootstrap de la precisión y el sesgo de un estimador {#cap2-boot}

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
como se  mostró en la Sección \@ref(cap1-implementacion).

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


\BeginKnitrBlock{example}\iffalse{-91-65-112-114-111-120-105-109-97-99-105-243-110-32-98-111-111-116-115-116-114-97-112-32-100-101-32-108-97-32-112-114-101-99-105-115-105-243-110-32-100-101-32-101-115-116-105-109-97-99-105-111-110-101-115-32-100-101-108-32-116-105-101-109-112-111-32-100-101-32-118-105-100-97-32-109-101-100-105-111-32-100-101-32-109-105-99-114-111-111-114-103-97-110-105-115-109-111-115-93-}\fi{}<div class="example"><span class="example" id="exm:estimacion-boot-precision"><strong>(\#exm:estimacion-boot-precision)  \iffalse (Aproximación bootstrap de la precisión de estimaciones del tiempo de vida medio de microorganismos) \fi{} </strong></span><br> \vspace{0.5cm}

Continuando con el ejemplo de los tiempos de vida de microorganismos,
supongamos que queremos estimar
la precisión de dos estimadores de su vida media: media muestral y
mediana muestral, a partir de los datos observados: 0.143, 0.182, 0.256, 0.260, 0.270,
0.437, 0.509, 0.611, 0.712, 1.04, 1.09, 1.15, 1.46, 1.88, 2.08.</div>\EndKnitrBlock{example}

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
## [1] 0.1570513
```

```r
desvmedianaboot
```

```
## [1] 0.2517709
```

```r
# Aproximaciones sesgo
sesgomediaboot <- mean(media) - mean(muestra)
sesgomedianaboot <- mean(mediana) - muestra[8]
sesgomediaboot
```

```
## [1] 0.0004837933
```

```r
sesgomedianaboot
```

```
## [1] 0.0482999
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
##      original      bias    std. error
## t1* 0.8053333 -0.00046904   0.1550047
## t2* 0.6110000  0.04636960   0.2498358
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
##      original      bias    std. error
## t1* 0.8053333 -0.00046904   0.1550047
## t2* 0.6110000  0.04636960   0.2498358
```

