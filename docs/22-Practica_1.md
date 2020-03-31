# Práctica 1: Correlación y regresión lineal {#practica1}




Esta práctica debe entregarse en formato pdf, 
incluyendo el código R utilizado, las correspondientes salidas 
y los comentarios (o interpretaciones de los resultados) pertinentes
(para ello se recomienda emplear RMarkdown, 
a partir de un fichero *.Rmd* o un fichero *.R* mediante spin).
Se debe establecer la semilla para la generación de las réplicas bootstrap
igual a los dos último dígitos del DNI mediante la función `set.seed()`.

En esta práctica se empleará el conjunto de datos `Prestige` del paquete `car`, 
considerando como variable respuesta `prestige` (puntuación de ocupaciones obtenidas 
a partir de una encuesta) y como variables explicativas: 
`income` (media de ingresos en la ocupación) y `education` (media de los años de educación).


```r
library(car)
data(Prestige)
# ?Prestige
```


## Inferencia sobre el coeficiente de correlación lineal

Supongamos que estamos interesados en estudiar la correlación lineal 
entre las variables `income` ($X$) y `prestige` ($Y$) 
del conjunto de datos `Prestige`. 
Para ello podemos considerar el coeficiente de correlación lineal de Pearson:
$$\rho =\frac{ Cov \left( X, Y \right) }
{ \sigma \left( X \right) \sigma \left( Y \right) }$$
Su estimador es el coeficiente de correlación muestral:
$$r=\frac{\sum_{i=1}^{n}(x_i-\overline{x})(y_i-\overline{y})}
{\sqrt{ \sum_{i=1}^{n}(x_i-\overline{x})^{2}} 
\sqrt{\sum_{i=1}^{n}(y_i-\overline{y})^{2}}},$$
que podemos calcular en `R` empleando la función `cor()`:

```r
# with(Prestige, cor(income, prestige))
cor(Prestige$income, Prestige$prestige)
```

```
## [1] 0.7149057
```
Para realizar inferencias sobre $\rho$ podemos emplear la función
`cor.test()`:

```r
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
Esta función^[Se puede obtener el código tecleando en la consola `stats:::cor.test.default`.],
además de realizar el contraste $H_0: \rho = 0$^[Empleando el estadístico 
$\frac{r\sqrt{n - 2}}{\sqrt{1 - r^2}} \underset{aprox}{\sim } \mathcal{t}_{n-2}$
bajo la hipótesis nula de que la verdadera correlación es cero.], 
calcula un intervalo de confianza a partir de la transformación $Z$ de Fisher:
$$Z = \frac{1}{2}\ln \left( \frac{1+r}{1-r} \right) = \operatorname{arctanh}(r),$$
que es una transformación (aprox.) normalizadora y estabilizadora de varianza.
Suponiendo que $(X, Y)$ es normal bivariante y que hay independencia entre las observaciones:
$$Z \sim \mathcal{N}\left( \frac{1}{2}\ln \left( \frac{1+\rho}{1-\rho} \right), \frac{1}{n-3} \right).$$
El intervalo de confianza asintótico se obtiene empleando esta transformación, 
la aproximación normal tradicional y la transformación inversa:
$$r = \frac{\exp(2Z)-1}{\exp(2Z)+1} = \operatorname{tanh}(Z)$$.

### Ejercicio (para entregar)

1.  Emplear bootstrap uniforme para aproximar el sesgo y la precisión de la transformación
    $Z$ de Fisher del coeficientes de correlación entre `income` y `prestige`:
    
    ```r
    data <- Prestige
    r <- cor(data$income, data$prestige)
    z <- atanh(r)
    z
    ```
    
    ```
    ## [1] 0.8971466
    ```
    En La aproximación implementada en `cor.test()` se supone que 
    el sesgo es nulo y que el error estándar es:
    
    ```r
    n <- nrow(data)
    sigma <- 1/sqrt(n - 3)
    sigma
    ```
    
    ```
    ## [1] 0.1005038
    ```
    Aunque realmente el sesgo de esta transformación no es nulo 
    (aprox. $r/(2(n-1))$ según Pearson y Hartley, 1954, p. 29).
    

```r
library(boot)

# set.seed(DNI)
# ...
```

### Ejemplo: Intervalo de confianza bootstrap en escala transformada

Empleando la función `boot.ci()` del paquete `boot` se pueden obtener intervalos de confianza
calculados en una escala transfomada del estadístico estableciendo los parámetros:

- `h`: función vectorial que define la transformación. 
  Los intervalos se calculan en la escala de $h(t)$ y se aplica la función inversa 
  (si se especifica) para transformarlos a la escala original.

- `hinv`: (opcional) función inversa de la transformación 
  (si no se especifica solo se calculan los intervalos en la escala transformada). 

- `hdot`: (opcional) función derivada de la transformación 
  (empleada por algunos métodos para aproximar la varianza en la escala transformada 
  mediante el método delta).

Por ejemplo, para considerar la transformación $Z$ de Fisher del coeficiente de correlación, 
se podría emplear el siguiente código:

```r
library(boot)

statistic <- function(data, i){
  remuestra <- data[i, ]
  r <- cor(remuestra$income, remuestra$prestige)
  r
}

set.seed(1)
res.boot <- boot(Prestige, statistic, R = 1000)
res.boot
```

```
## 
## ORDINARY NONPARAMETRIC BOOTSTRAP
## 
## 
## Call:
## boot(data = Prestige, statistic = statistic, R = 1000)
## 
## 
## Bootstrap Statistics :
##      original      bias    std. error
## t1* 0.7149057 0.006306905  0.04406473
```

```r
plot(res.boot)
```



\begin{center}\includegraphics[width=0.7\linewidth]{22-Practica_1_files/figure-latex/unnamed-chunk-8-1} \end{center}

```r
h <- function(t) atanh(t)
hdot <- function(t) 1/(1 - t^2)
hinv <- function(t) tanh(t)

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

```r
# La distribución bootstrap en la escala transformada
# se aproxima más a la normalidad...
ht <- h(res.boot$t)
hist(ht, freq = FALSE, breaks = "FD",
     main = "Distribución bootstrap en la escala transformada")
curve(dnorm(x, mean=mean(ht), sd=sd(ht)), lty = 2, add = TRUE)
```



\begin{center}\includegraphics[width=0.7\linewidth]{22-Practica_1_files/figure-latex/unnamed-chunk-8-2} \end{center}

## Inferencia sobre un modelo de regresión lineal {#bootmod}

Para ajustar un modelo de regresión lineal que trate de explicar `prestige` 
a partir de `income` y `education` podemos emplear el siguiente código:

```r
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

Como se muestra en Canty (2002), para realizar inferencias sobre modelos 
de regresión se podrían emplear dos algoritmos bootstrap (entre otros,
incluyendo el *Wild Bootstrap*, descrito en la Sección \@ref(wild-bootstrap) 
para modelos de regresión no paramétricos).
El primero consistiría en utilizar directamente bootstrap uniforme,
remuestreando las observaciones, y sería adecuado para el caso de 
diseño aleatorio. Por ejemplo, empleando el siguiente código
podríamos realizar inferencias sobre el coeficiente de determinación ajustado:

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

```r
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

La otra alternativa, más adecuada para el caso de diseño fijo, sería
lo que se conoce como *remuestreo residual*, *remuestreo basado en modelos* o *bootstrap semiparamétrico*.
En esta aproximación se mantienen fijos los valores de las variables
explicativas y se remuestrean los residuos:
$$\mathbf{r} = \mathbf{Y} - X\hat{\mathbf{\beta}} = \mathbf{Y} - \hat{\mathbf{Y}}$$
obteniéndose las réplicas bootstrap:
$$\mathbf{Y}^{\ast} = \hat{\mathbf{Y}} + \mathbf{r}^{\ast}$$
Por ejemplo, adaptando el código en Canty (2002) para este conjunto de 
datos, podríamos emplear:

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

Sin embargo, la variabilidad de los residuos no reproduce la de los
verdaderos errores, por lo que podría ser preferible (especialmente 
si el tamaño muestral es pequeño) emplear la modificación descrita
en Davison y Hinkley (1997, Alg. 6.3, p. 271).
La idea es remuestrear los residuos reescalados y centrados:
$$e_i = \frac{r_i}{\sqrt{1 - h_{ii} }} - \bar{e}$$
donde $h_{ii}$ es el valor de influencia o leverage, 
el elemento $i$-ésimo de la diagonal de la matriz de proyección
$H = X\left( X^{\prime }X\right)^{-1}X^{\prime }$.

En `R` podríamos obtener estos residuos mediante los comandos:

```r
pres.dat$sres <- residuals(modelo)/sqrt(1 - hatvalues(modelo))
pres.dat$sres <- pres.dat$sres - mean(pres.dat$sres)
```
Sim embargo puede ser más cómodo emplear la función `Boot()` del paquete `car`
(que internamente llama a la función `boot()`), 
como se describe en el apéndice "Bootstrapping Regression Models in R" del libro 
"An R Companion to Applied Regression" de Fox y Weisberg (2018), disponible
[aquí](https://socialsciences.mcmaster.ca/jfox/Books/Companion/appendices/Appendix-Bootstrapping.pdf).

Esta función es de la forma:

```r
Boot(object, f=coef, labels=names(f(object)), R=999, method=c("case", "residual"))
```
donde:

- `object`: es un objeto que contiene el ajuste de un modelo de regresión.

- `f`: es la función de estadísticos (utilizando el ajuste como argumento).

- `method`: especifíca el tipo de remuestreo: remuestreo de observaciones (`"case"`)
  o de residuos (`"residual"`), empleando la modificación descrita anteriormente.

### Ejercicio (para entregar)

2. Emplear la función `Boot()` del paquete `car` para hacer inferencia sobre 
  el coeficiente de determinación ajustado del modelo de regresión lineal 
  que explica `prestige` a partir de `income` y `education` 
  (obtener una estimación del sesgo y de la predicción,
  y opcionalmente una estimación por intervalo de confianza).


```r
library(car)

# set.seed(DNI)
# ...
```
