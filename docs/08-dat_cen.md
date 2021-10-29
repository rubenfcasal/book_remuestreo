# El Bootstrap con datos censurados {#bootcen}




En este capítulo se hace una introducción a los datos censurados y se
presentan diversos métodos de remuestreo para este contexto, analizando
la validez de los mismos.

## Introducción a los datos censurados

Considérese una variable de interés, $X$, no negativa que no siempre es
posible observar (por ejemplo un tiempo de vida) pues, en ocasiones,
ocurre otro fenómeno previo, cuyo tiempo hasta su ocurrencia, $C$, puede
ser anterior a la variable de interés (es decir, $C<X$). Cuando $X$ es
un tiempo de vida ante una enfermedad mortal, la variable $C$ suele
representar el tiempo hasta el fin del estudio, el tiempo hasta que el
individuo fallezca por otra causa o el tiempo hasta que se produce una
pérdida de seguimiento. Es habitual definir el indicador de no censura
$\delta
=\mathbf{1}_{\left\{ X\leq C\right\} }$. Si $C<X$ diremos que la
observación es censurada y sólo seremos capaces de observar $C$ junto
con el valor de $\delta$. Cuando $X\leq C$ entonces somos capaces de
observar la variable de interés y además el valor de $\delta$.

En resumen, en lugar de observar la muestra 
$\left( X_1, X_2, \ldots, X_n \right)$, sólo podemos observar 
$$\left( \left( T_1, \delta _1 \right), \left( T_2, \delta _2 \right),
\ldots ,\left( T_n, \delta_n \right) \right),$$
siendo $T_i=\min \left\{ X_i,C_i\right\}$
los tiempos de vida observados y
$$\delta _i=\mathbf{1}_{\left\{ X_i\leq
C_i\right\} }=\mathbf{1}_{\left\{ T_i=X_i\right\} }$$ 
los indicadores de censura, para $i=1,2,\ldots ,n$. 
En el modelo de censura aleatoria por la derecha, 
que es el más habitual, se supone que $X_i$
y $C_i$ ($i=1,2,\ldots ,n$) son independientes. 
Además $\left( X_1, X_2, \ldots ,X_n \right)$ son mutuamente independientes, 
como también lo son $\left( C_1, C_2,\ldots ,C_n \right)$.

Denotando por $F$ (respectivamente $G$ y $H$) la función de distribución
de la variable aleatoria $X$ (respectivamente $C$ y $T$), la condición
de independencia implica que 
$$1-H\left( t \right) =\left( 1-F\left(
t \right) \right) \left( 1-G\left( t \right) \right).$$

### Estimador de Kaplan-Meier

Puede verse fácilmente que, bajo censura aleatoria por la derecha, la
distribución empírica $F_n\left( t \right) =\frac{1}{n}
\sum_{i=1}^{n}\mathbf{1}_{\left\{ T_i\leq t\right\} }$ deja de ser
consistente. En este contexto el estimador no paramétrico de máxima
verosimilitud de la función de distribución es el estimador
límite-producto, propuesto por Kaplan y Meier (1958), obtenido
a partir de la función de supervivencia 
($S\left( t \right) = 1-F\left( t \right)$):
$$\hat{S}\left( t \right) = 1-\hat{F}\left( t \right) =
\prod_{T_{(i)}\leq t}\left( \frac{n-i}{n-i+1} \right)^{\delta _{(i)}},$$
siendo $\left( T_{(1)},T_{\left( 2 \right)},\ldots ,T_{\left(
n \right)} \right)$ la muestra de estadísticos ordenados de los tiempos
de vida observados y $\left( \delta _{(1)},\delta _{\left(2 \right)},
\ldots ,\delta _{(n)} \right)$ los correspondientes concomitantes.

\BeginKnitrBlock{example}\iffalse{-91-69-115-116-105-109-97-99-105-243-110-32-100-101-32-75-97-112-108-97-110-45-77-101-105-101-114-93-}\fi{}
<span class="example" id="exm:kaplan-meier"><strong>(\#exm:kaplan-meier)  \iffalse (Estimación de Kaplan-Meier) \fi{} </strong></span>
<br> \vspace{0.5cm}

Se observan los datos censurados: $\left( 2.1,0 \right)$,
$\left(3.2,1 \right)$, $\left( 1.2,1 \right)$, $\left( 4.3,0 \right)$, $\left(
1.8,1 \right)$, $\left( 3.9,1 \right)$, $\left( 2.7,0 \right)$, $\left(
2.5,1 \right)$. El estimador resulta:

$$\hat{F}\left( t \right) =\left\{ 
\begin{array}{ll}
0& \text{si } t<1.2 \\ 
0.125& \text{si } 1.2\leq t<1.8 \\ 
0.25& \text{si } 1.8\leq t<2.5 \\ 
0.4& \text{si } 2.5\leq t<3.2 \\ 
0.6& \text{si } 3.2\leq t<3.9 \\ 
0.8& \text{si } 3.9\leq t
\end{array}
\right.$$
\EndKnitrBlock{example}


En `R` se recomienda emplear el paquete `survival` para el análisis de
datos censurados. Podemos utilizar la función `survfit()` para obtener 
la estimación Kaplan-Meier de la función de supervivencia 
(y a partir de ella la de la distribución).
En este caso podríamos utilizar el siguiente código [Figura \@ref(fig:survival)]:

```r
datcen <- data.frame(t = c(2.1, 3.2, 1.2, 4.3, 1.8, 3.9, 2.7, 2.5), 
                 cen = c(0, 1, 1, 0, 1, 1, 0, 1))

library(survival)
fit <- survfit(Surv(t, cen)~1, data = datcen)
summary(fit)
```

```
## Call: survfit(formula = Surv(t, cen) ~ 1, data = datcen)
## 
##  time n.risk n.event survival std.err lower 95% CI upper 95% CI
##   1.2      8       1    0.875   0.117       0.6734            1
##   1.8      7       1    0.750   0.153       0.5027            1
##   2.5      5       1    0.600   0.182       0.3315            1
##   3.2      3       1    0.400   0.203       0.1477            1
##   3.9      2       1    0.200   0.174       0.0363            1
```

```r
old.par <- par(mfrow = c(1, 2))
plot(fit, main = "Método plot de un objeto 'survfit'")
legend("bottomleft",  c("supervivencia", "conf.int"), lty = 1:2)

with(fit, {
  plot(c(0, time), c(1, surv), type = "s", lty = 2,
       main = "Estimaciones funciones supervicencia y distribución", 
       xlab = "t", ylab = "", ylim = c(0, 1))
  lines(c(0, time), 1 - c(1, surv), type = "s")
  legend("bottomright",  c("supervivencia", "distribución"), lty = 2:1)
})
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{08-dat_cen_files/figure-latex/survival-1} 

}

\caption{Estimaciones Kaplan-Meier de la función de supervivencia y de la función de distribución.}(\#fig:survival)
\end{figure}

```r
par(old.par)
```


### Distribución asintótica del estimador de Kaplan-Meier

El estimador de Kaplan-Meier sólo otorga pesos positivos a los datos no
censurados, aunque la forma de distribuirse los datos censurados en
medio de los no censurados afecta a los pesos de estos últimos. Por otra
parte, en ausencia de censura (es decir $\delta _i=1$,
$i=1,2,\ldots ,n$), el estimador de Kaplan-Meier coincide con la
distribución empírica.

La obtención de la propiedades de sesgo y varianza asintóticos y
distribución límite del estimador de Kaplan-Meier es mucho más laboriosa
que en el caso de la distribución empírica, en un contexto sin censura.
Esto es así porque el estimador de Kaplan-Meier deja de ser una suma de
variables iid, como sí ocurre con la empírica.

Breslow y Crowley (1974) obtienen el siguiente resultado para la distribución
límite para el estimador de Kaplan-Meier:
$$\sqrt{n}\left( \hat{F}\left( t \right) -F\left( t \right) \right) 
\overset{d}{\longrightarrow} \mathcal{N}\left( 0,\sigma^2\left( t \right)
\right),$$
siendo
$$\begin{aligned}
\sigma^2\left( t \right) &= \left( 1-F\left( t \right) \right)
^2\int_{0}^{t}\frac{dH_1\left( u \right)}{\left( 1-H\left( u \right)
 \right)^2}\text{, }t\leq H^{-1}(1) , \\
H_1\left( u \right) &= P\left( X\leq u,X\leq C \right) =P\left( T\leq
u,\delta =1 \right).
\end{aligned}$$
También existen resultados de convergencia en distribución del 
proceso estocástico
$$\left\{ \sqrt{n}\left( \hat{F}\left( t \right) - F\left( t \right)
 \right) : t \in \left[ 0,H^{-1}(1) \right] \right\}$$
a un proceso gaussiano límite.


## Remuestreos Bootstrap en presencia de censura

Estos métodos tratan del mecanismo bootstrap para aproximar la
distribución de un estadístico, $R\left( \mathbf{T},
\boldsymbol{\delta} \right)$, siendo $\mathbf{T}=\left( T_1, T_2,
\ldots,T_n \right)$ y $\boldsymbol{\delta}=\left( \delta _1,\delta_2,
\ldots ,\delta _n \right)$. Los dos siguientes métodos de
remuestreo fueron propuestos por Efron (1981).

### El bootstrap simple

Procede de la siguiente forma:

1.  Construir la distribución empírica bidimensional,
    $F_n^{T,\delta }$, de la muestra
    $\left\{ \left( T_1,\delta _1 \right), \left( T_2,\delta _2 \right), \ldots, 
    \left( T_n,\delta _n \right) \right\}$.

2.  Arrojar remuestras $\left\{ \left( T_1^{\ast},\delta _1^{\ast} \right), 
    \left( T_2^{\ast},\delta _2^{\ast} \right), \ldots, 
    \left( T_n^{\ast},\delta _n^{\ast} \right) \right\}$
    a partir de dicha distribución empírica. Esto es tanto como decir
    que$$P^{\ast}\left( \left( T^{\ast},\delta^{\ast} \right) =\left( T_i,\delta
    _i \right) \right) =\frac{1}{n}\text{, para }i=1,2,\ldots ,n\text{.}$$

3.  Evaluar el estadístico de interés en el vector que contiene la
    remuestra bootstrap: $R^{\ast}=R\left( \mathbf{T}^{\ast},
    \boldsymbol{\delta}^{\ast} \right)$, con
    $\mathbf{T}^{\ast}
    =\left( T_1^{\ast},T_2^{\ast},\ldots ,T_n^{\ast} \right)$ y
   
    $\boldsymbol{\delta}^{\ast}=\left( \delta _1^{\ast},\delta
    _2^{\ast},\ldots ,\delta _n^{\ast} \right)$.

4.  Aproximar la distribución en el muestreo del estadístico
    $R\left( \mathbf{T}, \boldsymbol{\delta} \right)$ por la
    distribución en el remuestreo de $R\left( 
    \mathbf{T}^{\ast},\boldsymbol{\delta}^{\ast} \right)$.

Este método es de muy rápida implementación y ejecución.

### El bootstrap obvio {#bootcen-obvio}

Para detallar el método es necesario definir el estimador de
Kaplan-Meier, $\hat{G}\left( t \right)$, de la variable censurante, a
partir
de$$1-\hat{G}\left( t \right) =\prod_{T_{(i)}\leq t}\left( \frac{n-i
}{n-i+1} \right)^{1-\delta _{(i)}}.$$ Observemos que este
estimador es totalmente semejante al de Kaplan-Meier de la variable de
interés pero simplemente reemplazando cada valor
$\delta _{(i)}$ por $1-\delta _{(i)}$.

El mecanismo de remuestreo procede como sigue:

1.  Construir los estimadores de Kaplan-Meier de las distribuciones de
    la variable de interés, $\hat{F}\left( t \right)$, y de la variable
    censurante, $\hat{G}\left( t \right)$.

2.  Para cada índice $i=1,2,\ldots ,n$, arrojar observaciones bootstrap
    independientes, $X_i^{\ast}$ con distribución $\hat{F}\ $y
    $C_i^{\ast}$ con distribución $\hat{G}$.

3.  Definir $T_i^{\ast}=\min \left\{ X_i^{\ast},C_i^{\ast}\right\}$ y
    $\delta_i^{\ast}=\mathbf{1}_{\left\{ X_i^{\ast}\leq
    C_i^{\ast}\right\}}$, para $i = 1, 2, \ldots, n$, 
    y considerar la remuestra bootstrap
    $\left( \mathbf{T}^{\ast},\boldsymbol{\delta}^{\ast}\right)$, con
    $\mathbf{T}^{\ast}=\left( T_1^{\ast},T_2^{\ast}, \ldots,
    T_n^{\ast} \right)$ y $\boldsymbol{\delta}^{\ast} = \left(
    \delta_1^{\ast}, \delta_2^{\ast},\ldots ,\delta_n^{\ast} \right)$.

4.  Aproximar la distribución en el muestreo del estadístico
    $R\left( \mathbf{T},\boldsymbol{\delta} \right)$ por la
    distribución en el remuestreo de su análogo bootstrap, $R\left( 
    \mathbf{T}^{\ast},\boldsymbol{\delta}^{\ast} \right)$.

Obviamente, este método de remuestreo imita fielmente el modelo de datos
censurados por la derecha. Su ejecución es considerablemente más lenta
que la del método simple, pues necesita de la construcción de los
estimadores de Kaplan-Meier, de la obtención de remuestras a partir de
ellos y de algunos cálculos adicionales.


## Relaciones entre los métodos de remuestreo bajo censura

### Equivalencia entre el bootstrap simple y el obvio

Es fácil demostrar que el bootstrap simple y el obvio son planes de
remuestreo equivalentes (cuando se supone que en la muestra no existe
ninguna observación no censurada que esté empatada con otra censurada).
Esta equivalencia se establece en el sentido de que la distribución
bootstrap de $\left( T^{\ast},\delta^{\ast} \right)$ es la misma para
cualquiera de los dos métodos.

Así, si $\left( T^{\ast},\delta^{\ast} \right)$ se genera mediante el
método obvio, entonces

$$\begin{aligned}
P^{\ast}\left( T^{\ast}>t \right) &= P^{\ast}\left( X^{\ast}>t,C^{\ast
}>t \right) \\
&= P^{\ast}\left( X^{\ast}>t \right) P^{\ast}\left( C^{\ast}>t \right)
=\left( 1-\hat{F}\left( t \right) \right) \left( 1-\hat{G}\left( t \right)
 \right) \\
&= \left[ \prod_{T_{(i)}\leq t}\left( \frac{n-i}{n-i+1} \right)
^{\delta _{(i)}}\right] \left[ \prod_{T_{(i)}\leq
t}\left( \frac{n-i}{n-i+1} \right)^{1-\delta _{(i)}}\right] \\
&= \prod_{T_{(i)}\leq t}\frac{n-i}{n-i+1}=\prod_{i=1}^{\#\left
\{ T_{(j)}\leq t\right\} }\frac{n-i}{n-i+1} \\
&= \frac{n-1}{n}\cdot \frac{n-2}{n-1}\cdot \cdots \cdot \frac{n-\#\left\{
T_{(j)}\leq t\right\} }{n-\#\left\{ T_{(j)}\leq
t\right\} +1} \\
&= \frac{n-\#\left\{ T_{(j)}\leq t\right\} }{n}=1-H_n\left(
t \right) =\frac{\#\left\{ T_{(j)}>t\right\} }{n},
\end{aligned}$$

siendo $H_n\left( t \right)$ la distribución empírica de la muestra
$\left( T_1,T_2,\ldots ,T_n \right)$.

Esto demuestra que la distribución bootstrap marginal de $T^{\ast}$ es
la misma para ambos remuestreos. Sólo resta probar pues que la
distribución condicionada
$\left. \delta^{\ast}\right\vert _{T^{\ast}=T_i}$ es idéntica en
ambos casos. Pero esto es inmediato ya que, en los dos remuestreos esa
distribución condicionada es la degenerada en el valor
$\delta _i$.

### El bootstrap de Reid {#bootcen-reid}

Es otro método alternativo propuesto por Reid (1981). Consta de los
siguientes pasos:

1.  Construir el estimador de Kaplan-Meier, $\hat{F}\left( t \right)$,
    de la muestra original.

2.  Arrojar remuestras bootstrap (todas formadas por observaciones no
    censuradas, $T_i^{\ast}$, $i=1,2,\ldots ,n$) a partir de
    $\hat{F}\left(t \right)$.

3.  Aproximar la distribución en el muestreo de $R\left( 
    \mathbf{T},\boldsymbol{\delta} \right)$, por la
    distribución bootstrap de
    $R\left( \mathbf{T}^{\ast},\mathbf{1} \right)$, 
    siendo $\mathbf{1}$ el vector formado
    por $n$ unos.

### Validez de los planes de remuestreo

Akritas (1986) demuestra que los procesos bootstrap
$$\begin{aligned}
\sqrt{n}\left( \hat{F}^{\ast}_{Efron}\left( t \right) 
- \hat{F}\left( t \right) \right) \\ 
\sqrt{n}\left( \hat{F}^{\ast}_{Reid} \left( t \right) 
- \hat{F}\left( t \right) \right)
\end{aligned}$$ 
tienden a sendos procesos límite distintos. Aquí
$\hat{F}^{\ast}_{Efron}$ denota la versión bootstrap del estimador
de Kaplan-Meier bajo el remuestreo de Efron (cualquiera de ellos, ya que
el remuestreo simple y el obvio son equivalentes) y
$\hat{F}^{\ast}_{Reid}$ es la correspondiente versión bootstrap del
estimador de Kaplan-Meier bajo el remuestreo de Reid (una distribución
empírica, al fin y al cabo, porque en el remuestreo de Reid todas las
observaciones son no censuradas).

Además el proceso límite del estimador de Kaplan-Meier, $\sqrt{n}
\left( \hat{F}\left( t \right) -F\left( t \right) \right)$, es el mismo
que el del bootstrap de Efron. Como consecuencia el remuestreo de Efron
es consistente y el de Reid es inconsistente.

## Implementación en `R` (con los paquetes `boot` y `survival`)

La función `censboot()` del paquete `boot` implementa distintos métodos 
de remuestreo para datos censurados. Por defecto utiliza el bootstrap simple
(`sim = "ordinary"`) y su uso es prácticamente igual al del bootstrap uniforme 
con la función `boot()` (descrita en la Sección \@ref(intro-pkgboot)), 
la única diferencia es que la función `statistic` solo tiene los datos 
como único parámetro (aunque en este caso podríamos emplear también
la función `boot()`). 

### Bootstrap simple

Como ejemplo utilizaremos el conjunto de datos `channing` del paquete `boot`, 
que contiene la edad de entrada y de partida o muerte de las personas 
que pasaron por el centro de retiro 'Channing House' (Palo Alto, California), 
desde su apertura en 1964 hasta el 1 de julio de 1975
(ver Sección 3.5 y 'Practical 3.2', de Davison y Hinkley, 1997).
En primer lugar consideraremos únicamente la muestra de hombres:


```r
# Datos
library(boot)
data(channing)
# Calcular edad (de partida o muerte) en años
channing$age <- (channing$entry + channing$time)/12
# Seleccionar hombres (y de paso hacer que `index = c(1, 2)` para `censboot()`)
chan <- subset(channing, sex=="Male", c(age, cens))

# Estimación supervivencia
library(survival)
chan.F <- survfit(Surv(age, cens)~1, data = chan)
chan.F
```

```
## Call: survfit(formula = Surv(age, cens) ~ 1, data = chan)
## 
##       n  events  median 0.95LCL 0.95UCL 
##    97.0    46.0    87.0    85.8    90.4
```

```r
# plot(chan.F)
# Estimaciones de interés
with(chan.F, 
    c(s75 = max(surv[time > 75]), s85 = max(surv[time > 85]),
      p75 = min(time[surv <= 0.75]), p50 = min(time[surv <= 0.5])) 
)
```

```
##        s75        s85        p75        p50 
##  0.9160745  0.6347541 82.4166667 87.0000000
```

```r
# Bootstrap
# library(boot)
chan.stat <- function(data) {
    s <- survfit(Surv(age, cens)~1, data = data)
    with(s, c(s75 = max(surv[time > 75]), s85 = max(surv[time > 85]),
            p75 = min(time[surv <= 0.75]), p50 = min(time[surv <= 0.5])))
}
set.seed(1)
chan.boot <- censboot(chan, chan.stat, R = 199) # sim = "ordinary"
chan.boot
```

```
## 
## CASE RESAMPLING BOOTSTRAP FOR CENSORED DATA
## 
## 
## Call:
## censboot(data = chan, statistic = chan.stat, R = 199)
## 
## 
## Bootstrap Statistics :
##       original       bias    std. error
## t1*  0.9160745 -0.006995081  0.03133656
## t2*  0.6347541 -0.003538939  0.05595748
## t3* 82.4166667 -0.040619765  1.22517032
## t4* 87.0000000  0.195142379  1.07267425
```

### Otros métodos de remuestreo

La función `censboot()` implementa otros dos métodos de remuestreo, 
`sim = c("cond", "weird")`, aunque en ambos casos hay que establecer en el
parámetro `F.surv` la estimación de Kaplan-Meier de la supervivencia y,
si `sim = "cond"`, la correspondiente a la variable censurante en `G.surv`.
Se recomienda estimarlas con la función `survfit()` del paquete `survival`
(ver Figura \@ref(fig:survfit-f-g)):


```r
# Estimación supervivencia variable censurante
chan.G <- survfit(Surv(age, 1-cens)~1, data = chan)
# Representación
old.par <- par(mfrow = c(1, 2))
plot(chan.F, main = "Supervivencia (edad)", mark.time = TRUE, 
    xlim = c(60, 100))
plot(chan.G, main = "Supervivencia variable censurante (partida)", 
     mark.time = TRUE, xlim = c(60, 100))
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{08-dat_cen_files/figure-latex/survfit-f-g-1} 

}

\caption{Estimaciones de la supervivencia (izquierda; indicando los tiempos de las observaciones censuradas) y de la variable censurante (derecha; indicando los de las no censuradas).}(\#fig:survfit-f-g)
\end{figure}

```r
par(old.par)
```

En el *boostrap condicional* (`sim = "cond"`) se condiciona el muestreo al 
patrón de censura observado (en lugar de fijarlo a $\mathbf{1}$ como en el bootstrap de Reid; Sección \@ref(bootcen-reid)). 
El mecanismo es similar al del bootstrap obvio (Sección \@ref(bootcen-obvio)):

1.  Construir los estimadores de Kaplan-Meier de las distribuciones de
    la variable de interés, $\hat{F}\left( t \right)$, y de la variable
    censurante, $\hat{G}\left( t \right)$.

2.  Para cada índice $i=1,2,\ldots ,n$, generar $X_i^{\ast}$ independientes
    con distribución $\hat{F}$. 
    
3.  Si la $i$-ésima observación está censurada ($\delta_i=0$) se toma 
    $C_i^{\ast}=X_i$ y si no ($\delta_i=1$) se genera un valor de la estimación
    de la distribución de la variable censurante condicionada a $C > X_i$:
    $$\hat G \left(\left. t \ \right\vert_{\ t > X_i} \right) 
    = \frac{\hat G(t) - \hat G(X_i)}{1- \hat G(X_i)}.$$

3.  Definir $T_i^{\ast}=\min \left\{ X_i^{\ast},C_i^{\ast}\right\}$ y
    $\delta_i^{\ast}=\mathbf{1}_{\left\{ X_i^{\ast}\leq
    C_i^{\ast}\right\}}$, para $i = 1, 2, \ldots, n$, 
    y considerar la remuestra bootstrap
    $\left( \mathbf{T}^{\ast},\boldsymbol{\delta}^{\ast}\right)$, con
    $\mathbf{T}^{\ast}=\left( T_1^{\ast},T_2^{\ast}, \ldots,
    T_n^{\ast} \right)$ y $\boldsymbol{\delta}^{\ast} = \left(
    \delta_1^{\ast}, \delta_2^{\ast},\ldots ,\delta_n^{\ast} \right)$.    
    
El otro método (`sim = "weird"`) es el denominado *weird bootstrap* 
(Andersen et al., 1993) que emplea la estimación de Nelson-Aalen de la 
función de riesgo acumulada para generar los valores 
(e.g. Sección 3.5.2 de Davison y Hinkley, 1997).

El siguiente código muestra un ejemplo de la aplicación de ambos métodos:

```r
chan.boot2 <- censboot(chan, chan.stat, R = 199, F.surv = chan.F, 
                  G.surv = chan.G, sim = "cond")
chan.boot2
```

```
## 
## CONDITIONAL BOOTSTRAP FOR CENSORED DATA
## 
## 
## Call:
## censboot(data = chan, statistic = chan.stat, R = 199, F.surv = chan.F, 
##     G.surv = chan.G, sim = "cond")
## 
## 
## Bootstrap Statistics :
##       original       bias    std. error
## t1*  0.9160745  0.002340590  0.02796666
## t2*  0.6347541 -0.001057025  0.05382609
## t3* 82.4166667  0.019681742  1.21793756
## t4* 87.0000000  0.286013400  0.95041994
```

```r
chan.boot3 <- censboot(chan, chan.stat, R = 199, F.surv = chan.F, 
                  sim = "weird")
chan.boot3
```

```
## 
## WEIRD BOOTSTRAP FOR CENSORED DATA
## 
## 
## Call:
## censboot(data = chan, statistic = chan.stat, R = 199, F.surv = chan.F, 
##     sim = "weird")
## 
## 
## Bootstrap Statistics :
##       original        bias    std. error
## t1*  0.9160745 -4.082326e-05  0.02825475
## t2*  0.6347541  1.639197e-03  0.05086460
## t3* 82.4166667  2.303183e-02  1.14079354
## t4* 87.0000000  2.458124e-01  0.96511648
```

## Ejercicios

\BeginKnitrBlock{exercise}\iffalse{-91-66-111-111-116-115-116-114-97-112-32-99-101-110-115-117-114-97-100-111-32-112-111-114-32-101-115-116-114-97-116-111-115-93-}\fi{}
<span class="exercise" id="exr:censboot-strata-ej"><strong>(\#exr:censboot-strata-ej)  \iffalse (Bootstrap censurado por estratos) \fi{} </strong></span>
Analizar el conjunto de datos `channing` completo, teniendo en cuenta el sexo como estrato
(i.e. `Surv(age, cens) ~ sex` y `strata = chan$sex`)

\EndKnitrBlock{exercise}


```r
# Datos
data(channing)
# Calcular edad (de partida o muerte) en años
channing$age <- (channing$entry + channing$time)/12
# Seleccionar variables
chan <-channing[c("age", "cens", "sex")]

# Estimación supervivencia
library(survival)
chan.F <- survfit(Surv(age, cens) ~ sex, data = chan)
chan.F
```

```
## Call: survfit(formula = Surv(age, cens) ~ sex, data = chan)
## 
##              n events median 0.95LCL 0.95UCL
## sex=Female 365    130     88    86.7    89.5
## sex=Male    97     46     87    85.8    90.4
```

```r
plot(chan.F, lty = 1:2, xlim = c(60, 100))
```

\begin{figure}[!htb]

{\centering \includegraphics[width=0.7\linewidth]{08-dat_cen_files/figure-latex/surv-strata-1} 

}

\caption{Estimaciones de la supervivencia.}(\#fig:surv-strata)
\end{figure}

```r
res <- summary(chan.F)
# res
str(res)
```

```
## List of 19
##  $ n            : int [1:2] 365 97
##  $ time         : num [1:146] 67 68.5 69.2 70 70.4 ...
##  $ n.risk       : num [1:146] 364 359 355 353 352 346 344 340 335 334 ...
##  $ n.event      : num [1:146] 1 1 1 1 1 1 1 1 1 1 ...
##  $ n.censor     : num [1:146] 2 3 3 1 0 6 0 3 4 0 ...
##  $ surv         : num [1:146] 0.997 0.994 0.992 0.989 0.986 ...
##  $ std.err      : num [1:146] 0.00274 0.0039 0.00479 0.00554 0.00619 ...
##  $ cumhaz       : num [1:146] 0.00275 0.00553 0.00835 0.01118 0.01402 ...
##  $ std.chaz     : num [1:146] 0.00275 0.00391 0.00482 0.00559 0.00627 ...
##  $ strata       : Factor w/ 2 levels "sex=Female","sex=Male": 1 1 1 1 1 1 1 1 1 1 ...
##  $ type         : chr "right"
##  $ logse        : logi TRUE
##  $ conf.int     : num 0.95
##  $ conf.type    : chr "log"
##  $ lower        : num [1:146] 0.992 0.987 0.982 0.978 0.974 ...
##  $ upper        : num [1:146] 1 1 1 1 0.998 ...
##  $ call         : language survfit(formula = Surv(age, cens) ~ sex, data = chan)
##  $ table        : num [1:2, 1:9] 365 97 365 97 365 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:2] "sex=Female" "sex=Male"
##   .. ..$ : chr [1:9] "records" "n.max" "n.start" "events" ...
##  $ rmean.endtime: num [1:2] 101 101
##  - attr(*, "class")= chr "summary.survfit"
```

```r
# Estimaciones de interés
res$table[, c("*rmean", "median")]
```

```
##              *rmean median
## sex=Female 88.46153     88
## sex=Male   86.89935     87
```

```r
as.numeric(res$table[, c("*rmean", "median")])
```

```
## [1] 88.46153 86.89935 88.00000 87.00000
```


\BeginKnitrBlock{exercise}\iffalse{-91-66-111-111-116-115-116-114-97-112-32-99-101-110-115-117-114-97-100-111-32-99-111-110-32-114-105-101-115-103-111-32-112-114-111-112-111-114-99-105-111-110-97-108-32-100-101-32-67-111-120-93-}\fi{}
<span class="exercise" id="exr:censboot-cox-ej"><strong>(\#exr:censboot-cox-ej)  \iffalse (Bootstrap censurado con riesgo proporcional de Cox) \fi{} </strong></span>
Reproducir el ejemplo en Canty (2002, [Rnews_2002-3](http://cran.fhcrc.org/doc/Rnews/Rnews_2002-3.pdf)) del modelo de riesgo proporcional de Cox (Cox, 1972):

\EndKnitrBlock{exercise}


```r
# Datos
data(melanoma)
mel <- melanoma[melanoma$ulcer == 1, ]
mel$cens <- 1 * (mel$status == 1)
# Estimación supervivencia
library(survival)
# Modelo de riesgo proporcional de Cox
mel.cox <- coxph(Surv(time, cens) ~ thickness, data = mel)
mel.cox
```

```
## Call:
## coxph(formula = Surv(time, cens) ~ thickness, data = mel)
## 
##              coef exp(coef) se(coef)    z      p
## thickness 0.09968   1.10481  0.04052 2.46 0.0139
## 
## Likelihood ratio test=5  on 1 df, p=0.02541
## n= 90, number of events= 41
```

```r
# summary(mel.cox)
# Estadísticos de interés
mel.cox$coefficients
```

```
##  thickness 
## 0.09967665
```
