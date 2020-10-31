# Bootstrap y estimación no paramétrica de la densidad {#npden}





En este capítulo se presentarán diversos métodos bootstrap adecuados
para realizar inferencia en algunos problemas en el contexto de la
estimación no paramétrica de la densidad. Concretamente, se abordará el
problema de construcción de intervalos de confianza para la funciones de
densidad en un punto dado, así como la selección del parámetro de
suavizado para el estimador tipo núcleo de la función de densidad. A
continuación se incluye una breve introducción a estos métodos no
paramétricos de estimación de curvas.

## Estimación no paramétrica de la función de densidad

Como ya se introdujo en la Sección \@ref(modunif-boot-suav), si 
$\left( X_1, X_2, \ldots, X_n \right)$ es una muestra aleatoria simple
(m.a.s.),, de una población con función de distribución $F$, absolutamente
continua, y función de densidad $f$, el estimador tipo núcleo propuesto por
Parzen (1962) y Rosenblatt (1956) viene dado por
$$\hat{f}_{h}\left( x \right) =\frac{1}{nh}\sum_{i=1}^{n}K\left( \frac{x-X_i}{
h} \right) =\frac{1}{n}\sum_{i=1}^{n}K_{h}\left( x-X_i \right),$$
donde $K_{h}\left( u \right) =\frac{1}{h}K\left( \frac{u}{h} \right)$, 
$K$ es una función núcleo (normalmente una densidad simétrica en torno al cero)
y $h>0$ es una parámetro de suavizado, llamado ventana, que regula el
tamaño del entorno que se usa para llevar a cabo la estimación. Es
habitual exigir que la función núcleo $K$ sea no negativa y su integral
sea uno:
$$K\left( u \right) \geq 0,~\forall u,~\int_{-\infty }^{\infty }
K\left( u \right) du=1.$$
Además también es frecuente exigir que $K$ sea una
función simétrica ($K\left( -u \right) =K\left( u \right)$).

<!-- 
Incluir más detalles sobre la selección
de la ventana con `density()`
* `"nrd0"`: regla del pulgar (Silverman, 1986).
* `"nrd"`: regla de Scott (1992).
* `"ucv"`: validación cruzada insesgada.
* `"bcv"`: validación cruzada sesgada.
* `"SJ"`: métodos de Sheather y Jones (1991).
También se pueden obtener estas ventanas llamando a la función `bw.xxx()` correspondiente, donde `xxx` son los posibles valores del parámetro `bw`.
-->
 
## Sesgo, varianza y error cuadrático medio

### Sesgo

Mediante cálculos sencillos puede obtenerse el sesgo del estimador de
Parzen-Rosenblatt:
$$\begin{aligned}
Sesgo\left( \hat{f}_{h}\left( x \right) \right) &= E\left( \hat{f}_{h}\left(
x \right) \right) -f\left( x \right) =\int \frac{1}{h}K\left( \frac{x-y}{h}
 \right) f\left( y \right) dy-f\left( x \right) \\
&= \left( K_{h}\ast f \right) \left( x \right) -f\left( x \right),
\end{aligned}$$
siendo $\ast$ el operador convolución:
$$\left( f\ast g \right) \left( x \right) 
= \int f\left( x-y \right) g\left( y \right) dy.$$
A partir de la expresión del sesgo puede obtenerse otra
asintótica para el mismo:
$$E\left( \hat{f}_{h}\left( x \right) \right) -f\left( x \right) =\frac{d_{K}}{2}
h^2f^{\prime \prime }\left( x \right) +O\left( h^{4} \right),$$con
$d_{K}=\int t^2K\left( t \right) dt$.

### Varianza

La varianza puede tratarse análogamente:
$$\begin{aligned}
Var\left( \hat{f}_{h}\left( x \right) \right) &= \frac{1}{nh^2}Var\left(
K\left( \frac{x-X_1}{h} \right) \right) \\
&= \frac{1}{nh^2}\left[ \int K\left( \frac{x-y}{h} \right)^2f\left(
y \right) dy-\left( \int K\left( \frac{x-y}{h} \right) f\left( y \right)
dy \right)^2\right] \\
&= \frac{1}{n}\left[ \left( \left( K_{h} \right)^2\ast f \right) \left(
x \right) -\left( \left( K_{h}\ast f \right) \left( x \right) \right)^2
\right] \\
&= \frac{1}{nh}\left[ \left( K^2 \right) _{h}\ast f\right] \left( x \right) -
\frac{1}{n}\left[ \left( K_{h}\ast f \right) \left( x \right) \right]^2.\end{aligned}$$Su
expresión asintótica resulta:
$$Var\left( \hat{f}_{h}\left( x \right) \right) =\frac{c_{K}}{nh}f\left(
x \right) - \frac{1}{n}f\left( x \right)^2 + O\left( \frac{h}{n} \right),$$
con $c_{K}=\int K\left( t \right)^2dt$.

### Error cuadrático medio

Como consecuencia el error cuadrático medio del estimador es:
$$\begin{aligned}
MSE\left( \hat{f}_{h}\left( x \right) \right) =&\ E\left( \hat{f}_{h}\left(
x \right) -f\left( x \right) \right)^2=Sesgo\left( \hat{f}_{h}\left(
x \right) \right)^2+Var\left( \hat{f}_{h}\left( x \right) \right) \\
=&\ \left[ \left( K_{h}\ast f \right) \left( x \right) -f\left( x \right) \right]
^2+\frac{1}{nh}\left[ \left( K^2 \right) _{h}\ast f\right] \left( x \right) \\
&-\frac{1}{n}\left[ \left( K_{h}\ast f \right) \left( x \right) \right]^2.
\end{aligned}$$
Además, su expresión asintótica
es:$$MSE\left( \hat{f}_{h}\left( x \right) \right) =\frac{d_{K}^2}{4}
h^{4}f^{\prime \prime }\left( x \right)^2+\frac{c_{K}}{nh}f\left( x \right)
-\frac{1}{n}f\left( x \right)^2+O\left( h^{6} \right) +O\left( \frac{h}{n}
 \right).$$

### Error cuadrático medio integrado (MISE)

Una medida global (no para un $x$ particular) del error cometido por el
estimador es el error cuadrático medio integrado:
$$\begin{aligned}
& & MISE\left( \hat{f}_{h}\left( x \right) \right) =\int E\left[ \left( \hat{f}
_{h}\left( x \right) -f\left( x \right) \right)^2\right] dx=\int MSE\left( 
\hat{f}_{h}\left( x \right) \right) dx= \\
&&\int \left[ \left( K_{h}\ast f \right) \left( x \right) -f\left( x \right) 
\right]^2dx+\frac{c_{K}}{nh}-\frac{1}{n}\int \left[ \left( K_{h}\ast
f \right) \left( x \right) \right]^2dx.
\end{aligned}$$
Una expresión asintótica para el mismo es la siguiente:
$$\begin{aligned}
MISE\left( \hat{f}_{h}\left( x \right) \right) =&\ \frac{d_{K}^2}{4}h^4\int
f^{\prime \prime }\left( x \right)^2dx+\frac{c_{K}}{nh}-\frac{1}{n}\int
f\left( x \right)^2dx \\
&+O\left( h^{6} \right) +O \left( \frac{h}{n} \right).
\end{aligned}$$
En ella se puede ver el efecto negativo de tomar ventanas ($h$) demasiado
grandes o demasiado pequeñas.

## Aproximación Bootstrap de la distribución del estimador de Parzen-Rosenblatt {#aproximacion-bootstrap}

Antes de proceder a abordar el bootstrap en este contexto conviene
presentar la distribución asintótica del estimador y otras
aproximaciones posibles. Pueden encontrarse más detalles sobre estos
resultados en Cao (1990).

### Distribución asintótica del estimador de Parzen-Rosenblatt

Las condiciones mínimas necesarias para que el sesgo y la varianza del
estimador tiendan a cero cuando el tamaño muestral tiende a infinito son
$h\rightarrow 0$, $nh\rightarrow \infty$. En tales circunstancias se
tiene
$$\sqrt{nh}\left( \hat{f}_{h}\left( x \right) -f\left( x \right) \right) \overset
{d}{\rightarrow }\mathcal{N}\left( B,V \right) \text{.}$$Además, puede probarse que
el valor asintóticamente óptimo de $h$, en el sentido del $MSE$, es
$h=c_{0}n^{-1/5}$, con
$$c_{0}=\left( \frac{c_{K}f\left( x \right)}{d_{K}^2f^{\prime \prime }\left(
x \right)^2} \right)^{1/5}.$$

Con esa elección de $h$ los valores de media y varianza de la
distribución normal límite son$$\begin{aligned}
B &= \frac{1}{2}c_{0}^{5/2}d_{K}f^{\prime \prime }\left( x \right), \\
V &= c_{K}f\left( x \right).\end{aligned}$$

Para utilizar la distribución asintótica anterior en la construcción de
intervalos de confianza para $f\left( x \right)$ podemos

1.  Estimar $B$ y $V$ y utilizarlos en la correspondiente distribución
    normal (**metodo plug-in**).

2.  Diseñar un plan de remuestreo y utilizar el **método bootstrap**.

### Aproximación plug-in

Pasa por estimar $B$ y $V$ mediante$$\begin{aligned}
\hat{B} &= \frac{1}{2}\hat{c}_{0}^{5/2}d_{K}\hat{f}_{g}^{\prime \prime
}\left( x \right), \\
\hat{V} &= c_{K}\hat{f}_{h}\left( x \right),\end{aligned}$$siendo $g$
una ventana adecuada para estimar la derivada segunda de la función de
densidad. Utilizando la desigualdad de Berry-Esséen se obtiene el
siguiente orden de
convergencia:$$\sup_{z\in \boldsymbol{R}}\left\vert P\left[ \sqrt{nh}\left( \hat{f}
_{h}\left( x \right) -f\left( x \right) \right) \leq z\right] -\Phi \left( 
\frac{z-\hat{B}}{\hat{V}^{1/2}} \right) \right\vert =O_{P}\left(
n^{-1/5} \right),$$

que empeora la tasa teórica de la aproximación normal basada en la media
y varianza exactas ($B_n=E\left[ \sqrt{nh}\left( \hat{f}_{h}\left(
x \right) -f\left( x \right) \right) \right]$ y
$V_n=Var\left[ \sqrt{nh} \left( \hat{f}_{h}\left( x \right) -f\left( x \right) \right) \right]$):
$$\sup_{z\in \boldsymbol{R}}\left\vert P\left[ \sqrt{nh}\left( \hat{f}
_{h}\left( x \right) -f\left( x \right) \right) \leq z\right] -\Phi \left( 
\frac{z-B_n}{V_n^{1/2}} \right) \right\vert =O\left( n^{-2/5} \right),$$
aunque no la de la normal asintótica, $\mathcal{N}\left( B,V \right)$, cuya tasa es
igualmente de orden $O_{P}\left( n^{-1/5} \right)$.

### Aproximación bootstrap

Se procede según el siguiente plan de remuestreo.

1.  A partir de la muestra $\left( X_1,X_2,\ldots ,X_n \right)$ y
    utilizando una **ventana piloto** $g$, se calcula el estimador de
    Parzen-Rosenblatt $\hat{f}_{g}$.

2.  Se arrojan remuestras bootstrap $\left( X_1^{\ast},X_2^{\ast
    },\ldots ,X_n^{\ast} \right)$ a partir de la densidad
    $\hat{f}_{g}$.

3.  Se construye el análogo bootstrap del estimador de Parzen-Rosenblatt
    $$\hat{f}_{h}^{\ast}\left( x \right) =\frac{1}{nh}\sum_{i=1}^{n}K\left( \frac{
    x-X_i^{\ast}}{h} \right).$$

4.  Se aproxima la distribución en el muestreo de $\sqrt{nh}\left( 
    \hat{f}_{h}\left( x \right) -f\left( x \right) \right)$ por la
    distribución en el remuestreo de
    $\sqrt{nh}\left( \hat{f}_{h}^{\ast}\left( x \right) -
    \hat{f}_{g}\left( x \right) \right)$.

Si nuestro interés estuviese en el sesgo o la varianza de $\hat{f}
_{h}\left( x \right)$ entonces utilizaríamos, en el paso 4 del algoritmo
anterior, los análogos bootstrap del sesgo o la varianza:
$E^{\ast}\left( \hat{f}_{h}^{\ast
}\left( x \right) -\hat{f}_{g}\left( x \right) \right)$ o
$Var^{\ast}\left( 
\hat{f}_{h}^{\ast}\left( x \right) \right)$.

En el algoritmo anterior, la ventana $g$ ha de ser asintóticamente mayor
que $h$. De hecho, una elección razonable para $g$ es aquella que
minimiza $E\left[ \left( \hat{f}_{g}^{\prime \prime }\left( x \right)
-f^{\prime \prime }\left( x \right) \right)^2\right]$.
Asintóticamente esa ventana viene dada por
$$g\simeq \left( \frac{5f\left( x \right) \int K^{\prime \prime }\left( t
\right)^2dt}{d_{K}^2f^{\left( 4 \right)}\left( x \right)^2n} \right)^{1/9}.$$

El orden de convergencia de para la aproximación bootstrap viene dado por
$$\begin{aligned}
&\sup_{z\in \boldsymbol{R}}\left\vert P\left[ \sqrt{nh}\left( \hat{f}
_{h}\left( x \right) -f\left( x \right) \right) \leq z\right] -P^{\ast}\left[ 
\sqrt{nh}\left( \hat{f}_{h}^{\ast}\left( x \right) -\hat{f}_{g}\left(
x \right) \right) \leq z\right] \right\vert \\
&= O_{P}\left( n^{-2/9} \right),\end{aligned}$$que mejora los ofrecidos
por la aproximación normal teórica y el método plug-in.

## El Bootstrap en la selección del parámetro de suavizado.

### Expresión asintótica de la ventana óptima

El $MISE$ tiene una expresión asintótica que puede usarse como criterio
para obtener un valor óptimo del parámetro de
suavizado:$$MISE\left( h \right) =AMISE\left( h \right) +O\left( h^{6} \right) +O\left( 
\frac{h}{n} \right),$$con$$AMISE\left( h \right) =\frac{d_{K}^2}{4}h^{4}\int f^{\prime \prime }\left(
x \right)^2dx+\frac{c_{K}}{nh}-\frac{1}{n}\int f\left( x \right)^2dx.$$El
parámetro de suavizado que minimiza el $AMISE$
es$$h_{AMISE}=\left( \frac{c_{K}}{nd_{K}^2\int f^{\prime \prime }\left(
x \right)^2dx} \right)^{1/5}.$$

Existen multitud de métodos encaminados a dar respuesta al problema de
selección del parámetro de suavizado. Entre ellos destacamos los métodos
plug-in, los de validación cruzada (suavizada o no) y, desde luego, los
métodos bootstrap (ver, por ejemplo, Marron (1992)).

### Análogo bootstrap del $MISE$

La idea básica (Cao (1993)) consiste en diseñar un plan de remuestreo,
del tipo bootstrap suavizado, para estimar el $MISE$:

1.  A partir de la muestra $\left( X_1,X_2,\ldots ,X_n \right)$ y
    utilizando una ventana piloto $g$, se calcula el estimador de
    Parzen-Rosenblatt $\hat{f}_{g}$.

2.  Se arrojan remuestras bootstrap $\left( X_1^{\ast},X_2^{\ast
    },\ldots ,X_n^{\ast} \right)$ de la densidad $\hat{f}_{g}$.

3.  Para cada $h>0$, se obtiene el análogo bootstrap del estimador de
    Parzen-Rosenblatt
    $$\hat{f}_{h}^{\ast}\left( x \right) =\frac{1}{nh}\sum_{i=1}^{n}K\left( \frac{
    x-X_i^{\ast}}{h} \right).$$

4.  Se construye la versión bootstrap del
    $MISE$:$$MISE^{\ast}\left( h \right) =\int E^{\ast}\left[ \left( \hat{f}_{h}^{\ast
    }\left( x \right) -\hat{f}_{g}\left( x \right) \right)^2\right] dx.$$

5.  Se minimiza $MISE^{\ast}\left( h \right)$ en $h>0$ y se obtiene el
    selector bootstrap:
    $$h_{MISE}^{\ast}=\arg \min_{h>0}MISE^{\ast}\left( h \right)$$

### Expresión cerrada para $MISE^{\ast}$

A diferencia de lo que es habitual, en este contexto es posible obtener
una expresión cerrada para el análogo bootstrap del $MISE$:
$$\begin{aligned}
MISE^{\ast}\left( h \right) =&\ \int \left[ \left( K_{h}\ast 
\hat{f}_{g} \right) \left( x \right) -\hat{f}_{g}\left( x \right) \right]^2dx \\
&+\frac{c_{K}}{nh}-\frac{1}{n}\int \left[ \left( K_{h}\ast 
\hat{f}_{g} \right) \left( x \right) \right]^2dx \\
=&\ \frac{c_{K}}{nh}-\frac{1}{n^{3}}\sum_{i,j=1}^{n}\left[ \left( K_{h}\ast
K_{g} \right) \ast \left( K_{h}\ast K_{g} \right) \right] \left(
X_i-X_j \right) \\
&+\frac{1}{n^2}\sum_{i,j=1}^{n}\left[ \left( K_{h}\ast K_{g}-K_{g} \right)
\ast \left( K_{h}\ast K_{g}-K_{g} \right) \right] \left( X_i-X_j \right).\end{aligned}$$

### Elección de la ventana piloto

De nuevo ocurre que el problema de elección óptima de la ventana piloto,
$g$, viene ligado al de estimación óptima de la curvatura de la función
de densidad. Así, una buena elección de $g$ es la que
minimiza$$E\left[ \left( \int \hat{f}_{g}^{\prime \prime }\left( x \right)^2dx-\int
f^{\prime \prime }\left( x \right)^2dx \right)^2\right] .$$El valor
asintótico de dicha ventana $g$
es$$g_{0}=\left( \frac{\int K^{\prime \prime }\left( t \right)^2dt}{nd_{K}\int
f^{\left( 3 \right)}\left( x \right)^2dx} \right)^{1/7}.$$

### Resultados teóricos

Utilizando cualquier ventana piloto determinística que cumpla 
$\frac{g-g_{0}}{g_{0}}=O\left( n^{-1/14} \right)$, se tiene
$$\begin{aligned}
\frac{h_{MISE}^{\ast}-h_{MISE}}{h_{MISE}} &= O_{P}\left( n^{-5/14} \right),\\
\frac{MISE\left( h_{MISE}^{\ast} \right) -MISE\left( h_{MISE} \right)}{
MISE\left( h_{MISE} \right)} &= O_{P}\left( n^{-5/7} \right).
\end{aligned}$$

Mediante técnicas más sofisticadas que permiten que $g$ dependa de
$h$ pueden obtenerse tasas ligeramente
mejores:$$\frac{h_{MISE}^{\ast}-h_{MISE}}{h_{MISE}}=O_{P}\left( n^{-1/2} \right).$$

### Caso particular de núcleo gaussiano

Cuando el núcleo $K$ es la función de densidad de una
$\mathcal{N}\left( 0,1 \right)$:

-   $K_{h}$ es la densidad de una $\mathcal{N}\left( 0,h^2 \right)$

-   $K_{g}$ es la densidad de una $\mathcal{N}\left( 0,g^2 \right)$

-   $K_{h}\ast K_{g}$ es la densidad de una
    $\mathcal{N}\left( 0,h^2+g^2 \right)$

-   $\left( K_{h}\ast K_{g} \right) \ast \left( K_{h}\ast K_{g} \right)$
    es la densidad de una $\mathcal{N}\left( 0,2h^2+2g^2 \right)$

-   $\left( K_{h}\ast K_{g} \right) \ast K_{g}$ es la densidad de una
    $\mathcal{N}\left( 0,h^2+2g^2 \right)$

-   $K_{g}\ast K_{g}$ es la densidad de una $\mathcal{N}\left( 0,2g^2 \right)$

con lo cual
$$\begin{aligned}
MISE^{\ast}\left( h \right) =&\ \frac{c_{K}}{nh}-\frac{1}{n^{3}}
\sum_{i,j=1}^{n}K_{\sqrt{2h^2+2g^2}}\left( X_i-X_j \right) \\
&+\frac{1}{n^2}\sum_{i,j=1}^{n}\left[ K_{\sqrt{2h^2+2g^2}}\left(
X_i-X_j \right) \right. \\
&\left. -2K_{\sqrt{h^2+2g^2}}\left( X_i-X_j \right) +K_{\sqrt{2g^2}
}\left( X_i-X_j \right) \right] .
\end{aligned}$$

### Comparación con otros selectores

El método bootstrap presentado es muy semejante al de validación cruzada
suavizada (SCV) propuesto por Hall, Marron y Park (1992). En estudios de
simulación comparativos (ver Cao, Cuevas y González-Manteiga (1993))
puede verse como este método ofrece resultados muy competitivos con
otros métodos de selección del parámetro de suavizado. En general es el
que mejor comportamiento ofrece junto con el método plug-in tipo
solve-the-equation de Sheather y Jones (1991) y el método SCV.

Otros selectores bootstrap con mucho peor comportamiento son:

-   Hall (1990), en el que se remuestrea de la distribución empírica,
    con lo cual no se imita el sesgo.

-   Faraway y Jhun (1990), que eligen $g$ como la ventana de validación
    cruzada, que resulta ser demasiado pequeña.

-   Taylor (1989), que elige $g=h$ , con lo cual $MISE^{\ast}\left(
    h \right) \rightarrow 0$, cuando $h\rightarrow \infty$, lo cual
    produce un mínimo global de $MISE^{\ast}$ inconsistente con
    $h_{MISE}$.


## Estimación no paramétrica de la densidad en R {#npden-r}

Como ya se comentó en la Sección \@ref(modunif-boot-suav),
en `R` podemos emplear la función `density()` del paquete base para obtener
una estimación tipo núcleo de la densidad. 
Los principales parámetros (con los valores por defecto) son los siguientes:

```
density(x, bw = "nrd0", adjust = 1, kernel = "gaussian", n = 512, from, to)
```

- `bw`: ventana, puede ser un valor numérico o una cadena de texto que la determine
  (en ese caso llamará internamente a la función `bw.xxx()` donde `xxx` se corresponde
  con la cadena de texto). Las opciones son:

    - `"nrd0"`, `"nrd"`: Reglas del pulgar de Silverman (1986, page 48, eqn (3.31)) y 
      Scott (1992), respectivamente. Como es de esperar que la densidad objetivo 
      no sea tan suave como la normal, estos criterios tenderán a seleccionar 
      ventanas que producen un sobresuavizado de las observaciones.

    - `"ucv"`, `"bcv"`: Métodos de validación cruzada insesgada y sesgada, respectivamente.
    
    - `"sj"`, `"sj-ste"`, `"sj-dpi"`: Métodos de Sheather y Jones (1991), 
        "solve-the-equation" y "direct plug-in", respectivamente.
 
-   `adjust`: parameto para reescalado de la ventana, las estimaciones se calculan 
    con la ventana `adjust*bw`.

-   `kernel`: cadena de texto que determina la función núcleo, las opciones son: `"gaussian"`,
    `"epanechnikov"`, `"rectangular"`, `"triangular"`, `"biweight"`, `"cosine"` y `"optcosine"`.
    
-   `n`, `from`, `to`: permiten establecer la rejilla en la que se obtendrán las estimaciones
    (si $n>512$ se emplea `fft()` por lo que se recomienda establecer `n` a un múltiplo de 2;
    por defecto `from` y `to` se establecen como `cut = 3` veces la ventana desde los extremos 
    de las observaciones).

Utilizaremos como punto de partida el código empleado en la Sección \@ref(modunif-boot-suav).
Considerando el conjunto de datos `precip` (que contiene el promedio de precipitación, 
en pulgadas de lluvia, de 70 ciudades de Estados Unidos).


```r
x <- precip
h <- bw.SJ(x)
npden <- density(x, bw = h)
# npden <- density(x, bw = "SJ")

# plot(npden)
hist(x, freq = FALSE, main = "Kernel density estimation",
     xlab = paste("Bandwidth =", formatC(h)), lty = 2,
     border = "darkgray", xlim = c(0, 80), ylim = c(0, 0.04))
lines(npden, lwd = 2)
rug(x, col = "darkgray")
```



\begin{center}\includegraphics[width=0.7\linewidth]{06-npden_files/figure-latex/unnamed-chunk-2-1} \end{center}

Alternativamente podríamos emplear implementaciones en otros paquetes de `R`.
Uno de los más empleados es `ks` (Duong, 2019), que admite estimación 
incondicional y condicional multidimensional.
También se podrían emplear los paquetes `KernSmooth` (Wand y Ripley, 2019), 
`sm` (Bowman y Azzalini, 2019), `np` (Tristen y Jeffrey, 2019), 
`kedd` (Guidoum, 2019), `features` (Duong y Matt, 2019) y `npsp` (Fernández-Casal, 2019), 
entre otros.

## Ejemplos

En esta sección nos centraremos en el bootstrap
en la estimación tipo núcleo de la densidad 
para la aproximación de la precisión y el sesgo,
y también para el cálculo de intervalos de confianza.


### Bootstrap y estimación del sesgo

La idea sería aproximar la distribución del error de estimación 
$\hat f(x) - f(x)$ por la distribución bootstrap de
$\hat f^{\ast}(x) - \hat f(x)$ (bootstrap percentil básico).

Como se comentó en la Sección \@ref(aproximacion-bootstrap) 
la ventana $g$ ha de ser asintóticamente mayor que $h$ (de orden $n^{-1/5}$) 
y la recomendación sería emplear la ventana óptima para la estimación de 
$f^{\prime \prime }\left( x \right)$, de orden $n^{-1/9}$. 


```r
# Remuestreo
set.seed(1)
n <- length(x)
g <- h * n^(4/45) # h*n^(-1/9)/n^(-1/5)
range.x <- range(npden$x) # Para fijar las posiciones de estimación
B <- 1000
stat_den_boot <- matrix(nrow = length(npden$x), ncol = B)
for (k in 1:B) {
    # x_boot <- sample(x, n, replace = TRUE) + rnorm(n, 0, g)
    x_boot <- rnorm(n, sample(x, n, replace = TRUE), g)
    den_boot <- density(x_boot, bw = h, from = range.x[1], to = range.x[2])$y
    # Si se quiere tener en cuenta la variabilidad debida a la selección de
    # la ventana habría que emplear el mismo criterio en la función `density`.
    stat_den_boot[, k] <- den_boot - npden$y
}

# Calculo del sesgo y error estándar 
bias <- apply(stat_den_boot, 1, mean)
std.err <- apply(stat_den_boot, 1, sd)

# Representar estimación y corrección de sesgo bootstrap
plot(npden, type="l", ylim = c(0, 0.05), lwd = 2)
lines(npden$x, npden$y - bias)
```



\begin{center}\includegraphics[width=0.7\linewidth]{06-npden_files/figure-latex/unnamed-chunk-3-1} \end{center}

```r
# lines(npden$x, pmax(0, npden$y - bias))
```


### Estimación por intervalos de confianza {#npden-r-ic}

Empleando la aproximación descrita en la Sección \@ref(icboot-basic)
podemos cálcular de estimaciones por intervalo de confianza (puntuales)
por el método percentil (básico).


```r
alfa <- 0.05
pto_crit <- apply(stat_den_boot, 1, quantile, probs = c(alfa/2, 1 - alfa/2))
# ic_inf_boot <- npden$y - pto_crit[2, ]
ic_inf_boot <- pmax(0, npden$y - pto_crit[2, ])
ic_sup_boot <- npden$y - pto_crit[1, ]

plot(npden, type="l", ylim = c(0, 0.05), lwd = 2)
lines(npden$x, pmax(0, npden$y - bias))
lines(npden$x, ic_inf_boot, lty = 2)
lines(npden$x, ic_sup_boot, lty = 2)
```



\begin{center}\includegraphics[width=0.7\linewidth]{06-npden_files/figure-latex/unnamed-chunk-4-1} \end{center}


### Implementación con el paquete `boot`

Como también se comentó en la Sección \@ref(modunif-boot-suav), 
la recomendación es implementar el bootstrap suavizado como un bootstrap paramétrico:


```r
library(boot)

# Los objetos necesarios para el cálculo del estadístico
# hay que pasarlos a traves del argumento `data` de `boot`.
range.x <- range(npden$x)
data.precip <- list(x = x, h = h, range.x = range.x)

ran.gen.smooth <- function(data, mle) {
    # Función para generar muestras aleatorias mediante
    # bootstrap suavizado con función núcleo gaussiana,
    # mle contendrá la ventana
    n <- length(data$x)
    g <- mle
    xboot <- rnorm(n, sample(data$x, n, replace = TRUE), g)
    out <- list(x = xboot, h = data$h, range.x = data$range.x)
}

statistic <- function(data) 
                density(data$x, bw = data$h, from = range.x[1], to = range.x[2])$y

set.seed(1)
res.boot <- boot(data.precip, statistic, R = B, sim = "parametric",
                 ran.gen = ran.gen.smooth, mle = g)

# Calculo del sesgo y error estándar
bias <- with(res.boot, apply(t, 2, mean, na.rm = TRUE) -  t0)
std.err <- apply(res.boot$t, 2, sd, na.rm = TRUE)
```

Además, la función `boot.ci()` solo permite el cálculo del intervalo de 
confianza para cada valor de $x$ de forma independiente (parámetro `index`). 
Por lo que podría ser recomendable obtenerlo a partir de las réplicas 
bootstrap del estimador:


```r
# Método percentil básico calculado directamente 
# a partir de las réplicas bootstrap del estimador
alfa <- 0.05
pto_crit <- apply(res.boot$t, 2, quantile, probs = c(alfa/2, 1 - alfa/2))
ic_inf_boot <- pmax(0, 2*npden$y - pto_crit[2, ])
ic_sup_boot <- 2*npden$y - pto_crit[1, ]

plot(npden, ylim = c(0, 0.05), lwd = 2)
lines(npden$x, pmax(0, npden$y - bias))

lines(npden$x, ic_inf_boot, lty = 2)
lines(npden$x, ic_sup_boot, lty = 2)
```



\begin{center}\includegraphics[width=0.7\linewidth]{06-npden_files/figure-latex/unnamed-chunk-6-1} \end{center}

En la práctica, en muchas ocasiones se trabaja directamente con
las réplicas bootstrap del estimador. Por ejemplo, es habitual
generar envolventes como medida de la precisión de la estimación
(que se interpretan de forma similar a una banda de confianza):


```r
matplot(npden$x, t(res.boot$t), type = "l", col = "darkgray")
lines(npden, lwd = 2)
lines(npden$x, npden$y - bias)
```



\begin{center}\includegraphics[width=0.7\linewidth]{06-npden_files/figure-latex/unnamed-chunk-7-1} \end{center}

Pero la recomendación es emplear bootstrap básico (o percentil-*t*) en lugar
de bootstrap percentil (directo) en la presencia de sesgo:


```r
matplot(npden$x, 2*npden$y - t(res.boot$t), type = "l", col = "darkgray")
lines(npden, lwd = 2)
lines(npden$x, npden$y - bias)
```



\begin{center}\includegraphics[width=0.7\linewidth]{06-npden_files/figure-latex/unnamed-chunk-8-1} \end{center}

