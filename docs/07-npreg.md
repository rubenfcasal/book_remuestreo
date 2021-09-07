# Bootstrap y regresión no paramétrica {#npreg}




Por simplicidad nos centraremos en el caso bivariante, el caso multivariante sería análogo.
Supongamos que $\left\{ \left( X_1,Y_1 \right),\left( X_2,Y_2 \right), \ldots, \left( X_n,Y_n \right) \right\}$ es una m.a.s. de una población bidimensional $\left( X,Y \right)$, con $E\left( \left\vert Y\right\vert \right) <\infty$, y que el objetivo es realizar inferencias sobre la distribución condicional $\left. Y \right\vert_{X=x}$, principalmente estimar la función de regresión de $Y$ dada $X$: 
$$m\left( x \right) =E\left( \left. Y\right\vert_{X=x} \right).$$

En esta sección se introducirá la estimación no paramétrica de la función de regresión. 
En primer lugar se considerará el *estimador de Nadaraya-Watson* (en la Sección \@ref(locpol-r) se mostrarán generalizaciones de este estimador desde un punto de vista práctico en R). 
En la siguiente sección se introducirán distintos métodos de remuestreo, diseñados inicialmente para este estimador, y resultados para ellos.


## Estimador de Nadaraya-Watson {#nadaraya-watson}

La función de regresión $m\left( x \right) =E\left( \left. Y\right\vert_{X=x} \right)$ puede escribirse así:
$$\begin{aligned}
m\left( x \right) &= \int yf_{2|1}\left( \left. y\right\vert _{x} \right)
dy=\int y\frac{f\left( x,y \right)}{f_1\left( x \right)}dy=\frac{\int
yf\left( x,y \right) dy}{f_1\left( x \right)} \\
&= \frac{\int yf_{1|2}\left( \left. x\right\vert _{y} \right) f_2\left(
y \right) dy}{f_1\left( x \right)}=\frac{\Psi \left( x \right)}{f_1\left(
x \right)},
\end{aligned}$$
siendo $f_1\left( x \right)$ la función de densidad marginal de $X$ y
$$\Psi \left( x \right) =\int yf_{1|2}\left( \left. x\right\vert _{y} \right)
f_2\left( y \right) dy=E\left( Yf_{1|2}\left( \left. x\right\vert_{Y} \right)
\right).$$

Las funciones $\Psi \left( x \right)$ y $f_1\left( x \right)$ pueden
estimarse mediante el método núcleo:
$$\begin{aligned}
\hat{f}_{1,h}\left( x \right) &= \frac{1}{nh}\sum_{i=1}^{n}K\left( \frac{
x-X_i}{h} \right), \\
\hat{\Psi}_{h}\left( x \right) &= \frac{1}{nh}\sum_{i=1}^{n}K\left( \frac{
x-X_i}{h} \right) Y_i,
\end{aligned}$$
resultando así el estimador tipo núcleo de Nadaraya-Watson 
(ver Nadaraya (1964) y Watson (1964)):
$$\hat{m}_{h}\left( x \right) =\frac{\hat{\Psi}_{h}\left( x \right)}{\hat{f}
_{1,h}\left( x \right)}=\frac{\frac{1}{n}\sum_{i=1}^{n}K_{h}\left(
x-X_i \right) Y_i}{\frac{1}{n}\sum_{i=1}^{n}K_{h}\left( x-X_i \right)},$$
donde $K_{h}\left( x-X_i \right) =\frac{1}{h}K\left( \frac{x-X_i}{h} \right)$.

Para este estimador se pueden probar propiedades semejantes a las
mencionadas para el estimador de Parzen-Rosenblatt de la función de
densidad.

En esta sección se presentarán métodos de remuestreo bootstrap adecuados
para el contexto de la función de regresión. El objetivo es aproximar la
distribución en el muestreo del estimador de Nadaraya-Watson. Los
resultados reflejan el comportamiento de los métodos de remuestreo
bootstrap, tanto en un aspecto condicional a la muestra de la variable
explicativa como incondicionalmente.

### Distribución asintótica del estimador de Nadaraya-Watson

Antes de proceder a abordar el bootstrap en este contexto conviene
presentar la distribución asintótica del estimador de Nadaraya-Watson,
dado por
$$\hat{m}_{h}\left( x \right) =\frac{\frac{1}{n}\sum_{i=1}^{n}K_{h}\left(
x-X_i \right) Y_i}{\frac{1}{n}\sum_{i=1}^{n}K_{h}\left( x-X_i \right)}.$$

De forma semejante al caso de la densidad, puede comprobarse que las
condiciones mínimas necesarias para la consistencia del estimador, en
términos del parámetro de suavizado, son $h\rightarrow 0$,
$nh\rightarrow \infty$, cuando $n\rightarrow \infty$. En tales
circunstancias se tiene
$$\sqrt{nh}\left( \hat{m}_{h}\left( x \right) -m\left( x \right) \right) \overset
{d}{\rightarrow }\mathcal{N}\left( B,V \right) \text{.}$$

Además, puede probarse que el valor asintóticamente óptimo de $h$, en el
sentido del $MSE$, es de la forma $h=c_{0}n^{-1/5}$. En tal caso, los
valores de media y varianza de la distribución normal límite son
$$\begin{aligned}
B &= \frac{1}{2}c_{0}^{5/2}d_{K}\frac{m^{\prime \prime }\left( x \right)
f\left( x \right) +2m^{\prime}\left( x \right) f^{\prime}\left( x \right)}{
f\left( x \right)}, \\
V &= c_{K}\frac{\sigma^2\left( x \right)}{f\left( x \right)},
\end{aligned}$$
siendo $f\left( x \right)$ la función de densidad marginal de $X$ y
$\sigma^2\left( x \right) =Var\left( \left. Y\right\vert _{X=x} \right)$
la varianza condicional de $Y$ dado $X=x$.

Al igual que en el caso de la densidad, para utilizar la distribución
asintótica anterior en la construcción de intervalos de confianza para
$m\left( x \right)$ podemos

1.  Estimar $B$ y $V$ y utilizarlos en la correspondiente distribución
    normal (**metodo plug-in**).

2.  Diseñar un plan de remuestreo y utilizar el **método bootstrap**.

### Órdenes de convergencia de la distribución del estimador de Nadaraya-Watson a su distribución asintótica

Los órdenes de convergencia de la aproximación de la distribución
(condicional o incondicional) del estadístico a la distribución normal
límite vienen dados por:
$$\begin{aligned}
\sup_{z\in \boldsymbol{R}}\left\vert P^{\left. Y\right\vert _{X}}\left[ 
\sqrt{nh}\left( \hat{m}_{h}\left( x \right) -m\left( x \right) \right) \leq z
\right] -\Phi \left( \frac{z-B}{V^{1/2}} \right) \right\vert &= O_{P}\left(
n^{-1/5} \right), \\
\sup_{z\in \boldsymbol{R}}\left\vert P\left[ \sqrt{nh}\left( \hat{m}
_{h}\left( x \right) -m\left( x \right) \right) \leq z\right] -\Phi \left( 
\frac{z-B}{V^{1/2}} \right) \right\vert &= O\left( n^{-2/5} \right),
\end{aligned}$$
donde $P^{\left. Y\right\vert_{X}}\left( A \right)$ denota 
$P\left( \left. A \right\vert_{X_1,X_2,\ldots ,X_n} \right)$.

### Aproximación plug-in

Consiste en estimar $B$ y $V$ mediante estimadores apropiados de
$f\left(x \right)$, $f^{\prime}\left( x \right)$, $m\left( x \right)$,
$m^{\prime}\left( x \right)$, $m^{\prime \prime }\left( x \right)$ y 
$\sigma^2\left( x \right)$. Usando, para cada una de estas seis curvas,
selectores de los parámetros de suavizado encaminados a aproximar las
ventanas óptimas para cada una de ellas (proceso bastante laborioso),
pueden obtenerse estimadores del sesgo, $\hat{B}$, y la varianza,
$\hat{V}$, que cumplen $\hat{B}-B=O_{P}\left( n^{-2/9} \right)$ y
$\hat{V}-V=O_{P}\left( n^{-2/5} \right)$. 
Como consecuencia se tienen los siguientes órdenes de convergencia 
(condicional e incondicional) para la aproximación plug-in:
$$\begin{aligned}
\sup_{z\in \boldsymbol{R}}\left\vert P^{\left. Y\right\vert _{X}}\left[ 
\sqrt{nh}\left( \hat{m}_{h}\left( x \right) -m\left( x \right) \right) \leq z
\right] -\Phi \left( \frac{z-\hat{B}}{\hat{V}^{1/2}} \right) \right\vert
&= O_{P}\left( n^{-1/5} \right), \\
\sup_{z\in \boldsymbol{R}}\left\vert P\left[ \sqrt{nh}\left( \hat{m}
_{h}\left( x \right) -m\left( x \right) \right) \leq z\right] -\Phi \left( 
\frac{z-\hat{B}}{\hat{V}^{1/2}} \right) \right\vert &= O_{P}\left(
n^{-2/9} \right).
\end{aligned}$$
que iguala y empeora, respectivamente,
la tasa teórica de la aproximación normal límite (ver Cao (1991)).

## Métodos de remuestreo en regresión no paramétrica

En este caso se podría emplear también bootstrap uniforme (Sección \@ref(boot-unif-reg)) y el bootstrap residual, de forma totálmente análoga a como se mostró en la Sección \@ref(boot-residual).
Sin embargo, en el caso heterocedástico es habitual emplear *Wild bootstrap* y en el caso de diseño aleatorio podría ser recomendable emplear *bootstrap suavizado en la variable explicativa*.

### Wild bootstrap {#wild-bootstrap}

Este método de remuestreo bootstrap, propuesto por Wu (1986) y estudiado
por Härdle y Marron (1991), procede del siguiente modo:

1.  A partir del estimador de Nadaraya-Watson de $m\left( x \right)$ y
    tomando el parámetro ventana de partida, $h$, se construyen los
    residuos
    $r_i = Y_i - \hat{m}_{h}\left( X_i \right)$,
    $i=1, 2, \ldots, n$.

2.  Para cada índice $i=1,2,\ldots ,n$, se arroja, condicionalmente a la
    muestra observada, $\left\{ \left( X_1,Y_1 \right), \ \left(
    X_2,Y_2 \right),\right.$ 
    $\left.\ldots ,\ \left( X_n,Y_n \right) \right\}$, 
    un error bootstrap $\hat{\varepsilon}_i^{\ast}$ de una
    distribución de probabilidad que cumpla,
    $E^{\ast}\left( \hat{\varepsilon}_i^{\ast} \right) =0$,
    $E^{\ast}\left( \hat{\varepsilon}_i^{\ast 2} \right) =r_i^2$ y 
    $E^{\ast}\left( \hat{\varepsilon}_i^{\ast 3} \right) =r_i^{3}$. 
    Aunque la condición del momento de orden 3 no es estrictamente necesaria, 
    es útil para las demostraciones de validez del método.

3.  Usando una ventana piloto $g$, asintóticamente mayor que $h$ (i.e.
    $g/h\rightarrow \infty$), se arrojan análogos bootstrap de las
    observaciones de la variable respuesta:
    $Y_i^{\ast}=\hat{m}_{g}\left(X_i \right)
     +\hat{\varepsilon}_i^{\ast}$, $i=1,2,\ldots ,n$.

4.  A partir de la remuestra bootstrap $\left\{ \left( X_1,Y_1^{\ast
    } \right),\left( X_2,Y_2^{\ast} \right),\ldots ,\left(
    X_n,Y_n^{\ast} \right) \right\}$ se construye el análogo
    bootstrap del estimador de Nadaraya-Watson:
    $$\hat{m}_{h}^{\ast}\left( x \right) =\frac{\frac{1}{n}\sum_{i=1}^{n}K_{h}
    \left( x-X_i \right) Y_i^{\ast}}{\frac{1}{n}\sum_{i=1}^{n}K_{h}\left(
    x-X_i \right)}.$$

5.  Se aproxima la distribución en el muestreo de $\sqrt{nh}\left( 
    \hat{m}_{h}\left( x \right) -m\left( x \right) \right)$ por la
    distribución en el remuestreo de
    $\sqrt{nh}\left( \hat{m}_{h}^{\ast}\left( x \right) -
    \hat{m}_{g}\left( x \right) \right)$.

El paso 2 suele llevarse a cabo encontrando una variable aleatoria,
$V^{\ast}$, que cumpla $E^{\ast}\left( V^{\ast} \right) =0$, $E^{\ast}\left(
V^{\ast 2} \right) =1$ y $E^{\ast}\left( V^{\ast 3} \right) =1$,
arrojando una muestra de tamaño $n$ de la misma,
$\left( V_1^{\ast},V_n^{\ast},\ldots ,V_n^{\ast} \right)$, y luego 
definiendo $\hat{\varepsilon}_i^{\ast} = r_iV_i^{\ast}$ 
para $i=1, 2, \ldots, n$.

Una de las elecciones más habituales para la distribución de
$V^{\ast}$ es la distribución discreta con masa de probabilidad en dos
puntos ($P^{\ast}\left( V^{\ast}=a \right) =p$ y
$P^{\ast}\left( V^{\ast}=b \right) =1-p$) que es solución del sistema 
de tres ecuaciones dadas por los tres primeros momentos:
$$\begin{aligned}
ap+b\left( 1-p \right) &= 0, \\
a^2p+b^2\left( 1-p \right) &= 1, \\
a^{3}p+b^{3}\left( 1-p \right) &= 1.
\end{aligned}$$
Esto da lugar al llamado bootstrap de la sección aurea 
(golden section bootstrap), con
$a=\frac{1-\sqrt{5}}{2}$, $b=\frac{1+\sqrt{5}}{2}$, $p=\frac{
5+\sqrt{5}}{10}$, es decir
$$\begin{aligned}
P^{\ast}\left( V^{\ast}=\frac{1-\sqrt{5}}{2} \right) &= \frac{5+\sqrt{5}}{10} \\
P^{\ast}\left( V^{\ast}=\frac{1+\sqrt{5}}{2} \right) &= \frac{5-\sqrt{5}}{10}
\end{aligned}$$

La elección de la ventana $g$, que aparece en el paso 3, guarda relación
con la estimación de $m^{\prime \prime }\left( x \right)$, pues esa es
la cantidad crítica a la hora de estimar $B$ y $V$. Tomando una ventana
piloto de orden óptimo en ese sentido, $g_{0}\simeq d_{0}n^{-1/9}$, 
se obtienen las siguientes tasas de convergencia (condicionales e
incondicionales) para la aproximación dada por el wild bootstrap:
$$\begin{gathered}
\sup_{z\in \boldsymbol{R}} \left\vert P^{\left. Y\right\vert _{X}}\left[ 
\sqrt{nh}\left( \hat{m}_{h}\left( x \right) -m\left( x \right) \right) \leq z
\right] - P^{\ast}\left[ \sqrt{nh}\left( \hat{m}_{h}^{\ast}\left( x \right) -
\hat{m}_{g}\left( x \right) \right) \leq z\right] \right\vert = O_{P}\left( n^{-2/9} \right), 
\\
\sup_{z\in \boldsymbol{R}} \left\vert P\left[ \sqrt{nh}\left( \hat{m}
_{h}\left( x \right) -m\left( x \right) \right) \leq z\right] 
 - P^{\ast}\left[ \sqrt{nh}\left( \hat{m}_{h}^{\ast}\left( x \right) -\hat{m}_{g}\left(
x \right) \right) \leq z\right] \right\vert = O_{P}\left( n^{-1/5} \right).
\end{gathered}$$

### Bootstrap suavizado en la variable explicativa

La idea es tratar, por un lado, de considerar la variabilidad inherente
a la variable explicativa (en el wild bootstrap esa parte de la
remuestra se mantiene fija) y, por otro, que la distribución en el
remuestreo de $\left. Y^{\ast}\right\vert _{X^{\ast}=X_i}$ no sea 
degenerada (como sí lo sería en un bootstrap naïve bidimensional).

El plan de remuestreo, propuesto por Cao y González-Manteiga (1993)
consta de los siguientes pasos:

1.  Dada la muestra $\left\{ \left( X_1,Y_1 \right),\left(
    X_2,Y_2 \right),\ldots ,\left( X_n,Y_n \right) \right\}$, se
    construye un estimador (empírico en la variable respuesta y
    suavizado en la explicativa) de la distribución conjunta de
    $\left( X,Y \right)$:
    $$\hat{F}_{g}\left( x,y \right) =\frac{1}{n}\sum_{i=1}^{n}\mathbf{1}_{\left\{
    Y_i\leq y\right\} }\int_{-\infty }^{x}K_{g}\left( t-X_i \right) dt.$$

2.  Se arrojan remuestras bootstrap, $\left\{ \left( X_1^{\ast
    },Y_1^{\ast} \right),\left( X_2^{\ast},Y_2^{\ast} \right),\ldots
    ,\left( X_n^{\ast},Y_n^{\ast} \right) \right\}$, con
    distribución
    $\hat{F}_{g}\left( x,y \right)$.

3.  Se construye el análogo bootstrap del estimador de Nadaraya-Watson:
    $$\hat{m}_{h}^{\ast}\left( x \right) =\frac{\frac{1}{n}\sum_{i=1}^{n}K_{h}
    \left( x-X_i^{\ast} \right) Y_i^{\ast}}{\frac{1}{n}\sum_{i=1}^{n}K_{h}
    \left( x-X_i^{\ast} \right)}.$$

4.  Se utiliza la distribución en el remuestreo de $\sqrt{nh}\left( 
    \hat{m}_{h}^{\ast}\left( x \right) -\hat{m}_{g}\left( x \right) \right)$
    para aproximar la distribución del estadístico de interés:
    $\sqrt{nh}\left( \hat{m}_{h}\left( x \right) -m\left( x \right) \right)$.

La ventana piloto, $g$, óptima vuelve a ser de orden $n^{-1/9}$, es
decir asintóticamente mayor que $h$.

La distribución bidimensional de la que se remuestrea en el paso 2,
$\hat{F}_{g}\left( x,y \right)$, puede sustituirse por una distribución
suavizada en ambas variables:
$$\tilde{F}_{g}\left( x,y \right) =\frac{1}{n}\sum_{i=1}^{n}
\int_{-\infty}^{y}K_{g}\left( s-Y_i \right) ds \int_{-\infty }^{x}
K_{g}\left( t-X_i \right) dt.$$
Esto es lo mismo que remuestrear de la densidad bidimensional
$$\hat{f}_{g}\left( x,y \right) = \frac{1}{n}\sum_{i=1}^{n}
K_{g}\left( x-X_i \right) K_{g}\left( y-Y_i \right),$$
que es el estimador tipo
núcleo de Parzen-Rosenblatt de la variable bidimensional
$\left( X,Y \right)$.

Cálculos sencillos permiten demostrar que si $\left( X^{\ast},Y^{\ast} 
\right)$ tiene distribución $\hat{F}_{g}\left( x,y \right)$, entonces,

-   $X^{\ast}$ tiene densidad marginal bootstrap $\hat{f}_{g}\left(
    x \right)$.

-   La distribución marginal bootstrap de $Y^{\ast}$ es la empírica de
    las $Y_i$: $\hat{F}_n^{Y}\left( y \right) =\frac{1}{n}
    \sum_{i=1}^{n}\mathbf{1}_{\left\{ Y_i\leq y\right\} }$.

-   La función de regresión del plan de remuestreo bootstrap coincide
    con la estimación de Nadaraya-Watson con ventana $g$, es decir,
    $$E^{\ast}\left( \left. Y^{\ast}\right\vert _{X^{\ast}=x} \right)
    =\hat{m}_{g}\left( x \right).$$

-   De hecho, la distribución condicional $\left. Y^{\ast}\right\vert
    _{X^{\ast}=x}$ es
    $$\hat{F}_{g}\left( \left. y\right\vert _{x} \right) =\frac{\frac{1}{n}
    \sum_{i=1}^{n}K_{g}\left( x-X_i \right) \mathbf{1}_{\left\{ Y_i\leq
    y\right\} }}{\frac{1}{n}\sum_{i=1}^{n}K_{g}\left( x-X_i \right)},$$
    es decir, el estimador tipo núcleo Nadaraya-Watson de la distribución
    condicional.

Esta última observación da pie a diseñar un método que permita simular
valores de $\left( X^{\ast},Y^{\ast} \right)$, tal y como se requiere
en el paso 2 del plan de remuestreo. Para ello basta con simular
$X^{\ast}$ a partir del estimador de Parzen-Rosenblatt construído con
la muestra de la variable explicativa (es decir, el bootstrap suavizado
clásico) y luego simular $Y^{\ast}$ a partir de la distribución
discreta que da a cada dato $Y_i$ la probabilidad
$$w_i\left( X^{\ast} \right) =\frac{\frac{1}{n}K_{g}\left( X^{\ast}
-X_i \right)}{\frac{1}{n}\sum_{j=1}^{n}K_{g}\left( X^{\ast}-X_j \right)}
\text{, }i=1,2,\ldots ,n\text{.}$$

Las tasas de convergencia de la aproximación bootstrap proporcionadas
por este método resultan:
$$\begin{gathered}
\sup_{z\in \boldsymbol{R}}\left\vert P^{\left. Y\right\vert _{X}}\left[ 
\sqrt{nh}\left( \hat{m}_{h}\left( x \right) -m\left( x \right) \right) \leq z
\right]  -P^{\left. Y^{\ast}\right\vert _{X^{\ast}}}\left[ \sqrt{nh}\left( 
\hat{m}_{h}^{\ast}\left( x \right) -\hat{m}_{g}\left( x \right) \right) \leq z
\right] \right\vert \\
=O_{P^{\ast}}\left( n^{-2/9} \right) \text{, en probabilidad }P, \\
\sup_{z\in \boldsymbol{R}}\left\vert P\left[ \sqrt{nh}\left( \hat{m}
_{h}\left( x \right) -m\left( x \right) \right) \leq z\right] -P^{\ast}\left[ 
\sqrt{nh}\left( \hat{m}_{h}^{\ast}\left( x \right) -\hat{m}_{g}\left(
x \right) \right) \leq z\right] \right\vert \\
=O_{P}\left( n^{-2/9} \right).
\end{gathered}$$

### Resumen comparativo

La siguiente tabla recoge un resumen de las tasas de convergencia
obtenidas con cada una de las aproximaciones estudiadas:

+--------------------------+--------------------------------------+------------------------------+
| Aproximación             | condicional                          | incondicional                |
+==========================+======================================+==============================+
| Normal teórica           | $O_{P}\left( n^{-1/5}\right)$        | $O\left(n^{-2/5}\right)$     |
+--------------------------+--------------------------------------+------------------------------+
| Plug-in                  | $O_{P}\left( n^{-1/5}\right)$        | $O_{P}\left( n^{-2/9}\right)$|
+--------------------------+--------------------------------------+------------------------------+
| Wild bootstrap           | $O_{P}\left( n^{-2/9}\right)$        | $O_{P}\left(n^{-1/5}\right)$ |
+--------------------------+--------------------------------------+------------------------------+
| Bootstrap suavizado      | $O_{P^{\ast}}\left( n^{-2/9}\right)$ | $O_{P}\left(n^{-2/9}\right)$ |
| en la variable           | en probabilidad $P$                  |                              |
| explicativa              |                                      |                              |
+--------------------------+--------------------------------------+------------------------------+

Exceptuando las tasas de convergencia de la normal teórica (aproximación
inutilizable en la práctica) se observa que el método bootstrap
suavizado en la variable explicativa presenta órdenes que igualan o
mejoran al resto de los métodos, tanto en un aspecto condicional como
condicionalmente. Así, condicionalmente los dos remuestreos bootstrap
son los que ofrecen una mejor tasa de convergencia ($n^{-2/9}$, 
frente a $n^{-1/5}$ de la aproximación plug-in). En el
sentido incondicional el bootstrap suavizado en la variable explicativa
y la aproximación plug-in son los que presentan un mejor orden
($n^{-2/9}$, frente a $n^{-1/5}$ del wild bootstrap).


## Regresión polinómica local en R {#locpol-r}

El estimador de Nadaraya-Watson de $m(x)$ descrito en la Sección \@ref(nadaraya-watson)
es un caso particular de una clase más amplia de estimadores no paramétricos, 
denominados estimadores polinómicos locales. 


En el caso univariante, para cada $x_0$ se ajusta un polinomio:
$$\beta_0+\beta_{1}\left(x - x_0\right) + \cdots 
+ \beta_{p}\left( x-x_0\right)^{p}$$ 
por mínimos cuadrados ponderados, con pesos
$w_{i} = \frac{1}{h}K\left(\frac{x-x_0}{h}\right)$. 

-   La estimación en $x_0$ es $\hat{m}_{h}(x_0)=\hat{\beta}_0$.

-   Adicionalmente^[Se puede pensar que se están estimando los coeficientes de 
    un desarrollo de Taylor de $m(x_0)$.]: 
    $\widehat{m_{h}^{(r)}}(x_0) = r!\hat{\beta}_{r}$.

Por tanto, la estimación polinómica local de grado $p$, $\hat{m}_{h}(x)=\hat{\beta}_0$, se obtiene al minimizar:
$$\begin{aligned}
    \min_{\beta_0 ,\beta_1, \ldots, \beta_p}
    \sum_{i=1}^{n}\left\{ Y_{i} - \beta_0 
    -\beta_1(x - X_i) - \ldots \right. \nonumber \\
    \left. -\beta_p(x - X_i)^p \right\}^{2}
    K_{h}(x - X_i).
\end{aligned}$$ 

Explícitamente:
$$\hat{m}_{h}(x) = \mathbf{e}_{1}^{T} \left(
X_{x}^{T} {W}_{x} 
X_{x} \right)^{-1} X_{x}^{T} 
{W}_{x}\mathbf{Y} \equiv {s}_{x}^{T}\mathbf{Y},$$
donde $\mathbf{e}_{1} = \left( 1, \cdots, 0\right)^{T}$, $X_{x}$ 
es la matriz con $(1,x - X_i, \ldots, (x - X_i)^p)$ en la fila $i$,
$W_{x} = \mathtt{diag} \left( K_{h}(x_{1} - x), \ldots, K_{h}(x_{n} - x) \right)$
es la matriz de pesos, e $\mathbf{Y} = \left( Y_1, \cdots, Y_n\right)^{T}$ es el vector de observaciones de la respuesta.

Se puede pensar que se obtiene aplicando un suavizado polinómico a 
$(X_i, Y_i)$:
$$\hat{\boldsymbol{m}} = S\mathbf{Y},$$ 
siendo $S$ la matriz de suavizado con $\mathbf{s}_{X_{i}}^{T}$ en la fila $i$.

Habitualmente se considera:

-   $p=0$: Estimador Nadaraya-Watson.

-   $p=1$: Estimador lineal local.

Asintóticamente el estimador lineal local tiene un sesgo menor que el de 
Nadaraya-Watson (pero del mismo orden) y la misma varianza (e.g. Fan and Gijbels, 1996).
Sin embargo, su principal ventaja es que se ve menos afectado por el denominado
efecto frontera (*edge effect*).

Aunque el paquete base de `R` incluye herramientas para la estimación
tipo núcleo de la regresión (`lowess()`, `ksmooth()`), recomiendan
el uso del paquete `KernSmooth` (Wand y Ripley, 2019). 
Otros paquetes incluyen más funcionalidades: `sm` (Bowman y Azzalini, 2019), 
`np` (Tristen y Jeffrey, 2019), `npsp` (Fernández-Casal, 2019), entre otros.

Como ejemplo emplearemos el conjunto de datos `MASS::mcycle` que contiene mediciones 
de la aceleración de la cabeza en una simulación de un accidente de motocicleta, 
utilizado para probar cascos protectores.


```r
data(mcycle, package = "MASS")
x <- mcycle$times
y <- mcycle$accel  

library(KernSmooth)
h <- dpill(x, y) # Método plug-in de Ruppert, Sheather y Wand (1995)
fit <- locpoly(x, y, bandwidth = h) # Estimación lineal local
plot(x, y)
lines(fit)
```



\begin{center}\includegraphics[width=0.7\linewidth]{07-npreg_files/figure-latex/unnamed-chunk-2-1} \end{center}

Hay que tener en cuenta que el paquete `KernSmooth` no implementa los métodos
`predict()` y `residuals()`:


```r
est <- approx(fit, xout = x)$y # est <- predict(fit)
resid <- y - est # resid <- residuals(fit)
```

Tampoco calcula medidas de bondad de ajuste, aunque podríamos obtener fácilmente un (pseudo) R-cuadrado:


```r
r.squared <- 1 - sum(resid^2)/sum((y - mean(y))^2)
r.squared
```

```
## [1] 0.8023864
```


### Estimación de la varianza

En el caso heterocedástico, se puede obtener una estimación de la varianza 
$\sigma^2(x)$ mediante suavizado local de los residuos al cuadrado
(Fan y Yao, 1998). Mientras que en el caso homocedástico, se puede obtener 
una estimación de la varianza a partir de la suma de cuadrados residual y la 
matriz de suavizado:
$$\hat\sigma^2 = \frac{RSS}{df_e},$$
siendo $RSS=\Sigma_{i=1}^n \left( Y_i - \hat m(X_i) \right)^2$
y $df_e = tr(I - S)$ (de forma análoga al caso lineal), o alternativamente $df_e = tr \left( (I - S^{T})(I - S)\right)$, es una aproximación de los grados de libertad del error.

Adicionalmente:
$$\widehat{Var}\left(\hat{m}_{h}(x_i)\right) = \hat\sigma^2\sum_{j=1}^n s^2_{ij}$$

Uno de los pocos paquetes de `R` que implementan la estimación de la varianza
y el cálculo de intervalos de confianza es el paquete^[El paquete `np` calcula 
estimaciones similares aunque no documenta la aproximación que emplea. También 
implementa bootstrap uniforme y bootstrap por bloques.] `sm` (Bowman y Azzalini, 2019).


```r
library(sm)
hcv <- hcv(x, y) # Método de validación cruzada
fit.sm <- sm.regression(x, y, h = hcv, display = "se")
```



\begin{center}\includegraphics[width=0.7\linewidth]{07-npreg_files/figure-latex/unnamed-chunk-5-1} \end{center}

```r
fit.sm$sigma
```

```
## [1] 22.82508
```

Alternativamente se podría emplear bootstrap.


## Ejemplos

En esta sección nos centraremos en el bootstrap en la estimación tipo núcleo de la función de regresión, para la aproximación de la precisión y el sesgo, y también para el cálculo de intervalos de confianza y de predicción.

<!-- (en futuras versiones de los apuntes se incluirá también detalles sobre la construcción de bandas de confianza). -->


### Bootstrap residual

El modelo ajustado de regresión se puede emplear para estimar la respuesta media $m(x_0)$ cuando la variable explicativa toma un valor concreto $x_0$.
En este caso también podemos emplear el bootstrap residual (Sección \@ref(boot-residual)) para realizar inferencias acerca de la media.
La idea sería aproximar la distribución del error de estimación $\hat{m}(x_0) - m(x_0)$ por la distribución bootstrap de $\hat{m}^{\ast}(x_0) - \hat{m}(x_0)$.

Para reproducir adecuadamente el sesgo del estimador, la ventana $g$ ha de ser asintóticamente mayor que $h$ (de orden $n^{-1/5}$).
Análogamente al caso de la densidad, la recomendación es emplear la ventana óptima para la estimación de $m^{\prime \prime }\left( x_0 \right)$, de orden $n^{-1/9}$ (Sección \@ref(wild-bootstrap)). 


```r
n <- length(x)
g <- h * n^(4/45) # h*n^(-1/9)/n^(-1/5)
fit2 <- locpoly(x, y, bandwidth = g)
est2 <- approx(fit2, xout = x)$y # est2 <- predict(fit2)
# resid2 <- y - est2 # resid2 <- residuals(fit2)

# Remuestreo
set.seed(1)
B <- 1000
stat_fit_boot <- matrix(nrow = length(fit$x), ncol = B)
resid0 <- resid - mean(resid)
for (k in 1:B) {
    y_boot <- est2 + sample(resid0, replace = TRUE)
    fit_boot <- locpoly(x, y_boot, bandwidth = h)$y
    stat_fit_boot[ , k] <- fit_boot - fit$y
}

# Calculo del sesgo y error estándar 
bias <- apply(stat_fit_boot, 1, mean)
std.err <- apply(stat_fit_boot, 1, sd)

# Representar estimación y corrección de sesgo bootstrap
plot(x, y)
lines(fit, lwd = 2)
lines(fit$x, fit$y - bias)
```



\begin{center}\includegraphics[width=0.7\linewidth]{07-npreg_files/figure-latex/unnamed-chunk-6-1} \end{center}

NOTA: De forma análoga al caso lineal (Sección \@ref(boot-residual)), se podrían reescalar los residuos a partir de la matriz de suavizado (empleando los paquetes `sm` o `npsp`).


### Intervalos de confianza y predicción

De forma análoga al caso de la estimación de la densidad mostrado en la Sección \@ref(npden-r-ic), podemos cálcular de estimaciones por intervalo de confianza (puntuales) por el método percentil (básico):


```r
alfa <- 0.05
pto_crit <- apply(stat_fit_boot, 1, quantile, probs = c(alfa/2, 1 - alfa/2))
ic_inf_boot <- fit$y - pto_crit[2, ]
ic_sup_boot <- fit$y - pto_crit[1, ]

plot(x, y)
lines(fit, lwd = 2)
lines(fit$x, fit$y - bias)
lines(fit$x, ic_inf_boot, lty = 2)
lines(fit$x, ic_sup_boot, lty = 2)
```



\begin{center}\includegraphics[width=0.7\linewidth]{07-npreg_files/figure-latex/unnamed-chunk-7-1} \end{center}


El modelo ajustado también es empleado para predecir una nueva respuesta individual $Y(x_0)$ para un valor concreto $x_0$ de la variable explicativa. 
En el caso de errores independientes $\hat{Y}(x_0) = \hat{m}(x_0)$, pero si estamos interesados en realizar inferencias sobre el error de predicción $r(x_0) = Y(x_0) - \hat{Y}(x_0)$, a la variabilidad de $\hat{m}(x_0)$ debida a la muestra, se añade la variabilidad del error $\varepsilon(x_0)$.

La idea sería aproximar la distribución del error de predicción:
$$r(x_0) = Y(x_0) - \hat{Y}(x_0) = m(x_0) + \varepsilon(x_0) - \hat{m}(x_0)$$
por la distribución bootstrap de:
$$r^{\ast}(x_0) = Y^{\ast}(x_0) - \hat{Y}^{\ast}(x_0) = \hat{m}(x_0) + \varepsilon^{\ast}(x_0) - \hat{m}^{\ast}(x_0)$$


```r
# Remuestreo
set.seed(1)
n_pre <- length(fit$x)
stat_pred_boot <- matrix(nrow = n_pre, ncol = B)
for (k in 1:B) {
    y_boot <- est2 + sample(resid0, replace = TRUE)
    fit_boot <- locpoly(x, y_boot, bandwidth = h)$y
    pred_boot <- fit_boot + sample(resid0, n_pre, replace = TRUE)
    stat_pred_boot[ , k] <- pred_boot - fit_boot
}

# Cálculo de intervalos de predicción
# por el método percentil (básico)
alfa <- 0.05
pto_crit_pred <- apply(stat_pred_boot, 1, quantile, probs = c(alfa/2, 1 - alfa/2))
ip_inf_boot <- fit$y + pto_crit_pred[1, ]
ip_sup_boot <- fit$y + pto_crit_pred[2, ]

plot(x, y, ylim = c(-150, 75))
lines(fit, lwd = 2)
lines(fit$x, fit$y - bias)
lines(fit$x, ic_inf_boot, lty = 2)
lines(fit$x, ic_sup_boot, lty = 2)
lines(fit$x, ip_inf_boot, lty = 3)
lines(fit$x, ip_sup_boot, lty = 3)
```



\begin{center}\includegraphics[width=0.7\linewidth]{07-npreg_files/figure-latex/unnamed-chunk-8-1} \end{center}

En este caso puede no ser recomendable considerar errores i.i.d., sería de esperar heterocedásticidad (e incluso dependencia temporal).
El bootstrap residual se puede extender al caso heterocedástico y/o dependencia (e.g. Castillo-Páez *et al.*, 2019, 2020).


### Wild bootstrap {#r-wild-bootstrap}

Como se describe en la Sección \@ref(wild-bootstrap) en el caso heterocedástico se puede emplear wild bootstrap.


```r
# Remuestreo
set.seed(1)
B <- 1000
fit_boot <- matrix(nrow = n_pre, ncol = B)
for (k in 1:B) {
		rwild <- sample(c((1 - sqrt(5))/2, (1 + sqrt(5))/2), n, replace = TRUE, 
		                prob = c((5 + sqrt(5))/10, 1 - (5 + sqrt(5))/10))
    y_boot <- est2 + resid*rwild
    fit_boot[ , k] <- locpoly(x, y_boot, bandwidth = h)$y
    # OJO: bootstrap percetil directo
}
		

# Calculo del sesgo y error estándar
bias <- apply(fit_boot, 1, mean, na.rm = TRUE) -  fit$y
std.err <- apply(fit_boot, 1, sd, na.rm = TRUE)

# Representar estimación y corrección de sesgo bootstrap
plot(x, y)
lines(fit, lwd = 2)
lines(fit$x, fit$y - bias)
```



\begin{center}\includegraphics[width=0.7\linewidth]{07-npreg_files/figure-latex/unnamed-chunk-9-1} \end{center}


### Ejercicio

Siguiendo con el conjunto de datos `MASS::mcycle`, emplear wild bootstrap para 
obtener estimaciones por intervalo de confianza de la función de regresión
de `accel` a partir de `times` mediante bootstrap percentil básico.
Comparar los resultados con los obtenidos mediante bootstrap residual (comentar las diferencias y cuál de las aproximaciones sería más adecuada para este caso).
