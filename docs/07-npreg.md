# Bootstrap y regresión no paramétrica {#cap7}




En este capítulo se introducirá la estimación no paramétrica de la función de regresión y distintos métodos de remuestreo y resultados para ellos...


## Estimación no paramétrica de la función de regresión

Sea $\left\{ \left( X_1,Y_1 \right),\left( X_2,Y_2 \right), \ldots, \left( X_n,Y_n \right) \right\}$ una m.a.s. de una población
bidimensional $\left( X,Y \right)$, con $E\left( \left\vert Y\right\vert
 \right) <\infty$, para la cual queremos estimar la función de regresión
de $Y$ dada $X$: 
$$m\left( x \right) =E\left( \left. Y\right\vert_{X=x} \right).$$
La función de regresión puede escribirse así:
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
x-X_i \right) Y_i}{\frac{1}{n}\sum_{i=1}^{n}K_{h}\left( x-X_i \right)}.$$

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

## Distintos métodos de remuestreo y resultados para ellos

### Wild bootstrap {#wild-bootstrap}

Este método de remuestreo bootstrap, propuesto por Wu (1986) y estudiado
por Härdle y Marron (1991), procede del siguiente modo:

1.  A partir del estimador de Nadaraya-Watson de $m\left( x \right)$ y
    tomando el parámetro ventana de partida, $h$, se construyen los
    residuos
    $\hat{\varepsilon}_i = Y_i - \hat{m}_{h}\left( X_i \right)$,
    $i=1, 2, \ldots, n$.

2.  Para cada índice $i=1,2,\ldots ,n$, se arroja, condicionalmente a la
    muestra observada, $\left\{ \left( X_1,Y_1 \right), \ \left(
    X_2,Y_2 \right),\right.$ 
    $\left.\ldots ,\ \left( X_n,Y_n \right) \right\}$, 
    un residuo bootstrap $\hat{\varepsilon}_i^{\ast}$ de una
    distribución de probabilidad que cumpla,
    $E^{\ast}\left( \hat{\varepsilon}_i^{\ast} \right) =0$,
    $E^{\ast}\left( \hat{\varepsilon}_i^{\ast 2} \right) =\hat{
    \varepsilon}_i^2$ y $E^{\ast}\left( \hat{\varepsilon}_i^{\ast
    3} \right) =\hat{\varepsilon}_i^{3}$. Aunque la condición del
    momento de orden 3 no es estrictamente necesaria, es útil para las
    demostraciones de validez del método.

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
definiendo $\hat{\varepsilon}_i^{\ast}=\hat{\varepsilon}_iV_i^{\ast}$ 
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

