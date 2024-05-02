# 一般の関数の最小二乗法による近似

区間 $[-\pi, \pi]$ において任意の関数 $f(x)$ を指定された n 個の関数 ${\phi_{i}(x)}$ の線形結合で以下のように近似することを考える。

$$
f(x) \approx c_{1}\phi_{1}(x) + c_{2}\phi_{2}(x) + \dots + c_{n}\phi_{n}(x)
$$

このような近似を考えるとき、次の積分のような最小二乗法を考えればよい。

$$
J = \frac{1}{2}\int_{a}^{b}\left(f(x)-\sum_{k=1}^{n}c_{k}\phi_{k}(x)\right)\mathrm{d}x \to \min
$$

このように近似できる $c_{1}$, ..., $c_{n}$ を求めるには以下を解けばよい。

$$
\frac{\partial{J}}{\partial{c_{1}}} = 0,\quad \frac{\partial{J}}{\partial{c_{2}}} = 0,\quad \dots,\quad \frac{\partial{J}}{\partial{c_{n}}} = 0
$$

$i$ 番目の式の左辺について考えると、

$$
\begin{align*}
\frac{\partial{J}}{\partial{c_{i}}}
&= -\int_{a}^{b}\left(f(x)-\sum_{k=1}^{n}c_{k}\phi_{k}(x)\right)\phi_i(x)\mathrm{d}x \\
&= \sum_{k=1}^{n}c_{k}\int_{a}^{b}\phi_{k}(x)\phi_{i}(x)\mathrm{d}x - \int_{a}^{b}f(x)\phi_{i}(x)\mathrm{d}x
\end{align*}
$$

これが 0 と等しくなるとき、以下の式が成り立つことになる。

$$
\sum_{k=1}^{n}c_{k}\int_{a}^{b}\phi_{k}(x)\phi_{i}(x)\mathrm{d}x = \int_{a}^{b}f(x)\phi_{i}(x)\mathrm{d}x
$$

これが $i = 1, \dots, n$ に対して成り立つので、それらを行列で表記すると以下の正規方程式が得られる。

$$
\begin{pmatrix}
\int_{a}^{b}\phi_{1}(x)^2\mathrm{d}x & \int_{a}^{b}\phi_{1}(x)\phi_{2}(x)\mathrm{d}x & \dots & \int_{a}^{b}\phi_{1}(x)\phi_{n}(x)\mathrm{d}x \\
\int_{a}^{b}\phi_{2}(x)\phi_{1}(x)\mathrm{d}x & \int_{a}^{b}\phi_{2}(x)^2\mathrm{d}x & \dots & \int_{a}^{b}\phi_{2}(x)\phi_{n}(x)\mathrm{d}x \\
\vdots & \vdots & \ddots & \vdots \\
\int_{a}^{b}\phi_{n}(x)\phi_{1}(x)\mathrm{d}x & \int_{a}^{b}\phi_{n}(x)\phi_{2}(x)\mathrm{d}x & \dots & \int_{a}^{b}\phi_{n}(x)^2\mathrm{d}x
\end{pmatrix}
\begin{pmatrix}
c_{1} \\
c_{2} \\
\vdots \\
c_{n}
\end{pmatrix}=
\begin{pmatrix}
\int_{a}^{b}\phi_{1}(x)f(x)\mathrm{d}x \\
\int_{a}^{b}\phi_{2}(x)f(x)\mathrm{d}x \\
\vdots \\
\int_{a}^{b}\phi_{n}(x)f(x)\mathrm{d}x
\end{pmatrix}
$$

# 直交関数の最小二乗法による近似

## 直交とは

区間 $[a, b]$ 上の関数 $f(x), g(x)$ が $\int_{a}^{b}f(x)g(x)\mathrm{d}x = 0$ を満たすとき、 $f(x), g(x)$ は直交するという。
また、関数 $\phi_{0}(x), \phi_{1}(x), \dots, \phi_{n}(x)$ が $\int_{a}^{b}\phi_{i}(x)\phi_{j}(x)\mathrm{d}x = 0,\quad i \ne j$ を満たすとき、これらの関数は区間 $[a, b]$ 上の直交関数系であるという。

## 直交関数系で近似するときの正規方程式

ある関数 $f(x)$ を直交関数系 $\phi_{0}(x), \phi_{1}(x), \dots, \phi_{n}(x)$ で $f(x) \approx c_{1}\phi_{1}(x) + c_{2}\phi_{2}(x) + \dots + c_{n}\phi_{n}(x)$ のように近似するとき、 $\int_{a}^{b}\phi_{i}(x)\phi_{j}(x)\mathrm{d}x = 0,\quad i \ne j$ が成り立つので前述した正規方程式は以下のように変形できる。

$$
\begin{pmatrix}
\int_{a}^{b}\phi_{1}(x)^2\mathrm{d}x &  & &  \\
& \int_{a}^{b}\phi_{2}(x)^2\mathrm{d}x & & \\
& & \ddots & \\
& & & \int_{a}^{b}\phi_{n}(x)^2\mathrm{d}x
\end{pmatrix}
\begin{pmatrix}
c_{1} \\
c_{2} \\
\vdots \\
c_{n}
\end{pmatrix}=
\begin{pmatrix}
\int_{a}^{b}\phi_{1}(x)f(x)\mathrm{d}x \\
\int_{a}^{b}\phi_{2}(x)f(x)\mathrm{d}x \\
\vdots \\
\int_{a}^{b}\phi_{n}(x)f(x)\mathrm{d}x
\end{pmatrix}
$$

すなわち、近似するときの係数は以下の式でただちに求められる。

$$
c_{i} = \frac{\int_{a}^{b}\phi_{i}(x)f(x)\mathrm{d}x}{\int_{a}^{b}\phi_{i}(x)^2\mathrm{d}x},\quad i = 1, \dots, n
$$

# フーリエ級数とフーリエ係数

唐突ではあるが、ここで区間 $[-\pi, \pi]$ における $\frac{1}{2},\quad \cos{kx},\quad \sin{kx},\quad k = 1, 2, \dots$ について考えると、以下のようにしてこれらの関数が直交関数系であることを示すことができる。

$\cos{kx},\quad \sin{kx},\quad k = 1, 2, \dots$ は周期 $2\pi$ の周期関数なので、 $\int_{-\pi}^{\pi}\cos{kx}\mathrm{d}x$, $\int_{-\pi}^{\pi}\sin{kx}\mathrm{d}x$ は 0 である。このことから、

$$
\int_{-\pi}^{\pi}\frac{1}{2}\cos{kx}\mathrm{d}x = 0 \\
\int_{-\pi}^{\pi}\frac{1}{2}\sin{kx}\mathrm{d}x = 0
$$

また、

$$
\begin{align*}
\int_{-\pi}^{\pi}\cos{kx}\sin{lx}\mathrm{d}x
&= \int_{-\pi}^{\pi}\frac{1}{2}\left(\sin{(k+l)x} - \sin{(k-l)x}\right)\mathrm{d}x \\
&= \frac{1}{2}\int_{-\pi}^{\pi}\sin{(k+l)x}\mathrm{d}x - \frac{1}{2}\int_{-\pi}^{\pi}\sin{(k-l)x}\mathrm{d}x \\
&= 0
\end{align*}
$$

$$
\begin{align*}
\int_{-\pi}^{\pi}\cos{kx}\cos{lx}\mathrm{d}x
&= \int_{-\pi}^{\pi}\frac{1}{2}\left(\cos{(k+l)x} + \cos{(k-l)x}\right)\mathrm{d}x \\
&= \frac{1}{2}\int_{-\pi}^{\pi}\cos{(k+l)x}\mathrm{d}x + \frac{1}{2}\int_{-\pi}^{\pi}\cos{(k-l)x}\mathrm{d}x \\
&= 0
\end{align*}
$$

$$
\begin{align*}
\int_{-\pi}^{\pi}\sin{kx}\sin{lx}\mathrm{d}x
&= \int_{-\pi}^{\pi}\frac{1}{2}\left(\cos{(k+l)x} - \cos{(k-l)x}\right)\mathrm{d}x \\
&= \frac{1}{2}\int_{-\pi}^{\pi}\cos{(k+l)x}\mathrm{d}x - \frac{1}{2}\int_{-\pi}^{\pi}\cos{(k-l)x}\mathrm{d}x \\
&= 0
\end{align*}
$$

以上より、区間 $[-\pi, \pi]$ における $\frac{1}{2},\quad \cos{kx},\quad \sin{kx},\quad k = 1, 2, \dots$ が直交関数系であることが示された。

関数 $f(x)$ を区間 $[-\pi, \pi]$ 上で直交関数系 $\frac{1}{2},\quad \cos{kx},\quad \sin{kx},\quad k = 1, 2, \dots$ で以下のように近似する。

$$
f(x) \approx \frac{a_o}{2} + a_{1}\cos{x} + b_{1}\sin{x} + a_{2}\cos{2x} + b_{2}\sin{2x} + \dots
$$

これは直交関数系による近似なので、正規方程式より $a_{0}, a_{1}, b_{1}, \dots$ は以下のように計算できる。

$$
\begin{align*}
a_{0}
&= \frac{\int_{-\pi}^{\pi}\frac{1}{2}f(x)\mathrm{d}x}{\int_{-\pi}^{\pi}(\frac{1}{2})^2\mathrm{d}x} \\
&= \frac{\frac{1}{2}\int_{-\pi}^{\pi}f(x)\mathrm{d}x}{\frac{\pi}{2}} \\
&= \frac{1}{\pi}\int_{-\pi}^{\pi}f(x)\mathrm{d}x
\end{align*}
$$

$$
\begin{align*}
a_{k}
&= \frac{\int_{-\pi}^{\pi}\cos{kx}f(x)\mathrm{d}x}{\int_{-\pi}^{\pi}(\cos{kx})^2\mathrm{d}x} \\
&= \frac{\int_{-\pi}^{\pi}\cos{kx}f(x)\mathrm{d}x}{\frac{1}{2}\int_{-\pi}^{\pi}(1 + \cos{2kx})\mathrm{d}x} \\
&= \frac{\int_{-\pi}^{\pi}\cos{kx}f(x)\mathrm{d}x}{\frac{1}{2}2\pi} \\
&= \frac{1}{\pi}\int_{-\pi}^{\pi}\cos{kx}f(x)\mathrm{d}x
\end{align*}
$$

$$
\begin{align*}
b_{k}
&= \frac{\int_{-\pi}^{\pi}\sin{kx}f(x)\mathrm{d}x}{\int_{-\pi}^{\pi}(\sin{kx})^2\mathrm{d}x} \\
&= \frac{\int_{-\pi}^{\pi}\sin{kx}f(x)\mathrm{d}x}{\frac{1}{2}\int_{-\pi}^{\pi}(1 - \cos{2kx})\mathrm{d}x} \\
&= \frac{\int_{-\pi}^{\pi}\sin{kx}f(x)\mathrm{d}x}{\frac{1}{2}2\pi} \\
&= \frac{1}{\pi}\int_{-\pi}^{\pi}\sin{kx}f(x)\mathrm{d}x
\end{align*}
$$

すなわち、区間 $[-\pi, \pi]$ 上の連続関数 $f(x)$ は以下のように近似できる。

$$
f(x) = \frac{a_{0}}{2} + \sum_{k=1}^{\infty}{(a_{k}\cos{kx} + b_{k}\sin{kx})}
$$

$$
a_{0} = \frac{1}{\pi}\int_{-\pi}^{\pi}f(x)\mathrm{d}x,\quad
a_{k} = \frac{1}{\pi}\int_{-\pi}^{\pi}f(x)\cos{kx}\mathrm{d}x,\quad
b_{k} = \frac{1}{\pi}\int_{-\pi}^{\pi}f(x)\sin{kx}\mathrm{d}x
$$

この $\frac{1}{2},\quad \cos{kx},\quad \sin{kx},\quad k = 1, 2, \dots$ を用いた直交関数展開をフーリエ級数と呼び、 $a_{0}, a_{k}, b_{k}$ の係数をフーリエ係数と呼ぶ。

これまでは区間 $[-\pi, \pi]$ での $\frac{1}{2},\quad \cos{kx},\quad \sin{kx},\quad k = 1, 2, \dots$ を用いた直交関数展開を考えてきたが、ここで一般の幅 $T$ の区間 $[-\frac{T}{2}, \frac{T}{2}]$ でのフーリエ級数を考える。上式において $x = \frac{2\pi t}{T}$ と置くと、 $-\pi \leqq x \leqq \pi$ のとき $-\frac{T}{2} \leqq t \leqq \frac{T}{2}$、 $\mathrm{d}x = \frac{2\pi}{T}\mathrm{d}t$ となる。さらに範囲 $[-\frac{T}{2}, \frac{T}{2}]$ において $\omega_{o} = \frac{2\pi}{T}$ と置くと $x = \omega_{o}t,\quad dx = \omega_{o}dt$ となり、

$$
f(t) = \frac{a_{0}}{2} + \sum_{k=1}^{\infty}(a_{k}\cos{k\omega_{o}t + b_{k}\sin{k\omega_{o}t}})
$$

$$
a_{k} = \frac{2}{T}\int_{-\frac{T}{2}}^{\frac{T}{2}}f(t)\cos{k\omega_{o}t}\mathrm{d}t,\quad
b_{k} = \frac{2}{T}\int_{-\frac{T}{2}}^{\frac{T}{2}}f(t)\sin{k\omega_{o}t}\mathrm{d}t
$$

と表すことができる。この $\omega_{o}$ を基本周波数、 $a_{k}, b_{k}$ を実フーリエ係数と呼ぶ。

# オイラーの式と複素フーリエ係数
指数部が虚数の指数関数 $e^{i\theta}$ を以下のように定義する。

$$
e^{i\theta} = \cos{\theta} + i\sin{\theta}
$$

この式をオイラーの式と呼ぶ。

オイラーの式において $\theta = -\theta$ と置き換えると

$$
e^{-i\theta} = \cos{\theta} - i\sin{\theta}
$$

これら 2 式を足し引きすることで以下の関係をみちびける。

$$
\cos{\theta} = \frac{e^{i\theta} + e^{-i\theta}}{2},\quad
\sin{\theta} = \frac{e^{i\theta} - e^{-i\theta}}{2i}
$$

この式を用いてフーリエ級数を変形すると、

$$
\begin{align*}
f(t)
&= \frac{a_{0}}{2} + \sum_{k=1}^{\infty}\left(a_{k}\frac{e^{ik\omega_{o}t} + e^{-ik\omega_{o}t}}{2} + b_{k}\frac{e^{ik\omega_{o}t} - e^{-ik\omega_{o}t}}{2i}\right) \\
&= \frac{a_{0}}{2} + \frac{1}{2}\sum_{k=1}^{\infty}\left((a_{k} - ib_{k})e^{ik\omega_{o}t}+(a_{k} + ib_{k})e^{-ik\omega_{o}t}\right) \\
&= \sum_{-\infty}^{\infty}C_{k}e^{ik\omega_{o}t}
\end{align*}
$$

ただし、

$$
C_{k} =
\begin{cases}
\frac{(a_{k} - ib_{k})}{2} & k \gt 0 \\
\frac{a_{0}}{2} & k = 0 \\
\frac{(a_{-k} + ib_{-k})}{2} & k \lt 0
\end{cases}
$$

$C_{k}$ の各場合分けについて、 $a_{k} = \frac{2}{T}\int_{-\frac{T}{2}}^{\frac{T}{2}}f(t)\cos{k\omega_{o}t}\mathrm{d}t,\quad b_{k} = \frac{2}{T}\int_{-\frac{T}{2}}^{\frac{T}{2}}f(t)\sin{k\omega_{o}t}\mathrm{d}t$ を代入すると

$$
\begin{align*}
\frac{(a_{k} - ib_{k})}{2}
&= \frac{2}{T}\int_{-{\frac{T}{2}}}^{\frac{T}{2}}f(t)\frac{\cos{k\omega_{o}t} - i\sin{k\omega_{o}t}}{2}\mathrm{d}t \\
&= \frac{1}{T}\int_{-{\frac{T}{2}}}^{\frac{T}{2}}f(t)e^{-ik\omega_{o}t}\mathrm{d}t
\end{align*}
$$

$$
\frac{a_{0}}{2} = \frac{1}{T}\int_{-\frac{T}{2}}^{\frac{T}{2}}f(t)\mathrm{d}t
$$

$$
\begin{align*}
\frac{(a_{-k} + ib_{-k})}{2}
&= \frac{2}{T}\int_{-{\frac{T}{2}}}^{\frac{T}{2}}f(t)\frac{\cos{(-k\omega_{o}t)} + i\sin{(-k\omega_{o}t)}}{2}\mathrm{d}t \\
&= \frac{1}{T}\int_{-{\frac{T}{2}}}^{\frac{T}{2}}f(t)(\cos{k\omega_{o}t} - i\sin{k\omega_{o}t)}\mathrm{d}t \\
&= \frac{1}{T}\int_{-{\frac{T}{2}}}^{\frac{T}{2}}f(t)e^{-ik\omega_{o}t}\mathrm{d}t
\end{align*}
$$

すなわち、以下のようにまとめられる。

$$
f(t) = \sum_{k=-\infty}^{\infty}C_{k}e^{ik\omega_{o}t},\quad C_{k} = \frac{1}{T}\int_{-{\frac{T}{2}}}^{\frac{T}{2}}f(t)e^{-ik\omega_{o}t}\mathrm{d}t
$$

$C_{k}$ を複素フーリエ係数と呼ぶ。

# フーリエ変換とたたみこみ積分 (合成積)

複素フーリエ係数を用いたフーリエ級数への展開の式において、 $\omega_{k} = k\omega_{o},\quad k = 0, \pm1, \pm2, \dots$ と置き、

$$
\Delta\omega = \omega_{k+1} - \omega_{k} = \omega_{o} = \frac{2\pi}{T},\quad
C_{k} = \frac{F(\omega_{k})}{T}
$$

と書くと、以下のような式が得られる。

$$
\begin{align*}
f(t)
&= \sum_{k=-\infty}^{\infty}C_{k}e^{ik\omega_{o}t} \\
&= \sum_{k=-\infty}^{\infty}\frac{F(\omega_{k})}{T}e^{i\omega_{k}t} \\
&= \frac{1}{2\pi}\sum_{k=-\infty}^{\infty}F(\omega_{k})e^{i\omega_{k}t}\Delta\omega
\end{align*}
$$

$$
\begin{align*}
F(\omega_{k})
&= TC_{k} \\
&=\int_{-{\frac{T}{2}}}^{\frac{T}{2}}f(t)e^{-ik\omega_{o}t}\mathrm{d}t \\
&=\int_{-{\frac{T}{2}}}^{\frac{T}{2}}f(t)e^{-i\omega_{k}t}\mathrm{d}t
\end{align*}
$$

周期 $T \to \infty$ とすると $\Delta\omega \to 0$ となり、 $f(t)$ も積分で表すことができ、以下のようになる。

$$
f(t) = \frac{1}{2\pi}\int_{-\infty}^{\infty}F(\omega)e^{i\omega t}\mathrm{d}\omega,\quad
F(\omega) = \int_{-\infty}^{\infty}f(t)e^{-i\omega t}\mathrm{d}t
$$

第 2 式を信号 $f(t)$ のフーリエ変換、その逆変換である第 1 式を逆フーリエ変換と呼ぶ。逆フーリエ変換は信号 $f(t)$ をあらゆる周波数の振動の重ね合わせで表すものである。 $F(\omega)$ は周波数 $\omega$ の成分 $e^{i\omega t}$ の大きさを表し、 $f(t)$ のスペクトルと呼ばれる。

信号 $f(t), g(t)$ のたたみこみ積分 (合成積) $f(t) * g(t)$ は以下のように定義される。

$$
f(t) * g(t) = \int_{-\infty}^{\infty}f(s)g(t-s)\mathrm{d}s
$$

信号 $f(t), g(t)$ のフーリエ変換をそれぞれ $F(\omega), G(\omega)$ とするとき、たたみこみ積分 $f(t)*g(t)$ のフーリエ変換は $F(\omega)G(\omega)$ になる。これは以下のような式変形によって示すことができる。

$$
\begin{align*}
\int_{-\infty}^{\infty}f(t)*g(t)e^{-i\omega t}\mathrm{d}t
&= \int_{-\infty}^{\infty}\left(\int_{-\infty}^{\infty}f(s)g(t-s)\mathrm{d}s\right)e^{-i\omega t}\mathrm{d}t \\
&= \int_{-\infty}^{\infty}f(s)\left(\int_{-\infty}^{\infty}g(t-s)e^{-i\omega t}\mathrm{d}t\right)\mathrm{d}s \\
\end{align*}
$$

ここで $t' = t - s$ と置くと $t = t' + s,\quad dt = dt'$ なので

$$
\begin{align*}
\int_{-\infty}^{\infty}f(s)\left(\int_{-\infty}^{\infty}g(t-s)e^{-i\omega t}\mathrm{d}t\right)\mathrm{d}s
&= \int_{-\infty}^{\infty}f(s)\left(\int_{-\infty}^{\infty}g(t')e^{-i\omega(t'+s)}\mathrm{d}t'\right)\mathrm{d}s \\
&= \int_{-\infty}^{\infty}f(s)\left(\int_{-\infty}^{\infty}g(t')e^{-i\omega t'}\mathrm{d}t'\right)e^{-i\omega s}\mathrm{d}s \\
&= \left(\int_{-\infty}^{\infty}f(s)e^{-i\omega s}\mathrm{d}s\right)\left(\int_{-\infty}^{\infty}g(t')e^{-i\omega t'}\mathrm{d}t'\right) \\
&= F(\omega)G(\omega)
\end{align*}
$$