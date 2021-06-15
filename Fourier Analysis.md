# 傅里叶分析笔记

Yeah

[TOC]



## Chapter 1. The Genesis of Fourier Analysis

### 行波、驻波与波传导方程

首先需要引入的是用行波法和驻波法求解波传导方程。何为驻波？
$$
u(x,t)=\varphi(x)\psi(t)
$$
这一个定义的核心就是分离变量。我们这个波并不会乱变，给定一个由$\varphi(x)$和$t_0$规定的初相，之后的变化都由$\psi(t)$决定。

行波法代表着随着时间变化，波本身具有的局部形状不变，而体现出逐步位移的性质。
$$
u(x,t)=F(x-ct)
$$
之后随着物理推导得出波传导方程：
$$
\frac{1}{c^2}\frac{\part^2u}{\part t^2}=\frac{\part^2u}{\part x^2}
$$
之后用比例方所变换将方程化简为这样一个形式
$$
\frac{\part^2u}{\part t^2}=\frac{\part^2u}{\part x^2}\tag1\\on\,\,\, 0\leq x\leq \pi, with \,\,\,t\geq0.
$$

#### 行波法

行波法怎么解？观察到$F(x+t)$和$F(x-t)$都能解出这样一个方程。又根据线性，观察到：
$$
u(x,t)=F(x+t)+G(x-t)
$$
物理意义就是方向相反的两个行波交叉可以解这一方程。这一点可以用线性变换从(1)中积分积出来。之后的限制就在于，初值问题，还有定点问题。转换为条件即为：
$$
u(x,0)=f(x)\\
u(0,t)=u(\pi,t)=0 \quad\forall t\geq 0\quad0\leq x\leq \pi
$$
下一步是延拓函数，让f和u为奇函数，且以$2\pi$为周期。这样把x的区间拓宽了，不用困在小区间。

还需要一个作为速度的的处置条件。
$$
\frac{\part u}{\part t}(x,0)=g(x)
$$
联立，积分，得到d'Alembert's formula
$$
u(x,t)=\frac{1}{2}[f(x+t)+f(x-t)]+\frac{1}{2}\int_{x-t}^{x+t}g(y)dy
$$

#### 驻波法

先分离变量
$$
u(x,t)=\varphi(x)\psi(t)
$$
带入pde，分离变量得：
$$
\begin{equation}
\left\{
\begin{array}{lr}
\psi''(t)-\lambda\psi(t)=0\\
\varphi''(x)-\lambda\varphi(x)=0
\end{array}
\right.
\end{equation}
$$
俩简单ODE，解出来
$$
\begin{equation}
\left\{
\begin{array}{lr}
\psi(t)=Acosmt+Bsinmt\\
\varphi(x)=\tilde{A}cosmt+\tilde{B}sinmt
\end{array}
\right.
\end{equation}
$$
然后前后固定点，得到$\tilde{A}=0$，m必须是一个整数，所以对任意正整数m：
$$
u_m(x,t)=(A_mcosmt+B_msinmt)sinmx
$$
随后叠加。
$$
u(x,t)=\sum_{m=1}^{\infin} (A_mcosmt+B_msinmt)sinmx
$$
如果这玩意给出了所有解，那么初值带进去吧：
$$
\sum_{m=1}^{\infin} A_msinmx=f(x)
$$
哎，这就引出傅里叶级数了！看看怎么推出$A_m$:
$$
\int_0^\pi f(x)sinnxdx=\int_0^\pi(\sum_{m=1}^\infin A_msinmx)sinnxdx=\frac{\pi}{2}A_n
$$

$$
A_n=\frac{2}{\pi}\int_0^\pi f(x)sinnxdx
$$

真是太有趣了！$f(x)$这一个作为所有以$A_m$为系数的傅里叶级数的加总，反倒可以显身于每一个$A_n$的决定过程之中。

随后开始傅里叶级数的引入了！f(x)为奇函数，g(x)为偶函数。一个初等函数的结论是，一个函数能表示为奇函数和偶函数的和。
$$
F(x)=\sum_{m=1}^{\infin} A_msinmx+\sum_{m=0}^{\infin} A'_mcosmx
$$
用欧拉等式呢，如果可以将傅里叶系数设为复数，那么这么写就很简单了。
$$
F(x)=\sum_{m=-\infin}^{\infin} a_me^{imx}
$$
用类似的积分手法处理，可以得出傅里叶系数：
$$
a_n=\frac{1}{2\pi}\int_{-\pi}^\pi F(x)e^{-inx}dx
$$

### 热传导

$$
\frac{\sigma}{\kappa}\frac{\part u}{\part t}=\frac{\part^2 u}{\part x^2}+\frac{\part^2 u}{\part y^2}
$$

这是time-dependent heat equation。如果已经过去很久，再也无时变，那么就有steady-state heat equation.
$$
\frac{\part^2 u}{\part x^2}+\frac{\part^2 u}{\part y^2}=0\\
\Delta u=0
$$
这一个方程的解就是所谓的调和函数。

考察单位圆与其边界：
$$
D=\{(r,\theta):0\leq r < 1\}\quad and \quad C=\{(r,\theta):r=1\}
$$
所以这一个方程的条件就是：
$$
\Delta u=0
\\
u(1,\theta)=f(\theta)
$$
极坐标下呢：
$$
\Delta u=\frac{\part^2 u}{\part r^2}+\frac{1}{r}\frac{\part u}{\part r}+\frac{1}{r^2}\frac{\part^2 u}{\part \theta^2}=0
$$
然后分离变量！
$$
\frac{r^2F''(r)+rF'(r)}{F(r)}=-\frac{G''(\theta)}{G(\theta)}=\lambda
$$

$$
G(\theta)=Ae^{im\theta}+Be^{-im\theta}
$$

对于F呢，当$\lambda=m^2$,$F(r)=r^m$或$F(r)=r^{-m}$

并在一起得到$u_m$，再加总得到最后结果
$$
u_m(r,\theta)=r^{|m|}e^{im\theta},\quad m\in\mathbb{Z}
\\
u(r,\theta)=\sum_{m=-\infin}^{\infin}a_mr^{|m|}e^{im\theta}
$$
来一个初值条件：
$$
u(1,\theta)=\sum_{m=-\infin}^{\infin}a_me^{im\theta}=f(\theta)
$$


### Homework

Constructing......

## Chapter 2. Basic Properties of Fourier Series

### 傅里叶级数的正式介绍

For a $f$ defined on $[0,L]$, the Fourier coefficients are:
$$
a_n=\frac{1}{L}\int_0^Lf(x)e^{-2\pi inx/L}dx,\quad for \,n\in\mathbb{Z}
$$
可积性怎么说？首先，复值函数为主。另外，函数定义在圆上，也就是极坐标。假定本书所选用的函数基本都可积。

单位圆上的点都可以用$e^{i\theta}$ 表示。若$F$ 是定义在圆上的函数，则我们可以定义：
$$
f(\theta)=F(e^{i\theta})
$$
f周期为$2\pi$。可积，连续，可微均能传递。

正式写起来：
$$
\hat{f}(n)=\frac{1}{L}\int_a^bf(x)e^{-2\pi inx/L}dx,\quad for \,n\in\mathbb{Z}
$$
这是傅里叶系数

那么傅里叶级数：
$$
\sum_{n=-\infin}^{\infin} \hat{f}(n)e^{2\pi inx/L}
$$
则我们用这个标记表示右侧是f的傅里叶级数
$$
f~\sim\sum_{n=-\infin}^{\infin} \hat{f}(n)e^{2\pi inx/L}
$$
傅里叶级数是三角级数的一种。

我们记其Nth Partial Sum：
$$
S_N(f)(x)=\sum_{n=-N}^N\hat{f}(n)e^{2\pi inx/L}
$$
对称的。

比较重要的是Kernel。Nth Dirichlet kernel is defined as:
$$
D_N(x)=\sum_{n=-N}^Ne^{inx}
$$
其有闭式解：
$$
D_N(x)=\frac{sin((N+\frac{1}{2})x)}{sin(x/2)}
$$
Poisson kernel is:
$$
P_r(\theta)=\sum_{n=-\infin}^{\infin}r^{|n|}e^{in\theta}\\
=\frac{1-r^2}{1-2rcos\theta+r^2}
$$
问题，是否收敛于f？问题的提法：
$$
\lim_{N\to\infin}S_N(f)(\theta)=f(\theta)~~~~~~for ~every~\theta?
$$
首先，所有点一起收敛，不可能。因为该点不影响可积。若f连续周期，其实也不是。

### 傅里叶级数的唯一性

有限不同点的函数都有相同的傅里叶级数。

Theorem 2.1 若f在单位圆可积，且傅里叶系数为0，则$f(\theta_0)=0$若f在该点连续。

Proof:

首先
$$
\hat{f}(n)=0
\\which ~means~for~every~n
\\
\int_{-\pi}^\pi f(\theta)e^{i\theta}d\theta=0
\\
Using ~e^{ix}=cosx+isinx, and ~the ~addictivity ~of~intergrals:
\\
for ~every~n:
\int_{-\pi}^\pi f(\theta)cos(n\theta)d\theta=\int_{-\pi}^\pi f(\theta)sin(n\theta)d\theta=0
$$
而由于三角多项式对乘积操作封闭，则可以知道对于我们选定的
$$
p(\theta)=\epsilon+cos\theta
\\
p_k(\theta)=[p(\theta)]^k
$$
都是三角多项式，也就是三角函数的线性组合，也就可以解释说：
$$
\int_{-\pi}^\pi f(\theta)p_k(\theta)d\theta=0
$$
是必须满足的。

那如何估计？

B是函数绝对值上限。我们选一个$\delta$，让这个区间内的函数值都大于$f(0)/2$。反正$\delta$可以再小，小了以后区间外的函数就有$p(\theta)<1-\epsilon/2$，再取$\eta<\delta$让更小区间内有$p(\theta)\geq1+\epsilon/2$

有了这些小伎俩，我们可以估计了。首先，区间外积分的绝对值怎么看呢？
$$
\left|\int_{\delta\leq|\theta|} f(\theta)p_k(\theta)d\theta\right|\leq B\left|\int_{\delta\leq|\theta|} p_k(\theta)d\theta\right|\\
\leq B(1-\epsilon/2)^k\left|\int_{\delta\leq|\theta|} d\theta\right|\leq 2\pi B(1-\epsilon/2)^k
$$
又$\eta,\delta$之间的区间内，p和f都大于零，这一部分积分也大于0.

外面两边无所谓了。就看里面
$$
\int_{|\theta|<\eta} f(\theta)p_k(\theta)d\theta
$$
此时p_k很大，越来越大。eta，f也有取值，区间也有，就会越来越大，随着k不断往上。需要注意的是我们取区间不是根据p_k取，而是p取，这样有保证。

确实这么统摄，但要找到重点，才可以。

#### Corollary 2.3  傅里叶级数的一致收敛

假设f是一个单位圆上的连续函数，傅里叶系数绝对可加收敛，那么傅里叶级数一致收敛到f。
$$
\lim_{N\to\infin}S_N(f)(\theta)=f(\theta)~~~uniformly~in~\theta
$$
Proof:

傅里叶系数绝对可加收敛，则：
$$
g(\theta)=\lim_{N\to\infin}\sum_{n=-N}^N\hat{f}(n)e^{in\theta}
$$
连续，由连续函数序列一致收敛，则极限一致收敛得。

#### Corollary 2.4 傅里叶系数的收敛速度

假设两次连续可微，单位元，那么
$$
\hat{f}(n)=O(1/|n|^2)~~~as~|n|\to\infin
$$
证明有点长，记一下思路吧。

先由傅里叶系数的定义积分，然后对$e^{-in\theta}$分部积分，积出来那个项是0。然后再分部积分！那个项还是0。为何不能再动了呢？再动下去是不是i不好处理？并不是。越平滑，收敛越快。(Exercise 10)

顺便得到小钥匙：
$$
\hat{f}'(n)=in\hat{f}(n)
\\f\sim\sum a_ne^{in\theta}
\\f'\sim\sum a_nine^{in\theta}
$$
另外，Holder condition也可以用来证明其收敛。

### 卷积

卷起来！

给定两个以$2\pi$为周期的可积函数f和g，在R上，那么他们的卷积$f*g~on~[-\pi,\pi]$定义为：
$$
(f*g)(x)=\frac{1}{2\pi}\int_{-\pi}^\pi f(y)g(x-y)dy=\frac{1}{2\pi}\int_{-\pi}^\pi f(x-y)g(y)dy
$$
这是如何与傅里叶级数产生关系的？
$$
S_N(f)(x)=\sum_{n=-N}^N\hat{f}(n)e^{-inx}\\=\sum_{n=-N}^N\frac{1}{2\pi}\int_{-\pi}^\pi f(y)e^{iny}dye^{-inx}
\\=\frac{1}{2\pi}\int_{-\pi}^\pi f(y)(\sum_{n=-N}^N e^{in(x-y)})dy
\\
=(f*D_N)(x)
\\
where~ D_N(x)=\sum_{n=-N}^N e^{inx}
$$
傅里叶部分和可以表现为f和迪利克雷核的卷积结果

卷积的主要性质：
$$
f*(g+h)=(f*g)+(f*h)\\
(cf)*g=c(f*g)=f*(cg) ~for~any~c\in \mathbb{C}
\\
f*g=g*f
\\(f*g)*h=f*(g*h)
\\f*g~is~continuous\\
\widehat{f*g}(n)=\hat{f}(n)\hat{g}(n)
$$
最后一个最重要。f和g的卷积的傅里叶系数等于f的傅里叶系数乘以g的傅里叶系数，非常关键！

###  Good kernels

首先是一系列三角多项式。

Definition: A family of kernels $\{K_n(x)\}_{n=1}^\infin$ on the circle is said to be a family of good kernels if it satisfies the following properties:

(a) For all $n \geq1$,
$$
\frac{1}{2\pi}\int_{-\pi}^\pi K_n(x)dx=1
$$
(b) There exists M>0 such that for all n $\geq 1$,
$$
\int_{-\pi}^\pi |K_n(x)|dx\le M
$$
(c) For every $\delta>0$
$$
\int_{\delta\le|x|\le\pi}|K_n(x)|dx\to0,as~n\to\infin
$$
a就是，单位区间上质量综合为1.c就是越靠近0，质量越集中，当n变大。

好kernel的作用还是要做到当n趋向于无穷时，中间的那坨东西会让原函数与他的卷积收敛到原函数

#### Theorem 4.1

Let $\{K_n\}^\infin_{n=1}$ be a family of good kernels, and f an integrable function on the circle, then:
$$
\lim_{n\to\infin} (f*K_n)(x)=f(x)
$$
whenever f is continuous on x. If f is continuous everywhere, then the above limit is uniform.

加权平均嘛，在x点处逐渐有full mass。

Proof的思路：

控制$|f(x-y)-f(x)|<\epsilon$，于是$(f*K_n)(x)-f(x)$可以分段估计，小于$\delta$那一段肯定epsilon控制以内，大于$\delta$那一段，由good kernel性质三可以得到小于$\epsilon$
很可惜，迪利克雷核不是good kernel。所以出问题，傅里叶级数的收敛就很艰难。

### Cesaro和Abel求和

开始sum！

#### Cesari Sum:

$s_n=\sum_{k=0}^nc_k$

define:
$$
\sigma_N=\frac{s_0+s_1+...+s_{N-1}}{N}
$$
这玩意要用在Fejer Kernel上
$$
F_N(x)=\frac{D_0(x)+...+D_{N-1}(x)}{N}
$$
Theorem 5.2：

若f圆上可积，则f的傅里叶级数，cesaro 可加收敛于f对于任何一个f上连续点

若f圆上连续，则傅里叶级数一致cesaro可加于f

这可以用来证明一些引理。如，f的傅里叶级数的cesaro和，就是一个f的一致估计量。

#### Abel和

Definition: A series of complex numbers $\sum_{k=0}^\infin c_k$ is said to be Abel summable to s if for every 0 $\le$ r <1 , the series:
$$
A(r)=\sum_{k=0}^\infin c_kr^k
$$
converges, and 
$$
\lim_{r\to1 }A(r)=s
$$
若原级数收敛于s，则abel和收敛于s。cesaro收敛，则abel收敛，反之不然。

#### 单位圆上的Poisson Kernel和Dirichlet's problem

$$
f(\theta)\sim\sum_{n=-\infin}^\infin a_ne^{in\theta}\\
A_r(f)(\theta)=\sum_{n=-\infin}^\infin r^{|n|}a_ne^{in\theta}
$$

如何和poisson kernel产生联系？
$$
A_r(f)(\theta)=(f*P_r)(\theta)
\\
A_r(f)(\theta)=\sum_{n=-\infin}^\infin r^{|n|}a_ne^{in\theta}
\\
=\sum_{n=-\infin}^\infin r^{|n|}(\frac{1}{2\pi}\int_{-\pi}^\pi f(\varphi)e^{-in\varphi}d\varphi)e^{in\theta}
\\
=\frac{1}{2\pi}\int_{-\pi}^\pi f(\varphi)(\sum_{n=-\infin}^\infin r^{|n|}e^{-in(\varphi-\theta)})d\varphi
$$
Poisson kernel 是good kernel

若f圆上可积，则f的傅里叶级数，abel 可加收敛于f对于任何一个f上连续点

若f圆上连续，则傅里叶级数一致abel可加于f。

用途：热传导：
$$
u(r,\theta)=(f*P_r)(\theta)
$$
这玩意有啥特点呢？

1、两方面可导且导数连续，$\Delta u=0$

2、越接近r=1，越到f($\theta$)如果$\theta$是f的连续点。如果f都连续，则limit是uniform
$$
\lim_{r\to1}u(r,\theta)=f(\theta)
$$
3、如果f连续，则$u(r,\theta)$是稳态热传导唯一解。

## Chapter 3. Convergence of Fourier Series

### 向量空间与内积

vertor space defined over $\mathbb{R}$ :
$$
1.两元素和仍然在向量空间中
\\2.可以有标量乘它
$$
复数一样。

内积是什么？

定义一个$(X,Y)\to\mathbb{R}$

首先：镜像：$(X,Y)=(Y,X)$

其次，有线性性：
$$
(\alpha X+\beta Y,Z)=\alpha(X,Z)+\beta(Y,Z)
$$
然后，正定，也就是$(X,X)\ge0$ 

最后，定义范数：
$$
||X||=(X,X)^\frac{1}{2}
$$
若范数为零可以得出X为零，则严格正定。

对于映射为复数的范数，那么就有变化。

1、Hermitian $(X,Y)=\overline{(Y,X)}$

2、第一个变量线性，第二个变量共轭：
$$
(\alpha X+\beta Y,Z)=\alpha(X,Z)+\beta(Y,Z)\\
(X,\alpha X+\beta Z)=\overline{\alpha}(X,Z)+\overline\beta(Y,Z)
$$
正定照旧。

接下来是正交的定义

if $(X,Y)=0$，则X正交于Y

正交有三个性质：

1、毕达哥拉斯定理

若X，Y正交，则：
$$
||X+Y||^2=||X||^2+||Y||^2
$$
2、柯西施瓦兹不等式

对任意向量空间中的X和Y
$$
|(X,Y)|\le||X||~||Y||
$$
3、三角不等式：

对任意向量空间中的X和Y
$$
||X+Y||\le||X||+||Y||
$$
Proof:

1、
$$
||X+Y||^2=(X+Y,X+Y)=(X,X)+2(X,Y)+(Y,Y)=||X||^2+||Y||^2
$$
2、
$$
若||Y||=0,则要证明的是(X,Y)=0对所有x都成立
\\
用0\le||X+tY||^2=||X||^2+2tRe(X,Y)
\\由t的任意性可以知道Re(X,Y)必须满足，同理Im(X,Y)=0
\\
若||Y||\neq0,设c=\frac{(X,Y)}{(Y,Y)}，那么X-cY正交于Y，为啥？
 \\(X-cY,Y)=(X,Y)-\frac{(X,Y)}{(Y,Y)}(Y,Y)=0\\
 ez......
 \\X=X-cY+cY
 用毕达哥拉斯定理
 \\||X||^2=||X-cY||^2+||cY||^2\ge|c|^2||Y||^2
$$
3、
$$
||X+Y||^2\le||X||^2+||Y||^2+2||X||||Y||=(||X||+||Y||)^2
$$
重要的例子:

two infinite-dimensional vector spaces.

先要有一个绝对可加的无限复数序列：
$$
\sum_{n\in\mathbb{Z}}|a_n|^2<\infin
$$
就可以定义inner product：
$$
(A,B)=\sum_{n\in\mathbb{Z}}a_n\overline{b_n}
$$
范数为：
$$
||A||=(A,A)^{1/2}=(\sum_{n\in\mathbb{Z}}|a_n|^2)^\frac{1}{2}
$$
我们称这样一个vector space为$\ell^2(\mathbb{Z})$

问题，是否是向量空间？其实是。

这些向量空间满足一些性质：

1、严格正定

2、完备，意味着极限点都在该空间内。

这两个性质加起来的向量空间，称其为Hilbert空间。如果有一个不满足就叫pre-Hilbert空间。

$\mathcal{R}$所有复值黎曼可积函数，定义于$[0,2\pi]$上。这是向量空间！那么内积呢？
$$
(f,g)=\frac{1}{2\pi}\int_0^{2\pi}f(\theta)\overline{g(\theta)}d\theta
$$
这玩意吧，柯西施瓦兹不等式和三角不等式都满足的，但Hilbert空间的俩性质，都不满足！一方面，范数为0只能说连续点都是0，测度为0的点仍然可以乱取值，笑死。另一方面，不完备，比如我当然可以做一个接近于0的柯西列，使其极限不属于$\mathcal{R}$

### 均方可积

既然我们用了这么一个残废的内积空间，那还得用啊。用它来证傅里叶级数的收敛：
$$
||f-S_N(f)||\to 0 ~as~N\to~\infin
$$
我们先定义一组函数：
$$
e_n(\theta)=e^{in\theta}\\
(e_n,e_m)=1 ~if~n=m
\\=0~if~n~\neq~m 
$$
那么正好咱的傅里叶系数可以写成：
$$
(f,e_n)=a_n
$$
这样定义的$a_n$有很好的性质，也就是：
$$
(f-\sum_{|n|\le N}a_ne_n)正交于\sum_{|n|\le N}b_ne_n
\\Proof一样简单的，线性展开嘛。
$$
那么开始毕达哥拉斯定理！
$$
f=f-\sum_{|n|\le N}a_ne_n+\sum_{|n|\le N}a_ne_n
\\
||f||^2=||f-S_N(f)||^2+\sum_{|n|\le N}|a_n|^2
$$
第二个结果叫做Best approximation:
$$
||f-S_N(f)||\le||f-\sum_{|n|\le N}c_ne_n||
$$
也就是说傅里叶系数估计出来的误差的范数是所有用三角级数估计出的最小的。这一个范数用积分衡量，很坑。一个毕达哥拉斯就搞定了。

然后可以正式考察均方收敛了。

首先由三角级数估计可得每个点都可以用同一个足够高阶的三角级数估计出来，其次，随后，平方，积分，得到$||f-P||<\epsilon$，用best approximation$||f-S_N(f)||<\epsilon$。

如果仅仅是可积，那么先用积分类似的连续函数去逼近他，也就是其差的绝对值的积分够小，那么其范数就好办，也够小，最后用三角级数逼近g，再用个啥三角不等式，最后best approximation，结束了。

最后有Parseval's identity:
$$
||f||^2=||f-S_N(f)||^2+\sum_{|n|\le N}|a_n|^2
$$
由于最优估计，$N\to\infin$时，误差项趋向于0，则
$$
||f||^2=\sum_{n=-\infin}^\infin|a_n|^2
$$
这把两个空间的范数联系在一起了。

如果改变这个基，也就是$e_n$，那么如果其正交，$a_n=(f,e_n)$，那么有Bessel inequality:
$$
||f||^2\ge\sum_{n=-\infin}^\infin|a_n|^2
$$
好玩的是，对于$\ell^2(\mathbb{Z})$中的序列，可能没有可积函数F的傅里叶系数和这序列完全相等。

#### Riemann-Lebesgue lemma

如果f可积，那么n趋于无穷，傅里叶系数趋于0.

这看着还挺直接的，虽然也不知道怎么证，但无所谓了，总之反正理解一下我说的意思就行。

更加广义的Parseval等式吧，可以处理两个序列
$$
F,G可积\\
F~\sim~\sum a_ne^{in\theta}~~~and~~~G~\sim~\sum b_ne^{in\theta}
\\
Then~~~\frac{1}{2\pi}=\int_0^{2\pi}F(\theta)\overline{G(\theta)}d\theta=\sum_{n=-\infin}^\infin a_n\overline{b_n}
\\Proof:
\\(F,G)=\frac{1}{4}[||F+G||^2-||F-G||^2+i(||F+iG||^2-||F-iG||^2)]
$$

### 点收敛

若方程点可微，那么点傅里叶级数收敛。有连续，但不收敛的例子。这个傅里叶级数的收敛是个大坑啊。

这结论，很有趣。你看首先，傅里叶级数需要对整个区间的f积分，但是，在某点处的收敛，牵涉到整个区间，却只和该点的可微条件有关系。这是为什么？

### 连续函数，有扩散傅里叶系数

Symmetric breaking!

E.g. $S_N$,Fejer, Dirichlet,Poisson kernels are symmetric.如果我们按照n的正负分开来，那么问题就不一样了。

回头再看吧。

## Chapter 4. Some Application of Fourier Series

### The isoperimetric inequality

### Weyl's equidistribution theorem

### A continuous but nowhere differentiable function

## Chapter 5. The Fourier Transform on $\mathbb{R}$

## Chapter 6. The Fourier Transform on $\mathbb{R}^d$

## Chapter 7. Finite Fourier Analysis



