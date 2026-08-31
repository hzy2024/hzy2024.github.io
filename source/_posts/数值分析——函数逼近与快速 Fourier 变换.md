---
title: 数值分析——函数逼近与快速 Fourier 变换
date: 2026-07-31 09:41:42
tags: [数值分析, 函数逼近与拟合, FFT]
categories: [学习笔记]
---


{% note default no-icon %}

本文内容参考了李庆扬《数值分析（第 6 版）》及徐利治《函数逼近的理论与方法》的内容，并对某些证明进行了方法上的改进. 阅读本文事先假定读者已经具备关于高等代数或泛函分析中有关线性空间、范数及内积、Gram-Schmidt 正交化方法、点列收敛等基本知识，如需参考可以直接查看 {% post_link "泛函分析课程笔记" 泛函分析课程笔记 %}或{% post_link "高等代数中的一些定义及结论" 高等代数中的一些定义及结论 %}，文中就不再赘述了. 同时，由于本文内容涉及理论相关已经很多，所以有些例子这里不再说明了，若无例子请读者自行查找或询问 AI.

{% endnote %}

首先我们来谈一谈本文的内容是在干什么. 当我们对真实的自然现象建立数学模型，其中的真实函数太复杂无法直接使用或者精确计算的成本太高时，如何找到一个简单、可计算的函数来近似描述这种关系，就是**函数逼近**所要解决的问题. 我们对于简单的多项式函数和三角函数的性质已经研究得比较深入了，而对非线性的函数，在某些假设条件下，一些定理保证了我们可以用足够多的多项式或三角函数作为基函数去逼近原始的复杂函数，并且收敛到原始函数.
当我们有很多观测的数据点而不知道真实的函数长什么样时，我们希望用自己的先验猜测这些数据点之间的函数关系，并且用我们所认为的函数来描述这种关系，这就是**曲线拟合**所要做的事. 注意，直接采用观测数据点插值的做法在很多时候并不是好的选择，因为这很容易导致过拟合.
函数逼近和曲线拟合通常被放在一起来说，但这两种想法所产生的动机正如上述所言是不同的. 但从应用场景和操作方法上来讲两者具有很多相似性. 机器学习本质上就是通过一个可训练复杂函数来逼近未知的真实函数，而目前还无法精确求解的 N-S 方程则可以通过离散化使用拟合的方式来得到近似解. 可以说，两种方法是人类认识世界的内在规律的重要手段.
在正文开始前，首先说明比较常用的符号释义：记所有次数不超过 $ n $ 的一元多项式构成的集合为 $ \mathscr{P}_n $，所有一元代数多项式构成的集合为 $ \mathscr{P} $ ，在区间 $ [a,b] $ 上的所有连续函数构成的集合记为 $ C[a,b] $，而记 $ C^n[a,b] $ 为区间 $ [a,b] $ 上有 $ n $ 阶连续导函数的实（复）值函数集合，$ L^2[a,b] $ 是区间 $ [a,b] $ 上平方可积函数的集合. 本文内容的安排主要按照多项式逼近与曲线最小二乘拟合、有理逼近及三角多项式逼近三块进行.

## 一般函数逼近与多项式逼近

### 引入

在数学分析的课程中，有关函数逼近的一个非常经典的结论是 **Weierstrass 一致逼近定理**：设 $ f(x)\in C[a,b] $，则对任何 $ \varepsilon>0 $，存在一个代数多项式 $ p(x)\in\mathscr{P} $，使得

$$
\begin{equation*}
    \max_{a\le x\le b}|f(x)-p(x)|<\varepsilon.
\end{equation*}
$$

该定理的构造性证明是由 Bernstein 给出的，并且在定理证明中给出了著名的 Bernstein 多项式. 具体证明从略，可以参考数学分析相关教材. 但由于这一多项式逼近的速度太慢，所以实际很少使用.
一般而言，我们更希望使用一组在 $ C[a,b] $ 上线性无关的函数集合 $ \lbrace\varphi_i(x) \rbrace_{i=0}^n $ 中元素的线性组合来逼近 $ f(x)\in C[a,b] $. 若记 $ \varphi(x)\in\Phi=\mathrm{span}\lbrace\varphi_0(x),\varphi_1(x),\cdots,\varphi_n(x) \rbrace\subset C[a,b] $，那么逼近函数 $ \varphi(x) $ 就有形式 $ \varphi(x)=a_0\varphi_0(x)+a_1\varphi_1(x)+\cdots+a_n\varphi_n(x) $. 现在我们就可以正式说明函数逼近的定义：**函数逼近问题**就是对任何 $ f(x)\in C[a,b] $（连续函数已经可以解决大部分问题了），在 $ \Phi $ 中找一个元素 $ \varphi^{\ast} $ 使得 $ f(x)-\varphi^{\ast}(x) $ 在某种度量下达到最小（$ \varphi^{\ast} $ 称为**最佳逼近函数**）. 如果我们把逼近函数限制在 $ n $ 次多项式空间里，并且 $ p^{\ast}(x)\in\mathscr{P}_n $ 使得在给定范数下的逼近误差达到最小，即

$$
\begin{equation*}
    ||f(x)-p^{\ast}(x)||=\min_{p\in\mathscr{P}_n}||f(x)-p(x)||,
\end{equation*}
$$

则称 $ p^{\ast}(x) $ 是 $ f(x) $ 在 $ [a,b] $ 上的**最佳逼近多项式**. 若取范数为无穷范数 $ ||\cdot||_{\infty} $，则称 $ p^{\ast}(x) $ 是 $ f(x) $ 在 $ [a,b] $ 上的**最佳一致逼近多项式**，若取 2-范数，那么 $ p^{\ast}(x) $ 是**最佳平方逼近多项式**.

### 正交多项式

这一小节的内容中，我们将会看到在某种内积下的正交多项式在函数逼近中会起到重要作用. 在高等代数或泛函分析的课程中，一般讲到的在欧氏空间或 $ C[a,b] $ 空间中的内积都是所谓的**标准内积**. 实际上我们可以给坐标分量或积分带上权重，得到带权内积. 在 $ \mathbb{R}^n $ 中，给定**正实数**序列 $ \lbrace\omega_i\rbrace_{i=1}^n $，对任意 $ \boldsymbol{x}=(x_1,\cdots,x_n)',\boldsymbol{y}=(y_1,\cdots,y_n)' $，定义带权内积

$$
\begin{equation*}
    (\boldsymbol{x},\boldsymbol{y})=\sum_{i=1}^n\omega_ix_iy_i,
\end{equation*}
$$

可以验证其确实满足内积的所有性质，并由此可以导出对应的范数. 复向量空间 $ \mathbb{C}^n $ 上也可以定义类似的内积，只要把 $ \boldsymbol{y} $ 进行共轭运算即可. 在 $ C[a,b] $ 上我们可以定义类似的内积，为此首先给出权函数的定义：

{% note info no-icon %}

设 $ [a,b] $ 是有限或无限区间，在 $ [a,b] $ 上的非负函数 $ \rho(x) $ 满足条件：

(1) $ \displaystyle{\int_{a}^bx^k\rho(x)\mathrm{d}x} $ 存在且为有限值 $ (k=0,1,\cdots) $；
(2) 对 $ [a,b] $ 上的非负连续函数 $ g(x) $，如果 $ \displaystyle{\int_{a}^bg(x)\rho(x)\mathrm{d}x=0} $，则 $ g(x)\equiv0 $.

则称 $ \rho(x) $ 为 $ [a,b] $ 上的一个**权函数**.

{% endnote %}

有了权函数，对 $ f(x),g(x)\in C[a,b] $，我们可以定义带权内积

$$
\begin{equation*}
    (f(x),g(x))=\int_{a}^b\rho(x)f(x)g(x)\mathrm{d}x.
\end{equation*}
$$

同样容易验证其满足内积的定义，由此可以导出对应的范数. 当然 $ \rho(x)\equiv1 $ 的情形还是最常用的，复值函数相应的定义仍只需要把第二变元作共轭即可. 现在我们就可以把正交、正交函数族、标准正交函数族等概念推广至带权情形：

{% note info no-icon %}

设 $ f(x),g(x)\in C[a,b],\rho(x) $ 是 $ [a,b] $ 上的权函数. 若

$$
\begin{equation*}
    (f(x),g(x))=\int_{a}^b\rho(x)f(x)g(x)\mathrm{d}x=0,
\end{equation*}
$$

则称 $ f(x) $ 与 $ g(x) $ 在 $ [a,b] $ 上带权 $ \rho(x) $ **正交**. 若函数族 $ \varphi_0(x),\varphi_1(x),\cdots,\varphi_n(x),\cdots $ 满足关系

$$
\begin{equation*}
    (\varphi_j,\varphi_k)=\int_{a}^b\rho(x)\varphi_j(x)\varphi_k(x)\mathrm{d}x=\begin{cases}
        0, & j\ne k,\\
        A_k>0, & j=k.
    \end{cases}
\end{equation*}
$$

则称 $ \lbrace\varphi_k(x)\rbrace_{k=0}^{\infty} $ 是 $ [a,b] $ 上的带权 $ \rho(x) $ 的**正交函数族**；若 $ A_k\equiv1 $，则称之为**标准正交函数族**. 特别，若 $ \varphi_n(x) $ 是 $ [a,b] $ 上首项系数 $ a_n\ne0 $ 的 $ n $ 次多项式，则称 $ \varphi_n(x) $ 为 $ [a,b] $ 上带权 $ \rho(x) $ 的 **$ n $ 次正交多项式**.

{% endnote %}

在给定区间 $ [a,b] $ 和权函数 $ \rho(x) $ 后，只要使用幂函数族 $ \lbrace 1,x,\cdots,x^n,\cdots\rbrace $ 通过 Gram-Schmidt 正交化方法就能构造正交多项式序列 $ \lbrace \varphi_n(x)\rbrace_{n=0}^{\infty} $：

$$
\begin{equation*}
    \begin{cases}
        \varphi_0(x)=1,\\
        \varphi_n(x)=x^n-\displaystyle{\sum_{j=0}^{n-1}}\dfrac{(x^n,\varphi_j(x))}{(\varphi_j(x),\varphi_j(x))}\varphi_j(x),~~n=1,2,\cdots
    \end{cases}
\end{equation*}
$$

仍然可以证明 $ \lbrace \varphi_n(x)\rbrace_{n=0}^{\infty} $ 中的任意有限项线性无关，并且由于构造特点可知 $ \varphi_n(x) $ **都是首一多项式，不妨称为首一正交多项式（族）**. 事实上这样构造出来的正交多项式就是 $ \mathscr{P} $ 的另一组基，因此对于任意 $ p(x)\in\mathscr{P}\_n $，取前 $ n+1 $ 个正交化得到的正交多项式就可以将其线性表示，并且由上述线性表示很容易知道 $ \varphi_n(x) $ 与任意 $ p(x)\in\mathscr{P}\_{n-1} $ 正交.
除了上面提到的两个性质，我们还有关于正交多项式的两个重要性质：

{% note primary no-icon %}

设 $ \lbrace\varphi_n(x) \rbrace_{n=0}^{\infty} $ 是 $ [a,b] $ 上带权 $ \rho(x) $ 的正交多项式，则

1. 对 $ n\ge0 $ 成立递推关系

$$
\begin{equation*}
    \varphi_{n+1}(x)=(x-\alpha_n)\varphi_n(x)-\beta_n\varphi_{n-1}(x),~~n=0,1,2,\cdots,
\end{equation*}
$$

其中 $ \varphi_0(x)=1,\varphi_{-1}(x)=0 $，

$$
\begin{equation*}
    \alpha_n=\frac{(x\varphi_n(x),\varphi_n(x))}{(\varphi_n(x),\varphi_n(x))},~~\beta_n=\frac{(\varphi_n(x),\varphi_n(x))}{(\varphi_{n-1}(x),\varphi_{n-1}(x))},~~n=0,1,2,\cdots
\end{equation*}
$$

2. $ \varphi_n(x) $ 在区间 $ (a,b) $ 有 $ n $ 个不同的零点.

{% endnote %}

**证明**

1. 我们最主要需要说明为什么新的项只与旧的最后两个多项式有关. 注意 $ \deg x\varphi_n(x)=n+1 $，因此其可以写成 $ \varphi_0(x),\cdots,\varphi_{n+1}(x) $ 的线性组合的形式

$$
\begin{equation*}
    x\varphi_n(x)=a_{n+1}\varphi_{n+1}(x)+\sum_{i=0}^na_i\varphi_i(x).
\end{equation*}
$$

而注意**从幂函数基正交化得到的多项式都是首一的**，因此 $ a_{n+1}=1 $. 上式可以改写为

$$
\begin{equation*}
    \varphi_{n+1}(x)=(x-a_n)\varphi_n(x)+\sum_{i=0}^{n-1}\varphi_i(x).
\end{equation*}
$$

对于 $ k\le n-2 $，注意多项式的次数 $ \deg x\varphi_k(x)\le n-1 $，因此其和 $ \varphi_n(x) $ 正交，从两者内积的表达式可知有 $ (x\varphi_n(x),\varphi_k(x))=(\varphi_n(x),x\varphi_k(x))=0 $，在上述展开式两侧和 $ \varphi_k(x) $ 作内积就知道 $ a_k=0,k\le n-2 $. 因此最终结果仅与 $ \varphi_n(x),\varphi_{n-1}(x) $ 有关. 接下来两侧和 $ \varphi_n(x) $ 作内积即得 $ \alpha_n $，和 $ \varphi_{n-1}(x) $ 作内积即得 $ \beta_n $. 注意此时得到的式子不是我们想要的，实际上再次利用递推公式有

$$
\begin{align*}
    \beta_n&=\frac{(x\varphi_n(x),\varphi_{n-1}(x))}{(\varphi_{n-1}(x),\varphi_{n-1}(x))}=\frac{(\varphi_n(x),x\varphi_{n-1}(x))}{(\varphi_{n-1}(x),\varphi_{n-1}(x))}\\
    &=\frac{(\varphi_n(x),\varphi_n(x)+\alpha_{n-1}\varphi_{n-1}(x)+\beta_{n-1}\varphi_{n-2}(x))}{(\varphi_{n-1}(x),\varphi_{n-1}(x))}\\&=\frac{(\varphi_n(x),\varphi_{n}(x))}{(\varphi_{n-1}(x),\varphi_{n-1}(x))}.
\end{align*}
$$

由此得证. 后面将会看到这一推导过程中的恒等式在 Legendre 多项式性质的推导中会起到重要作用.

2. 设 $ \varphi_n(x) $ 在 $ (a,b) $ 内的零点都是偶数重的, 则 $ \varphi_n(x) $ 在 $ (a,b) $ 上不变号，这与 $ (\varphi_n,\varphi_0)=0 $ 矛盾. 现假设 $ x_i(i=1,2,\cdots,l) $ 是其在 $ (a,b) $ 内从小到大排列的奇数重零点，则 $ \varphi_n(x) $ 在 $ x_i $ 处变号，若令

$$
\begin{equation*}
    q(x)=(x-x_1)(x-x_2)\cdots(x-x_l),
\end{equation*}
$$

那么 $ \varphi_n(x)q(x) $ 不变号，从而 $ (\varphi_n(x),q(x))\ne0 $. 若 $ l<n $，则由 $ \varphi_n $ 与 $ q $ 正交知 $ (\varphi_n(x),q(x))=0 $，矛盾. 故 $ l\ge n $. 又 $ \varphi_n(x) $ 只有 $ n $ 个零点，因此只能是 $ l=n $，即 $ n $ 个零点是互异的且都落在区间 $ (a,b) $ 内.

{% note warning no-icon %}

值得注意的是，虽然我们这里的递推式是针对首一正交多项式而言，但是对一般的正交多项式，**由于正交性仍然保持不变，所以我们可以假设一般递推式形式为**

$$
\begin{equation*}
    P_{n+1}(x)=c_{n}(x-a_n)P_n(x)+b_nP_{n-1}(x),
\end{equation*}
$$

**然后通过具体的内积形式来确定每一个系数**. 这一方法将在后续的 Legendre 多项式递推公式的建立时用到.

{% endnote %}

接下来我们给出一些常用的正交多项式和它们的性质.

#### Legendre 多项式

若取区间为 $ [-1,1] $，权函数为 $ \rho(x)=1 $，由 $ \lbrace1,x,\cdots,x^n,\cdots \rbrace $ 正交化得到的多项式称为 **Legendre 多项式**，并使用 $ P_0(x),P_1(x),\cdots,P_n(x),\cdots $ 来表示. Rodrigul 给出了 Legendre 多项式的简单表达式

$$
\begin{equation*}
    P_0(x)=1,~~P_n(x)=\frac{1}{2^n\cdot n!}\frac{\mathrm{d}^n}{\mathrm{d}x^n}(x^2-1)^n,~~n=1,2,\cdots.
\end{equation*}
$$

注意这里 $ P_n(x) $ 并不是首一多项式，注意其形式是由一个 $ 2n $ 次多项式求 $ n $ 次导数得到的，可知 $ P_n(x) $ 的最高次项系数 $ a_n=\dfrac{(2n)!}{2^n(n!)^2} $，而对应的首一多项式则为

$$
\begin{equation*}
    \widetilde{P}_n(x)=\frac{n!}{(2n)!}\frac{\mathrm{d}^n}{\mathrm{d}x^n}(x^2-1)^n.
\end{equation*}
$$

Legendre 多项式有如下几个重要性质：

{% note primary no-icon %}

1. （正交性）

$$
\begin{equation*}
    \int_{-1}^1P_n(x)P_m(x)\mathrm{d}x=\begin{cases}
        0, & m\ne n,\\
        \dfrac{2}{2n+1}, & m=n.
    \end{cases}
\end{equation*}
$$

2. （奇偶性）

$$
\begin{equation*}
    P_n(-x)=(-1)^nP_n(x).
\end{equation*}
$$

3. （递推关系）

$$
\begin{equation*}
    (n+1)P_{n+1}(x)=(2n+1)xP_n(x)-nP_{n-1}(x),~~,n=1,2,\cdots.
\end{equation*}
$$

4. $ P_n(x) $ 在区间 $ (-1,1) $ 内有 $ n $ 个不同的实零点.

{% endnote %}

**证明** 我们来逐一说明. 令 $ \varphi(x)=(x^2-1)^n=(x+1)^n(x-1)^n $，利用 Newton-Lebniz 公式易知 $ \varphi^{(k)}(\pm1)=0(k=0,1,\cdots,n-1) $. 若 $ Q(x)\in C^n[-1,1] $，那么由分部积分法知道

$$
\begin{align*}
    &\int_{-1}^1P_n(x)Q(x)\mathrm{d}x=\frac{1}{2^n\cdot n!}\int_{-1}^1Q(x)\varphi^{(n)}(x)\mathrm{d}x\\
    &=\frac{-1}{2^n\cdot n!}\int_{-1}^1Q'(x)\varphi^{(n-1)}(x)\mathrm{d}x=\cdots=\frac{(-1)^n}{2^n\cdot n!}\int_{-1}^1Q^{(n)}(x)\varphi(x)\mathrm{d}x.
\end{align*}
$$

若 $ Q(x) $ 是次数小于 $ n $ 的多项式，马上得到 $ Q^{(n)}(x)=0 $，从而结合交换性得到 $ (P_n,P_m)=0,m\ne n $. 若 $ Q(x)=P_n(x) $，那么利用先前 $ P_n(x) $ 的首项系数，易知求导 $ n $ 次后有 $ Q^{(n)}(x)=P_n^{(n)}(x)=\dfrac{(2n)!}{2^n\cdot n!} $. 从而

$$
\begin{align*}
    \int_{-1}^1P_n^2(x)\mathrm{d}x&=\frac{(-1)^n(2n)!}{2^{2n}(n!)^2}\int_{-1}^1(x^2-1)^n\mathrm{d}x\\
    &=\frac{(2n)!}{2^{2n-1}(n!)^2}\int_{0}^{\frac{\pi}{2}}\cos^{2n+1}x\mathrm{d}x\\
    &=\frac{2\cdot(2n)!}{[(2n)!!]^2}\cdot\frac{(2n)!!}{(2n+1)!!}=\frac{2}{2n+1}.
\end{align*}
$$

而利用先前关于首一正交多项式的性质，我们有

$$
\begin{align*}
    (xP_n(x),P_{n-1}(x))&=\frac{(2n)!}{2^n(n!)^2}\frac{(2n-2)!}{2^{n-1}((n-1)!)^2}(x\widetilde{P}_n(x),\widetilde{P}_{n-1}(x))\\
    &=\frac{(2n)!}{2^n(n!)^2}\frac{(2n-2)!}{2^{n-1}((n-1)!)^2}(\widetilde{P}_n(x),\widetilde{P}_{n}(x))\\
    &=\frac{(2n-2)!}{2^{n-1}((n-1)!)^2}\frac{2^n(n!)^2}{(2n)!}\cdot\frac{2}{2n+1}=\frac{2n}{4n^2-1}.
\end{align*}
$$

接下来讨论奇偶性. 由于 $ \varphi(x) $ 为偶函数，因此求奇数阶导得到的函数为奇函数，求偶数阶导得到的函数为偶函数，因此当 $ n $ 为偶数时 $ P_n $ 为偶函数，否则为奇函数，因此奇偶性式子成立. 从奇偶性出发，我们马上知道 $ x\widetilde{P}_n^2(x) $ 为奇函数，因此在构建递推式时，利用首一正交多项式的递推式表达式我们马上知道 $ \alpha_n=0 $，故可设递推表达式为

$$
\begin{equation*}
    P_{n+1}(x)=c_1xP_n(x)+c_2P_{n-1}(x),
\end{equation*}  
$$

结合前面得到的几个计算结果，对上式两侧分别和 $ P_n(x),P_{n+1}(x) $ 做内积很容易可以确定

$$
\begin{equation*}
    c_1=\frac{2n+1}{n+1},c_2=-\frac{n}{n+1}.
\end{equation*}
$$

从而有递推式

$$
\begin{equation*}
    (n+1)P_{n+1}(x)=(2n+1)xP_n(x)-nP_{n-1}(x),~~n=1,2,\cdots.
\end{equation*}
$$

最后的零点性质直接由首一正交多项式的性质即得.
通过递推公式，可以确定前几个 Legendre 多项式为

$$
\begin{equation*}
    P_1(x)=x,P_2(x)=\frac{1}{2}(3x^2-1),P_3(x)=\frac{1}{2}(5x^3-3x),P_4(x)=\frac{1}{8}(35x^4-30x^2+3),\cdots
\end{equation*}
$$

#### Chebyshev 多项式

取区间为 $ [-1,1] $，权函数为 $ \dfrac{1}{\sqrt{1-x^2}} $，由幂函数基正交化得到的正交多项式就是 **Chebyshev 多项式**，可以表示为

$$
\begin{equation*}
    T_n(x)=\cos(n\arccos x),~~|x|\le 1.
\end{equation*}
$$

令 $ x=\cos\theta,\theta\in[0,\pi] $，则 $ T_n(x)=\cos n\theta $. **注意，该定义下多项式仍然不是首一的**. 由该形式可以导出很多 Chebyshev 多项式的重要性质：

{% note primary no-icon %}

1. （递推关系）

$$
\begin{equation*}
    \begin{cases}
        T_0(x)=1,~~T_1(x)=x,\\
        T_{n+1}(x)=2xT_n(x)-T_{n-1}(x),~~n=1,2,\cdots
    \end{cases}
\end{equation*}
$$

进而有 $ T_{2k}(x) $ 只含 $ x $ 的偶次幂， $ T_{2k+1} $ 只含 $ x $ 的奇次幂，且 $ T_n(x) $ 的首项 $ x^n $ 的系数为 $ 2^{n-1} $.

2. （正交性） $ \lbrace T_k(x) \rbrace $ 在区间 $ [-1,1] $ 上带权 $ \rho(x)=1/\sqrt{1-x^2} $ 正交，且

$$
\begin{equation*}
    \int_{-1}^{1}\frac{T_n(x)T_m(x)}{\sqrt{1-x^2}}=\begin{cases}
        0, & n\ne m,\\
        \dfrac{\pi}{2}, & n=m\ne 0,\\
        \pi, & n=m=0.
    \end{cases}
\end{equation*}
$$

3. $ T_n(x) $ 在区间 $ (-1,1) $ 上有 $ n $ 个不同的实零点

$$
\begin{equation*}
    x_k=\cos\frac{2k-1}{2n}\pi,~~k=1,2,\cdots,n.
\end{equation*}
$$

4. 令 $ \widetilde{T}_0(x)=1,\widetilde{T}_n(x)=\dfrac{1}{2^{n-1}}T_n(x) $ 为首一的 Chebyshev 多项式， $ \widetilde{\mathscr{P}}_n $ 为所有次数小于等于 $ n $ 的首一多项式集合，那么有

$$
\begin{equation*}
    \max_{-1\le x\le 1}|\widetilde{T}_n(x)|\le\max_{-1\le x\le 1}|p(x)|, ~~ \forall p(x)\in\widetilde{\mathscr{P}}_n,
\end{equation*}
$$

且 $ \max\limits_{-1\le x\le 1}|\widetilde{T}_n(x)|=\dfrac{1}{2^{n-1}} $.

{% endnote %}

性质 1-3 的证明非常简单，这里不展开说明了. 递推关系只要由三角恒等式 $ \cos(n+1)\theta=2\cos\theta\cos n\theta-\cos(n-1)\theta $ 就可得到. 我们使用反证法来证明最后一点. 易知 $ \widetilde{T}\_n(x) $ 在区间 $ [-1,1] $ 内（包括端点）有 $ n+1 $ 个极值点 $ x_k=\cos\dfrac{k\pi}{n},k=0,1,\cdots,n $，对应的极值为 $ (-1)^k\cdot2^{1-n} $. 假设另有一首一多项式 $ q(x) $ 满足 $ \max\limits_{-1\le x\le 1}|q(x)|<2^{1-n}=\max\limits_{-1\le x\le 1}|\widetilde{T}\_n(x)| $，令 $ R_n(x)=\widetilde{T}\_n(x)-q(x) $，则 $ \deg R_n(x)\le n-1 $. 在点 $ x_k $ 上，由于 $ |q(x)|<\dfrac{1}{2^{n-1}} $，因此通过简单的大小关系比较可知 $ \widetilde{T}\_n(x_k) $ 和 $ R_n(x_k) $ 同号且 $ R_n(x_k) $ 均不为零. 因而由连续函数的介值定理可知，在区间 $ (x_k,x_{k+1})(k=0,\cdots,n-1) $ 上 $ R_n(x) $ 各至少存在一个零点，即 $ R_n(x) $ 至少有 $ n $ 个互异零点，这与 $ R_n(x) $ 的次数为 $ n-1 $ 相矛盾，定理得证.
我们现在构建的 Chebyshev 多项式是在区间 $ [-1,1] $ 上的，如果需要构造一般区间 $ [a,b] $ 上的 Chebyshev 多项式，只需要使用变换

$$
\begin{equation*}
    x=\frac{1}{2}[(b-a)t+a+b]\Leftrightarrow t=\frac{1}{b-a}(2x-a-b)
\end{equation*}
$$

即可将 $ x\in[a,b] $ 和 $ t\in[-1,1] $ 一一对应起来，通过这样的变换不会改变 Chebyshev 多项式上述的性质，仅是在对应数值上会产生等比例的放缩.
Chebyshev 多项式在区间 $ [-1,1] $ 上的零点和极值点这两组点集称为 **Chebyshev 点**，其为单位圆周上等距分布点的横坐标，在接近区间端点处比较密集. 利用 Chebyshev 点做插值可使得插值区间的最大误差最小化. 为了说明这一点，设 $ f(x)\in C^n[-1,1] $，插值节点为 $ x_0,x_1,\cdots,x_n\in[-1,1] $， $ L_n(x) $ 为相应的 Lagrange 插值多项式，那么插值余项

$$
\begin{equation*}
    R_n(x)=f(x)-L_n(x)=\frac{f^{(n+1)}(\xi)}{(n+1)!}\omega_n(x)=\frac{f^{(n+1)}(\xi)}{(n+1)!}(x-x_0)(x-x_1)\cdots(x-x_n),
\end{equation*}
$$

若记 $ M_{n+1}=||f^{(n+1)}(x)||\_{\infty}=\max\limits_{-1\le x\le 1}|f^{(n+1)}(x)| $，并选择插值节点为 $ T_{n+1}(x) $ 的零点，那么由上面的第 4 条性质我们知道

$$
\begin{equation*}
    \max_{-1\le x\le 1}|f(x)-L_n(x)|\le\frac{1}{2^n(n+1)!}||f^{(n+1)}(x)||_{\infty}.
\end{equation*}
$$

当然对于一般区间 $ [a,b] $ ，我们仍然可以使用先前提到的变换得到变换后的插值节点为

$$
\begin{equation*}
    x_k=\frac{b-a}{2}\cos\frac{2k+1}{2(n+1)}\pi+\frac{a+b}{2},~~k=0,1,\cdots,n.
\end{equation*}
$$

对于高次插值多项式，其不一定能够逼近真实函数 $ f(x) $，误差可能会随着次数的升高而不断增大，这种现象称为**龙格（Runge）现象**. 龙格给出的例子是使用等距节点进行插值得到的，但若用 Chebyshev 多项式的零点进行插值，却能避免龙格现象，保证整个区间上的收敛性.

#### 其他常用正交多项式

如果我们选取不同的区间 $ [a,b] $ 和权函数 $ \rho(x) $，那么得到的正交多项式也不同. 这里我们仅列举一些常用的正交多项式和它们的性质.

**第二类 Chebyshev 多项式**

若选取区间 $ [-1,1] $ 和权函数 $ \rho(x)=\sqrt{1-x^2} $，得到的正交多项式称为第二类 Chebyshev 多项式，表达式为

$$
\begin{equation*}
    U_n(x)=\frac{\sin[(n+1)\arccos x]}{\sqrt{1-x^2}},
\end{equation*}
$$

其带权内积为

$$
\begin{equation*}
    (U_n(x),U_m(x))=\int_{0}^\pi\sin(n+1)\theta\sin(m+1)\theta\mathrm{d}\theta=\begin{cases}
        0, & m\ne n,\\
        \dfrac{\pi}{2}, & m=n,
    \end{cases}
\end{equation*}
$$

递推关系式为

$$
\begin{equation*}
    \begin{cases}
        U_0(x)=1,~~U_1(x)=2x,\\
        U_{n+1}(x)=2xU_n(x)-U_{n-1}(x),~~n=1,2,\cdots.
    \end{cases}
\end{equation*}
$$

**Laguerre 多项式**

若选取区间 $ [0,+\infty) $ 和权函数 $ \rho(x)=\mathrm{e}^{-x} $，得到的正交多项式称为 Laguerre 多项式，表达式为

$$
\begin{equation*}
    L_n(x)=\mathrm{e}^{x}\frac{\mathrm{d}^n}{\mathrm{d}x^n}(x^n\mathrm{e}^{-x}),
\end{equation*}
$$

其带权内积为

$$
\begin{equation*}
    (L_n(x),L_m(x))=\int_{0}^{+\infty}\mathrm{e}^{-x}L_n(x)L_m(x)\mathrm{d}x=\begin{cases}
        0, & m\ne n,\\
        (n!)^2, & m=n,
    \end{cases}
\end{equation*}
$$

递推关系式为

$$
\begin{equation*}
    \begin{cases}
        L_0(x)=1,~~L_1(x)=1-x,\\
        L_{n+1}(x)=(1+2n-x)L_n(x)-n^2L_{n-1}(x),~~n=1,2,\cdots.
    \end{cases}
\end{equation*}
$$

**Hermite 多项式**

若选取区间 $ (-\infty,+\infty) $ 和权函数 $ \rho(x)=\mathrm{e}^{-x^2} $，得到的正交多项式称为 Hermite 多项式，表达式为

$$
\begin{equation*}
    H_n(x)=(-1)^n\mathrm{e}^{x^2}\frac{\mathrm{d}^n}{\mathrm{d}x^n}(\mathrm{e}^{-x^2}),
\end{equation*}
$$

其带权内积为

$$
\begin{equation*}
    (H_n(x),H_m(x))=\int_{-\infty}^{+\infty}\mathrm{e}^{-x^2}L_n(x)L_m(x)\mathrm{d}x=\begin{cases}
        0, & m\ne n,\\
        2^n\cdot n!\sqrt{\pi}, & m=n,
    \end{cases}
\end{equation*}
$$

递推关系式为

$$
\begin{equation*}
    \begin{cases}
        H_0(x)=1,~~H_1(x)=2x,\\
        H_{n+1}(x)=2xH_n(x)-2nH_{n-1}(x),~~n=1,2,\cdots.
    \end{cases}
\end{equation*}
$$

### 最佳平方逼近

我们在关于 Chebyshev 多项式的讨论中得到了其一致逼近任一连续函数的最佳性，接下来我们将目光转向本章提及的 2-范数意义下的最佳逼近问题——最佳平方逼近. 设 $ f(x)\in C[a,b],\varphi=\mathrm{span}\lbrace\varphi_0(x),\varphi_1(x),\cdots,\varphi_n(x) \rbrace\subset C[a,b] $，若存在 $ S^\ast(x)\in\varphi $ 使得

$$
\begin{equation*}
    ||f(x)-S^\ast(x)||_2^2=\min_{S(x)\in\varphi}||f(x)-S(x)||_2^2=\min_{S(x)\in\varphi}\int_{a}^b\rho(x)[f(x)-S(x)]^2\mathrm{d}x,
\end{equation*}
$$

则称 $ S^\ast(x) $ 是 $ f(x) $ 在子空间 $ \varphi $ 中的**最佳平方逼近函数**. 上述问题等价于求多元函数

$$
\begin{equation*}
    I(a_0,a_1,\cdots,a_n)=\int_{a}^b\rho(x)\left[\sum_{j=0}^na_j\varphi_j(x)-f(x)\right]^2\mathrm{d}x
\end{equation*}
$$

的最小值. 利用多元函数求极值的必要条件，求解 $ I $ 的驻点可得到驻点方程

$$
\begin{equation*}
    \sum_{j=0}^n(\varphi_k(x),\varphi_j(x))a_j=(f(x),\varphi_k(x)),~~k=0,1,\cdots,n.
\end{equation*}
$$

具体过程比较简单，这里略去. 上述 $ n+1 $ 个关于 $ a_0,a_1,\cdots,a_n $ 的线性方程组称为**法方程**，我们还可以将其写成矩阵形式

$$
\begin{equation*}
    \begin{pmatrix}
        (\varphi_0,\varphi_0) & (\varphi_0,\varphi_1) & \cdots &(\varphi_0,\varphi_n)\\
        (\varphi_1,\varphi_0) & (\varphi_1,\varphi_1) & \cdots & (\varphi_1,\varphi_n)\\
        \vdots & \vdots & & \vdots\\
        (\varphi_n,\varphi_1) & (\varphi_n,\varphi_2) & \cdots & (\varphi_n,\varphi_n)
    \end{pmatrix}\begin{pmatrix}
        a_0 \\
        a_1 \\
        \vdots\\
        a_n
    \end{pmatrix}=\begin{pmatrix}
        (f,\varphi_0) \\
        (f,\varphi_1) \\
        \vdots\\
        (f,\varphi_n)
    \end{pmatrix},
\end{equation*}
$$

简记为 $ \boldsymbol{Ga}=\boldsymbol{d} $，其中 $ \boldsymbol{G}=\boldsymbol{G}(\varphi_0,\varphi_1,\cdots,\varphi_n) $ 为 **Gram 矩阵**，所有内积都是给定权函数 $ \rho(x) $ 的带权内积. $ \boldsymbol{G} $ 非奇异等价于 $ \varphi_0,\varphi_1,\cdots,\varphi_n $ 线性无关，这一证明是高等代数有关内积空间的熟知结论，不再赘述. 因此我们只要选取好基函数，就能唯一确定上述线性方程组的解 $ \boldsymbol{a}=\boldsymbol{a}^\ast=(a_0^\ast,a_1^\ast,\cdots,a_n^\ast) $. 由于 $ I $ 关于其变元 $ \boldsymbol{a} $ 是严格凸函数，实际上可以注意到

$$
\begin{equation*}
    I(\boldsymbol{a})=\boldsymbol{a}'\boldsymbol{Ga}-2\boldsymbol{d}'\boldsymbol{a}+||f||_2^2,
\end{equation*}
$$

而 $ \boldsymbol{G}\succ\boldsymbol{O} $，因此其为严格凸二次函数，其极小值点就是全局的唯一最小值点，所以由法方程确定的 $ S^\ast(x) $ 就是最佳平方逼近函数. 而利用法方程我们还可以直接得到最佳平方逼近的误差为

$$
\begin{equation*}
    ||\delta(x)||_2^2=||f(x)-S^\ast(x)||_2^2=I(\boldsymbol{a}^\ast)=||f(x)||_2^2-\boldsymbol{d}'\boldsymbol{a}^{\ast}=||f(x)||_2^2-\sum_{k=0}^na_k^\ast(f(x),\varphi_k(x)).
\end{equation*}
$$

若我们取 $ \varphi_k(x)=x^k,\rho(x)\equiv1,f(x)\in C[0,1] $，那么就是要在子空间 $ \mathscr{P}_n $ 中求 $ n $ 次最佳平方逼近多项式，此时得到的 Gram 矩阵称为 **Hilbert 矩阵**，不难计算得到

$$
\begin{equation*}
    \boldsymbol{H}=\begin{pmatrix}
        1 & 1/2 & \cdots & 1/(n+1)\\
        1/2 & 1/3 & \cdots & 1/(n+2)\\
        \vdots & \vdots & & \vdots\\
        1/(n+1) & 1/(n+2) & \cdots & 1/(2n+1)
    \end{pmatrix}.
\end{equation*}
$$

根据具体的函数 $ f $ 的形式就可以解出相应的最佳平方逼近多项式. **但是，Hilbert 矩阵在 $ n $ 较大时是高度病态的，因此直接求解法方程相当困难，还会带来数值上的严重不稳定**. 为此，我们采用正交多项式作为基. 这样的好处是得到的 Gram 矩阵是对角阵且对角元大于 $ 0 $，此时方程组 $ \boldsymbol{Ga}=\boldsymbol{d} $ 的解就是

$$
\begin{equation*}
    a_k^\ast=\frac{(f(x),\varphi_k(x))}{(\varphi_k(x),\varphi_k(x))},~~k=0,1,\cdots,n,
\end{equation*}
$$

平方逼近误差可写作

$$
\begin{equation*}
    ||\delta_n(x)||_2=\left(||f(x)||_2^2-\sum_{k=0}^n\left[\frac{(f(x),\varphi_k(x))}{||\varphi_k(x)||_2}\right]^2\right)^{\frac{1}{2}}.
\end{equation*}
$$

从 $ ||\delta_n(x)||_2\ge0 $ 还可推得 **Bessel 不等式**

$$
\begin{equation*}
    \sum_{k=0}^n(a_k^\ast||\varphi_k(x)||_2)^2\le||f(x)||_2^2.
\end{equation*}
$$

若现在给定的基函数族中的基有可列个，按上述最优系数 $ a_k^\ast $ 将 $ f(x) $ 展开，则称级数

$$
\begin{equation*}
    \sum_{k=0}^{\infty}a_k^\ast\varphi_k(x)
\end{equation*}
$$

为 $ f $ 的**广义 Fourier 级数**，并称 $ a_k^\ast $ 为**广义 Fourier 系数**.
**如果基函数族 $ \lbrace\varphi_k(x)\rbrace_{k=0}^{\infty} $ 是由 $ \lbrace 1,x,\cdots,x^n,\cdots \rbrace $ 正交化得到的正交多项式族**，实际上由实变函数/泛函分析的知识有 $ \mathscr{P}[a,b] $ 在 $ L^2[a,b] $ 中稠密，自然也在 $ C[a,b] $ 中稠密. 如果对每个 $ \varphi_k(x) $ 标准化为 $ \tilde{\varphi}_k(x) $，即乘上非零常数使其范数为 $ 1 $，结合 $ \mathscr{P}=\mathrm{span}\lbrace\varphi_0(x),\varphi_1(x),\cdots,\varphi_n(x),\cdots\rbrace=\mathrm{span}\lbrace\tilde{\varphi}_0(x),\tilde{\varphi}_1(x),\cdots,\tilde{\varphi}_n(x),\cdots\rbrace $ 可知在 $ \mathscr{P}[a,b] $ 上的标准正交基 $ \lbrace\tilde{\varphi}_k(x)\rbrace\_{k=0}^{\infty} $ 是完全标准正交基，因此在 $ C[a,b] $ 上成立 Parseval 等式

$$
\begin{equation*}
    ||f(x)||_2^2=\sum_{k=0}^{\infty}b_k^2=\sum_{k=0}^{\infty}(a_k^*||\varphi_k(x)||_2)^2.
\end{equation*}
$$

其中 $ b_k $ 是针对标准正交基得到的最佳系数. 由此还可推得定理

$$
\begin{equation*}
    \lim_{n\to\infty}||f(x)-S_n^\ast(x)||_2=0.
\end{equation*}
$$

上述说明中提及的定理可以参考{% post_link "泛函分析课程笔记" 泛函分析课程笔记 %}中有关希尔伯特空间的内容. 因此，如果我们选取 Legendre 多项式作为基函数，那么就能在 $ [-1,1] $ 上逼近真实函数 $ f(x)\in C[a,b] $，相应的展开式和逼近误差这里就不再说明. 如果 $ f\in C^2[-1,1] $，那么还有结论：任给 $ \varepsilon>0 $，当 $ n $ 充分大时，对任意 $ x\in[-1,1] $ 都有

$$
\begin{equation*}
    |f(x)-S_n^*(x)|\le\dfrac{\varepsilon}{\sqrt{n}}.
\end{equation*}
$$

此外，若将 Legendre 多项式首项归一化（即先前提到的 $ \widetilde{P}_n(x) $），那么其在 $ [-1,1] $ 上与零的平方逼近误差最小. 证明的想法是对任意首一 $ n $ 次多项式按 Legendre 多项式展开，证明其范数不小于 $ \widetilde{P}_n(x) $ 即可.
一般区间上仍然可以使用 $ x=\dfrac{b-a}{2}t+\dfrac{b+a}{2} $ 得到对应的多项式展开. 由于我们采用正交多项式和幂函数作为基生成的都是同样的空间，所以求解得到的最佳平方多项式是完全一致的，但是在计算上使用正交多项式就具有巨大优势了.
当然，除了 Legendre 多项式，我们也可以使用 Chebyshev 多项式在 $ [-1,1] $ 上将 $ f(x)\in C[-1,1] $ 上展成广义 Fourier 级数，只要 $ f' $ 是分段光滑的，那么上述级数能够一致收敛到 $ f $，其取有限项可作为近似最佳一致逼近多项式，这部分内容就不用再展开说明了.

### 曲线拟合的最小二乘法

**曲线拟合**就是要在事先观测到的一组实验数据 $ \lbrace(x_i,y_i),i=0,1,\cdots,m\rbrace $ 的基础上，寻找一个函数 $ y=S^{\ast}(x) $ 使得其与所给数据尽可能接近. 我们将平方逼近中的误差函数离散化，采用完全相同的记号，如果采用误差平方和的形式作为度量，那么有

$$
\begin{equation*}
    ||\boldsymbol{\delta}||_2^2=\sum_{i=0}^m[S^\ast(x_i)-y_i]^2=\min_{S(x)\in\varphi}\sum_{i=0}^m[S(x_i)-y_i]^2,
\end{equation*}
$$

其中 $ S(x)=a_0\varphi_0(x)+\cdots+a_n\varphi_n(x),n<m $. 这是一般的**最小二乘逼近**，几何语言就是曲线拟合的**最小二乘法**. 一般我们都考虑加权平方和，这是说对点 $ (x_i,y_i) $ 可以有重复的观测，从而有对应的权函数 $ \omega(x) $ . 最小化误差平方和仍然可以得到一个法方程 $ \boldsymbol{Ga}=\boldsymbol{d} $，但这里的元素含义会因为离散化发生变化. 若记 $ \boldsymbol{y}=(y_0,y_1,\cdots,y_m),\boldsymbol{\varphi}_j=(\varphi_j(x_0),\varphi_j(x_1),\cdots,\varphi_j(x_m)),j=0,1,\cdots,n $，采用标准内积

$$
\begin{align*}
    (\boldsymbol{\varphi}_j,\boldsymbol{\varphi}_k)&=\sum_{i=0}^m\omega(x_i)\varphi_j(x_i)\varphi_k(x_i),\\
    (\boldsymbol{y},\boldsymbol{\varphi}_k)&=\sum_{i=0}^m\omega(x_i)y_i\varphi_k(x_i)\equiv d_k,~~k=0,1,\cdots,n,
\end{align*}
$$

那么就可以得到离散情况下的法方程. **注意，基函数的无关性此时并不能推出 Gram 矩阵非奇异**，这是和连续场合下有区别的地方. 但是 **Haar 条件**却能够保证 $ \boldsymbol{G} $ 非奇异：

{% note success no-icon %}

（**Haar 条件**）设 $ \varphi_0(x),\varphi_1(x),\cdots,\varphi_n(x)\in C[a,b] $ 的任意线性组合在点集 $ \lbrace x_i\rbrace_{i=0}^m(m\ge n) $ 上至多只有 $ n $ 个不同的零点，则称 $ \varphi_0(x),\varphi_1(x),\cdots,\varphi_n(x) $ 在点集 $ \lbrace x_i\rbrace_{i=0}^m(m\ge n) $ 上满足 Haar 条件.

{% endnote %}

显然 $ 1,x,\cdots,x^n $ 在任意 $ m(m\ge n) $ 个点上满足 Haar 条件. 只要基函数满足这个条件，那么得到的系数解唯一，且极小化平方误差，和连续场合类似. 当然采用幂函数基得到的矩阵仍然是病态的，这时我们仍可以采用带权 $ \omega(x_i)(i=0,1,\cdots,m) $ 的正交函数族作为基来直接得到相应系数. 这时我们仍然有和连续场合下**完全一致的递推式形式，只要把内积形式修改一下即可**. 不能像连续场合那样给出多种通式的原因是我们的观测点不同，以及选取的权函数随实际变化比较大，所以每次都需要重新构造一族正交多项式.
最小二乘法的内容也是较为简单且广泛使用的，案例也很常见，这里也就不需要花太多篇幅介绍了.

## 有理逼近

多项式逼近闭区间上的连续函数计算简便且能够保证收敛性，然而如果函数在某点（非无穷远）无界或在无穷远处收敛时，有限维多项式空间中的任何多项式都不具备这样的性质，这是采用多项式逼近的效果就并不理想. 但这时如果采用如 $ \dfrac{ax+b}{x-c} $ 这样的分式函数可能就能取得良好效果，因为其反映出函数在 $ x=c $ 附近无界而当 $ x\to\infty $ 时趋于定值 $ a $ 的特性. 更一般地，我们可以用形如

$$
\begin{equation*}
    R_{nm}(x)=\frac{p_n(x)}{q_m(x)}=\frac{\sum\limits_{k=0}^na_kx^k}{\sum\limits_{k=0}^mb_kx^k}
\end{equation*}
$$

的函数逼近 $ f(x) $，其中 $ p_n(x),q_m(x) $ 互素. 这种函数逼近称为**有理逼近**. 类似于多项式逼近，若取 $ ||f(x)-R_{nm}(x)||\_{\infty} $ 达到最小，即得**最佳有理一致逼近**，取 $ ||f(x)-R_{nm}(x)||\_{2} $ 最小则得到**最佳有理平方逼近**. 最佳有理逼近这里不做讨论，本节主要还是关注 Padé 逼近.

### 连分式法有理逼近

我们先从连分式开始来讨论有理逼近. 我们知道 $ \ln(1+x) $ 在收敛域上有 Taylor 展开式

$$
\begin{equation*}
    \ln(1+x)=\sum_{k=1}^{\infty}(-1)^{k-1}\frac{x^k}{k},~~x\in[-1,1],
\end{equation*}
$$

其部分和（取前 $ n $ 项）记为 $ S_n(x) $. 根据连分式理论对 Taylor 展开式做连分式展开，可得

$$
\begin{equation*}
    \ln(1+x)=\frac{x}{1+\dfrac{1^2\cdot x}{2+\dfrac{1^2\cdot x}{3+\dfrac{2^2\cdot x}{4+\dfrac{2^2\cdot x}{5+\cdots}}}}},
\end{equation*}
$$

我们并不打算对连分式理论展开说明，只说明其对函数的逼近作用. 若取上述展开的前 $ 2,4,6,8 $ 项（即看每一阶分子的项数，然后把分母加号后面的部分舍去）可得如下几个有理逼近

$$
\begin{equation*}
    \begin{cases}
        R_{11}(x)=\dfrac{2x}{2+x},~~R_{22}(x)=\dfrac{6x+3x^2}{6+6x+x^2},\\
        R_{33}(x)=\dfrac{60x+60x^2+11x^3}{60+90x+36x^2+3x^3},\\
        R_{44}(x)=\dfrac{420x+630x^2+260x^3+25x^4}{420+840x+540x^2+120x^3+6x^4}.
    \end{cases}
\end{equation*}
$$

若取同阶的 Taylor 展开式 $ S_{2n}(x) $ 和上面的连分式进行对比，对 $ \ln 2 $ 的值做近似，会发现在计算量差不多的情况下，有理逼近结果明显比 Taylor 展开快非常多，下标说明 $ R_{44}(1) $ 的精度比 $ S_8(1) $ 要高出近十万倍.

|   $ n $    |   $ S_{2n}(1) $    |  $ \varepsilon_S=\vert\ln 2-S_{2n}(1)\vert $    |  $ R_{nn}(1) $ |    $ \varepsilon_R=\vert\ln 2-R_{nn}(1)\vert $    |
|:-------:|:------:|:------:|:------:|:------:|
|   $ 1 $     |   $ 0.50  $  |  $ 0.19  $   |  $ 0.667 $    |  $ 0.026  $   |
|   $ 2  $   |   $ 0.58  $  | $  0.11  $   |  $  0.69231 $   |  $ 0.00084  $   |
|  $  3  $   |  $  0.617   $ | $   0.076 $   | $ 0.693122    $  |  $ 0.000025 $    |
|  $  4  $   | $   0.634 $   |  $  0.058 $   | $  0.69314642 $    |  $ 0.00000076   $  |

一般而言，利用连分式作为有理分式逼近函数的好处是在计算量相当的情况下取得比 Taylor 展开好得多的效果，这部分内容不作为重点说明.

### Padé 逼近

假设 $ f(x) $ 在 $ x=0 $ 处可以进行 Taylor 展开，展开式记为

$$
\begin{equation*}
    f(x)=\sum_{k=0}^N\frac{f^{(k)}(0)}{k!}x^k+\frac{f^{(N+1)(\xi)}}{(N+1)!}x^{N+1}:=\sum_{k=0}^Nc_kx^k+\frac{f^{(N+1)(\xi)}}{(N+1)!}x^{N+1},
\end{equation*}
$$

前面的部分和记作 $ p(x) $，现在我们想找一个有理函数来逼近 $ p(x) $，以便保留 $ f $ 可能有的无界性质并且又在一个邻域内和原始函数接近，这就引出了 Padé 逼近.

{% note info no-icon %}

设 $ f(x)\in C^{N+1}(-a,a),N=n+m $，如果有理函数

$$
\begin{equation*}
    R_{nm}(x)=\frac{a_0+a_1x+\cdots+a_nx^n}{1+b_1x+\cdots+b_mx^m}=\frac{p_n(x)}{q_m(x)},
\end{equation*}
$$

且 $ p_n(x),q_m(x) $ 互素，且满足 $ R_{nm}^{(k)}(0)=f^{(k)}(0)(k=0,1,\cdots,N) $，则称 $ R_{nm}(x) $ 为函数 $ f(x) $ 在 $ x=0 $ 处的 $ (n,m) $ 阶 **Padé 逼近**，记作 $ R(n,m) $，简称 $ R(n,m) $ 的 Padé 逼近. 上述条件也等价于

$$
\begin{equation*}
    f(x)-R_{nm}(x)=O(x^{N+1}).
\end{equation*}
$$

{% endnote %}

$ R_{nm}(x) $ 共有 $ n+m+1 $ 个系数需要确定，和要求条件中的方程个数一致，因此可以通过解线性方程组确定. 若令

$$
\begin{equation*}
    h(x)=p(x)q_m(x)-p_n(x)=q_m(x)(p(x)-R_{nm}(x)),
\end{equation*}
$$

那么通过 Lebniz 求导公式和 $ q_m(0)=1 $ 即知条件等价于 $ h^{(k)}(0)=0,k=0,1,\cdots,N $. 由于 $ p_n^{(k)}(0)=k!a_k $，再应用 Lebniz 求导公式就知道

$$
\begin{equation*}
    h^{(k)}(0)=\sum_{j=0}^{k}\binom{k}{j}f^{(j)}(0)(k-j)!b_{k-j}-k!a_k=k!\sum_{j=0}^{k}c_jb_{k-j}-k!a_k=0.
\end{equation*}
$$

这里 $ c_j $ 含义同本小节开始，而且 $ b_0=1,b_j=0~(j>m) $. 于是经过化简我们有

$$
\begin{align*}
    a_k&=\sum_{j=0}^{k-1}c_jb_{k-j}+c_k,\quad k=0,1,\cdots,n,\\
    &-\sum_{j=0}^{k-1}c_jb_{k-j}=c_k,\quad k=n+1,\cdots,n+m.
\end{align*}
$$

注意 $ b_j=0,j>m $，因此上面第二个式子还可以化简. 根据求和式的下标特点，第二个方程组可以写成

$$
\begin{equation*}
    \begin{cases}
        -c_{n-m+1}b_m-\cdots-c_{n-1}b_2-c_nb_1=c_{n+1},\\
        -c_{n-m+2}b_m-\cdots-c_{n}b_2-c_{n+1}b_1=c_{n+2},\\
        \cdots~\cdots\\
        -c_{n}b_m-\cdots-c_{n+m-2}b_2-c_{n+m-1}b_1=c_{n+m}.
    \end{cases}
\end{equation*}
$$

其中 $ c_j=0,j<0 $. 上述方程组也可以将其写成矩阵的形式，这里不再赘述. 由于 $ c_j $ 都已知，因此从第二个方程组解出 $ b_1,\cdots,b_m $，代回第一个方程组就可得到所有系数. **当然，我们这里还没有说明为什么上述方程组一定有解.** 事实上，任一在 $ 0 $ 附近能进行 Taylor 展开的函数都存在 Padé 逼近，**但是结果可能是退化的**. 这是因为，我们在定义中假定了 $ q_m(x) $ 的常数项不为零（**即标准化条件**），这可能会导致上面的第二个方程组无解. 如果我们修改定义，不作限制，那么可能会出现新的方程组有无穷解的情况，并且此时 $ p_n(x),q_m(x) $ 不是互素的，从而相应的阶数降低了，而对应的低阶逼近的方程是存在唯一解的. 因此，**如果我们允许退化解的存在，我们可以认为任意函数的 Padé 逼近总是存在且唯一的.** 这个定理我们不做证明，但我们可以讨论在标准化条件下存在解时逼近的唯一性：

{% note primary no-icon %}

对于任意形式的幂级数 $ f(x) $，若其 Padé 逼近 $ R(m,n) $ 存在，则必唯一.

{% endnote %}

**证明** 设有两个这样的逼近 $ R_{nm}(x)=\dfrac{p_n(x)}{q_m(x)},R_{nm}'(x)=\dfrac{p_n'(x)}{q_m'(x)} $，那么由 Padé 逼近定义中的等价条件可知

$$
\begin{equation*}
    R_{nm}(x)-R_{nm}'(x)=O(x^{N+1}).
\end{equation*}
$$

上式同乘 $ q_m(x)q_m'(x) $，则左端次数小于等于 $ N $，但要求和右端一样和 $ x^{N+1} $ 同阶，于是只能是 $ p_n(x)q_m'(x)-p_n'(x)q_m(x)\equiv0 $. 由于两个分母不恒为零，再由标准化条件以及 $ p_n,q_m $ 与 $ p_n',q_m' $ 互素可知只能是 $ p_n(x)=p_n'(x),q_m(x)=q_m'(x) $，唯一性得证.
事实上，根据上面两个方程，可以直接写出 Padé 逼近的结果为下面两个行列式相除的形式：

$$
\begin{equation*}
    R_{nm}(x)=\frac{\begin{vmatrix}
        c_{n-m+1} & \cdots & c_n & c_{n+1} \\
        c_{n-m+2} & \cdots & c_{n+1} & c_{n+2} \\
        \vdots & \vdots & \vdots & \vdots \\
        c_{n} & \cdots & c_{n+m-1} & c_{n+m} \\
        \sum_{i=0}^{m-n}c_ix^{m+i} & \cdots & \sum_{i=0}^{m-1}c_ix^{1+i} & \sum_{i=0}^{m}c_ix^{i}
    \end{vmatrix}}{\begin{vmatrix}
        c_{n-m+1} & \cdots & c_n & c_{n+1} \\
        c_{n-m+2} & \cdots & c_{n+1} & c_{n+2} \\
        \vdots & \vdots & \vdots & \vdots \\
        c_{n} & \cdots & c_{n+m-1} & c_{n+m} \\
        x^m & \cdots & x & 1
    \end{vmatrix}}.
\end{equation*}
$$

其中求和上标为负时结果为 $ 0 $，且 $ c_j=0,j<0 $. 这个结果是由 Jacobi 得到的，实际上这一结果的得到并不很困难，是完全根据上面的矩阵形式和行列式展开得到的，这里就不做证明了.
实践中常把前几个 $ R_{n,m} $ 列在一张二维表——“ Padé 表”中，其形式为

$$
\begin{equation*}
    \begin{matrix}
        R(0,0) & R(0,1) & R(0,2) & \cdots\\
        R(1,0) & R(1,1) & R(1,2) & \cdots\\
        R(2,0) & R(2,1) & R(2,2) & \cdots\\
        \vdots & \vdots & \vdots & \ddots
    \end{matrix}
\end{equation*}
$$

大量计算表明在两阶数和相同时，效果较好的 Padé 逼近集中在对角线附近，因此我们总是选择满足 $ |m-n|\le1 $ 的阶数来对函数进行逼近. 此外，上一小节的连分式逼近实际上和 Padé 表的对角元有一一对应关系.

## 三角多项式逼近与快速 Fourier 变换

自然界中的复杂振动现象常常由不同频率及振幅的波叠加得到，而复杂的波还可分解为一系列谐波. 在模型数据具有周期性时，使用正弦函数和余弦函数作为基函数拟合较为合适，能较好地反映周期性. 在数学分析的课程中有介绍使用三角多项式作为基函数表示或逼近任意函数的方法，称为 **Fourier 变换**，在计算机分析中主要用到了三角函数逼近给定样本函数的最小二乘和插值，主要是用于解决当函数只能用有限个离散点上的测量值来表示时，如何像傅里叶分析那样提取出其中的频率成分的问题. 这种方法称为**离散 Fourier 变换（DFT）**. 然而由于 DFT 的计算量很大，一时间难以应用，直到**快速 Fourier 变换（FFT）**的出现，对原始算法进行了运算量级的大幅改进，DFT 才得以大规模运用. 我们从离散 Fourier 系数开始，逐步深入说明算法.

### 最佳三角平方逼近与三角插值

我们知道周期为 $ 2\pi $ 的函数 $ f(x)\in L^2[a,b] $ 的 Fourier 级数为

$$
\begin{equation*}
    f(x)\sim\frac{1}{2}a_0+\sum_{k=1}^{\infty}(a_k\cos kx+b_k\sin kx),
\end{equation*}
$$

其中

$$
\begin{equation*}
    \begin{cases}
        a_k=\dfrac{1}{\pi}\displaystyle{\int_{0}^{2\pi}}f(x)\cos kx\mathrm{d}x, & k=0,1,\cdots,\\
        b_k=\dfrac{1}{\pi}\displaystyle{\int_{0}^{2\pi}}f(x)\sin kx\mathrm{d}x, & k=1,2,\cdots,\\
    \end{cases}
\end{equation*}
$$

并且取部分和函数 $ S_n(x) $（即将上述级数求和从无穷改为有限的 $ n $）就是 $ f $ 的**最佳平方三角逼近**多项式. 若 $ f(x) $ 在 $ [0,2\pi] $ 上是分段光滑的，那么可将上述级数展开中的 $ \sim $ 换成等号，级数是一致收敛到 $ f(x) $ 的. 即使没有这一性质，我们仍有对应的 Bessel 不等式

$$
\begin{equation*}
    \frac{1}{2}a_0^2+\sum_{k=1}^n(a_k^2+b_k^2)\le\frac{1}{\pi}\int_{0}^{2\pi}f^2(x)\mathrm{d}x,
\end{equation*}
$$

同时由 $ f\in L^2[0,2\pi] $ 可知 $ \lim\limits_{n\to\infty}a_n=\lim\limits_{n\to\infty}b_n=0 $.
但上述连续情况在实际中却不是常见情况，因为我们对数据的采样很多时候是离散的. 但我们依然可以在离散点集上定义正交性和相应的离散 Fourier 系数，利用 Euler 公式和积化和差公式不难证明，如果我们在 $ [0,2\pi] $ 上等距地取 $ N=2m+1 $ 个点（这时取奇数个点的原因将在后面说明），即

$$
\begin{equation*}
    x_j=\frac{2\pi j}{2m+1},~~j=0,1,\cdots,2m,
\end{equation*}
$$

那么函数族 $ \lbrace1,\cos x,\sin x,\cdots,\cos mx,\sin mx\rbrace $ 在 $ \left\lbrace x_j=\dfrac{2\pi j}{2m+1}\right\rbrace $ 上正交. 若令 $ f_j=f(x_j)(j=0,1,\cdots,2m) $，那么 $ f $ 的最小二乘三角逼近是

$$
\begin{equation*}
    S_n(x)=\frac{1}{2}a_0+\sum_{k=0}^{n}(a_k\cos kx+b_k\sin kx), ~~ n<m,
\end{equation*}
$$

其中

$$
\begin{align*}
    a_k&=\frac{2}{2m+1}\sum_{j=0}^{2m}f_j\cos\frac{2\pi jk}{2m+1}, \quad k=0,1,\cdots,n,\\
    b_k&=\frac{2}{2m+1}\sum_{j=0}^{2m}f_j\sin\frac{2\pi jk}{2m+1}, \quad k=1,2,\cdots,n.
\end{align*}
$$

最小二乘的计算过程这里就不再说明了. 当 $ n=m $ 时，可以证明

$$
\begin{equation*}
    S_m(x_j)=f_j,~~j=0,1,\cdots,2m,
\end{equation*}
$$

于是上述 $ S_n(x)=S_m(x) $ 就是**三角插值多项式**. 一般地，不对 $ N $ 的奇偶性做限制，那么只要 $ m\le\left\lfloor\dfrac{N-1}{2}\right\rfloor $，就能保证三角函数族的正交性，否则会因为频率的重叠破坏正交性，即所谓的**混叠现象**. 要解释这一点，我们直接将情形推广到复数域上，假设现在 $ f $ 是 $ [0,2\pi] $ 上的复值函数，并且在 $ N $ 个等分点 $ x_j=\dfrac{2\pi j}{N}(j=0,1,\cdots,n-1) $ 上的值 $ f_j=f(x_j) $，接下来我们简要说明函数族 $ \lbrace 1,\mathrm{e}^{\mathrm{i}x},\cdots,\mathrm{e}^{\mathrm{i}(N-1)x}\rbrace $ 在 $ [0,2\pi] $ 上是正交的. 记函数 $ \mathrm{e}^{\mathrm{i}jx} $ 在上述等分点上的值组成的向量记作

$$
\begin{equation*}
    \boldsymbol{\phi}_j=\left(1,\mathrm{e}^{\mathrm{i}j\frac{2\pi}{N}},\cdots,\mathrm{e}^{\mathrm{i}j\frac{2\pi}{N}(N-1)}\right)'
\end{equation*}
$$

那么使用等比数列求和公式很容易得到

$$
\begin{equation*}
    (\boldsymbol{\phi}_l,\boldsymbol{\phi}_s)=\sum_{k=0}^{N-1}\mathrm{e}^{\mathrm{i}(l-s)\frac{2\pi}{N}k}=\begin{cases}
        0, & l\ne s,\\
        N, & l =s.
    \end{cases}
\end{equation*}
$$

这就说明了正交性. **注意，此时我们选取的基是正余弦函数成对在实部虚部搭配出现的，而我们有恒等式 $ \cos x=\dfrac{\mathrm{e}^{\mathrm{i}x}+\mathrm{e}^{-\mathrm{i}x}}{2},\sin x=\dfrac{\mathrm{e}^{\mathrm{i}x}-\mathrm{e}^{-\mathrm{i}x}}{2\mathrm{i}} $**，这对实数域上的 $ x $ 仍然成立. 在我们采样的离散点 $ x_j $ 上，经过简单的计算就可知 $ \mathrm{e}^{\mathrm{i}(N-j)x_k}=\mathrm{e}^{\mathrm{i}jx_k},\forall j,k\in\lbrace0,1,\cdots,N-1\rbrace $，如果三角函数族里的元素族太多，那么根据这个恒等式就知道 $ \cos(N-j)x=\cos jx,\sin(N-j)x=-\sin jx $ 在选取的离散点上恒成立，产生了线性相关性，从而破坏了正交性. 因此，在选取基函数时，必须有 $ 2m<N $，这里的 $ m $ 就是在实数域上考虑的那个 $ m $，从而得到之前提到的 $ m\le\left\lfloor\dfrac{N-1}{2}\right\rfloor $ 这一条件.
从上述结论就可导出复数域上的 $ f $ 在 $ N $ 个点 $ \left\lbrace x_j=\dfrac{2\pi}{N}j,j=0,1,\cdots,N-1\right\rbrace $ 上的最小二乘 Fourier 逼近为

$$
\begin{equation*}
    S(x)=\sum_{k=0}^{n-1}c_k\mathrm{e}^{\mathrm{i}kx},~~n\le N,~~c_k=\frac{1}{N}\sum_{j=0}^{N-1}f_j\mathrm{e}^{-\mathrm{i}kj\frac{2\pi}{N}},k=0,1,\cdots,n-1.
\end{equation*}
$$

若 $ n=N $，那么 $ S(x) $ 就是 $ f(x) $ 在点 $ x_j(j=0,1,\cdots,N-1) $ 上的插值函数，从而由 $ S(x) $ 的表达式马上可知

$$
\begin{equation*}
    S(x_j)=f(x_j)=f_j=\sum_{k=0}^{N-1}c_k\mathrm{e}^{\mathrm{i}k\frac{2\pi}{N}j},~~j=0,1,\cdots,N-1.
\end{equation*}
$$

由 $ \lbrace f_j \rbrace $ 确定 $ \lbrace c_k\rbrace $ 的过程称为 $ f(x) $ 的**离散 Fourier 变换**，而由插值表达式得到的由 $ \lbrace c_k\rbrace $ 确定 $ \lbrace f_j\rbrace $ 的过程称为**反变换**. 这里要特别注意 $ c_k $ 的指数上有负号，这是因为复数域上的共轭操作导致的，求解过程这里就不写出了. 两个过程是计算机进行 Fourier 分析的主要方法，在数字信号处理、全息技术、光谱和声谱分析等很多领域有广泛应用.

### $ N $ 点 DFT 和 FFT 算法

在上一节中提到的关于离散 Fourier 变换和逆变换的过程，都可归结为计算

$$
\begin{equation*}
    c_j=\sum_{k=0}^{N-1}x_k\omega_N^{kj},~~j=0,1,\cdots,N-1,
\end{equation*}
$$

这里 $ c,f $ 的记号实际可以看作两个地位等价的东西，只是根据实际需要来解释其含义. 其中 $ \lbrace x_k\rbrace_{k=0}^{N-1} $ 是输入数据，$ \omega_N=\mathrm{e}^{\frac{2\pi}{N}\mathrm{i}} $ 是 $ N $ 次单位根之一. 上述归结式被称为 $ N $ 点 DFT，直观地看，计算 $ c_j $ 只需要进行 $ N $ 次复数乘法和 $ N $ 次加法，统称为 $ N $ 次操作，计算所有结果则需 $ N^2 $ 次操作，时间复杂度即 $ O(N^2) $，这在 $ N $ 很大时计算时间开销是难以接受的. 在 DFT 算法基础上的改进产生了 FFT 算法，让 DFT 得到了广泛应用. 其思想是减少乘法次数减少耗时，这就需要用到单位根的性质，根据 $ \omega_N $ 的表达式，下述等式是比较显然的：

$$
\begin{align*}
    \omega_N^j\omega_N^k=\omega_N^{j+k},\quad\omega_N^{Nj+k}=\omega_N^k~(周期性),\quad\omega_{jN}^{jk}=\omega_{N}^{k}.
\end{align*}
$$

由正余弦函数的周期性，$ \omega_N^{jk}(j,k=0,1,\cdots,N-1) $ 最多只能取 $ N $ 个不同的值，特别地，当 $ N $ 为偶数时，

$$
\begin{equation*}
    \omega_{N}^{N/2}=-1,~~\omega_{N}^{jk+N/2}=-\omega_{N}^{jk}~(对称性)
\end{equation*}
$$

从而 $ \omega_{N}^{jk} $ 只有 $ N/2 $ 个相对意义上不同的值，即成对出现. 利用这一性质，我们可以将表达式改写：

$$
\begin{equation*}
    c_j=\sum_{k=0}^{N/2-1}x_k\omega_{N}^{jk}+\sum_{k=0}^{N/2-1}x_{N/2+k}\omega_{N}^{j(N/2+k)}=\sum_{k=0}^{N/2-1}\left[x_k+(-1)^jx_{N/2+k}\right]\omega_{N}^{jk}.
\end{equation*}
$$

再按下标奇偶性考察，易知

$$
\begin{align*}
    c_{2j}&=\sum_{k=0}^{N/2-1}\left[x_k+x_{N/2+k}\right]\omega_{N/2}^{jk}:=\sum_{k=0}^{N/2-1}y_k\omega_{N/2}^{jk},\\
    c_{2j+1}&=\sum_{k=0}^{N/2-1}\left[(x_k-x_{N/2+k}\right)\omega_{N}^k]\omega_{N/2}^{jk}:=\sum_{k=0}^{N/2-1}y_k\omega_{N/2}^{jk}.
\end{align*}
$$

于是我们将 $ N $ 点 DFT 归结为了两个 $ N/2 $ 点 DFT. 如果 $ N/2 $ 为偶数，那么上述步骤可重复进，如此反复进行二分手续即可得 FFT 算法. 因此只要取 $ N $ 为 $ 2 $ 的幂次，这样就可以不断进行二分的步骤，如果不是 $ 2 $ 的幂，则可以通过增加若干零值点来满足条件. 补零会导致 Fourier 逼近的结果改变，但在实际问题中我们并不关心这一点，由于本人并非工科类专业学生，对这方面内容并不了解，因此主要还是聚焦于具体的 FFT 算法如何进行.
下面去取 $ N=2^3 $ 来说明 FFT 算法是如何使用的. 此时 $ k,j=0,1,\cdots,7 $ 记 $ \omega_N=\omega_8=\omega $，现在我们要计算

$$
\begin{equation*}
    c_j=\sum_{k=0}^{7}x_k\omega^{jk}, ~~j=0,1,\cdots,7,
\end{equation*}
$$

我们将 $ k,j $ 使用二进制表示为

$$
\begin{equation*}
    k=k_22^2+k_12^1+k_02^0:=(k_2k_1k_0),\quad j=j_22^2+j_12^1+j_02^0:=(j_2j_1j_0),
\end{equation*}
$$

那么可以记 $ c_j=c(j_2j_1j_0),x_k=x(k_2k_1k_0) $，则

$$
\begin{align*}
    c(j_2j_1j_0)&=\sum_{k_0=0}^1\sum_{k_1=0}^1\sum_{k_2=0}^1x(k_2k_1k_0)\omega^{(k_2k_1k_0)(j_22^2+j_12^1+j_02^0)}\\
    &=\sum_{k_0=0}^1\left\lbrace\sum_{k_1=0}^1\left[\sum_{k_2=0}^1x(k_2k_1k_0)\omega^{j_0(k_2k_1k_0)}\right]\omega^{j_1(k_1k_00)}\right\rbrace\omega^{j_2(k_000)}.
\end{align*}
$$

第二个等式的原因实际上是二进制表达式的位移导致的，由于 $ \omega^8=1 $，因此可以约去多余部分. 如果再引进符号

$$
\begin{align*}
    A_0(k_2k_1k_0)&=x(k_2k_1k_0),\\
    A_1(k_1k_0j_0)&=\sum_{k_2=0}^1A_0(k_2k_1k_0)\omega^{j_0(k_2k_1k_0)},\\
    A_2(k_0j_1j_0)&=\sum_{k_1=0}^1A_1(k_1k_0j_0)\omega^{j_1(k_1k_00)},\\
    A_3(j_2j_1j_0)&=\sum_{k_0=0}^1A_2(k_0j_1j_0)\omega^{j_2(k_000)},
\end{align*}
$$

那么有 $ c(j_2j_1j_0)=A_3(j_2j_1j_0) $. 注意 $ \omega^{j_02^{3-1}}=\omega^{j_0N/2}=(-1)^{j_0} $，那么对于 $ A_1 $，其可以化简为

$$
\begin{align*}
    A_1(k_1k_0j_0)&=A_0(0k_1k_0)\omega^{j_0(0k_1k_0)}+A_0(1k_1k_0)\omega^{j_0(1k_1k_0)}\\
    &=[A_0(0k_1k_0)+(-1)^{j_0}A_0(1k_1k_0)]\omega^{j_0(0k_1k_0)},\\
    A_1(k_1k_00)&=A_0(0k_1k_0)+A_0(1k_1k_0),\\
    A_1(k_1k_01)&=[A_0(0k_1k_0)-A_0(1k_1k_0)]\omega^{(0k_1k_0)}.
\end{align*}
$$

将二进制表示全部还原为十进制，令 $ k=(0k_1k_0)=k_12^1+k_02^0\Rightarrow k=0,1,2,3 $，那么还可以将上述规则写成十进制表示下的递推式

$$
\begin{equation*}
    \begin{cases}
        A_1(2k)=A_0(k)+A_0(k+2^2),\\
        A_1(2k+1)=[A_0(k)-A_0(k+2^2)]\omega^k,
    \end{cases}
    k=0,1,2,3.
\end{equation*}
$$

因此我们可以确定所有 $ A_1 $ 的值，完全类似地，我们有

$$
\begin{align*}
    A_2(k_00j_0)&=A_1(0k_0j_0)+A_1(1k_0j_0),\\
    A_2(k_01j_0)&=[A_1(0k_0j_0)-A_1(1k_0j_0)]\omega^{(0k_00)},
\end{align*}
$$

对应的十进制表示为

$$
\begin{equation*}
    \begin{cases}
        A_2(k2^2+j)=A_1(2k+j)+A_1(2k+j+2^2),\\
        A_2(k2^2+j+2)=[A_1(2k+j)-A_1(2k+j+2^2)]\omega^{2k},
    \end{cases}
    k,j=0,1.
\end{equation*}
$$

再次类似可得

$$
\begin{equation*}
    \begin{cases}
        A_3(j)=A_2(j)+A_2(j+2^2),\\
        A_3(j+2^2)=A_2(j)-A_2(j+2^2),
    \end{cases}
    j=0,1,2,3.
\end{equation*}
$$

根据上述式子可以递推计算出所有 $ c_j $. 在计算机中，我们只需要预先知道 $ \omega^{k},k=0,1,2,3 $ 的值以及所有 $ x_k=A_0(k) $ 并将它们存放在数组里，再开辟新的二维数组保存所有递推计算结果，就可以得到最终值. 对于上述 $ 8 $ 个点采样的例子，递推过程中只用到了 $ 8 $ 次乘法运算和 $ 24 $ 次加法运算.
我们可以将结果推广到 $ N=2^p $ 的情形，一般的 FFT 计算公式为

$$
\begin{equation*}
    \begin{cases}
        A_q(k2^q+j)=A_{q-1}(k2^{q-1}+j)+A_{q-1}(k2^{q-1}+j+2^{p-1}),\\
        A_q(k2^q+j+2^{q-1})=[A_{q-1}(k2^{q-1}+j)-A_{q-1}(k2^{q-1}+j+2^{p-1})]\omega^{k2^{q-1}},\\
        q=1,2,\cdots,p;k=0,1,\cdots,2^{p-q}-1;j=0,1,\cdots,2^{q-1}-1.
    \end{cases}
\end{equation*}
$$

一组 $ A_q $ 占用 $ N $ 个复数单元，只要从 $ A_0(m)(m=0,1,\cdots,N-1) $ 出发， $ q $ 从 $ 1 $ 到 $ p $ 算到 $ A_p(j)=c_j $ 就可以求得离散频谱. 计算可分为两重循环，外层是关于 $ q=1,2,\cdots,p $ 的，内层是关于 $ k=0,1,\cdots,2^{p-q}-1,j=0,1,\cdots,2^{q-1}-1 $ 的，算一个 $ A_q $ 只需进行 $ 2^{p-q}\cdot2^{q-1}=N/2 $ 次复数乘法，最后一步计算 $ A_p $ 时甚至不需要复数乘法，因为 $ \omega^{k2^{p-1}}=1 $，因此总共需要 $ (p-1)N/2 $ 次乘法，因此整体计算的时间复杂度是 $ O(N\log N) $. 上述递推计算式被称为**改进 FFT 算法**，一般的程序步骤如下：

{% note info no-icon %}

**Step 1** 初始化数组 $ A_1(N),A_2(N),\omega(N/2),N=2^p $；
**Step 2** 将已知的记录复数数组 $ \lbrace x_k\rbrace $ 输入单元 $ A_1(k)(k=0,1,\cdots,N-1) $；
**Step 3** 计算 $ \omega(m)=\exp\left(-\mathrm{i}\dfrac{2\pi}{N}m\right)(m=0,1,\cdots,N/2-1) $；
**Step 4** $ q $ 循环从 $ 1 $ 到 $ p $，若 $ q $ 为奇数进行 Step 5，否则进行 Step 6；
**Step 5** $ k $ 循环从 $ 0 $ 到 $ 2^{p-q}-1 $， $ j $ 循环从 $ 0 $ 到 $ 2^{q-1}-1 $，计算

$$
\begin{equation*}
    \begin{cases}
        A_2(k2^q+j)=A_{1}(k2^{q-1}+j)+A_{1}(k2^{q-1}+j+2^{p-1}),\\
        A_2(k2^q+j+2^{q-1})=[A_{1}(k2^{q-1}+j)-A_{1}(k2^{q-1}+j+2^{p-1})]\omega(k2^{q-1})
    \end{cases}
\end{equation*}
$$

并转 Step 7；

**Step 6** $ k $ 循环从 $ 0 $ 到 $ 2^{p-q}-1 $， $ j $ 循环从 $ 0 $ 到 $ 2^{q-1}-1 $，计算

$$
\begin{equation*}
    \begin{cases}
        A_1(k2^q+j)=A_{2}(k2^{q-1}+j)+A_{2}(k2^{q-1}+j+2^{p-1}),\\
        A_1(k2^q+j+2^{q-1})=[A_{2}(k2^{q-1}+j)-A_{2}(k2^{q-1}+j+2^{p-1})]\omega(k2^{q-1})
    \end{cases}
\end{equation*}
$$

循环结束进行下一步；

**Step 7** 若 $ q=p $ 则转 Step 8，否则 $ q\leftarrow q+1 $，转 Step 4；
**Step 8** $ q $ 循环结束，若 $ p $ 为偶数，则 $ A_2(j)\leftarrow A_1(j) $，结果为 $ c_j=A_2(j)(j=0,1,\cdots,N-1) $.

{% endnote %}

上述算法对内存开销进行了优化，因为 $ A_q $ 的计算仅依赖于 $ A_{q-1} $，因此可将前面的记录覆盖，最后一步也可以不再进行赋值，根据 $ p $ 的奇偶性选择对应数组.
更多时候我们还是在实数域上考虑问题. 如果给出 $ N $ 个点就需要 $ N $ 个基函数来构建插值多项式. 当点的个数为 $ N=2m+1 $ 个时，假定的插值函数形式就是

$$
\begin{equation*}
    S_m(x)=\frac{a_0}{2}+\sum_{k=1}^m(a_k\cos kx+b_k\sin kx),
\end{equation*}
$$

对应的系数可以参考先前得到的关于奇数个点的最小二乘三角逼近的系数. 事实上这一系数对偶数点也成立，统一地可以写作

$$
\begin{align*}
    a_k&=\frac{2}{N}\sum_{j=0}^{N-1}f_j\cos\frac{2\pi jk}{N},~~k=0,1,\cdots,\left\lceil\frac{N-1}{2}\right\rceil,\\
    b_k&=\frac{2}{N}\sum_{j=0}^{N-1}f_j\sin\frac{2\pi jk}{N},~~k=0,1,\cdots,\left\lceil\frac{N-1}{2}\right\rceil.
\end{align*}
$$

当 $ k=0 $ 时总有 $ b_0=0 $，这和展开式的结构相对应. 当 $ N=2m $ 时，上述插值函数的末尾会多一个余弦项，如果把多的对应正弦部分加上，可以发现对应的 $ b_m=0 $，后面马上会说明这一观察的用处. 我们可以直接用复数形式的 $ c_j $ 确定对应的 $ a_j,b_j $，实际上根据两者的结构，不难列出下面的恒等式

$$
\begin{equation*}
    a_j+\mathrm{i}b_j=\frac{2}{N}c_j,\quad j=0,1,\cdots,\left\lceil\frac{N-1}{2}\right\rceil,
\end{equation*}
$$

比较实部和虚部就可得到每个值. 而且我们不需要担心有多余插值项的问题，因为偶数情形，$ c_{N/2} $ 的虚部恰好为零，和 $ b_{N/2} $ 正好对应，奇数项则恰好匹配.
