---
title: EM 算法
date: 2026-07-10 21:24:41
tags: [数理统计, 极大似然估计, 统计机器学习]
categories: [文献阅读, 学习笔记]
---

{% note info no-icon %}
本文是对经典文章 {% link "**Maximum Likelihood from Incomplete Data via the *EM* Algorithm**" https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/j.2517-6161.1977.tb01600.x %} 的回顾，并且参考文章 {% link **On The Convergence Properties of EM Algorithm** https://www.jstor.org/stable/2240463 %} 对原始文献中关于算法错误收敛性的证明进行了修正，篇幅上进行了删减.
{% endnote %}

**EM 算法**是用来求解不完全观测数据的极大似然估计的有效算法. 关于 EM 算法的讨论在二十世纪六七十年代就逐步开展了. 本篇内容所基于的文章实际上结合了很多有关作者的研究.
首先我们给出**不完全数据**的一般表达：

{% note primary no-icon %}
设有两个样本空间 $ \mathcal{X},\mathcal{Y} $ 并且有一个从 $ \mathcal{X} $ 到 $ \mathcal{Y} $ 的**多对一映射** $ \boldsymbol{x}\to\boldsymbol{y}(\boldsymbol{x}) $，其中 $ \boldsymbol{y} $ 是样本空间 $ \mathcal{Y} $ 的一个实现，而对应的 $ \boldsymbol{x} $ 无法直接观测到，其信息蕴含在 $ \mathcal{X}(\boldsymbol{y}):=\lbrace \boldsymbol{x}\in\mathcal{X}:\boldsymbol{y}=\boldsymbol{y}(\boldsymbol{x})\rbrace $ 这一集合中，其中 $ \boldsymbol{y} $ 是一条观测到的数据. $ \mathcal{X}(\boldsymbol{y}) $ 实际上可以看作 $ \mathcal{Y} $ 到 $ \mathcal{X} $ 的逆映射. 我们称 $ \boldsymbol{x} $ 为**完全数据（complete data）**，相应地， $ \boldsymbol{y} $ 为**不完全数据（incomplete data）**.
**注意**： $ \boldsymbol{x} $ 甚至可以是传统意义上被称作参数的数据.
{% endnote %}

假设我们有一族依赖于参数 $ \boldsymbol{\phi} $ 的完全数据的样本密度函数 $ f(\boldsymbol{x}|\boldsymbol{\phi}) $，对应的不完全数据的样本密度函数为 $ g(\boldsymbol{y}|\boldsymbol{\phi}) $，那么由上面的定义可以知道有如下关系

$$
\begin{equation*}
    g(\boldsymbol{y}|\boldsymbol{\phi})=\int_{\mathcal{X}(\boldsymbol{y})}f(\boldsymbol{x}|\boldsymbol{\phi})\mathrm{d}\boldsymbol{x}.
\end{equation*}
$$

EM 算法是为了找到使得一个参数 $ \boldsymbol{\phi} $ 使得 $ g(\boldsymbol{y}|\boldsymbol{\phi}) $ 在给定 $ \boldsymbol{y} $ 时最大化，但这也可以通过对完全数据空间的密度函数族的必要使用来实现. 有些时候我们可以对原始的密度函数族做出假定且这些假定比较直观，但有些时候需要通过某些方式去定义 $ f(\cdot|\cdot) $. EM 算法具体可以分为两个步骤，即 E 步（期望步）和 M 步（最大化步），我们用一个数值案例来说明其是怎么操作的，给出一些直觉.

## 引入：以多项分布为例

Rao 在一篇论文中给出了将 197 种动物分成四类的例子，其服从多项分布. 观测数据为

$$
\begin{equation*}
    \boldsymbol{y}=(y_1,y_2,y_3,y_4)=(125,18,20,34),
\end{equation*}
$$

经典的对群体建立的基因模型给出了每类出现的概率分布为

$$
\begin{equation*}
    \left(\frac{1}{2}+\frac{1}{4}\pi,\frac{1}{4}(1-\pi),\frac{1}{4}(1-\pi),\frac{1}{4}\pi\right), \quad 0\le\pi\le 1.
\end{equation*}
$$

于是观测数据的概率分布为

$$
\begin{equation*}
    g(\boldsymbol{y}|\pi)=\frac{(y_1+y_2+y_3+y_4)!}{y_1!y_2!y_3!y_4!}\left(\frac{1}{2}+\frac{1}{4}\pi\right)^{y_1}\left(\frac{1}{4}-\frac{1}{4}\pi\right)^{y_2}\left(\frac{1}{4}-\frac{1}{4}\pi\right)^{y_3}\left(\frac{1}{4}\pi\right)^{y_4}.
\end{equation*}
$$

我们的目标是求出 $ \pi $ 的极大似然估计. 然而，虽然我们可以直接得到对数似然函数的驻点并求出 $ \hat{\pi} $，但是通过一些计算可以看到，求解实际上需要解一个三次方程，直接计算极为麻烦. 为此，我们将观测得到的 $ \boldsymbol{y} $ 进行拆分，**改写成一个有五类的多项分布**的总体. 这个步骤**实际上就指定了映射 $ \boldsymbol{y}(\boldsymbol{x}) $**. 具体而言，我们将 $ y_1 $ 拆成 $ y_1=x_1+x_2 $ 并且 $ x_1 $ 出现概率指定为 $ \dfrac{1}{2} $， $ x_2 $ 出现概率指定为 $ \dfrac{\pi}{4} $，其他保持不变，即 $ y_2=x_3,y_3=x_4,y_4=x_5 $ 并且出现概率保持不变. 那么完全数据的概率质量函数就是

$$
\begin{equation*}
    f(\boldsymbol{x}|\pi)=\frac{(x_1+x_2+x_3+x_4+x_5)!}{x_1!x_2!x_3!x_4!x_5!}\left(\frac{1}{2}\right)^{x_1}\left(\frac{1}{4}\pi\right)^{x_2}\left(\frac{1}{4}-\frac{1}{4}\pi\right)^{x_3}\left(\frac{1}{4}-\frac{1}{4}\pi\right)^{x_4}\left(\frac{1}{4}\pi\right)^{x_5}.
\end{equation*}
$$

注意这里 $ x_1,x_2 $ 是随机变量，但满足 $ x_1+x_2=y_1=125 $. EM 算法的操作步骤具体体现为：首先在给定 $ \boldsymbol{y} $ 时得到完全数据 $ \boldsymbol{x} $ 的充分统计量的估计值，然后用极大似然方法在估计得到的完全数据下求出新的参数估计值. 在本例中，我们从 $ \pi^{(0)} $ 开始，重复上述步骤得到 $ \pi^{(p)},p=1,2,\cdots $. 为了说明具体的操作，首先注意 $ x_3,x_4,x_5 $ 都是固定的，实际上我们只要得到 $ x_1,x_2 $ 在第 $ p $ 步时的期望. 我们首先通过 Bayes 定理求出 $ x_1^{(p)}(x_2^{(p)}) $ 的分布，有

$$
\begin{align*}
    \mathbb{P}(x_1^{(p)}|\boldsymbol{y},\pi^{(p)})&=\frac{f(x_1^{(p)},y_1-x_1^{(p)},y_2,y_3,y_4|\pi^{(p)})}{g(\boldsymbol{y}|\pi^{(p)})}\\
    &=\binom{y_1}{x_1^{(p)}}\left(\dfrac{2}{2+\pi}\right)^{x_1^{(p)}}\left(\frac{\pi}{2+\pi}\right)^{y_1-x_1^{(p)}}\sim b\left(y_1,\frac{2}{2+\pi^{(p)}}\right).
\end{align*}
$$

于是对 $ x_1,x_2 $ 求期望分别有

$$
\begin{equation*}
    \mathbb{E}[x_1^{(p)}|\boldsymbol{y},\pi^{(p)}]=y_1\cdot\dfrac{2}{2+\pi^{(p)}},~~\mathbb{E}[x_2^{(p)}|\boldsymbol{y},\pi^{(p)}]=y_1\cdot\dfrac{\pi^{(p)}}{2+\pi^{(p)}}.
\end{equation*}
$$

这就是 EM 算法中的 E 步. 接下来使用估计得到的完全数据 $ \boldsymbol{x}^{(p)}=(x_1^{(p)},x_2^{(p)},x_3,x_4,x_5) $ 代入 $ f $，得到关于 $ \pi $ 的似然函数. 为了计算简便，这里我们直接使用对数似然，表达式可写作

$$
\begin{equation*}
    \ell(\pi|\boldsymbol{x}^{(p)})=(x_2^{(p)}+x_5)\log\pi+(x_3+x_4)\log(1-\pi)+C,
\end{equation*}
$$

其中 $ C $ 是与 $ \pi $ 无关的常量. 求解驻点易得更新等式

$$
\begin{equation*}
    \pi^{(p+1)}=\frac{x_2^{(p)}+x_5}{x_2^{(p)}+x_5+x_3+x_4}.
\end{equation*}
$$

反复迭代至指定收敛精度就能得到 $ \hat{\pi} $. 原始文章通过对上述迭代式解不动点给出了具体的数值 $ \hat{\pi}\approx 0.6268214980 $. 我们可以通过如下 `python` 代码来实现具体求解过程：

```python
import numpy as np

def EM(y0, pi0, pi_true, tol=1e-6, max_iter=1000):
    y1, pi_old = y0[0], pi0
    x3, x4, x5 = y0[1], y0[2], y0[3]
    for i in range(max_iter):
        # E 步
        # x1 = y1 * 2 / (2 + pi_old) 这一步其实可以省略
        x2 = y1 * pi_old / (2 + pi_old)
        # M 步
        pi_new = (x2 + x5) / (x2 + x5 + x3 + x4)
        # 检查收敛性
        if(np.abs(pi_new - pi_old) < tol):
            return pi_new
        print(f"Iteration {i+1}: pi = {pi_new:.9f}, pi - pi_true = {pi_new - pi_true:.9f}, rate = {(pi_new - pi_true) / (pi_old - pi_true):.4f}")
        pi_old = pi_new # 更新参数

y0, pi0, pi_true = np.array([125, 18, 20, 34]), np.array(0.5), np.array(0.6268214980)
res = EM(y0, pi0, pi_true)

# 输出
Iteration 1: pi = 0.608247423, pi - pi_true = -0.018574075, rate = 0.1465
Iteration 2: pi = 0.624321050, pi - pi_true = -0.002500448, rate = 0.1346
Iteration 3: pi = 0.626488879, pi - pi_true = -0.000332619, rate = 0.1330
Iteration 4: pi = 0.626777322, pi - pi_true = -0.000044176, rate = 0.1328
Iteration 5: pi = 0.626815632, pi - pi_true = -0.000005866, rate = 0.1328
Iteration 6: pi = 0.626820719, pi - pi_true = -0.000000779, rate = 0.1328
```

经测试结果和论文中得到的是完全一致的. 这样的迭代计算降低了原来求解三次方程的复杂度，并且能够保证收敛到真值. 关于收敛性我们后文会展开进一步讨论.

## EM 算法的一般形式

### 指数族形式

我们先在对完全数据的密度函数族 $ f(\boldsymbol{x}|\boldsymbol{\phi}) $ 限制较严格的基础上来谈 EM 算法. 我们首先假设 $ f $ 是**正则指数族**，其具有如下形式

$$
\begin{equation*}
    f(\boldsymbol{x}|\boldsymbol{\phi})=b(\boldsymbol{x})\exp[\boldsymbol{\phi}'\boldsymbol{T}(\boldsymbol{x})]/a(\boldsymbol{\phi}).
\end{equation*}
$$

其中 $ \boldsymbol{\phi}\in\Omega\subset\mathbb{R}^{r\times 1} $ 为参数向量，其中 $ \Omega $ 为凸集. $ \boldsymbol{T}(\boldsymbol{x})\in\mathbb{R}^{r\times1} $ 为完全数据的充分统计量. 上述关于参数的正则性条件保证了在不同充分统计量的选择下，通过一个非奇异 $ r\times r $ 线性变换和不同充分统计量所联系的参数向量是唯一确定的. 满足上述正则条件的参数被称为**自然参数**. 注意，自然参数和传统意义上的分布参数会有不同，可能和一个非线性变换联系在一起. 比如对二项分布 $ b(n,p) $（其中 $ n $ 是确定的），写成指数族可以发现 $ \phi=\log\left(\dfrac{p}{1-p}\right) $. **现在我们都将 $ \boldsymbol{\phi} $ 视作自然参数**. 此外我们要求，对所有的 $ \boldsymbol{\phi}\in\Omega $，密度函数 $ f(\cdot|\boldsymbol{\phi}) $ 都具有共同的支撑.
现在我们可以阐述**对于正则指数族的 EM 算法**. 假设 $ \boldsymbol{\phi}^{(p)} $ 是 $ p $ 步迭代后得到的参数估计，那么接下来两步可以归纳为：
{% note default no-icon %}

1. E-步：通过

$$
\begin{equation*}
    \boldsymbol{T}^{(p)}=\mathbb{E}[\boldsymbol{T}(\boldsymbol{X})|\boldsymbol{y},\boldsymbol{\phi}^{(p)}]
\end{equation*}
$$

得到关于完全数据充分统计量的估计；

2. M-步：通过方程

$$
\begin{equation*}
    \mathbb{E}[\boldsymbol{T}(\boldsymbol{X})|\boldsymbol{\phi}]=\boldsymbol{T}^{(p)}
\end{equation*}
$$

确定更新后的参数 $ \boldsymbol{\phi}^{(p+1)} $.
{% endnote %}

这里的 M 步可能有些让人匪夷所思，似乎和极大化似然函数没有关系. 但实际上，注意我们现在有指数族的形式，我们的等价目标是极大化 $ \ell(\boldsymbol{\phi})=\log f(\boldsymbol{x}|\boldsymbol{\phi}) $. 通过对参数 $ \boldsymbol{\phi} $ 求偏导数就能得到新的参数 $ \boldsymbol{\phi}^{(p+1)} $ 需要满足

$$
\begin{equation*}
    \frac{\partial\log a(\boldsymbol{\phi})}{\partial\boldsymbol{\phi}}=\boldsymbol{T}(\boldsymbol{x})=\boldsymbol{T}^{(p)}.
\end{equation*}
$$

而概率密度的积分必须为 $ 1 $，因此我们有如下的等式：

$$
\begin{align*}
    \frac{\partial}{\partial\boldsymbol{\phi}}\int f(\boldsymbol{x}|\boldsymbol{\phi})\mathrm{d}\boldsymbol{x}&=\frac{1}{a^2(\boldsymbol{\phi})}\int\left[a(\boldsymbol{\phi})\boldsymbol{T}(\boldsymbol{x})-\frac{\partial a(\boldsymbol{\phi})}{\partial\boldsymbol{\phi}}\right]b(\boldsymbol{x})\exp(\boldsymbol{\phi}'\boldsymbol{T}(\boldsymbol{x}))\mathrm{d}\boldsymbol{x}\\
    &=\mathbb{E}[\boldsymbol{T}(\boldsymbol{X})|\boldsymbol{\phi}]-\frac{\partial\log a(\boldsymbol{\phi})}{\partial\boldsymbol{\phi}}=0.
\end{align*}
$$

于是我们确定了对数配分函数的偏导数. 注意我们现在更新得到的充分统计量是 $ \boldsymbol{T}^{(p)} $，因此就得到了 M-步需要满足的等式. 可以看到，从形式上来看上我们不需要知道 $ \boldsymbol{x}\in\mathcal{X} $ 的具体数值，只要知道统计量的数值即可. 当然在对统计量的实际计算时，还是可能需要确定 $ \boldsymbol{x} $ 的值来计算. 对于 M-步而言，一个困难是对应的方程可能不可解，此时极大值点落在 $ \Omega $ 的边界上，这种情况我们留到更一般的情形讨论. 但如果可解，那么由凸函数及正则指数族的性质可知极大值点是唯一的.
在进行一般情形的讨论前，首先我们说明 EM 算法确实能够让参数估计收敛到 $ L(\boldsymbol{\phi})=\log g(\boldsymbol{y}|\boldsymbol{\phi}) $ 的 MLE $ \boldsymbol{\phi}^* $. 记给定 $ \boldsymbol{y},\boldsymbol{\phi} $ 时 $ \boldsymbol{x} $ 的条件密度为

$$
\begin{equation*}
    k(\boldsymbol{x}|\boldsymbol{y},\boldsymbol{\phi})=\frac{f(\boldsymbol{x}|\boldsymbol{\phi})}{g(\boldsymbol{y}|\boldsymbol{\phi})}\Rightarrow L(\boldsymbol{\phi})=\log f(\boldsymbol{x}|\boldsymbol{\phi})-\log k(\boldsymbol{x}|\boldsymbol{y},\boldsymbol{\phi})
\end{equation*}
$$

对于指数族，利用 $ f $ 和 $ g $ 的积分联系，我们可以将 $ k $ 写作

$$
\begin{equation*}
    k(\boldsymbol{x}|\boldsymbol{y},\boldsymbol{\phi})=b(\boldsymbol{x})\exp[\boldsymbol{\phi}'\boldsymbol{T}(\boldsymbol{x})]/a(\boldsymbol{\phi}|\boldsymbol{y}),
\end{equation*}
$$

其中

$$
\begin{equation*}
    a(\boldsymbol{\phi}|\boldsymbol{y})=\int_{\mathcal{X}(\boldsymbol{y})}b(\boldsymbol{x})\exp[\boldsymbol{\phi}'\boldsymbol{T}(\boldsymbol{x})]\mathrm{d}\boldsymbol{x}.
\end{equation*}
$$

利用概率密度积分恒等式马上可以得到

$$
\begin{equation*}
    a(\boldsymbol{\phi})=\int_{\mathcal{X}}b(\boldsymbol{x})\exp[\boldsymbol{\phi}'\boldsymbol{T}(\boldsymbol{x})]\mathrm{d}\boldsymbol{x},
\end{equation*}
$$

于是

$$
\begin{equation*}
    L(\boldsymbol{\phi})=-\log a(\boldsymbol{\phi})+\log a(\boldsymbol{\phi}|\boldsymbol{y}).
\end{equation*}
$$

利用对数配分函数偏导数的结果，简记 $ \boldsymbol{T}=\boldsymbol{T}(\boldsymbol{X}) $， 我们可以得到如下极为重要的形式：

$$
\begin{equation*}
    \frac{\partial L(\boldsymbol{\phi})}{\partial{\boldsymbol{\phi}}}:=\mathbf{D}L(\boldsymbol{\phi})=-\mathbb{E}[\boldsymbol{T}|\boldsymbol{\phi}]+\mathbb{E}[\boldsymbol{T}|\boldsymbol{y},\boldsymbol{\phi}].
\end{equation*}
$$

**该式是理解 EM 算法的关键.** 如果算法收敛到某个参数 $ \boldsymbol{\phi}^{\ast} $，那么有 $ \boldsymbol{\phi}^{(p)}=\boldsymbol{\phi}^{(p+1)}=\boldsymbol{\phi}^{\ast},p\to\infty $. 结合之前提及的关于指数族的 EM 算法更新步骤可知，有 $ \mathbb{E}[\boldsymbol{T}|\boldsymbol{\phi}^{\ast}]=\mathbb{E}[\boldsymbol{T}|\boldsymbol{y},\boldsymbol{\phi}^{\ast}]\Leftrightarrow\mathbf{D}L(\boldsymbol{\phi})|_{\boldsymbol{\phi}=\boldsymbol{\phi^{\ast}}}=0 $，这说明 $ \boldsymbol{\phi^{\ast}} $ 就是 MLE.
Sundberg（1974）通过对上面的配分函数进行不断求偏导，得到了两组比较显然但有重要意义的式子：

$$
\begin{equation*}
    \begin{cases}
        \mathbf{D}^ka(\boldsymbol{\phi})=a(\boldsymbol{\phi})\mathbb{E}[\boldsymbol{T}^k|\boldsymbol{\phi}]\\
        \mathbf{D}^ka(\boldsymbol{\phi}|\boldsymbol{y})=a(\boldsymbol{\phi}|\boldsymbol{y})\mathbb{E}[\boldsymbol{T}^k|\boldsymbol{y},\boldsymbol{\phi}]
    \end{cases}\quad
     \begin{cases}
        \mathbf{D}^k\log a(\boldsymbol{\phi})=\mathbf{K}^k(\boldsymbol{T}|\boldsymbol{\phi})\\
        \mathbf{D}^k\log a(\boldsymbol{\phi}|\boldsymbol{y})=\mathbf{K}^k(\boldsymbol{T}|\boldsymbol{y},\boldsymbol{\phi})
    \end{cases}
\end{equation*}
$$

注意这里 $ \mathbf{D}^k $ 表示的是 $ k $ 阶的所有偏导数构成的**张量**，即一个 $ k $ 维数组. 相应地 $ \boldsymbol{T}^k=\boldsymbol{T}\otimes\boldsymbol{T}\otimes\cdots\otimes\boldsymbol{T} $ 表示由 $ \boldsymbol{T} $ 产生的 $ k $ 阶**张量**. 通过积分号下的求导我们可以得到充分统计量的 $ k $ 阶原点矩. $ \mathbf{K} $ 是叫做**累积量**的一个东西，类似于高阶特征. 当 $ k=1 $ 时就是期望， $ k=2 $ 是就是协方差阵，依次类推. 利用累积量，我们还可以直接导出

$$
\begin{equation*}
    \mathbf{D}^k L(\boldsymbol{\phi})=-\mathbf{K}^k(\boldsymbol{T}|\boldsymbol{\phi})+\mathbf{K}^k(\boldsymbol{T}|\boldsymbol{y},\boldsymbol{\phi}).
\end{equation*}
$$

该式是说，似然函数的任意阶偏导数都可以通过将充分统计量的条件累积量和无条件累积量作差得到.
现在我们假定完全数据的密度函数族不再是正则指数族，而是**曲指数族**. 这是说， $ \boldsymbol{\phi} $ 包含在 $ r $ 维凸区域 $ \Omega $ 的一个曲子流形 $ \Omega_0 $ 中. 比如一维正态分布，我们假定位置参数和尺度参数为 $ (\mu,\mu^2) $，此时参数空间就从凸区域 $ \mathbb{R}\times\mathbb{R}_+ $ 退化到一条不含原点的抛物线. 这时我们需要将 M-步改写为如下形式：

{% note default no-icon %}

2. M-步：令 $ \boldsymbol{\phi}^{(p+1)} $ 为使得 $ -\log a(\boldsymbol{\phi})+\boldsymbol{\phi}'\boldsymbol{T}^{(p)} $ 在 $ \Omega_0 $ 上极大化的点 $ \boldsymbol{\phi} $.

{% endnote %}

### 一般形式

接下来我们去掉指数族的条件，定义

$$
\begin{equation*}
    Q(\boldsymbol{\phi}'|\boldsymbol{\phi})=\mathbb{E}[\log f(\boldsymbol{X}|\boldsymbol{\phi}')|\boldsymbol{y},\boldsymbol{\phi}],
\end{equation*}
$$

并假设对任意 $ \boldsymbol{\phi}\in\Omega $ 和几乎处处的 $ \boldsymbol{x}\in\mathcal{X} $ 均有 $ f(\boldsymbol{x}|\boldsymbol{\phi})>0 $. EM 算法可以阐述如下：

{% note primary no-icon %}

1. E-步：计算 $ Q(\boldsymbol{\phi}|\boldsymbol{\phi}^{(p)}) $；

2. M-步：选择 $ \boldsymbol{\phi}^{(p+1)} $ 为使得 $ Q(\boldsymbol{\phi}|\boldsymbol{\phi}^{(p)}) $ 在参数空间 $ \Omega $ 上极大化的 $ \boldsymbol{\phi} $.

**注意，这里的 $ \boldsymbol{\phi}^{(p+1)} $ 的选择可以是不唯一的，从极大化的过程定义了一个点到集合的映射. 后续在算法收敛性讨论时我们会进一步说明.**
{% endnote %}

这样做的经验想法是，我们要找 MLE $ \boldsymbol{\phi}^* $，但是我们并不清楚完全数据的密度函数，因此需要对在给定条件 $ \boldsymbol{y},\boldsymbol{\phi}^{(p)} $ 下的期望似然做极大化. 容易观察到指数族的 EM 算法只是上述一般化 EM 算法步骤的一个特例，其计算更为简单方便. 对于上述一般情形，我们需要对全空间的 $ \boldsymbol{\phi} $ 计算 $ Q(\boldsymbol{\phi}|\boldsymbol{\phi^{(p)}}) $ 来寻找下一个极大值点，而对于指数族，由于我们已经导出了简单形式，只要计算充分统计量的期望值 $ \boldsymbol{T}^{(p)} $ 就可以直接求解出下一个参数点.

### 后验形式

EM 算法的思想也可以直接迁移到贝叶斯统计中. 我们知道贝叶斯定理的形式是

$$
\begin{equation*}
    p(\boldsymbol{\phi}|\boldsymbol{x})=\frac{p(\boldsymbol{x}|\boldsymbol{\phi})\cdot p(\boldsymbol{\phi})}{p(\boldsymbol{x})}\propto p(\boldsymbol{x}|\boldsymbol{\phi})\cdot p(\boldsymbol{\phi}),
\end{equation*}
$$

这里为方便起见对概率密度函数的符号不做区分，释义是后验概率正比于对参数的先验乘上对事件的先验认知. 将其直接迁移到 EM 算法中，我们想要极大化后验概率，等价于极大化右边的分子，而由于我们并不清楚完全数据的概率分布，但知道给定不完全数据和参数下的条件分布，于是问题等价于极大化 $ Q(\boldsymbol{\phi}|\boldsymbol{\phi}^{(p)})+\log p(\boldsymbol{\phi})=Q(\boldsymbol{\phi}|\boldsymbol{\phi}^{(p)})+G(\boldsymbol{\phi}) $. 如果 $ G(\boldsymbol{\phi}) $ 是共轭先验族中的对数密度，那么前面的加和形式会和 $ Q(\boldsymbol{\phi}|\boldsymbol{\phi}^{(p)}) $ 有相似的函数形式，这方便了最大化问题的求解.

## EM 算法的一般理论性质

首先引入函数

$$
\begin{equation*}
    H(\boldsymbol{\phi}'|\boldsymbol{\phi})=\mathbb{E}[\log k(\boldsymbol{X}|\boldsymbol{y},\boldsymbol{\phi}')|\boldsymbol{y},\boldsymbol{\phi}],
\end{equation*}
$$

那么根据前述定义，取条件期望就有下述等式：

$$
\begin{equation*}
    Q(\boldsymbol{\phi}'|\boldsymbol{\phi})=L(\boldsymbol{\phi}')+H(\boldsymbol{\phi}'|\boldsymbol{\phi}).
\end{equation*}
$$

### 单调性相关讨论

{% note info no-icon %}

**引理 1** 对于任意的 $ (\boldsymbol{\phi}',\boldsymbol{\phi})\in\Omega\times\Omega $, 都有

$$
\begin{equation*}
    H(\boldsymbol{\phi}'|\boldsymbol{\phi})\le H(\boldsymbol{\phi}|\boldsymbol{\phi}),
\end{equation*}
$$

等号当且仅当 $ k(\boldsymbol{x}|\boldsymbol{y},\boldsymbol{\phi}')=k(\boldsymbol{x}|\boldsymbol{y},\boldsymbol{\phi}) $ 几乎处处成立.

{% endnote %}

这个不等式的证明实际上直接依据 K-L 散度的非负性即可，本质上是由 Jensen 不等式和对数函数的凹性导出的，这里不再赘述.
为了契合“迭代算法”这一框架，我们可以通过一个 $ \Omega\to\Omega $ 的映射法则 $ \mathbf{M} $ 来定义参数每一步的变化，即 $ \boldsymbol{\phi}^{p}\to\boldsymbol{\phi}^{(p+1)}\Leftrightarrow\boldsymbol{\phi}^{(p+1)}=\mathbf{M}(\boldsymbol{\phi}^{(p)}) $. 由此我们可以定义**广义 EM 算法**：

{% note default no-icon %}

**定义** 一个配备了映射 $ \mathbf{M}(\boldsymbol{\phi}) $ 的迭代算法被称为**广义 EM 算法（GEM 算法）**，如果

$$
\begin{equation*}
    Q(\mathbf{M}(\boldsymbol{\phi})|\boldsymbol{\phi})\ge Q(\boldsymbol{\phi}|\boldsymbol{\phi})
\end{equation*}
$$

对所有 $ \boldsymbol{\phi}\in\Omega $ 都成立.

{% endnote %}

可以注意到，先前的 EM 算法和 GEM 算法的定义是相容的，因为我们要在旧参数的基础上找到使得 $ Q $ 全局极大化的新参数. 我们有如下定理：

{% note info no-icon %}

**定理 1** 对所有的 GEM 算法，均有

$$
\begin{equation*}
    L(\mathbf{M}(\boldsymbol{\phi}))\ge L(\boldsymbol{\phi}),\quad\forall\boldsymbol{\phi}\in\Omega.
\end{equation*}
$$

等号成立当且仅当

$$
\begin{equation*}
    Q(\mathbf{M}(\boldsymbol{\phi})|\boldsymbol{\phi})=Q(\boldsymbol{\phi}|\boldsymbol{\phi}),\quad k(\boldsymbol{x}|\boldsymbol{y},\mathbf{M}(\boldsymbol{\phi}))=k(\boldsymbol{x}|\boldsymbol{y},\boldsymbol{\phi})
\end{equation*}
$$

同时几乎处处成立.

{% endnote %}

**证明** 注意等式

$$
\begin{equation*}
    L(\mathbf{M}(\boldsymbol{\phi}))-L(\boldsymbol{\phi})=\lbrace Q(\mathbf{M}(\boldsymbol{\phi})|\boldsymbol{\phi})-Q(\boldsymbol{\phi}|\boldsymbol{\phi})\rbrace+\lbrace H(\boldsymbol{\phi}|\boldsymbol{\phi})-H(\mathbf{M}(\boldsymbol{\phi})|\boldsymbol{\phi})\rbrace,
\end{equation*}
$$

根据引理 1 和 EM 算法的极大化要求即得其非负，取等条件结合引理 1 即得.

{% note warning no-icon %}

**实际上原文作者在定义映射 $ \mathbf{M} $ 时犯了一个错误，在寻找极大值点时，极大值点不一定是唯一的. 在先前正则指数族与参数空间为凸集的情况下，确实能够保证映射找到的点时唯一的，然而在更一般的情形下并不是这样.** 事实上，我们需要定义 $ \mathbf{M} $ 为**点到集合的映射**，具体而言，$ \mathbf{M} $ 将一个空间 $ \mathcal{X} $ 中的点映射成该空间中的一个集合. **后文若没有特殊说明，都采用这种理解方式.** 不过在修改定义后，对定理 1 和下面两个推论是没有影响的，结论仍然成立.

{% endnote %}

由定理 1 我们可以马上导出如下两个推论：

{% note info no-icon %}

**推论 1** 如果对某些 $ \boldsymbol{\phi}^{\ast}\in\Omega,L(\boldsymbol{\phi}^{\ast})\ge L(\boldsymbol{\phi}) $ 对任意 $ \boldsymbol{\phi}\in\Omega $ 成立. 那么对于任意的 GEM 算法，有

(a) $ L(\mathbf{M}(\boldsymbol{\phi}^{\ast}))=L(\boldsymbol{\phi}^{\ast}) $；
(b) $ Q(\mathbf{M}(\boldsymbol{\phi}^{\ast})|\boldsymbol{\phi}^{\ast})=Q(\boldsymbol{\phi}^{\ast}|\boldsymbol{\phi}^{\ast})$ 几乎处处成立；
(c) $ k(\boldsymbol{x}|\boldsymbol{y},\mathbf{M}(\boldsymbol{\phi}^{\ast}))=k(\boldsymbol{x}|\boldsymbol{y},\boldsymbol{\phi}^{\ast}) $ 几乎处处成立.

**推论 2** 如果对某个 $ \boldsymbol{\phi}^{\ast}\in\Omega, L(\boldsymbol{\phi}^{\ast})>L(\boldsymbol{\phi}) $ 对任意的 $ \boldsymbol{\phi}\ne\boldsymbol{\phi}^{\ast},\boldsymbol{\phi}\in\Omega $ 均成立，那么对任何 GEM 算法，均有 $ \mathbf{M}(\boldsymbol{\phi}^*)=\boldsymbol{\phi}^{\ast} $.

{% endnote %}

上述两个推论实际类似于极大似然估计的不变性和在严格极大值下的唯一性，推论 2 更强的条件保证了映射不动点唯一.

### 收敛性相关讨论

**接下来的定理 2 原文中的证明是有问题的，这里本人对照 Wu（1983）的文章进行了修正，在逻辑上保证通顺**. 这些定理主要是为了说明 EM 算法的收敛性质. 方便起见，我们直接说明需要增补的一些条件：

{% note default no-icon %}

I. $ \Omega_{\boldsymbol{\phi}_0}=\lbrace\boldsymbol{\phi}\in\Omega:L(\boldsymbol{\phi})\ge L(\boldsymbol{\phi}_0) \rbrace $ 对任意的 $ L(\boldsymbol{\phi}_0)>-\infty $ 是紧的；

II. $ L $ 在 $ \Omega $ 上连续，且在 $ \Omega $ 的内部可微；

III. 对于 $ \boldsymbol{\phi}_{0}\in\Omega $， $ \Omega\_{\boldsymbol{\phi}_0} $ 包含在 $ \Omega $ 的内部.

{% endnote %}

由于 Fermat 引理求得的驻点可以是局部极大值、鞍点以及全局最大值，因此我们需要上述更强的条件来支撑证明. 在正式说明 EM 算法的收敛性之前，我们不加证明地给出如下定理（详情见{% link Zangwill 1969 https://pubsonline.informs.org/doi/pdf/10.1287/mnsc.16.1.1?download=true %}）

{% note info no-icon %}

**全局收敛定理** 令序列 $ \lbrace x_{k}\rbrace_{k=0}^{\infty} $ 是由 $ x_{k+1}\in\mathbf{M}(x_k) $ 产生的，其中 $ \mathbf{M} $ 是 $ X $ 上的点到集合映射. 令 $ \Gamma\subset X $ 是给定的解集合，并假设 (i) 所有的点 $ x_k $ 都含于紧集 $ S\subset X $；(ii) $ \mathbf{M} $ 在 $ \Gamma $ 的补集上是闭的；(iii) 存在 $ X $ 上的连续函数 $ \alpha $ 使得 (a) 如果 $ x\notin\Gamma $，那么 $ \alpha(y)>\alpha(x),\forall y\in\mathbf{M}(x) $，以及 (b) 如果 $ x\in\Gamma $，那么 $ \alpha(y)\ge\alpha(x),\forall y\in\mathbf{M}(x) $.
那么所有的 $ \lbrace x_k\rbrace $ 的极限点都属于 $ \Gamma $，并且对于某些 $ x\in\Gamma $， $ \alpha(x_k) $ 单调收敛到 $ \alpha(x) $.
{% endnote %}

这里简单解释一下所需要的几个条件的含义：条件 (i) 保证了参数能够有收敛的极限点而不跑到无穷远处，条件 (ii) 保证了在迭代进行足够多步时，参数的位置不发生振荡或漂移，条件 (iii) 保证了迭代不会陷入死循环，保证算法始终能朝着最优值靠近.
现在我们令 $ \mathbf{M} $ 是 GEM 中的点到集合映射，并且 $ \alpha $ 是似然函数 $ L $，并选取集合 $ \Gamma $ 为 $ \mathscr{M} $（ $ \Omega $ 内部的所有局部极大值点） 或 $ \mathscr{I} $（ $ \Omega $ 内部的所有驻点），那么我们有如下定理：

{% note info no-icon %}

**定理 2.1** 令 $ \lbrace \boldsymbol{\phi}_{p}\rbrace $ 是由 $ \boldsymbol{\phi}\_{p+1}\in\mathbf{M}(\boldsymbol{\phi}\_p) $ 产生的 GEM 序列，并假设 (i) $ \mathbf{M} $ 是一个在 $ \mathscr{I}(\mathscr{M}) $ 的补集上闭的点到集合映射，(ii) 对任意的 $ \boldsymbol{\phi}\_{p}\notin\mathscr{I}(\mathscr{M}) $，都有 $ L(\boldsymbol{\phi}\_{p+1})>L(\boldsymbol{\phi}\_p) $.
那么序列 $ \lbrace\boldsymbol{\phi}\_p\rbrace $ 所有可能的极限点都是 $ L $ 的驻点（或局部极大值点），并且对于某些 $ \boldsymbol{\phi}^{\ast}\in\mathscr{I}(\mathscr{M}),L(\boldsymbol{\phi}\_p) $ 能够单调收敛到 $ L^*=L(\boldsymbol{\phi}^{\ast}) $.

{% endnote %}

对于 EM 算法，一个使 $ \mathbf{M} $ 是闭映射的简单的充分条件是

$$
\begin{equation*}
    Q(\boldsymbol{\phi}|\boldsymbol{\psi})\text{ 关于 }\boldsymbol{\phi},\boldsymbol{\psi}\text{ 均连续}.
\end{equation*}
$$

由此，我们有如下定理：

{% note info no-icon %}

**定理 2.2** 假设 $ Q $ 满足上述连续性假设，那么由 EM 算法生成的任何迭代序列 $ \lbrace\boldsymbol{\phi}\_p\rbrace $ 的所有极限点都是 $ L $ 的驻点，并且对于某些驻点 $ \boldsymbol{\phi}^{\ast},L(\boldsymbol{\phi}\_p) $ 能够单调收敛到 $ L^*=L(\boldsymbol{\phi}^{\ast}) $.

{% endnote %}

定理证明从略. 上述两个定理只能说明函数值能够收敛到驻点的函数值，而还无法说明参数点列的收敛性. Wu 的文章中指出原文的证明没有很直接的方式补救，相应地给出了能使得 $ \boldsymbol{\phi}_p\to\boldsymbol{\phi}^* $ 的一些条件，并由此出发推导了一些定理. 这里我们不加证明地给出最实用条件最容易验证的一个推论，具体相关定理及内容仍见 Wu 的文章：

{% note info no-icon %}

如果 $ L(\boldsymbol{\phi}) $ 在 $ \Omega $ 上是单峰的，仅有 $ \boldsymbol{\phi}^* $ 为其驻点，并且 $ \mathbf{D}^{10}Q(\boldsymbol{\phi}'|\boldsymbol{\phi}) $ 关于 $ \boldsymbol{\phi} $ 和 $ \boldsymbol{\phi}' $ 连续，那么任何 EM 序列 $ \lbrace\boldsymbol{\phi}_p\rbrace $ 都收敛到 $ L(\boldsymbol{\phi}) $ 的唯一极大值点 $ \boldsymbol{\phi}^* $.

{% endnote %}

**注意定理中 $ \mathbf{D}^{10} $ 符号表示对第一个变量求梯度，而对第二个变量不做操作，后文中的 $ \mathbf{D}^{20} $ 就表示第一个变量的 Hessian 矩阵，其他微分符号以此类推.** 我们在后文都对 $ L,Q,H $ 做出一阶及二阶偏导数连续的假定进行证明. 在收敛性条件的保证下，我们可以导出如下中值定理以及收敛速率的结论.

{% note info no-icon %}

**定理 3** 假设 $ \lbrace\boldsymbol{\phi}^{(p)}\rbrace $ 是一 GEM 序列且满足 $ \mathbf{D}^{10}Q(\boldsymbol{\phi}^{(p+1)}|\boldsymbol{\phi}^{(p)})=\boldsymbol{0} $，那么有如下中值定理成立：

$$
\begin{equation*}
    Q(\boldsymbol{\phi}^{(p+1)}|\boldsymbol{\phi}^{(p)})-Q(\boldsymbol{\phi}^{(p)}|\boldsymbol{\phi}^{(p)})=-(\boldsymbol{\phi}^{(p+1)}-\boldsymbol{\phi}^{(p)})'\mathbf{D}^{20}Q(\boldsymbol{\phi}_0^{(p+1)}|\boldsymbol{\phi}^{(p)})(\boldsymbol{\phi}^{(p+1)}-\boldsymbol{\phi}^{(p)}),
\end{equation*}
$$

其中 $ \boldsymbol{\phi}_0^{(p+1)} $ 位于连接 $ \boldsymbol{\phi}^{(p)} $ 和 $ \boldsymbol{\phi}^{(p+1)} $ 的线段上.

**定理 4** 假设 $ \lbrace\boldsymbol{\phi}^{(p)}\rbrace $ 是一 GEM 序列，使得
(1) $ \boldsymbol{\phi}^{(p)}\to\boldsymbol{\phi}^{\ast}\in\overline{\Omega},p\to\infty $；
(2) $ \mathbf{D}^{10}Q(\boldsymbol{\phi}^{(p+1)}|\boldsymbol{\phi}^{(p)})=\boldsymbol{0} $，

那么有

$$
\begin{equation*}
    \mathbf{D}L(\boldsymbol{\phi}^{\ast})=\boldsymbol{0},\mathbf{D}\mathbf{M}(\boldsymbol{\phi}^{\ast})=\mathbf{D}^{20}H(\boldsymbol{\phi}^{\ast}|\boldsymbol{\phi}^{\ast})\left[\mathbf{D}^{20}Q(\boldsymbol{\phi}^{\ast}|\boldsymbol{\phi}^{\ast})\right]^{-1}.
\end{equation*}
$$

{% endnote %}

这里需要做一些符号上的解释. $ \mathbf{DM}(\boldsymbol{\phi}^{\ast}) $ 代表的是映射 $ \mathbf{M} $ 在 $ \boldsymbol{\phi}^{\ast} $ 处的 Jacobi 矩阵，如果其所有特征值都在 $ [0,1) $ 内，其最大特征值越小，意味着算法的收敛速度越快，越靠近 $ 1 $ 则意味着收敛越慢. 定理 3 的推导非常容易，定理 4 的推导使用了得分函数及 Fisher 信息量相关的经典结论以及 GEM 算法的额外假定，并不十分困难，这里我们略去证明. 此外，上述两个定理本人删去了涉及原文错误收敛性条件的结论. 定理 4 其实反映了统计学中非常重要的一个原理：**缺失信息原理**. 当我们拥有的信息多而缺失部分少时，相对应地算法收敛会更快，因为我们不需要花费太多时间采用某种手段补全缺失的信息；否则，我们很难准确推断出缺失信息的真实情况，在算法上表现为无法收敛.
针对本文开头给出的例子，我们可以给出算法的收敛速度. 参考 E 步后给出的密度函数以及对数似然函数，经过一些代数变形有

$$
\begin{align*}
    \mathbf{D}^{20}Q(\pi'|\pi)&=-\lbrace\mathbb{E}[X_2|\boldsymbol{y},\pi]+y_4\rbrace/\pi'^2-(y_2+y_3)/(1-\pi')^2,\\
    \mathbf{D}^{20}H(\pi'|\pi)&=-\mathbb{E}[X_2|\boldsymbol{y},\pi]/\pi'^2+y_1/(2+\pi')^2.
\end{align*}
$$

将 $ \pi,\pi' $ 换为之前求出的 $ \hat{\pi} $，经过计算可知 $ \mathbf{DM}(\hat{\pi})\approx0.132778 $，先前的代码输出结果也印证了这一点. 原文还给出了让 EM 算法更快收敛或者运算更高效的方法，这里不展开叙述.

## 应用

EM 算法作为一个一般性算法，在不同的统计场景下有极为丰富的应用. 原始文献大多以文字描述，因此阅读起来并不方便，并且缺乏应用性，所以本文并不打算对原文进行复述. 这里我们通过更加具体的例子来展示 EM 算法到底怎么使用，并给出相应的 `python` 代码.

### 高斯混合模型（GMM）

假设我们现在有 $ K $ 类服从多元正态分布的总体，维度为 $ D $，每类服从正态分布 $ \mathcal{N}(\boldsymbol{\mu}_k,\boldsymbol{\Sigma}_k)=(2\pi)^{-D/2}|\boldsymbol{\Sigma}_k|^{-1/2}\exp\left[-\dfrac{1}{2}(\boldsymbol{x}-\boldsymbol{\mu}_k)'\boldsymbol{\Sigma}^{-1}_k(\boldsymbol{x}-\boldsymbol{\mu}_k)\right] $ （假设协方差阵均正定）以及观测到的数据集 $ \boldsymbol{x}_1,\cdots,\boldsymbol{x}_N $，但是我们不知道这些数据到底属于哪一类，即数据服从一个混合正态分布

$$
\begin{equation*}
    g(\boldsymbol{x}|\boldsymbol{\theta})=\sum_{k=1}^K\pi_k\cdot\mathcal{N}(\boldsymbol{x}|\boldsymbol{\mu}_k,\boldsymbol{\Sigma}_k),
\end{equation*}
$$

可以注意到参数 $ \boldsymbol{\theta}=(\pi_k,\boldsymbol{\mu}\_k,\boldsymbol{\Sigma}\_k),k=1,\cdots,K,\sum_{k=1}^K\pi_k=1 $. 直接写出数据的对数似然函数会发现，由于求和号的存在，求解变得异常困难甚至不可能. 但是如果我们知道每个数据具体属于哪一类，那么问题将会简单很多. 这时候就要构建一个完全数据，表示观测数据到底属于哪一类. 我们引入一个离散型随机变量 $ z_n $，其取 $ k $ 的概率为 $ \pi_k $，以及对应的服从两点分布的一个示性函数 $ \mathbb{I}(z_n=k):=z_{nk} $，那么这时候完全数据的似然函数可以写为

$$
\begin{equation*}
    L(\boldsymbol{\theta}|\boldsymbol{X},\boldsymbol{z})=\prod_{n=1}^Nf(\boldsymbol{x}_n,z_n|\boldsymbol{\theta})=\prod_{n=1}^N\prod_{k=1}^K[\pi_k\cdot\mathcal{N}(\boldsymbol{x}_n|\boldsymbol{\mu}_k,\boldsymbol{\Sigma}_k)]^{z_{nk}}.
\end{equation*}
$$

这时对数似然函数就有非常简洁的形式. 我们首先进行 E 步的操作：

$$
\begin{align*}
    Q(\boldsymbol{\theta}|\boldsymbol{\theta}^{(t)})&=\mathbb{E}_{\boldsymbol{Z}|\boldsymbol{X},\boldsymbol{\theta}^{(t)}}[\log L(\boldsymbol{\theta}|\boldsymbol{X},\boldsymbol{Z})|\boldsymbol{X},\boldsymbol{\theta}^{(t)}]\\
    &=\sum_{n=1}^N\sum_{k=1}^K\mathbb{E}[z_{nk}|\boldsymbol{x}_{n},\boldsymbol{\theta}^{(t)}]\cdot[\log\pi_k+\log\mathcal{N}(\boldsymbol{x}_n|\boldsymbol{\mu}_k,\boldsymbol{\Sigma}_k)]\\
    &=\sum_{n=1}^N\sum_{k=1}^K\mathbb{P}(z_{n}=k|\boldsymbol{x}_{n},\boldsymbol{\theta}^{(t)})\cdot[\log\pi_k+\log\mathcal{N}(\boldsymbol{x}_n|\boldsymbol{\mu}_k,\boldsymbol{\Sigma}_k)]\\
    &=\sum_{n=1}^N\sum_{k=1}^K\frac{f(\boldsymbol{x}_n,k|\boldsymbol{\theta}^{(t)})}{g(\boldsymbol{x}_n|\boldsymbol{\theta}^{(t)})}\cdot[\log\pi_k+\log\mathcal{N}(\boldsymbol{x}_n|\boldsymbol{\mu}_k,\boldsymbol{\Sigma}_k)]\\
    &:=\sum_{n=1}^N\sum_{k=1}^K\gamma_{nk}^{(t)}\cdot[\log\pi_k+\log\mathcal{N}(\boldsymbol{x}_n|\boldsymbol{\mu}_k,\boldsymbol{\Sigma}_k)].
\end{align*}
$$

**这里一定要分清楚参数的不同含义，求期望时用的是已有参数，注意符号的上下标.** 接下来进行 M 步的操作，现在我们有一个约束条件 $ \sum_{k=1}^K\pi_k=1 $，因此问题是带约束的，需要使用 Lagrange 乘子法. 写出 Lagrange 函数为

$$
\begin{equation*}
    F(\boldsymbol{\theta})=Q(\boldsymbol{\theta}|\boldsymbol{\theta}^{(t)})+\lambda\left(\sum_{k=1}^K\pi_k-1\right),
\end{equation*}
$$

对每类均值、方差、出现概率分别求导，可以得到

$$
\begin{align*}
    \frac{\partial F}{\partial \boldsymbol{\mu_k}}&=-\boldsymbol{\Sigma}_k^{-1}\left[\sum_{n=1}^N\gamma_{nk}^{(t)}(\boldsymbol{x}_n-\boldsymbol{\mu}_k)\right],\\
    \frac{\partial F}{\partial \boldsymbol{\Sigma}_k^{-1}}&=\frac{1}{2}\sum_{n=1}^N\gamma_{nk}^{(t)}\left[\boldsymbol{\Sigma}_k-(\boldsymbol{x}_n-\boldsymbol{\mu}_k)(\boldsymbol{x}_n-\boldsymbol{\mu}_k)'\right],\\
    \frac{\partial F}{\partial\pi_k}&=\sum_{n=1}^N\gamma_{nk}^{(t)}\pi_k^{-1}+\lambda,\\
    \dfrac{\partial F}{\partial\lambda}&=\sum_{k=1}^K\pi_k-1.
\end{align*}
$$

令偏导数全为零，可以解得

$$
\begin{align*}
    \boldsymbol{\mu}_k^{(t+1)}&=\frac{\sum_{n=1}^N\gamma_{nk}^{(t)}\boldsymbol{x}_n}{\sum_{n=1}^N\gamma_{nk}^{(t)}},\\
    \boldsymbol{\Sigma}_{k}^{(t+1)}&=\frac{\sum_{n=1}^N\gamma_{nk}^{(t)}(\boldsymbol{x}_n-\boldsymbol{\mu}_k^{(t+1)})(\boldsymbol{x}_n-\boldsymbol{\mu}_k^{(t+1)})'}{\sum_{n=1}^N\gamma_{nk}^{(t)}},\\
    \pi_k^{(t+1)}&=\frac{\sum_{n=1}^N\gamma_{nk}^{(t)}}{N}.
\end{align*}
$$

迭代至各类别稳定即可. GMM 模型实际上属于机器学习中的**无监督学习**，是用于聚类的. 实际中，只要给定输入数据和想要的聚类类别，算法最后就能输出对应类的均值、协方差阵以及出现概率. 下面给出模拟数据及 EM 算法聚类的 `python` 代码示例.

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal

# 生成模拟数据
np.random.seed(123)
N = 3000       # 数据总数
K_true = 3     # 真实的类别数
D = 2          # 数据维度 (2D 方便画图)

# 定义 3 个不同高斯分布的真实参数 (均值, 协方差, 混合权重)
means_true = np.array([
    [0, 4], [2 * 3**0.5, -2], [-2 * 3**0.5,-2]
])
covs_true = np.array([
    [[1, 0], [0, 1]], 
    [[4, 1.5], [1.5, 1]], 
    [[1, -0.8], [-0.8, 2]],
])
pi_true = np.array([0.2,0.5,0.3]) # 各自类别的概率

# 生成隐变量 z (真实标签) 和观测数据 X
z_true = np.random.choice(K_true, size=N, p=pi_true)
X = np.zeros((N, D))
for k in range(K_true):
    idx = (z_true == k)
    X[idx] = np.random.multivariate_normal(means_true[k], covs_true[k], size=np.sum(idx))


# EM 算法实现
K = 3 # 假设我们知道要聚成 3 类

# 1. 初始化参数
# 权重：平均分配
pi = np.ones(K) / K
# 均值：从数据中随机挑 K 个点作为初始中心
random_idx = np.random.choice(N, K, replace=False)
mu = X[random_idx].copy()
# 协方差：初始化为整个数据集的协方差矩阵
covariances = np.array([np.cov(X.T) for _ in range(K)])

max_iter = 200
tol = 1e-6
log_likelihoods = []

for it in range(max_iter):
    # E-步
    # weighted_pdfs 矩阵存储: pi_k * N(x_n | mu_k, Sigma_k)
    weighted_pdfs = np.zeros((N, K))
    for k in range(K):
        # 计算第 k 个高斯分布在所有点上的概率密度，并乘以其先验权重
        weighted_pdfs[:, k] = pi[k] * multivariate_normal.pdf(X, mean=mu[k], cov=covariances[k])
    
    # 每一行的和就是边缘概率 p(x_n)
    sum_pdfs = np.sum(weighted_pdfs, axis=1, keepdims=True)
    
    # 计算当前的对数似然函数并记录
    log_likelihood = np.sum(np.log(sum_pdfs))
    log_likelihoods.append(log_likelihood)
    
    # 判断是否收敛（如果似然函数增加极小，则停止）
    if it > 0 and (log_likelihood - log_likelihoods[-2]) < tol:
        print(f"算法在第 {it} 次迭代收敛！")
        break
        
    # 计算责任度 gamma: (N, K) 维矩阵
    gamma = weighted_pdfs / sum_pdfs

    # M-步
    # 计算第 k 类的“有效数据量” N_k
    N_k = np.sum(gamma, axis=0) # 形状: (K,)
    
    # 1. 更新混合权重 pi
    pi = N_k / N
    # 2. 更新均值 mu (加权平均)
    mu = (gamma.T @ X) / N_k[:, np.newaxis]
    # 3. 更新协方差阵
    for k in range(K):
        # 计算所有点到当前类均值的向量差: (N, D)
        diff = X - mu[k]
        # 对差值乘以责任度权重: (N, D)
        weighted_diff = diff * gamma[:, k:k+1]
        # 计算加权协方差矩阵，除以有效数据量
        covariances[k] = (weighted_diff.T @ diff) / N_k[k]
        
        # 在对角线上加一个极小的数，防止协方差矩阵由于数据太少变成奇异矩阵（不可逆）
        covariances[k] += np.eye(D) * 1e-6

# 打印算法估计得到的各参数
print("\n" + "="*45)
print("EM 算法估计的最终 GMM 参数报告：")
print("="*45)
for k in range(K):
    print(f"【类别 {k}】")
    print(f"出现概率 (pi) : {pi[k]:.4f}")
    print(f"样本均值 (mu) : [{mu[k, 0]:.4f}, {mu[k, 1]:.4f}]")
    print(f"协方差矩阵    :\n{covariances[k]}")
    print("-" * 45)

# 结果可视化
# 根据最终的责任度，把点分配给概率最大的那一类
z_pred = np.argmax(gamma, axis=1)

plt.figure(figsize=(15, 5))

# 图1：真实的分类情况
plt.subplot(1, 3, 1)
plt.scatter(X[:, 0], X[:, 1], c=z_true, cmap='viridis', s=10, alpha=0.7)
plt.title("True GMM Clusters")
plt.xlabel("X1"); plt.ylabel("X2")

# 图2：EM算法聚类的结果
plt.subplot(1, 3, 2)
plt.scatter(X[:, 0], X[:, 1], c=z_pred, cmap='viridis', s=10, alpha=0.7)
plt.scatter(mu[:, 0], mu[:, 1], c='red', marker='X', s=100, label='Est. Centers')
plt.title("EM Algorithm Predicted Clusters")
plt.xlabel("X1"); plt.ylabel("X2")
plt.legend()

# 图3：对数似然函数的单调递增曲线 (证明 Theorem 1)
plt.subplot(1, 3, 3)
plt.plot(log_likelihoods, marker='o', c='blue')
plt.title("Log-Likelihood Convergence")
plt.xlabel("Iteration")
plt.ylabel("Log-Likelihood")
plt.grid(True)

plt.tight_layout()
plt.show()

# 文字输出结果
算法在第 24 次迭代收敛！

=============================================
EM 算法估计的最终 GMM 参数报告：
=============================================
【类别 0】
出现概率 (pi) : 0.4944
样本均值 (mu) : [3.3666, -2.0111]
协方差矩阵    :
[[3.98001108 1.43869209]
 [1.43869209 0.95601346]]
---------------------------------------------
【类别 1】
出现概率 (pi) : 0.3016
样本均值 (mu) : [-3.5158, -2.0075]
协方差矩阵    :
[[ 1.04799398 -0.88332531]
 [-0.88332531  2.13381526]]
---------------------------------------------
【类别 2】
出现概率 (pi) : 0.2040
样本均值 (mu) : [-0.0149, 4.0394]
协方差矩阵    :
[[0.92494905 0.05520281]
 [0.05520281 0.99308669]]
---------------------------------------------
```

可以看到，EM 算法准确地将每一类进行了划分，并且给出了准确地出现概率. 由于是无监督学习，因此类别和原始设定的类会有所差异，但是我们并不关心类具体的顺序，只关心能否正确识别每一类.
GMM 模型还可以完全迁移到其他模型上，一个典型例子就是总共有 $ K $ 枚硬币，选择其中某枚的概率为 $ \pi_k,k=1,\cdots,K $，并且这些硬币间质地是不均的，现在进行 $ N $ 次实验，第 $ n $ 次实验抛掷 $ m_n $ 次选择的硬币，正面朝上的次数为 $ x_n $，我们可以用 EM 算法来估计出抛掷硬币类型被选中的概率和这枚硬币正面朝上的概率，具体操作仅仅需要将 GMM 模型中的正态分布换为 $ K $ 个二项分布再进行求解即可.

### 含有缺失数据的多元正态分布的恢复

假设我们有一组服从多元正态分布 $ \mathcal{N}(\boldsymbol{\mu},\boldsymbol{\Sigma}) $ 的独立同分布数据 $ \lbrace y_n\rbrace_{n=1}^N,y_n\in\mathbb{R}^D $，但是我们仅仅观测到了其中部分维度的数据（不妨认为就是前 $ d $ 维），记为 $ \boldsymbol{y}\_{n,o} $，而剩下的维度未观测到，记为 $ \boldsymbol{y}\_{n,m} $，那么潜变量就是未观测到的维度的数据. 现在我们想要估计完全的 $ \boldsymbol{\mu},\boldsymbol{\Sigma} $.
当 $ \lbrace y_n\rbrace $ 被完全观测到时，显然有对数似然函数

$$
\begin{equation*}
    \ell(\boldsymbol{\mu},\boldsymbol{\Sigma}|\boldsymbol{Y})=-\frac{N}{2}\log|\boldsymbol{\Sigma}|-\frac{1}{2}\sum_{n=1}^N(\boldsymbol{y}_n-\boldsymbol{\mu})'\boldsymbol{\Sigma}^{-1}(\boldsymbol{y}_n-\boldsymbol{\mu})+C.
\end{equation*}
$$

现在我们来构建 EM 算法的步骤. 对于 E 步，对原始数据每一条做同样的划分，即

$$
\begin{equation*}
    \boldsymbol{y}_n=\begin{pmatrix}
        \boldsymbol{y}_{n,o}\\
        \boldsymbol{y}_{n,m}
    \end{pmatrix},\boldsymbol{\mu}=\begin{pmatrix}
        \boldsymbol{\mu}_{o}\\
        \boldsymbol{\mu}_{m}
    \end{pmatrix},
    \boldsymbol{\Sigma}=\begin{pmatrix}
        \boldsymbol{\Sigma}_{oo}& \boldsymbol{\Sigma}_{om} \\
        \boldsymbol{\Sigma}_{mo} & \boldsymbol{\Sigma}_{mm}
    \end{pmatrix}.
\end{equation*}
$$

那么根据多元正态分布的条件期望公式，我们有

$$
\begin{align*}
    \mathbb{E}[\boldsymbol{y}_{n,m}|\boldsymbol{y}_{n,o}]&=\boldsymbol{\mu}_m+\boldsymbol{\Sigma}_{mo}\boldsymbol{\Sigma}_{oo}^{-1}(\boldsymbol{y}_{n,o}-\boldsymbol{\mu}_o),\\
    \text{Cov}(\boldsymbol{y}_{n,m}|\boldsymbol{y}_{n,o})&=\boldsymbol{\Sigma}_{mm}-\boldsymbol{\Sigma}_{mo}\boldsymbol{\Sigma}_{oo}^{-1}\boldsymbol{\Sigma}_{om}.
\end{align*}
$$

我们知道多元正态分布的充分统计量是样本均值和样本方差，与之相关联的两个重要统计量是 $ \boldsymbol{y}_n,\boldsymbol{y}_n\boldsymbol{y}_n' $，其能构建出对应的充分统计量. 利用协方差阵的恒等式可得

$$
\begin{equation*}
    \hat{\boldsymbol{y}}_n=\mathbb{E}[\boldsymbol{y}_{n,m}|\boldsymbol{y}_{n,o}],~~\mathbb{E}[\boldsymbol{y}_n\boldsymbol{y}_n'|\boldsymbol{y}_{n,o}]=\hat{\boldsymbol{y}}_n\hat{\boldsymbol{y}}_n'+\begin{pmatrix}
        \boldsymbol{O} & \boldsymbol{O} \\
        \boldsymbol{O} & \text{Cov}(\boldsymbol{y}_{n,m}|\boldsymbol{y}_{n,o})
    \end{pmatrix}.
\end{equation*}
$$

对于 M 步，我们要极大化对数似然函数，根据熟知的极大化相关结论（详见{% post_link 多元正态分布 "多元正态分布" %}一文），定义

$$
\begin{equation*}
    \tilde{S}_1=\sum_{n=1}^N\mathbb{E}[\boldsymbol{y}_n^{(t)}|\boldsymbol{y}_{n,o}],~~\tilde{S}_2=\sum_{n=1}^N\mathbb{E}[\boldsymbol{y}_n^{(t)}\boldsymbol{y}_{n}^{(t)}{'}|\boldsymbol{y}_{n,o}],
\end{equation*}
$$

那么更新步骤为

$$
\begin{align*}
    \boldsymbol{\mu}^{(t+1)}=\frac{1}{N}\tilde{S}_1,~~\boldsymbol{\Sigma}^{(t+1)}=\frac{1}{N}\tilde{S}_2-\boldsymbol{\mu}^{(t+1)}\boldsymbol{\mu}^{(t+1)}{'}.
\end{align*}
$$

注意协方差阵的计算是经过交叉项的化简的.**在设定迭代初始值时，注意观测部分的均值和协方差阵是可以直接计算得到的，其余部分可以进行随机初始化（也可以选择其他填充方式），经过迭代仍会收敛**. 具体代码这里不再给出.

{% note warning no-icon %}
说到这里，本人其实就产生了一个疑问：**EM 算法针对缺失数据的数字特征估计能做到如此精确吗？** 问了问 AI，本人认为这几点是比较对的，在满足如下条件时，EM 算法能够较精确地刻画缺失部分的数字特征：

1. 特征之间高度相关. 如果特征间几乎无关联，那么对于缺失数据会倾向于用均值填充，此时不能获取缺失数据足够有用的信息；
2. 缺失机制必须是“完全随机”或“条件随机”. 这是由著名统计学家，也是 EM 算法的提出者之一，Donald Rubin 所提出的. 这两个术语的意思实际上是对**缺失值产生的原因进行了刻画**. 大致意思是，数据的缺失如果是不慎遗失或者可以通过其他变量进行大致推断的，那么 EM 算法是非常有效的. 如果数据的缺失和缺失值的大小位置等有相关关系，那么 EM 算法就会失效.
3. 数据本身要符合多元正态分布. 如果数据间存在强烈的非线性关系，或者某一维度的数据本身有严重的偏倚，那么 EM 算法在建模时就已经失效了，更不用说算法本身.

{% endnote %}

### 隐马尔可夫模型（HMM）与 Baum-Welch 算法

隐马尔可夫模型是统计模型中经典的**概率图模型**. 用通俗的语言说，就是存在一条“明线”和“暗线”，观测就是明线，表示的是在每个时间点的一些已知的状态，而暗线就是潜变量，我们无法得知具体情况. 一个实际例子是，我知道我的朋友每天具体干了什么活动，我想要用这个序列来推测朋友所在城市每天的天气是怎么样的.
假设我们已经观测到了一序列数据 $ \boldsymbol{X}=\lbrace x_1,\cdots,x_T\rbrace $，现在我们要推断隐藏态 $ \boldsymbol{Z}=\lbrace z_1,\cdots,z_T\rbrace $ 在每一时刻的概率分布. HMM 模型可以抽象为一个五元组 $ (Q,V,\boldsymbol{\pi},\boldsymbol{A},\boldsymbol{B}) $，其中 $ Q=\lbrace q_1,\cdots,q_K\rbrace $ 为所有隐藏过程的可能状态， $ V=\lbrace v_1,\cdots,v_L \rbrace $. 简单起见，我们令 $ q_i=i,v_j=j,i=1,\cdots,K,j=1,\cdots,L. \boldsymbol{\pi}\in\mathbb{R}^K $ 表示初始状态概率向量，而 $ \boldsymbol{A}\in\mathbb{R}^{K\times K},\boldsymbol{B}\in\mathbb{R}^{K\times L} $ 分别是（隐藏状态的）转移概率矩阵和（从隐藏态到观测态的）发射矩阵. 这两个矩阵的定义是基于两个假设，即**隐藏状态的马尔可夫性**

$$
\begin{equation*}
    \mathbb{P}(z_t|z_{t-1},x_{t-1},\cdots,z_{1},x_{1})=\mathbb{P}(z_t|z_{t-1})
\end{equation*}
$$

和**观测独立性假设**

$$
\begin{equation*}
    \mathbb{P}(x_t|z_{T},x_{T},\cdots,z_t,\cdots,z_{1},x_{1})=\mathbb{P}(x_t|z_{t})
\end{equation*}
$$

只有这两个假设成立才能建立 HMM，并且建立后文的推导. 现在我们就需要对参数组 $ \boldsymbol{\theta}=(\boldsymbol{\pi},\boldsymbol{A},\boldsymbol{B}) $ 中的元素进行估计. 定义

$$
\begin{equation*}
    \pi_j=\mathbb{P}(z_1=j),A_{ij}=\mathbb{P}(z_{t+1}=j|z_t=i),B_j(x)=\mathbb{P}(x_{t}=x|z_t=j)
\end{equation*}
$$

来表示参数空间中的每一元，根据我们的假设以及条件概率公式，完全数据的似然函数及相应的对数似然可以写为

$$
\begin{align*}
    f(\boldsymbol{X},\boldsymbol{Z}|\boldsymbol{\theta})&=p(z_1)\cdot p(x_1|z_1)\cdot p(z_2|z_1,x_1)\cdots p(x_{T}|\boldsymbol{X}_{-T},\boldsymbol{Z})\\
    &=p(z_1)\cdot p(x_1|z_1)\cdot p(z_2|z_1)\cdots p(x_{T}|z_T)\\
    &=p(z_1)\prod_{t=1}^{T-1}p(z_{t+1}|z_t)\prod_{t=1}^Tp(x_t|z_t),\\
    \ell(\boldsymbol{\theta}|\boldsymbol{X},\boldsymbol{Z})&=\log\pi_{z_1}+\sum_{t=1}^{T-1}\log A_{z_{t},z_{t+1}}+\sum_{t=1}^T\log B_{z_t}(x_t).
\end{align*}
$$

这里为了式子的简便对符号不加以区分，但只要清楚这个序列是怎么生成的，就不会混淆符号的释义. 对数似然直接把乘号变为加号即可.
接下来我们需要进行 E 步. 观察似然函数的形式，后两个随机变量的期望需要计算两个后验概率

$$
\begin{align*}
    \gamma_t(i)&=\mathbb{P}(z_t=i|\boldsymbol{X},\boldsymbol{\theta}^{(n)}),\\
    \xi_t(i,j)&=\mathbb{P}(z_t=i,z_{t+1}=j|\boldsymbol{X},\boldsymbol{\theta}^{(n)}),
\end{align*}
$$

首先来看 $ \gamma_t(i) $. 由 Bayes 公式，我们有

$$
\begin{equation*}
    \gamma_t(i)=\frac{\mathbb{P}(z_{t}=i,\boldsymbol{X}|\boldsymbol{\theta}^{(n)})}{\mathbb{P}(\boldsymbol{X}|\boldsymbol{\theta}^{(n)})}=\frac{\mathbb{P}(z_{t}=i,\boldsymbol{X}|\boldsymbol{\theta}^{(n)})}{\sum_{k=1}^K\mathbb{P}(z_{t}=k,\boldsymbol{X}|\boldsymbol{\theta}^{(n)})}
\end{equation*}
$$

将分子的概率进行拆分，可以得到

$$
\begin{align*}
    \mathbb{P}(z_{t}=i,\boldsymbol{X}|\boldsymbol{\theta}^{(n)})&=\mathbb{P}(z_t=i,x_1,\cdots,x_t|\boldsymbol{\theta}^{(n)})\cdot\mathbb{P}(x_{t+1},\cdots,x_{T}|z_{t}=i,x_1,\cdots,x_t,\boldsymbol{\theta}^{(n)})\\
    &=\mathbb{P}(z_t=i,x_1,\cdots,x_t|\boldsymbol{\theta}^{(n)})\cdot\mathbb{P}(x_{t+1},\cdots,x_{T}|z_{t}=i,\boldsymbol{\theta}^{(n)}):=\alpha_t(i)\cdot\beta_t(i).
\end{align*}
$$

这把概率拆成了“过去”和未来两个部分. 利用 HMM 的假定，我们分别对两者建立递推式，称为前向与后向传播.

$$
\begin{align*}
    \alpha_{t+1}(j)&=\sum_{i=1}^K\mathbb{P}(z_{t+1}=j,z_t=i,x_1,\cdots,x_{t+1}|\boldsymbol{\theta}^{(n)})\\
    &=\sum_{i=1}^K\mathbb{P}(x_{t+1}|z_{t+1}=j,\boldsymbol{\theta}^{(n)})\cdot\mathbb{P}(z_{t+1}=j|z_t=i,\boldsymbol{\theta}^{(n)})\cdot\mathbb{P}(z_t=i,x_1,\cdots,x_{t}|\boldsymbol{\theta}^{(n)})\\
    &=\left[\sum_{i=1}^K \alpha_t(i)A_{ij}\right]B_{j}(x_{t+1}),\quad(\alpha_1(j)=\pi_j\cdot B_j(x_1))\\
    \beta_t(i)&=\sum_{j=1}^K\mathbb{P}(z_{t+1}=j,x_{t+1},\cdots,x_{T}|z_t=i,\boldsymbol{\theta}^{(n)})\\
    &=\sum_{j=1}^K\mathbb{P}(z_{t+1}=j|z_t=i,\boldsymbol{\theta}^{(n)})\cdot\mathbb{P}(x_{t+1}|z_{t+1}=j,\boldsymbol{\theta}^{(n)})\cdot\mathbb{P}(x_{t+2},\cdots,x_{T}|z_{t+1}=j,\boldsymbol{\theta}^{(n)})\\
    &=\sum_{j=1}^KA_{ij}B_{j}(x_{t+1})\beta_{t+1}(j).\quad(\beta_T(i)=1)
\end{align*}
$$

利用上述递推式，可得 $ \gamma_t(i)=\dfrac{\alpha_t(i)\beta_t(i)}{\sum_{i=1}^K\alpha_t(i)\beta_t(i)} $. 利用和 $ \gamma_t(i) $ 完全相同的方法，可以导出

$$
\begin{equation*}
    \xi_t(i,j)=\frac{\alpha_t(i)A_{ij}B_j(x_{t+1})\beta_{t+1}(j)}{\sum_{i',j'=1}^{K}\alpha_t(i')A_{i'j'}B_j'(x_{t+1})\beta_{t+1}(j')}.
\end{equation*}
$$

**注意上面的取期望过程都是在参数 $ \boldsymbol{\theta}^{(n)} $ 下进行的，因此实际上所有的更新量都需要带一个上标 $ (n) $，为了符号的简洁这里并未写出.**
算出了概率分布，接下来的 M 步更新就比较简单了. 我们要求解的是多个带约束的优化问题，类似于 GMM 的求法，通过 Lagrange 乘数法可以得到更新公式

$$
\begin{align*}
    \pi_k^{(n+1)}&=\gamma_1^{(n)}(k),\\
    A_{ij}^{(n+1)}&=\frac{\sum_{t=1}^{T-1}\xi_t^{(n)}(i,j)}{\sum_{t=1}^{T-1}\gamma_t^{(n)}(i)},\\
    B_{i}^{(n+1)}(v)&=\frac{\sum_{t=1}^T\gamma_t^{(n)}(i)\mathbb{I}(x_t=v)}{\sum_{t=1}^T\gamma_t^{(n)}(i)}.
\end{align*}
$$

上述算法通常被称作**Baum-Welch 算法**，其本质就是 EM 算法. 接下来给出相应的 `python` 代码，产生模拟数据并且测试算法的效果. **注意，由于序列较长，因此会导致概率的指数级别的衰减，在代码中有一个对 $ \alpha,\beta $ 的缩放技巧，且并不影响最终值.** 最后给出训练的结果对比以及似然函数值的变化.

```python
import numpy as np

class HMM:
    def __init__(self, K, M, seed=123):
        self.K = K  # 隐藏状态数
        self.M = M  # 观测符号数
        np.random.seed(seed)
        
        # 随机初始化参数，并确保行和为1
        self.pi = np.random.rand(K)
        self.pi /= self.pi.sum()
        
        self.A = np.random.rand(K, K)
        self.A /= self.A.sum(axis=1, keepdims=True)
        
        self.B = np.random.rand(K, M)
        self.B /= self.B.sum(axis=1, keepdims=True)

    def _forward(self, obs):
        T = len(obs)
        alpha = np.zeros((T, self.K))
        scales = np.zeros(T)
        
        # t = 0
        alpha[0] = self.pi * self.B[:, obs[0]]
        scales[0] = 1.0 / (np.sum(alpha[0]) + 1e-16)
        alpha[0] *= scales[0]
        
        # t = 1...T-1
        for t in range(1, T):
            alpha[t] = (alpha[t-1] @ self.A) * self.B[:, obs[t]]
            scales[t] = 1.0 / (np.sum(alpha[t]) + 1e-16)
            alpha[t] *= scales[t]
            
        return alpha, scales

    def _backward(self, obs, scales):
        T = len(obs)
        beta = np.zeros((T, self.K))
        
        # t = T-1
        beta[T-1] = np.ones(self.K) * scales[T-1]
        
        # t = T-2...0
        for t in range(T-2, -1, -1):
            beta[t] = (self.A @ (self.B[:, obs[t+1]] * beta[t+1])) * scales[t]
            
        return beta

    def fit(self, obs, max_iter=500, tol=1e-6):
        T = len(obs)
        ll_history = []
        
        for i in range(max_iter):
            # E-step
            alpha, scales = self._forward(obs)
            beta = self._backward(obs, scales)
            
            # 计算 gamma
            # 由于 alpha 和 beta 都已经缩放过，gamma 直接相乘归一化即可
            gamma = alpha * beta
            gamma /= gamma.sum(axis=1, keepdims=True)
            
            # 计算 xi: P(z_t=i, z_t+1=j | O)
            xi = np.zeros((T-1, self.K, self.K))
            for t in range(T-1):
                # 这里的计算需抵消 scales
                numer = alpha[t].reshape(-1, 1) * self.A * self.B[:, obs[t+1]] * beta[t+1]
                xi[t] = numer / (np.sum(numer) + 1e-16)
            
            # M-step
            # 更新 pi
            self.pi = gamma[0] / (np.sum(gamma[0]) + 1e-16)
            
            # 更新 A
            self.A = np.sum(xi, axis=0) / (np.sum(gamma[:-1], axis=0).reshape(-1, 1) + 1e-16)
            
            # 更新 B
            for v in range(self.M):
                mask = (obs == v)
                self.B[:, v] = np.sum(gamma[mask], axis=0) / (np.sum(gamma, axis=0) + 1e-16)
            
            # 计算当前对数似然
            log_likelihood = -np.sum(np.log(scales))
            ll_history.append(log_likelihood)
            
            if i > 0 and abs(ll_history[-1] - ll_history[-2]) < tol:
                print(f"Converged at iteration {i}")
                break
                
        return ll_history

# --- 测试部分 ---

# 1. 生成模拟数据
def generate_data(T, pi, A, B):
    K, M = B.shape
    states = [np.random.choice(K, p=pi)]
    obs = [np.random.choice(M, p=B[states[0]])]
    
    for _ in range(1, T):
        states.append(np.random.choice(K, p=A[states[-1]]))
        obs.append(np.random.choice(M, p=B[states[-1]]))
    return np.array(obs), np.array(states)

# 固定一个真实的参数 (3隐藏态, 4观测态)
true_pi = np.array([0.1, 0.6, 0.3])
true_A = np.array([
    [0.7, 0.2, 0.1],
    [0.1, 0.8, 0.1],
    [0.2, 0.3, 0.5]
])
true_B = np.array([
    [0.4, 0.3, 0.2, 0.1],
    [0.1, 0.1, 0.4, 0.4],
    [0.2, 0.5, 0.1, 0.2]
])

# 生成模拟序列
obs_seq, _ = generate_data(5000, true_pi, true_A, true_B)

# 2. 训练模型
model = HMM(K=3, M=4, seed=123)
history = model.fit(obs_seq, max_iter=500)

# 3. 报告结果
print("\n--- Final Results ---")
print("Log-Likelihood trend (first 5 and last 5):")
print(history[:5], "...", history[-5:])

np.set_printoptions(precision=3, suppress=True)
print("\nLearned Transition Matrix A:")
print(model.A)
print("\nTrue Transition Matrix A:")
print(true_A)

print("\nLearned Emission Matrix B:")
print(model.B)
print("\nTrue Emission Matrix B:")
print(true_B)

# 输出结果
--- Final Results ---
Log-Likelihood trend (first 5 and last 5):
[np.float64(-6992.780878169078), np.float64(-6870.557357936086), np.float64(-6866.9806435341), np.float64(-6862.821552993244), np.float64(-6858.055159476695)] ... 
[np.float64(-6807.073814597403), np.float64(-6807.073761859025), np.float64(-6807.073709180981), np.float64(-6807.073656561857), np.float64(-6807.07360400037)]

Learned Transition Matrix A:
[[0.463 0.287 0.25 ]
 [0.391 0.52  0.089]
 [0.05  0.145 0.804]]

True Transition Matrix A:
[[0.7 0.2 0.1]
 [0.1 0.8 0.1]
 [0.2 0.3 0.5]]

Learned Emission Matrix B:
[[0.164 0.029 0.254 0.553]
 [0.047 0.183 0.505 0.265]
 [0.344 0.362 0.171 0.123]]

True Emission Matrix B:
[[0.4 0.3 0.2 0.1]
 [0.1 0.1 0.4 0.4]
 [0.2 0.5 0.1 0.2]]
```

从结果上来看，训练的效果似乎并不好. **但要注意，这本质上还是无监督学习问题，最终的转移矩阵实际上可能是通过真实转移概率矩阵的合同变换得到的，转移矩阵的数值受到初始概率的影响会比较大. 因此在实际操作时，可以多次随机初始化，选择使得对数似然最大的那个解作为最终答案.** 对于本模拟例子，效果不好的另一个原因是发射概率矩阵在不同行为间的区分度不够明显，这就导致了多种可能性，进而导致似然函数可能有多个驻点，求解时很有可能是被吸引到局部极大值点上.

### 带有错分类标签的逻辑回归

考虑一个逻辑回归任务，我们的观测是带噪声的，即观测到的标签并不一定是真实标签. 于是真实标签就是潜变量. 逻辑回归要求标签只能在 $ \lbrace 0,1\rbrace $ 中取值，现在我们的任务是估计出真实的逻辑回归系数向量 $ \boldsymbol{w} $ 以及两个错分类概率 $ \alpha,\beta $，即参数向量 $ \boldsymbol{\theta}=(\boldsymbol{w},\alpha,\beta) $. 其中在完全数据下应该有

$$
\begin{align*}
    \mathbb{P}(t_n=1|\boldsymbol{x}_n,\boldsymbol{w})=\sigma(\boldsymbol{x}_n'\boldsymbol{w}),\quad\mathbb{P}(t_n=0|\boldsymbol{x}_n,\boldsymbol{w})=1-\sigma(\boldsymbol{x}_n'\boldsymbol{w}),
\end{align*}
$$

其中 $ \sigma(x)=\dfrac{1}{1+\mathrm{e}^{-x}} $ 为 Sigmoid 函数，以及错分概率

$$
\begin{align*}
    \mathbb{P}(y_n=0|t_n=1)=\alpha,~~\mathbb{P}(y_n=1|t_n=1)=1-\alpha,\\
    \mathbb{P}(y_n=1|t_n=0)=\beta,~~\mathbb{P}(y_n=0|t_n=0)=1-\beta.
\end{align*}
$$

其中 $ \lbrace y_n\rbrace_{n=1}^N,\lbrace t_n\rbrace_{n=1}^N $ 分别是观测标签和真实标签. 那么根据条件概率公式，完全数据的似然函数可以写为

$$
\begin{align*}
    \ell(\boldsymbol{w},\alpha,\beta|\boldsymbol{X},\boldsymbol{t},\boldsymbol{y})&=\log p(\boldsymbol{t},\boldsymbol{y}|\boldsymbol{X},\boldsymbol{w},\alpha,\beta)\\
    &=\log p(\boldsymbol{t}|\boldsymbol{X},\boldsymbol{w})+\log p(\boldsymbol{y}|\boldsymbol{t},\alpha,\beta)\\
    &=\sum_{n=1}^N[\log p(t_n|\boldsymbol{x}_n,\boldsymbol{w})+\log p(y_n|t_n,\alpha,\beta)].
\end{align*}
$$

注意，这个条件概率的拆分源于条件独立性. 在给定 $ \boldsymbol{t} $ 时， $ \boldsymbol{y} $ 就不能与原始数据 $ \boldsymbol{X} $ 产生联系了，就是在传导路径 $ \boldsymbol{X}\to\boldsymbol{t}\to\boldsymbol{y} $ 过程中被切断了联系. 现在我们就可以方便地写出两个密度函数的表达式：

$$
\begin{align*}
    \log p(t_n|\boldsymbol{x}_n,\boldsymbol{w})&=t_n\log\sigma(\eta_n)+(1-t_n)\log(1-\sigma(\eta_n))\\
    &=t_n\eta_n-\log(1+\mathrm{e}^{\eta_n}),~~(\eta_n=\boldsymbol{x}_n'\boldsymbol{w})\\
    \log p(y_n|t_n,\alpha,\beta)&=y_n[t_n\log(1-\alpha)+(1-t_n)\log\beta]\\&+(1-y_n)[t_n\log\alpha+(1-t_n)\log(1-\beta)].
\end{align*}
$$

现在我们需要 $ \boldsymbol{t} $ 的后验概率 $ r_n=p(t_n=1|y_n,\boldsymbol{x}_n,\boldsymbol{w}^{(t)},\alpha^{(t)},\beta^{(t)}) $ 来进行 E 步的计算. 首先考虑 $ y_n=1 $ 的情况. 我们有

$$
\begin{align*}
    p(t_n=1|y_n=1,\boldsymbol{x}_n)\propto p(y_n=1|t_n=1,x_n)\cdot p(t_n=1|x_n)=(1-\alpha^{(t)})\pi_n^{(t)},\\
    p(t_n=0|y_n=1,\boldsymbol{x}_n)\propto p(y_n=1|t_n=0,x_n)\cdot p(t_n=0|x_n)=\beta^{(t)}(1-\pi_n^{(t)}).
\end{align*}
$$

于是当 $ y_n=1 $ 时，有

$$
\begin{equation*}
    r_{n1}^{(t)}=\frac{(1-\alpha^{(t)})\pi_n^{(t)}}{(1-\alpha^{(t)})\pi_n^{(t)}+\beta^{(t)}(1-\pi_n^{(t)})},
\end{equation*}
$$

同理可得当 $ y_n=0 $ 时，

$$
\begin{equation*}
    r_{n0}^{(t)}=\frac{\alpha^{(t)}\pi_n^{(t)}}{\alpha^{(t)}\pi_n^{(t)}+(1-\beta^{(t)})(1-\pi_n^{(t)})}.
\end{equation*}
$$

我们不必把 $ Q(\boldsymbol{\theta}|\boldsymbol{\theta}^{(t)})=\mathbb{E}\_{\boldsymbol{t}}[\ell(\boldsymbol{w},\alpha,\beta|\boldsymbol{X},\boldsymbol{t},\boldsymbol{y})|\boldsymbol{X},\boldsymbol{y},\boldsymbol{\theta}^{(t)}] $ 的整体表达式写出来，后面将会看到可以对三个参数分别优化. 从上面的结果可以看到，实际上对 $ t_n $ 的期望可以写成 $ r_n^{(t)}=y_nr_{n1}^{(t)}+(1-y_n)r_{n0}^{(t)} $，或者 $ y_n $ 的指数幂的形式，**在工程实现中，我们不需要直接写 `if-else` 语句，而是可以直接利用向量选择的方式进行优化，从而实现硬件级的计算加速，具体见之后的代码**.
我们首先看关于 $ \boldsymbol{w} $ 的项. 我们有

$$
\begin{equation*}
    \frac{\partial Q}{\partial\boldsymbol{w}}=\sum_{n=1}^N[r_n^{(t)}-\sigma(\eta_n)]\boldsymbol{x}_n,
\end{equation*}
$$

这个驻点没有显式解，我们可以通过 Newton-Ralphson 方法来得到最优解. 接下来对于 $ \alpha $ 和 $ \beta $，分别有

$$
\begin{align*}
    \frac{\partial Q}{\partial\alpha}&=\sum_{n=1}^N\left[\frac{r_n^{(t)}(1-y_n)}{\alpha}-\frac{r_n^{(t)} y_n}{1-\alpha}\right],\\
    \frac{\partial Q}{\partial\beta}&=\sum_{n=1}^N\left[\frac{(1-r_n^{(t)})y_n}{\beta}-\frac{(1-r_n^{(t)})(1-y_n)}{1-\beta}\right].
\end{align*}
$$

经过化简合并可以得到

$$
\begin{align*}
    \alpha^{(t+1)}&=\frac{\sum_{n=1}^Nr_n^{(t)}(1-y_n)}{\sum_{n=1}^Nr_n^{(t)}},\\
    \beta^{(t+1)}&=\frac{\sum_{n=1}^N(1-r_n^{(t)})y_n}{\sum_{n=1}^N(1-r_n^{(t)})}.
\end{align*}
$$

接下来给出代码. 代码中的基线模型是传统的逻辑回归，仍然配备了一个数值模拟的案例并给出对比结果.

```python
import numpy as np
from scipy.special import expit
import time

class LabelNoiseEM:
    def __init__(self, tol=1e-6, max_iter=100, nr_iter=10):
        self.tol = tol
        self.max_iter = max_iter
        self.nr_iter = nr_iter
        self.w = None
        self.alpha = None
        self.beta = None

    def _compute_log_likelihood(self, X, y):
        pi = expit(X @ self.w)
        # 使用统一形式计算 P(y|x)
        p_y_cond_t1 = y * (1 - self.alpha) + (1 - y) * self.alpha
        p_y_cond_t0 = y * self.beta + (1 - y) * (1 - self.beta)
        p_y_cond_x = p_y_cond_t1 * pi + p_y_cond_t0 * (1 - pi)
        
        return np.sum(np.log(np.clip(p_y_cond_x, 1e-15, 1 - 1e-15)))

    def fit(self, X, y):
        N, D = X.shape
        self.w = np.zeros(D)
        self.alpha = 0.1
        self.beta = 0.1
        ll_history = []

        for i in range(self.max_iter):
            # --- E 步 ---
            pi = expit(X @ self.w)
            
            # 高效实现 r^(t) 的计算
            p_y_cond_t1 = y * (1 - self.alpha) + (1 - y) * self.alpha
            p_y_cond_t0 = y * self.beta + (1 - y) * (1 - self.beta)
            
            num = p_y_cond_t1 * pi
            den = num + p_y_cond_t0 * (1 - pi)
            r = num / (den + 1e-16)

            # --- M 步 ---
            # 1. 更新 w (Newton-Raphson)
            for _ in range(self.nr_iter):
                p_w = expit(X @ self.w)
                grad_w = X.T @ (r - p_w)
                W = p_w * (1 - p_w)
                H = - (X.T * W) @ X - 1e-6 * np.eye(D)
                self.w -= np.linalg.solve(H, grad_w)

            # 2. 更新 alpha 和 beta
            sum_r = np.sum(r)
            self.alpha = np.sum(r * (1 - y)) / (sum_r + 1e-16)
            self.beta = np.sum((1 - r) * y) / (np.sum(1 - r) + 1e-16)

            ll = self._compute_log_likelihood(X, y)
            ll_history.append(ll)
            if i > 0 and abs(ll_history[-1] - ll_history[-2]) < self.tol:
                break
        
        return ll_history

# --- 数值模拟与报告 ---

def run_simulation():
    np.random.seed(123)
    N, D = 50000, 5  # 增加数据量以观察加速效果
    
    # 真实参数
    true_w_with_bias = np.array([0.5, 1.5, -2.0, 0.8, 4.0, -3.5])
    true_alpha, true_beta = 0.12, 0.08
    
    # 生成数据
    X = np.hstack([np.ones((N, 1)), np.random.randn(N, D)])
    t = (np.random.rand(N) < expit(X @ true_w_with_bias)).astype(float)
    
    # 注入噪声
    y = np.copy(t)
    pos_mask = (t == 1)
    neg_mask = (t == 0)
    y[pos_mask] = np.where(np.random.rand(np.sum(pos_mask)) < true_alpha, 0, 1)
    y[neg_mask] = np.where(np.random.rand(np.sum(neg_mask)) < true_beta, 1, 0)
            
    # 训练 EM 模型
    model = LabelNoiseEM()
    start_time = time.time()
    history = model.fit(X, y)
    end_time = time.time()
    
    # 训练普通 LR 作为基准
    w_naive = np.zeros(X.shape[1])
    for _ in range(10):
        p = expit(X @ w_naive)
        H = -(X.T * (p * (1 - p))) @ X - 1e-6 * np.eye(X.shape[1])
        w_naive -= np.linalg.solve(H, X.T @ (y - p))

    # 报告结果
    print(f"训练完成，耗时: {end_time - start_time:.4f} 秒")
    print("\n" + "="*70)
    print(f"{'参数项':<15} | {'真实值':<12} | {'EM 估计值':<12} | {'普通 LR (忽略噪声)'}")
    print("-" * 70)
    print(f"{'alpha':<15} | {true_alpha:<12.4f} | {model.alpha:<12.4f} | {'N/A'}")
    print(f"{'beta':<15} | {true_beta:<12.4f} | {model.beta:<12.4f} | {'N/A'}")
    for j, val in enumerate(true_w_with_bias):
        print(f"{'w_' + str(j):<15} | {val:<12.4f} | {model.w[j]:<12.4f} | {w_naive[j]:<12.4f}")
    print("="*70)

if __name__ == "__main__":
    run_simulation()

# 输出结果
训练完成，耗时: 1.3146 秒

======================================================================
参数项             | 真实值          | EM 估计值       | 普通 LR (忽略噪声)
----------------------------------------------------------------------
alpha           | 0.1200       | 0.1176       | N/A
beta            | 0.0800       | 0.0811       | N/A
w_0             | 0.5000       | 0.5116       | 0.0515      
w_1             | 1.5000       | 1.5352       | 0.5057      
w_2             | -2.0000      | -2.0453      | -0.7086     
w_3             | 0.8000       | 0.8399       | 0.2836      
w_4             | 4.0000       | 4.0203       | 1.3734      
w_5             | -3.5000      | -3.4866      | -1.1952     
======================================================================
```

代码中有些工程实现的细节，比如防止数值下溢、矩阵奇异，以及 $ r_n^{(t)} $ 的高效计算等等. 从结果可以看到，带噪声的情况下，普通的逻辑回归对参数向量的估计产生了严重的偏差，而通过 EM 算法去噪可以得到准确的噪声存在情况及原始参数向量.