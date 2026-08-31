---
title: 2026年度 MIT 积分大赛 Regular Season 解答
date: 2026-02-11 16:24:30
tags: [数学]
categories: [杂项]
---

*基本是自己做的，有两三道题的思路其实自己已经想到了但莫名其妙做不出来，于是翻了一下官方的录像回放...*

### Problem 1

$$
\begin{align*}
    \int_{0}^{2026}\left(\sum_{k=0}^{\infty}\frac{x^{2026k+2025}}{k!}\right)\mathrm{d}x&=\int_{0}^{2026}x^{2025}\left(\sum_{k=0}^\infty\frac{(x^{2026})^k}{k!}\right)\mathrm{d}x\\
    &=\int_{0}^{2026}x^{2025}\mathrm{e}^{x^2026}\mathrm{d}x\\
    &=\left.\frac{\mathrm{e}^{x^{2026}}}{2026}\right|_{0}^{2026}=\frac{\mathrm{e}^{2026^{2026}}-1}{2026}.
\end{align*}
$$

### Problem 2

$$
\begin{align*}
    \int_{0}^{\log(34)}\left(2+\frac{1}{2\mathrm{e}^x-1}+\frac{1}{2\mathrm{e}^x-1}\right)\mathrm{d}x&=\int_{0}^{\log(34)}\frac{2\mathrm{e}^x}{2\mathrm{e}^x-1}+\frac{3\mathrm{e}^x}{3\mathrm{e}^x-1}\mathrm{d}x\\
    &=\int_{1}^{34}\frac{2}{2t-1}+\frac{3}{3t-1}\mathrm{d}t\\
    &=\left.\log[(2t-1)(3t-1)]\right|_{1}^{34}=\log\left(\frac{6767}{2}\right).
\end{align*}
$$

### Problem 3

$$
\begin{align*}
    &\int_{0}^{\frac{1}{2}}\left(\cos(\pi x)-\pi\left(\frac{1}{4}-x^2\right)\left(\frac{5}{4}-x^2\right)\right)\mathrm{d}x\\&=\left.\frac{1}{\pi}\sin(\pi x)\right|_{0}^{\frac{1}{2}}-\pi\int_{0}^{\frac{1}{2}}\left(x^4-\frac{3}{2}x^2+\frac{5}{16}\right)\mathrm{d}x\\
    &=\frac{1}{\pi}-\pi\left(\frac{1}{5}\cdot\frac{1}{32}-\frac{1}{2}\cdot\frac{1}{8}+\frac{5}{32}\right)=\frac{1}{\pi}-\frac{\pi}{10}.
\end{align*}
$$

### Problem 4

$$
\begin{align*}
    \int_{\int_{\int_{\int_{0}^1x\mathrm{d}x}^{1}\mathrm{d}x}^{1}\mathrm{d}x}^{\int_{0}^{\int_{0}^{\int_{0}^{1}x\mathrm{d}x}x\mathrm{d}x}x\mathrm{d}x}x\mathrm{d}x
    &=\int_{\int_{\int_{1/2}^{1}x\mathrm{d}x}^{1}\mathrm{d}x}^{\int_{0}^{\int_{0}^{1/2}x\mathrm{d}x}x\mathrm{d}x}x\mathrm{d}x\\
    &=\int_{\int_{3/8}^{1}x\mathrm{d}x}^{\int_{0}^{1/8}x\mathrm{d}x}x\mathrm{d}x\\
    &=\int_{55/128}^{1/128}x\mathrm{d}x=-\frac{189}{2048}.
\end{align*}
$$

### Problem 5

$$
\begin{align*}
    \int\frac{4x^2-1}{x^2(x^2-1)}\mathrm{d}x
    &=\int\frac{3}{x^2-1}+\frac{1}{x^2}\mathrm{d}x\\
    &=-\frac{1}{x}+\frac{3}{2}\log\left|\frac{x-1}{x+1}\right|+C.
\end{align*}
$$

### Problem 6

$$
\begin{align*}
    \int\frac{\mathrm{d}x}{27x-x^{-1/3}}&=\int\frac{x^{1/3}\mathrm{d}x}{27x^{4/3}-1}
    =\int\frac{3}{4}\cdot\frac{\mathrm{d}t}{27t-1}\quad(t=x^{4/3})\\&=\frac{1}{36}\log\left|27x^{\frac{4}{3}}-1\right|+C.
\end{align*}
$$

### Problem 7

$$
\begin{align*}
    \int\frac{2\sqrt{x}+1}{\sqrt{x^2+x\sqrt{x}}}\mathrm{d}x&=\int\frac{2+1/\sqrt{x}}{\sqrt{x+\sqrt{x}}}\mathrm{d}x\\
    &=2\int\frac{1}{\sqrt{x+\sqrt{x}}}\mathrm{d}(x+\sqrt{x})\\&=4\sqrt{x+\sqrt{x}}+C.
\end{align*}
$$

### Problem 8

$$
\begin{align*}
    \int\frac{\mathrm{e}^x}{\mathrm{e}^{2x}-\mathrm{e}^{-2x}}\mathrm{d}x&=\int\frac{t^2}{t^4-1}\mathrm{d}t\quad(t=\mathrm{e}^x)\\
    &=\frac{1}{2}\left(\int\frac{1}{t^2+1}\mathrm{d}t+\int\frac{1}{t^2-1}\mathrm{d}t\right)\\
    &=\frac{1}{2}\arctan(\mathrm{e}^x)+\frac{1}{4}\log\left|\frac{\mathrm{e}^x-1}{\mathrm{e}^x+1}\right|+C.
\end{align*}
$$

### Problem 9

$$
\begin{align*}
   &\int_{0}^\pi\frac{|\sin(x)+\sin(2x)+\sin(3x)+\sin(4x)|}{\sin(x)+\sin(2x)+\sin(3x)+\sin(4x)}\mathrm{d}x\\
    &=\int_{0}^{\pi}\text{sgn}(\sin(x)+\sin(2x)+\sin(3x)+\sin(4x))\mathrm{d}x\\
    &=\int_{0}^{\pi}\text{sgn}\left[\sin\left(\frac{5}{2}x\right)\cos\left(\frac{3}{2}x\right)+\sin\left(\frac{5}{2}x\right)\cos\left(\frac{1}{2}x\right)\right]\mathrm{d}x\\
    &=\int_{0}^\pi\text{sgn}\left[\sin\left(\frac{5}{2}x\right)\cos\left(x\right)\cos\left(\frac{1}{2}x\right)\right]\mathrm{d}x\\
    &=\frac{2\pi}{5}-\frac{\pi}{10}+\frac{3\pi}{10}-\frac{\pi}{5}=\frac{2\pi}{5}.
\end{align*}
$$

### Problem 10

$$
\begin{align*}
    \int_{0}^1\frac{1-x}{\sqrt{\mathrm{e}^{2x}-x^2}}\mathrm{d}x&=\int_{0}^{1}\frac{(1-x)\mathrm{e}^{-x}}{1-(x\mathrm{e}^{-x})^2}\mathrm{d}x\\
    &=\arcsin(x\mathrm{e}^{-x})|_0^1=\arcsin\left(\frac{1}{\mathrm{e}}\right).
\end{align*}
$$

### Problem 11

$$
\begin{align*}
    &\int_{0}^{2026}\left\{x+\frac{1}{2}\left\lfloor\frac{x}{2}\right\rfloor+\frac{1}{3}\left\lfloor\frac{x}{3}\right\rfloor+\frac{1}{4}\left\lfloor\frac{x}{4}\right\rfloor+\cdots\right\}\mathrm{d}x\\
    &=\sum_{n=0}^{2025}\int_{n}^{n+1}\left\{x+\frac{1}{2}\left\lfloor\frac{x}{2}\right\rfloor+\frac{1}{3}\left\lfloor\frac{x}{3}\right\rfloor+\frac{1}{4}\left\lfloor\frac{x}{4}\right\rfloor+\cdots\right\}\mathrm{d}x\\
    &=2026\int_{0}^1x\mathrm{d}x=1013.
\end{align*}
$$

**注**  本题需要注意到在 $ [n,n+1] $ 这个区间上，后面的取整项的值并不会发生改变，所以我们可以对每个长度为 $ 1 $ 的区间平移积分区域，由于 $ f(x)=\{x\} $ 是一个周期为 $ 1 $ 的周期函数，因此每一块积分都等于在 $ [0,1] $ 上的积分。

### Problem 12

$$
\begin{align*}
    \int_{1/2}^2\frac{x^8}{x^8-x^6+x^4-x^2+1}\mathrm{d}x&=\int_{1/2}^2\frac{1/x^2}{x^8-x^6+x^4-x^2+1}\mathrm{d}x\\
    &=\frac{1}{2}\int_{1/2}^2\frac{1}{x^2}\frac{x^{10}+1}{x^8-x^6+x^4-x^2+1}\mathrm{d}x\\
    &=\frac{1}{2}\int_{1/2}^2\frac{x^2+1}{x^2}\mathrm{d}x=\frac{3}{2}.
\end{align*}
$$

### Problem 13

$$
\begin{align*}
    \int(\text{arcsinh}(x))^2\mathrm{d}x&=x(\text{arcsinh}(x))^2-\int\frac{2x\text{arcsinh}x}{\sqrt{1+x^2}}\mathrm{d}x\\
    &=x(\text{arcsinh}(x))^2-2\text{arcsinh}(x)\sqrt{1+x^2}+2x+C.
\end{align*}
$$

**注**  本题需要注意到 $ (\text{arcsinh}(x))^{\prime}=(\log(\sqrt{1+x^{2}}+x))^{\prime}=1/(\sqrt{1+x^{2}}) $。

### Problem 14

$$
\begin{align*}
    \int\cos^2(2026x)\cos(1013x)\mathrm{d}x&=\frac{1}{2}\int[\cos(1013x)+\cos(4052x)\cos(1013x)]\mathrm{d}x\\
    &=\frac{1}{2026}\sin(1013x)+\frac{1}{4}\int[\cos(5065x)+\cos(3039x)]\mathrm{d}x\\
    &=\frac{\sin(1013x)}{2026}+\frac{\sin(3039x)}{12156}+\frac{\sin(5065x)}{20260}+C.
\end{align*}
$$

### Problem 15

$$
\begin{align*}
    \int\mathrm{e}^x\arcsin(\tanh(x))\mathrm{d}x&=\mathrm{e}^x\arcsin(\tanh(x))-\int\frac{\mathrm{e}^x/\cosh^2(x)}{\sqrt{1-\tanh^2(x)}}\mathrm{d}x\\
    &=\mathrm{e}^x\arcsin(\tanh(x))-\int\frac{2\mathrm{e}^x}{\mathrm{e}^x+\mathrm{e}^{-x}}\mathrm{d}x\\
    &=\mathrm{e}^x\arcsin(\tanh(x))-\log(\mathrm{e}^{2x}+1)+C.
\end{align*}
$$

### Problem 16

$$
\begin{align*}
    \int\frac{x((1-x^2)\cos(x)+2x\sin(x))}{(1+x^2)^2}\mathrm{d}x&=-\frac{1}{2}\int((1-x^2)\cos(x)+2x\sin(x))\mathrm{d}\left(\frac{1}{1+x^2}\right)\\
    &=-\frac{(1-x^2)\cos(x)+2x\sin(x)}{2(1+x^2)}+\frac{1}{2}\int\sin(x)\mathrm{d}x\\
    &=-\frac{\cos(x)+x\sin(x)}{1+x^2}+C.
\end{align*}
$$

### Problem 17

$$
\begin{align*}
    &\int_{0}^1\sin((x-1)(5x-1))\sin(x(3x-2))\mathrm{d}x\\&=\frac{1}{2}\int_{0}^1\cos(2x^2-4x+1)-\cos(8x^2-8x+1)\mathrm{d}x\\
    &=\int_{0}^{1/2}\cos(8x^2-8x+1)\mathrm{d}x-\frac{1}{2}\int_{0}^{1}\cos(8x^2-8x+1)\mathrm{d}x\\&=0.
\end{align*}
$$

**注**  第三步对前一个积分用到了换元 $ x=2t $，需要对式子结构很熟悉，最后等于 $ 0 $ 则是因为里面的二次函数对称轴为 $ 1/2 $，因此可使用区间再现。

### Problem 18

$$
\begin{align*}
    \int_{1/4}^{\sqrt{2}}\frac{\mathrm{d}x}{x(x^{x^{x^\cdots}}\log(x)-1)}&=\int_{1/2}^2\frac{\mathrm{d}((\log t)/t)}{\log t-1}\quad\left(t=x^{x^{x^{\cdots}}}\right)\\
    &=-\int_{1/2}^{2}\frac{1}{t^2}\mathrm{d}t=-\frac{3}{2}.
\end{align*}
$$

**注**   换元时积分上下限的确定使用了等式 $ \log x={\log t}/{t} $，以及 $ f(t)=\log t/t $ 在 $(0,2)$ 上的单调性。此外，无穷指数塔 $ x^{x^{x^{\cdots}}} $ 在 $ 0<x<\mathrm{e}^{1/\mathrm{e}} $ 时是收敛的，因此不用担心积分区域不存在的问题。

### Problem 19

$$
\begin{align*}
    &\int_0^1\frac{1+((2026+2x-x^2)^2-2)(2026+4x-4x^2)^2}{(2025+2x-x^2)(2027+2x-x^2)(2025+4x-4x^2)(2027+4x-4x^2)}\mathrm{d}x
\end{align*}
$$

这题需要注意到两个恒等式：

$$
\begin{align*}
    a(x)=(2025+2x-x^2)(2027+2x-x^2)&=(2026+2x-x^2)^2-1,\\
    b(x)=(2025+4x-4x^2)(2027+4x-4x^2)&=(2026+4x-4x^2)^2-1,
\end{align*}
$$

原积分就变为 

$$
\begin{equation*}
    \int_{0}^{1}\frac{1+(a(x)-1)(b(x)+1)}{a(x)b(x)}\mathrm{d}x=1+\int_{0}^{1}\left(\frac{1}{b(x)}-\frac{1}{a(x)}\right)\mathrm{d}x.
\end{equation*}
$$

再次利用第 $ 17 $ 题中的方法可以得到 

$$
\begin{align*}
    \int_{0}^{1}\frac{1}{b(x)}\mathrm{d}x=\frac{1}{2}\int_{0}^{2}\frac{1}{a(x)}\mathrm{d}x=\int_{0}^{1}\frac{1}{a(x)}\mathrm{d}x.
\end{align*}
$$

因此后面的积分值相互抵消，原积分值就是 $ 1 $ 。

### Problem 20

$$
\begin{align*}
    &\int_{1/2}^{1}4^{x-1}\left(1+4^{4^{x-1}-1}\left(1+4^{4^{4^{x-1}-1}-1}(1+\cdots)\right)\right)\mathrm{d}x=I\\
    &=\frac{1}{2\log4}+\frac{1}{\log4}\int_{1/2}^14^{t-1}\left(1+4^{4^{t-1}-1}\left(1+\cdots\right)\right)\mathrm{d}t
    \\&=\frac{1}{2\log4}+\frac{1}{\log4}I~~\Rightarrow~~I=\frac{1}{2(\log4-1)}.
\end{align*}
$$

**注**  该题的关键点在于使用 $ t=4^{x-1} $ 换元后，积分的上下限始终不改变。
