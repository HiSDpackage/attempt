---
layout: single
permalink: /
title: "iFEM: an Integrated Finite Element Methods Package in MATLAB"
excerpt: "A quick start guide."
sidebar:
    nav: docs
toc: true
toc_sticky: true
mathjax: true

---

# HiOSD算法简介
作者：肖子翔  

## 1. 什么是 $k$ 阶鞍点

给定一个在实希尔伯特空间 $\mathcal{H}$ 上定义的二次Fréchet可微的能量泛函 $E(\boldsymbol{x})$ ，其内积为 $\langle \cdot,\cdot \rangle$ ，我们令 $\boldsymbol{F}(\boldsymbol{x})=-\nabla E(\boldsymbol{x})$ 表示其自然力， $\mathbb{G}(\boldsymbol{x})=\nabla^2E(\boldsymbol{x})$ 表示其Hessian矩阵。

- 如果 $\lVert \boldsymbol{F}(\boldsymbol{\hat{x}}) \rVert =0$，则称 $\boldsymbol{\hat{x}}\in\mathcal{H}$ 为 $E(\boldsymbol{x})$ 的一个临界点。
- 一个不是局部极值的 $E(\boldsymbol{x})$ 的临界点被称为 $E(\boldsymbol{x})$ 的一个鞍点。特别的，有些时候我们也把局部极小值称为 $0$ 阶鞍点，局部极大值称为 $d$（系统维度）阶鞍点。
- 如果 $\mathbb{G}(\boldsymbol{\hat{x}})$ 有一个有界逆，则称临界点 $\boldsymbol{\hat{x}}$ 为非退化的。
- 根据Morse理论，非退化临界点 $\boldsymbol{\hat{x}}$ 的指数（Morse指数）是使得 $\mathbb{G}(\boldsymbol{\hat{x}})$ 为负定的最大子空间 $\mathcal{K}$ 的维数。我们的目标是找到势能面上的指数为 $k$ 的鞍点，简称为 $k$-saddle 或 $k$ 阶鞍点。
- 为简便起见，我们假设 $\mathcal{H}$ 的维数为 $d$，并将内积 $\langle \boldsymbol{x},\boldsymbol{y} \rangle$ 写为 $\boldsymbol{x}^{\top}\boldsymbol{y}$。

下面先从特征值和特征向量的角度给一个鞍点阶数的直观理解，便于之后问题的开展。

我们注意到Hessian矩阵是实对称矩阵（用到前面的二次Fréchet可微），其正交相似于一个对角矩阵，即存在正交矩阵 $\mathbb{T}$ 和对角矩阵 $\mathbb{D}$ 使得 

$$
\mathbb{T}^{-1}\mathbb{G}(\boldsymbol{\hat{x}})\mathbb{T}=\mathbb{D}.
$$

事实上这个对角矩阵的对角线上各个元素即为该实对称矩阵的特征值 $\hat{\lambda}_1,\hat{\lambda}_2,\ldots,\hat{\lambda}_d$，且不妨 $\hat{\lambda}_1\leq\hat{\lambda}_2\leq\ldots\leq\hat{\lambda}_d$，并设此时 $\mathbb{T}$ 的列向量为 $\boldsymbol{\hat{v}}_1,\boldsymbol{\hat{v}}_2,\ldots,\boldsymbol{\hat{v}}_d$。

可以证明，对于非退化临界点 $\hat{x}$，若 $\hat{\lambda}_1\leq\ldots\leq\hat{\lambda}\_k<0<\hat{\lambda}\_{k+1}\leq\ldots\leq\hat{\lambda}\_d$ ，则 $k$ 即为该鞍点的阶数。一方面， $\mathbb{G}(\boldsymbol{\hat{x}})$ 在 $\boldsymbol{\hat{v}}\_1,\boldsymbol{\hat{v}}\_2,\ldots,\boldsymbol{\hat{v}}_k$ 生成的子空间 $\hat{\mathcal{V}}$ 上是负定的，故 $\boldsymbol{\hat{x}}$ 至少为 $k$ 阶；另一方面，对于 $\mathcal{H}$ 的任意一个 $k+1$ 维子空间 $\mathcal{K'}$，有 $\boldsymbol{\hat{v}}\_{k+1},\ldots,\boldsymbol{\hat{v}}\_d$ 生成的子空间 $\hat{\mathcal{V}}^{\perp}$ 与 $\mathcal{K'}$ 交非零（否则 $\mathcal{K'}$ 中添上 $\boldsymbol{\hat{v}}\_{k+1},\ldots,\boldsymbol{\hat{v}}_d$ 生成 $d+1$ 维空间，矛盾！）。取交集中一个非零元素

$$
\boldsymbol{w}=\sum_{i=k+1}^{d} a_i\boldsymbol{\hat{v}}_i
$$

则有

$$
\boldsymbol{w^{\top}\mathbb{G}(\hat{x})w}=\sum_{i=k+1}^{d} a_i\boldsymbol{\hat{v}}_i^{\top}\sum_{i=k+1}^{d} \hat{\lambda}_ia_i\boldsymbol{\hat{v}}_i=\sum\_{i=k+1}^{d} \hat{\lambda}_ia_i^2>0
$$

故 $\mathbb{G}(\boldsymbol{\hat{x}})$ 在 $\mathcal{K'}$ 上不是负定的，即 $\boldsymbol{\hat{x}}$ 的阶数为 $k$。

## 将寻找$k$阶鞍点转化为优化问题
沿用上面的符号，注意到 $\mathbb{G}(\boldsymbol{\hat{x}})$ 在 $\hat{\mathcal{V}}$ 上是负定的，在其正交补 $\hat{\mathcal{V}}^{\perp}$ 上是正定的，这意味着 $\boldsymbol{\hat{x}}$ 是线性流形 $\boldsymbol{\hat{x}} + \hat{\mathcal{V}}$ 上的局部极大值，同时也是线性流形 $\boldsymbol{\hat{x}} + \hat{\mathcal{V}}^{\perp}$ 上的局部极小值。

考虑 $\boldsymbol{\hat{x}}\_{\hat{\mathcal{V}}}, \boldsymbol{\hat{x}}\_{\hat{\mathcal{V}}^{\perp}}$ 分别为 $\boldsymbol{\hat{x}}$ 在 $\hat{\mathcal{V}}, \hat{\mathcal{V}}^{\perp}$ 上的投影，则 $(\boldsymbol{v}, \boldsymbol{w}) = (\boldsymbol{\hat{x}}\_{\hat{\mathcal{V}}}, \boldsymbol{\hat{x}}\_{\hat{\mathcal{V}}^{\perp}})$ 是 minimax 问题

$$
\min\{\boldsymbol{w}\in\hat{\mathcal{V}}^{\perp}} \max\{\boldsymbol{v}\in\hat{\mathcal{V}}} E(\boldsymbol{v} + \boldsymbol{w})
$$

的一个解。但是，这并不是一个经典的 minimax 问题，因为空间 $\hat{\mathcal{V}}$ 是未知的，所以在求解优化问题的过程中我们的迭代法应该包括两个部分：一个是更新 $\boldsymbol{v}$ 和 $\boldsymbol{w}$（在这个问题中也就是更新 $\boldsymbol{x} = \boldsymbol{v} + \boldsymbol{w}$），还有一个就是要更新空间 $\mathcal{V}$（$\mathcal{V}$ 用于近似 $\hat{\mathcal{V}}$，一般用当前 $\boldsymbol{x}$ 处 Hessian 矩阵的最小 $k$ 个特征值对应的特征向量张成的子空间来描述）。

## $\boldsymbol{x}$的动力学
更新 $\boldsymbol{x}$ 直观上看是让 $\boldsymbol{\dot{x}}$ 在空间 $\mathcal{V}$ 上的投影 $\mathcal{P}_{\mathcal{V}}\boldsymbol{\dot{x}}$ 为能量函数 $E(\boldsymbol{x})$ 的上升方向，而在其补空间 $\mathcal{V}^{\perp}$ 上的投影 $\mathcal{P}\_{\mathcal{V}^{\perp}}\boldsymbol{\dot{x}}$ 为下降方向。

特别地，注意到自然力 $\boldsymbol{F}(\boldsymbol{x})=-\nabla E(\boldsymbol{x})$ 为最速下降方向，故可以考虑令 $\mathcal{P}_{\mathcal{V}}\boldsymbol{\dot{x}}=-\mathcal{P}\_{\mathcal{V}}\boldsymbol{F}(\boldsymbol{x})$ 以及

$$
\mathcal{P}_{\mathcal{V}^{\perp}}\boldsymbol{\dot{x}}=\mathcal{P}_{\mathcal{V}^{\perp}}\boldsymbol{F}(\boldsymbol{x})=\boldsymbol{F}(\boldsymbol{x})-\mathcal{P}\_{\mathcal{V}}\boldsymbol{F}(\boldsymbol{x})
$$

再取两个正的松弛常数 $\beta_{\mathcal{V}}$ 和 $\beta_{\mathcal{V}^{\perp}}$ 即可以给出 $\boldsymbol{x}$ 的动力学

$$
\boldsymbol{\dot{x}}=\beta_{\mathcal{V}}(-\mathcal{P}\_{\mathcal{V}}\boldsymbol{F}(\boldsymbol{x}))+\beta_{\mathcal{V}^{\perp}}(\boldsymbol{F}(\boldsymbol{x})-\mathcal{P}_{\mathcal{V}}\boldsymbol{F}(\boldsymbol{x}))
$$

进一步，如果简单地令 $\beta_{\mathcal{V}}=\beta_{\mathcal{V}^{\perp}}=\beta$ 则上式化为

$$
\beta^{-1}\boldsymbol{\dot{x}}=\boldsymbol{F}(\boldsymbol{x})-2\mathcal{P}_{\mathcal{V}}\boldsymbol{F}(\boldsymbol{x})
\label{the dynamics of x easy vesion}
$$

特别地，如果给出空间 $\mathcal{V}$ 的一组标准正交基 $\boldsymbol{v_1},\boldsymbol{v_2},\ldots,\boldsymbol{v_k}$ 则有投影变换 $\mathcal{P}_{\mathcal{V}}=\sum\_{i=1}^{k}\boldsymbol{v}_i\boldsymbol{v}^{\top}\_i$，从而公式

$$
\beta^{-1}\boldsymbol{\dot{x}}=\left(\mathbb{I}-2\sum_{i=1}^{k}\boldsymbol{v}_i\boldsymbol{v}^{\top}_i\right)\boldsymbol{F}(\boldsymbol{x})
\label{the dynamics of x easy vesion 2}
$$

其中 $\mathbb{I}$ 为单位矩阵。


## Installation

You can download the repository at [https://github.com/lyc102/ifem](https://github.com/lyc102/ifem), or alternatively you can get iFEM by using the following commands:

```bash
git clone https://github.com/lyc102/ifem.git
```

Then use MATLAB to run the `setpath.m` script in the root folder to add all the sub-folders to your MATLAB path. 

<!-- [Octave](www.octave.org) version is also available at [https://github.com/lyc102/ifemOctave](https://github.com/lyc102/ifemOctave). -->



## Help, Guides, and Contributing

This documentation website will be constantly updated. If you have any questions, please feel free to [contact us](mailto:lyc102@gmail.com). If you like to suggest an additional equation to be implemented in iFEM, please go to the [GitHub repo submit an issue](https://github.com/lyc102/ifem/issues). If you like to contribute to the development of iFEM, please see our 【xxx】.

### Use MATLAB help/doc
1. `help funexample` displays a description of and syntax for the function `funexample`. For example, `help mg` will show basic usage for `mg` function in the plain text.  
2. `ifem funexampledoc` show detailed description. For example, `ifem mgdoc` will explain the `mg` function step by step in `html` format. But not every function has a html documentation. Contribution to a documentation page is also welcome.

## License and References

iFEM can be freely distributed under GPL 3.0. If you feel iFEM is helpful for your research, please acknowledge your use by citing:

> L. Chen. iFEM: an integrated finite element method package in MATLAB. Technical Report, University of California at Irvine, 2009.

```bibtex
@techreport{Chen:2008ifem,
author = {Long Chen},
journal = {Technical Report, University of California at Irvine},
title = {$i$FEM: an integrated finite element methods package in {MATLAB}},
url = {https://github.com/lyc102/ifem},
year = {2009}}
```



## Acknowledgement

The author would like to thank [Professor Michael Holst](http://cam.ucsd.edu/~mholst/) in University of California at San Diego and [Professor Ludmil Zikatanov](http://www.personal.psu.edu/ltz1/) in Pennsylvania State University for many insightful discussion, and also Professor [Chensong Zhang](http://lsec.cc.ac.cn/~zhangcs/) in Chinese Academy of Sciences for the effort in the development of AFEM@matlab, an early version of iFEM.

The author thanks students or postdocs: [Shuhao Cao](https://scaomath.github.io/), [Huayi Wei](https://weihuayi.github.io), Ming Wang, Lin Zhong, and Jie Zhou for their contribution to iFEM in one way or another. Detailed credits can be found in the M-lint of several `.m` files.

The documentation website https://lyc102.github.io/ifem/ is set up by [Shuhao Cao](https://scaomath.github.io/) and his effort is greatly appreciated.

The author is also grateful to the NSF for the partial support over the years. 



<div style="width:400px" onclick="myhref('http://faculty.bicmr.pku.edu.cn/~zhanglei/');"><hr/>
Lei Zhang
<br>
Boya Distinguished Professor                
<br>
Beijing International Center for Mathematical Research
<br>
Center for Quantitative Biology
<br>
Center for Machine Learning Research
<br>
Peking University
<br>
http://faculty.bicmr.pku.edu.cn/~zhanglei/
</div>

<div style="width:400px" onclick="myhref('https://liuonly1121.github.io/');"><hr/>
Yuyang Liu
<br>
Ph.D. Candidate               
<br>
School of Mathematical Science
<br>
Peking University
<br>
https://liuonly1121.github.io/
</div>

<script type="text/javascript">
    function myhref(web){
      window.location.href = web;}
</script>
