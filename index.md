---
layout: single
permalink: /
title: "HiSD: an Algorithm Package for Computing High-index Saddle Points in Python"
excerpt: "A quick start guide."
sidebar:
    nav: docs
toc: true
toc_sticky: true
mathjax: true

---




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
