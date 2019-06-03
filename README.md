# QXmath

《用C语言开发的气象常用参数和物理量计算函数库》原文代码以及Python封装

## 安装方法
```
python setup.py build_ext --inplace
```
需要C编译器和`cython`模块

## 使用方法
在此目录新建Python解释器，即可import编译好的库。

### 使用原始函数
```python
>>> import pyqxmath as pq
>>> pq.showalter_index(16.6, 0.6, -15.9)
1.099999999999886
```

### 使用向量化函数
使用`numpy.vectorize`类进行函数的向量化，因此可以支持`numpy.ndarray`的直接运算。
```python
>>> import pyqxmath as pq
>>> import numpy as np
>>> t8 = np.array([16, 17, 18])
>>> td8 = np.array([2, 1, 0])
>>> t5 = np.array([-16, -17, -18])
>>> pq.vect.showalter_index(t8, td8, t5)
array([ 0.56  , -0.48  , -1.5456])
```
向量化函数和原始函数的名字是相同的。

可供Python使用的函数在`pyqxmath.pyx`中定义，由于文献中的函数名基本由拼音缩写而来，我在包装的过程中一般都会将其重命名成英文单词的组合，以方便记忆。函数的具体用法请参见`paper`文件夹里面的pdf文献。

## 未来计划实现的功能

1. 支持原C库中所有计算函数
2. 支持numpy array的运算
3. 完善函数的注释

## 备注

由于原文献扫描质量较差，可能会有（很多）排版造成的代码录入bug，如遇到计算不合理的地方，请参考`paper`文件夹里的pdf进行修改并重新编译，
欢迎提pull request。