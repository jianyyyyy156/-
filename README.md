# -
1.3、根据执行时间来看，稳定的数值算法比直接数值积分更快，但是结果的误差相对来说较大一些。直接积分在计算量较大时耗时可能会进一步增长，如果稳定的数值算法能够保证误差在可接受范围内，从效率角度考虑，稳定的数值算法方法会更好。

1.4.1、 循环5000次，其他参数基本一致下FTCS和BTCS比较； Dt：FTCS执行时间: 12.6617 秒 2Dt：FTCS执行时间: 9.9067 秒 Dt：BTCS执行时间: 10.5074 秒 2Dt：BTCS执行时间: 8.0995 秒 BTCS两端没有温度梯度，温度始终不变，而FTCS在计算时更新了边界值发生了变化；2Dt的计算速度更快，但是FTCS计算时Dt过大会计算出错，相同循环次数Dt越大扩散越快且计算越快。 BTCS的计算速度更快，可能是因为代码中除了计算部分之外还有其他处理过程可能会影响结果，这里的执行时间结果可能存在误差影响了结果；

1.5、 蛙越格式更稳定，迎风格式在2倍dt时计算就开始出错了；

1.6、 放大dt后BTCS计算更快，其他条件相同情况下BTCS需要100步循环就可以达到FTCS运行700步循环的结果；此外，BTCS在边界处变化与FTCS格式有所不同，这可能是因为边界的处理计算方式差异导致；

2.6、非线性模型的波动相对线性模型更加复杂；造成这一区别可能是因为非线性模型的结构更为复杂，考虑了更多的因素的影响。即使初始条件相同，随时间变化小的扰动可能会被放大并且互相作用，表现出更为复杂的波动。
