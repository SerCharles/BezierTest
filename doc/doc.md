1.Bezier曲线转多项式曲线
$$
公式推导：
Bezier曲线方程C(t)=\sum_{i=0}^{n}P_iB_{i,n}(t)，其中P_i为第i个控制顶点坐标。\\
多项式曲线方程D(t)=\sum_{j=0}^{n}Q_jt^j，其中Q_j为第j次项系数。\\
Bernstein基函数B_{i, n}(t)=\frac{n!}{i!(n-i)!}t^i(1-t)^{n-i}\\
设其第j项系数为B_{i,n,j},根据二项式定理，
B_{i,n,j}=\left\{
\begin{array}{}
0 & {j < i}\\
(-1)^{j-i}\frac{n!}{i!(n-i)!}\dbinom{n-i}{j-i}&{j \geq i}
\end{array} \right.\\
回到多项式曲线方程，有第j次项系数Q_j=\sum_{i=0}^{n}P_iB_{i,n,j}=\sum_{i=0}^{j}P_i(-1)^{j-i}\dbinom{n-i}{j-i}\frac{n!}{i!(n-i)!}\\
根据组合数性质，\dbinom{n-i}{j-i}=\frac{(n-i)!}{(n-j)!(j-i)!}，代入上个式子，有
Q_j=\sum_{i=0}^{j}(-1)^{j-i}P_i\frac{n!}{i!(n-j)!(j-i)!}\\
算法设计：
求解这个问题需要频繁进行从0!到n!的阶乘运算，因此考虑进行预处理，提前计算0到n的阶乘并且存储。\\
之后每次计算Q_j需要计算j+1项，每项的计算复杂度为O(1)，总共需要计算\sum_{j=0}^{n}j+1项。
综上所述，算法复杂度为O(n^2)。
$$

2.多项式曲线转Bezier曲线
$$
公式推导：函数逼近的Weierstrass定理：\\
设给定函数f(x)\in C[a, b],则对\forall \epsilon>0,存在一多项式p_{\epsilon}(x),使得|f(x)-p_{\epsilon}(x)|<\epsilon对\forall x \in [a, b]一致成立。\\
若函数f在[0, 1]上连续，则Bernstein多项式B_n(f, x)=\sum_{i=0}^{n}f(\frac{i}{n})\frac{n!}{i!(n-i)!}x^i(1-x)^{n-i}满足在[0, 1]上一致收敛于f。\\
因此，取P_i=f(\frac{i}{n})=\sum_{j=0}^{n}Q_j\frac{i^j}{n^j}即可。
算法设计：对于每个j，\frac{Q_j}{n^j}的函数值可以事先求出来存储下来，需要O(n)次运算。\\
这样对于每个P_i，只需要计算i^0到i^n的n+1个数值即可，每个P_i需要O(n)次运算，综上所述，算法复杂度为O(n^2)。
$$
