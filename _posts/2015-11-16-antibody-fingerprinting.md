---
layout: post
title:  "Antibody Fingerprinting"
author: Paul L. Maurizio
date: 2015-11-16
categories: antibody-footprinting DENV
use_math: true
---

In this [Science paper][science-paper] by Georgiev et al. (2013), the authors use a panel of antibody clusters against a panel of virus strains to generate a reference panel of _antibody fingerprints_ that are then used to characterize multiclonal sera. The notation is as follows. 

## Antibody Clustering
Stendarting with the antibody clustering analysis (on page 4 of the supplemental), let:

\begin{equation}
N\_j = \\{n\_{1j}, n\_{2j}, \cdots, n\_{ij}, \cdots, n\_{vj}\\} ,
\end{equation}

where $N_j$ is the neutralization fingerprint for antibody $j$, and where $n\_{ij}$ is the neutralization potency, for which antibody $j$ neutralizes virus strain $i$, and $v$ represents the total number of virus strains. Then, they transform $N_j$ into a rank vector $N\_j^R$, where the rank for strain $i$ is replaced by the rank, obtained by sorting the $n_{ij}$ values by potency (highest potency $=1$). For any two given antibodies $j_1$ and $j_2$, the Spearman rank correlation is calculated for the two rank vectors $N\_{j\_1}^R$ and $N\_{j\_2}^R$. To copy their example in `R`:

{% highlight R %}
na1 <- c(0.5, 0.3, 20.7, 11.3, 13.4)
na2 <- c(0.1, 11.4, 0.1, 0.7, 3.1)
nar1 <- rank(na1)
nar2 <- rank(na2)
cor(nar1, nar2, method="spearman")
[1] -0.4616903
{% endhighlight %}

They then use hierarchical clustering in Mathematica.  You can cluster using Manhattan distance in `R` using the following:

{% highlight R %}
d <- dist(as.matrix(mtcars), method="manhattan")
hc <- hclust(d)
plot(hc)
{% endhighlight %}

## Serum Specificity Delineation
Let $S$ be the set of all virus strains, and let $K$ be the set of antibody clusters:

\begin{equation}
S = \\{\mathsf{6101.10}, ~~\mathsf{Bal.01}, ~~\cdots, ~~\mathsf{ZM55.28a}\\}
\end{equation}

\begin{equation}
K = \\{\mathsf{VRC01-like}, ~~\mathsf{b12-like}, ~~\cdots, \mathsf{10E8-like}\\}
\end{equation}

For a given cluster $k \in K$, the neutralization fingerprint $R_k$ is defined as:

\begin{equation}
R_k = \\{\mu_{ik} \left\| i \in S \right . \\},
\end{equation}

where $\mu\_{ik}$ is the median of the ranks for strain $i$ with all antibodies $j$ within antibody cluster $k$. Let $\boldsymbol{R}$ be the matrix of representative neutralization fingerprints for all antibody clusters: $\boldsymbol{R} = \\{R\_k \left\| k \in K \right. \\}$, where

\begin{equation}
\left[R\_{k=1}\right] = \left[\mu\_{i=1, k=1} ~~ \mu\_{i=2, k=1} ~~ \cdots ~~ \mu\_{i=S\_T, k=1}\right],
\end{equation}

and where $\boldsymbol{R}$ is the reference set of epitope-specific neutralization fingerprints. This is the same as the data in the file `neut-rank.csv`, and saved as the object `abMatrixRanks`, where the first column of the `.csv` gives the strain names $S$ and the first row gives the antibody clusters $K$. Let $v$ be the total number of strains, and $K\_T$ be the total number of clusters.

$$
\boldsymbol{R} = \left( \begin{array}{cccc}
\mu_{i=1, k=1} & \mu_{i=1, k=2} & \cdots & \mu_{i=1, k=K_T} \\
\mu_{i=2, k=1} & \mu_{i=1, k=2} & \cdots & \mu_{i=2, k=K_T} \\
\vdots & \vdots & \ddots & \vdots \\
\mu_{i=v, k=1} & \mu_{i=v, k=2} & \cdots & \mu_{i=v, k=K_T} \end{array} \right)_{v \times K_T}.
$$

Otherwise represented as:

\begin{equation}
\boldsymbol{R} = \left( \begin{array}{cccc}
\left[R\_{k=1}\right]^T & \left[R\_{k=2}\right]^T & \cdots & \left[R\_{k=K_T}\right]^T \end{array} \right).
\end{equation}

Analagous to the definition ($N\_j$) of monoclonal antibody neutralization fingerprinting, let $N\_m=\\{n\_{im} \left\| i \in S \right.\\}$ be the neutralization pattern for serum $m \in M$, where $M$ represents the set of all sera being tested. You can then transform the serum neutralization pattern $N\_m$ into a rank vector $\boldsymbol{D\_m}$, based on neutralization / binding potency. Basically, at this point, they want to minimize the residual error between $\boldsymbol{D\_m}$ and the neutralization fingerprints $\boldsymbol{R}$ multiplied by some antibody cluster coefficients $\boldsymbol{C\_m} = \\{c\_k^m \left\| k \in K \right\. \\}$.

We want to minimize:

\begin{equation}
\parallel \boldsymbol{D\_m}-\boldsymbol{R} \cdot \boldsymbol{C\_m} \parallel
\end{equation}

where we have constrained $\sum_{k \in K} c_k^m =1$, and $0 \leq c_k^m \leq 1$. In our case $R$ is a $\\{v \times K\_T\\}$ matrix, $C\_m$ is a $\\{K\_T \times 1\\}$ column vector, and $D\_m$ is a $\\{1 \times v\\}$ row vector, and $\parallel x\parallel$ is the Euclidean or $\ell^2$ norm of $x$: $\parallel x\parallel = \sqrt{x_1^2 + \cdots + x_n^2}$

For instance, performing the calculation on the first row (first serum) would give us:

$$\parallel \boldsymbol{D_{m=1}}-\boldsymbol{R} \cdot \boldsymbol{C_{m=1}} \parallel,$$

or

{% highlight R %}
> seraFileRanks[,1]
   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16
17.0 18.0  2.5  7.0 15.0  2.5 20.0 21.0  2.5  5.0  8.0 19.0 13.0  6.0  2.5  9.0
  17   18   19   20   21
16.0 14.0 11.0 10.0 12.0
{% endhighlight %}


[science-paper]: http://www.ncbi.nlm.nih.gov/pubmed/23661761