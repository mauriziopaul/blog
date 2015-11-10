---
layout: post
title:  "Toward the Estimation of Treatment Effects in Allele-Specific RNA-Seq"
author: Paul L. Maurizio
date: 2015-11-09
categories: RNA-seq allelic-imbalance
use_math: true
---

![Founders](/blog/images/founderSNPs.png)

As a starting point for learning to analyze allele-specific RNA-seq reads, we will use existing data sets from this [Nature Genetics][natgen-ai-paper] article on allelic imbalance in crosses of wild-derived mice. The sequence read archive (sra) files from this study are available [here][ai-paper-sra-files]. We will start with WSB x PWK females, of which there are 6 files (SRR1924540, SRR1924539, SRR1924538, SRR1924535, SRR1924534, SRR1924529), and PWK x WSB females, of which there are 5 files (SRR1924516, SRR1924511, SRR1924510, SRR1924509, SRR1924508).

Unlike the original study, we will attempt to analyze allelic imbalance in the absence of any inbred strains. We will also be interested in treatment effects, but will have to resort to other data sets to take a look at modeling those effects (i.e. gene induction).

To begin, in order to access the sra files easily, download the [NCBI SRA Toolkit][ncbi-sra-tools].

After installing (decompressing) the archive and adding the package's */bin/ directory to the path, the sra files can be downloaded from the command line using:

  {% highlight bash %}
  fastq-dump -X 5 -Z SRX964397
  {% endhighlight %}

The files will be stored here (on a Mac OSX machine):

  {% highlight bash %}
  /Users/[username]/ncbi/public
  {% endhighlight %}

The raw data is hosted at this gene expression omnibus (GEO) [site][geo-site].

## Gene Expression: Total and Allele-Specific Read Counts (TReCASE)

Read counts are modeled for both autosomal and sex chromosomes. Some genes may only have a total read count (TReC) and not an allele-specific read count (ASReC / ASE), because there is a lack of SNPs in these transcripts between founder chromosomes, or because there is not sufficient depth to calculate allele-specificity.

For a given sample, $i$, and parental strain $B$, we can start by defining a binomial distribution to model the probability of the allele-specific read counts $ Y\_i $ being $ n\_{iB} $, given the total read counts $ Y\_i $. We do this for each gene, so I am omitting any gene-specific subscript. The sample-specific probability of an allele-specific read count for a given gene is $ \pi\_{iB}$, with $E(Y_i) = n_i \cdot \pi\_{iB}$, and variance $Var(Y_i)=n\_i \cdot \pi\_{iB} (1-\pi\_{iB})$:

$$Y_i \sim Bin(n_i, \pi_{iB})$$
$$Y_i \in \{0, 1, 2, ..., n_i\}$$
$$Pr(Y_i=n_{iB}; n_i, \pi_{iB}) = {n_i \choose n_{iB}} (\pi_{iB})^{n_{iB}}(1-\pi_{iB})^{n_i-n_{iB}}$$

We can then model the probability $\pi\_{iB}$ of an allele-specific read count as a beta-distribution:

$\pi\_{iB} \sim Beta(\alpha, \beta)$.

And we use a link function to estimate the genetic/POE effects on expression. Since my samples are all female, and we have no inbred mice, the strain effect can only be modeled in the reciprocals, and there is no sex-specific effect being modeled:

$$\log(\pi_{iB})=\beta_{B} + \beta_{B_{mat}}x_i$$

where $x\_i=0$ when the allele $B$ is paternal, whereas $x\_i=1$ when it is maternal.

However, we would like to adjust for the real-world scenario where the read counts are overdispersed, using the parameter ($\phi$):

$$Y_i \sim Bin(n_i, \pi_{iB}, \phi)$$
$$Pr(Y_i=n_{iB}; n_i, \pi_{iB}) = {n_i \choose n_{iB}} (\pi_{iB})^{n_{iB}}(1-\pi_{iB})^{n_i-n_{iB}}$$

$$Pr(N_B=n_B) = {n_T \choose n_B} \pi_B^{n_B}(1-\pi_B)^{n_T-n_B}$$

$$f_{BB}(n_{iB};n_i, \pi_i, \phi)={n_i \choose n_{iB}}
\frac{\prod_{k=0}^{n_{iB}-1}(\pi_i+k\phi)\pi_{k=0}^{n_i-n_{iB}-1}(1-\pi_i+k\phi)}
{\prod_{k=1}^{n_i-1}(1+k\phi)}$$

## Questions:

1. Should I simulate my data, and how should I do it?
2. How does restricting the search for DE / DAE to known imprinted genes influence my detction power?
3. Should I include a technical replicate, or replicates?
4. What about standard [spike in][spike-in] controls? 

Note on imprinted genes: Previous lists are available at [Gene Imprint][geneimprint], [U Otago][uotago], and [Mousebook][mousebook]. A new list from the Nat Gen study is available at [this site][new-imprinted].

[spike-in]: https://www.thermofisher.com/order/catalog/product/4456740
[ai-paper-sra-files]: http://www.ncbi.nlm.nih.gov/sra?term=SRP056236
[natgen-ai-paper]: http://www.nature.com/ng/journal/v47/n4/full/ng.3222.html
[geo-site]: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44555
[ncbi-sra-tools]: https://github.com/ncbi/sra-tools
[geneimprint]: http://www.geneimprint.com/
[uotago]: http://igc.otago.ac.nz/
[mousebook]: http://www.mousebook.org/catalog.php?catalog=imprinting
[new-imprinted]: http://www.nature.com/ng/journal/v47/n4/extref/ng.3222-S2.xlsx