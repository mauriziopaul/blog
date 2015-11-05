---
layout: post
title:  "Getting Started With Allele-Specific RNA-Seq"
categories: RNA-seq allelic-imbalance
---

1. To begin, in order to access the sequence read archives easily,
download the [NCBI SRA Toolkit][ncbi-sra-tools].

2. The sra files can be downloaded from the command line using:

{% highlight bash %}
fastq-dump -X 5 -Z SRR390728
{% endhighlight %}

and they will be stored here (on a Mac OSX machine):

{% highlight bash %}
/Users/[username]/ncbi/public
{% endhighlight %}

3. The raw data is hosted at this gene expression omnibus (GEO) [site][geo-site].

[geo-site]: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44555
[ncbi-sra-tools]: https://github.com/ncbi/sra-tools