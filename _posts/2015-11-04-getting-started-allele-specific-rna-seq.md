---
layout: post
title:  "Getting Started With Allele-Specific RNA-Seq"
author: Paul L. Maurizio
date: 2015-11-04
categories: RNA-seq allelic-imbalance
use_math: true
---

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

The model is given here:

$$p(n;n_i, \pi_i, \phi)$$

Note on imprinted genes: Previous lists are available at [Gene Imprint][geneimprint], [U Otago][uotago], and [Mousebook][mousebook]. A new list from the Nat Gen study is available at [this site][new-imprinted].

[ai-paper-sra-files]: http://www.ncbi.nlm.nih.gov/sra?term=SRP056236
[natgen-ai-paper]: http://www.nature.com/ng/journal/v47/n4/full/ng.3222.html
[geo-site]: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44555
[ncbi-sra-tools]: https://github.com/ncbi/sra-tools
[geneimprint]: http://www.geneimprint.com/
[uotago]: http://igc.otago.ac.nz/
[mousebook]: http://www.mousebook.org/catalog.php?catalog=imprinting
[new-imprinted]: http://www.nature.com/ng/journal/v47/n4/extref/ng.3222-S2.xlsx