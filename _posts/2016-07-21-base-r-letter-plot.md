---
layout: post
title:  "Plotting Block Letters Using Base R Graphics"
date: 2016-07-21
categories: coding-in-r 
---

This is a code snippet I came up with that can be used to make a stylized 
letter `L` using base R.

```
rvert <- runif(1000, min=-1, max=1)
rhoriz <- runif(500, min=0, max=0.4)

pdf("TheLetterL.pdf", width=8, height=8)
plot(x=rvert, y=seq(from=0, to=1, length.out=1000), xlim=c(-2, 4), 
	pch=16, col="black", yaxt='n', xaxt='n', ylab="", xlab="")
points(x=rvert, y=seq(from=0, to=1, length.out=1000), pch=16, 
	cex=0.8, col="white")
points(y=rhoriz, x=seq(from=1, to=3, length.out=500), pch=16, 
	col="black")
points(y=rhoriz, x=seq(from=1, to=3, length.out=500), pch=16, 
	cex=0.8, col="white")
dev.off()
```
