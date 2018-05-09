---
title: "Bioinformatics Class 11"
author: "Barry Grant"
date: "5/9/2018"
output: 
  html_document: 
    keep_md: yes
---



## PDB statistics

Import our PDB statistics CSV file and calculate precent structures by experimental method.


```r
p <- read.csv("Data Export Summary.csv", row.names=1)
```


```r
percent <- (p$Total / sum(p$Total)) * 100
names(percent) <- row.names(p)
percent
```

```
##               X-Ray                 NMR Electron Microscopy 
##         89.51673340          8.71321614          1.51239392 
##               Other        Multi Method 
##          0.16986775          0.08778879
```

