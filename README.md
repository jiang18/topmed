# topmed
```R
library(data.table)
d = fread("../funct.snp_info.csv", sep=",", head=T)
d = d[,!grepl("extend", colnames(d))]
d[is.na(d)] = 0
fcnt = apply(d[,2:ncol(d)], 2, sum)
scnt = apply(d[,2:ncol(d)], 1, sum)


```
