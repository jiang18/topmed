# topmed
```R
library(data.table)
d = fread("../funct.snp_info.csv", sep=",", head=T)
d = d[,!grepl("extend", colnames(d)), with=F]
d[is.na(d)] = 0
fcnt = apply(d[,2:ncol(d)], 2, sum)
scnt = apply(d[,2:ncol(d)], 1, sum)
snp = d[scnt==0, 1]
write.table(snp, file = "rest.extract.txt", quote=F, col.names=F, row.names=F)

```
