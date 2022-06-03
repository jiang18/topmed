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

######

library(data.table)
d = fread("../funct.snp_info.csv", sep=",", head=T)
d = d[,!grepl("extend", colnames(d)), with=F]
d[is.na(d)] = 0
scnt = apply(d[,-1], 1, sum)
d = cbind(0, d[,-1])
d[scnt==0, 1] = 1

sw = as.matrix(d)

# MPH output
mq = read.csv("LDL_ADJ.norm.mq.vc.csv")
mq = mq[-nrow(mq),]

# heritability estimate
lhs = t(sw) %*% sw
var = lhs %*% as.matrix(mq[,9:37]) %*% t(lhs)
h.est = lhs %*% mq$enrichment / sum(mq$m)
h.se = sqrt(diag(var)) / sum(mq$m) 

# enrichment estimate
e.est = h.est / (mq$m/mq$m[1])
e.se = h.se / (mq$m/mq$m[1])

cbind(h.est, h.se, e.est, e.se)


```
