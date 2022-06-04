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
nsnp = nrow(sw)

# trait name
trait = "HDL.norm"

# MPH output
mq = read.csv(paste(trait,"mq.vc.csv",sep="."))
mq = mq[-nrow(mq),]

# heritability estimate
lhs = t(sw) %*% sw
var = lhs %*% as.matrix(mq[,9:37]) %*% t(lhs)
h.est = lhs %*% mq$enrichment / sum(mq$m)
h.se = sqrt(diag(var)) / sum(mq$m) 

# enrichment estimate
e.est = h.est / (mq$m/nsnp)
e.se = h.se / (mq$m/nsnp)

cbind(h.est, h.se, e.est, e.se)
ee = cbind(e.est[-1], e.se[-1])
df = data.frame(FunctionalAnnotation=rownames(ee), Enrichment=ee[,1], SE=ee[,2])

library(ggplot2)
p<- ggplot(df, aes(x=FunctionalAnnotation, y=Enrichment)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Enrichment-SE, ymax=Enrichment+SE), width=.2,
                 position=position_dodge(.9)) 
p<- p + theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
p<- p+ geom_hline(yintercept=1, linetype="dashed", color = "red") + xlab("") + ggtitle(trait)
print(p)
dev.off()

z = (e.est-1)/e.se
pnorm(z, lower=F)

```
