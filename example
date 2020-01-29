##############################Simulações###################################
##########################################################################
setwd ("/home/anderson/sim_nadiv_high/C11/repli10")
p100 <- NULL
for (i in 1:dim(QTL100)[2]){
  p100[i] <- (as.numeric (table(QTL100[,i])[3]) + as.numeric (table (QTL100[,i])[2])/2)/as.numeric(dim(QTL100)[1])
}

summary(p100)
p100[which(is.na(p100))]<- 0
hist (p100)

dim (QTL100)
summary (p100)

dim (X1)
dim (qtlmap100)

head (rownames(QTL100))
head (rownames(X1))
head (phen1_c1$Progeny)

table (phen_b$G)

temp1 <- phen_b[which(phen_b$G>0&phen_b$G<5),]
indtemp <- sample (1:dim(temp1)[1], 3500)
temp1 <- temp1 [indtemp,]
temp2 <- phen_b[which(phen_b$G==5),]
indtemp2 <- sample (1:dim(temp2)[1], 500)
temp2 <- temp2 [indtemp2,]
phen <- rbind (temp1, temp2)
table (phen$G)



#Mantendo apenas animais selecionados
sample.ok<-match(phen$Progeny,rownames(X1))
which(is.na(sample.ok))
X <-X1[sample.ok,]
dim(X)
QTL_0.1 <- QTL100[sample.ok,]

head (rownames (X))
head (phen$Progeny)
head (rownames (QTL_0.1))
tail (rownames (X))
tail (phen$Progeny)
tail (rownames (QTL_0.1))

dim (QTL_0.1)

head (qtlmap100$ID)
head (colnames (QTL_0.1))

qtlef <- list("gamma", 0.42,1.66)
domef <- list("normal",0,1)

summary (p100)
var_a <- c(0.1,0.1, 0.1,0.3,0.3,0.3)
var_d <- c(0,0.05,0.1,0,0.15,0.3)
var_res <- c(0.9,0.85,0.8,0.7,0.55,0.4)

head (colnames (QTL_0.1))
head (qtlmap100$ID)

outall <- NULL

for (i in 1:10){
  out <- simdom(QTL=QTL_0.1,qtlmap=qtlmap100, var_a=var_a[5], var_d =var_d[5],var_eps = 0,  var_res = var_res[5],
                qtlef=qtlef,domef=domef,nDQTL = 300, p1=p100, npairs=1000, eps_ef = F)
  outall[i] <- list (out) 
}

outc1 <- outall[[1]]
outc3 <- outall[[2]]
outc5 <- outall[[3]]
outc7 <- outall[[4]]
outc9 <- outall[[5]]
outc11 <- outall[[6]]

phen1_c1 <- cbind (phen, outc1[[1]])
phen1_c3 <- cbind (phen, outc3[[1]])
phen1_c5 <- cbind (phen, outc5[[1]])
phen1_c7 <- cbind (phen, outc7[[1]])
phen1_c9 <- cbind (phen, outc9[[1]])
phen1_c11 <- cbind (phen, outc11[[1]])

qtlmap1_c1 <-  outc1[[2]]
qtlmap1_c3 <-  outc3[[2]]
qtlmap1_c5 <- outc5[[2]]
qtlmap1_c7 <-  outc7[[2]]
qtlmap1_c9 <-  outc9[[2]]
qtlmap1_c11 <-  outc11[[2]]
hist(qtlmap1_c5$m, breaks=15)

setwd ("/home/anderson/sim_nadiv_high/C1/repli10")
getwd()

write.table(phen1_c1, file = "phen1_c1.txt", sep= " ", col.names = T, row.names = F)
write.table(qtlmap1_c1, file = "qtlmap1_c1.txt", sep= " ", col.names = T, row.names = F)
write.table(snpmap10k, file = "snpmap1_c1.txt", sep= " ", col.names = T, row.names = F)
write.table(X, file = "X1_c1.txt", sep= " ", col.names = F, row.names = F)
write.table(p100, file = "p.txt", sep= " ", col.names = F, row.names = F)

