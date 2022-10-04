rep1_prefix <- "k562_pBC_exp_cBC_exp_rep1_collapsed"
rep2_prefix <- "k562_pBC_exp_cBC_exp_rep2_collapsed"
op_prefix <- "k562_pBC_exp_cBC_exp_rep12_collapsed"

np1 <- read.table(paste(rep1_prefix, "_numplasmid.tsv", sep = ""), head = T)
np2 <- read.table(paste(rep2_prefix, "_numplasmid.tsv", sep = ""), head = T)
np12 <- merge(np1, np2, by = c("pBC"))
print(nrow(np12))
print(ncol(np12))
write.table(np12, file = paste(op_prefix, "_numplasmid.tsv", sep = ""), row.names = F, quote = F, sep = "\t")

de1 <- read.table(paste(rep1_prefix, "_directexp.tsv", sep = ""), head = T)
de2 <- read.table(paste(rep1_prefix, "_directexp.tsv", sep = ""), head = T)
de12 <- merge(de1, de2, by = c("pBC"))
print(nrow(de12))
print(ncol(de12))
write.table(de12, file = paste(op_prefix, "_directexp.tsv", sep = ""), row.names = F, quote = F, sep = "\t")

ca1 <- read.table(paste(rep1_prefix, "_colannot.tsv", sep = ""), head = T)
ca2 <- read.table(paste(rep1_prefix, "_colannot.tsv", sep = ""), head = T)
ca12 <- rbind(ca1, ca2)
print(nrow(ca12))
print(ncol(ca12))
write.table(ca12, file = paste(op_prefix, "_colannot.tsv", sep = ""), row.names = F, quote = F, sep = "\t")
