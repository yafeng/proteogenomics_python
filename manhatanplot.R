library(qqman)
library(RColorBrewer)

setwd("") # set the folder where the input file is

# the input file should have four columns: "SNP","CHR","BP","ProbabilityScore"
df=read.table("ALL_hyperdiplo_lund_vardb_6rf_novpep.manplot.input.txt",header=T,sep="\t")
head(df)

pal <- brewer.pal(7,"Dark2")
pal
max(-log10(df$minPEP))
pdf("ALL_vardb_6rf_471novpeps_manplot.v2.pdf",width=10,height=7,useDingbats = FALSE)
manhattan(df,col=pal,
          suggestiveline = F,
          genomewideline = F,
          p = "minPEP", # the column name of probability score
          chrlabs=c(1:22, "X"),ylim=c(0,max(-log10(df$minPEP))+2),
          main="ALL dataset novel peptides (N=471)",
          ylab="-log10(PEP)")
dev.off()

pdf("ALL_novpep_manplot.23chrs.pdf",width=10,height=7,useDingbats = FALSE)
par(mfrow=c(2,3))
for (i in c(1:23)){
  df.chr=df[df$CHR==i,]
  manhattan(df.chr,col="black",
            suggestiveline = F,
            genomewideline = F,
            p = "minPEP",
            chrlabs=c(i),ylim=c(0,max(-log10(df.chr$minPEP))+2),
            ylab="-log10(PEP)")
}
dev.off()
