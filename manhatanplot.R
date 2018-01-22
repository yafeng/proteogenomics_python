library(qqman)
library(RColorBrewer)

setwd("") # set the folder where the input file is

# the input file should have the three columns: "chromosome","start","ProbabilityScore"
df=read.table("input.peptide.txt",header=T,sep="\t")
colnames(df) = c("CHR","BP","Prob.score")
# add additional column SNP, set it as row numbers
df$SNP  = rownames(df)
# reorder the data frame
df = df[,c(4,1,2,3)]

pal <- brewer.pal(7,"Dark2")
pdf("novpeps_manplot.all.pdf",width=10,height=7,useDingbats = FALSE)
manhattan(df,col=pal,
          suggestiveline = F,
          genomewideline = F,
          p = "Prob.score", # the column name of probability score
          chrlabs=df$CHR,ylim=c(0,max(-log10(df$Prob.score))+2),
          main=sprintf("novel peptides (N=%d)" % nrow(df)),
          ylab=("-log10(score)")

dev.off()

pdf("novpeps_manplot.individual.chromosomes.pdf",width=10,height=7,useDingbats = FALSE)
par(mfrow=c(2,2))
for (chr in df$CHR){
  df.chr = df[df$CHR==chr,]
  manhattan(df.chr,col="black",
            suggestiveline = F,
            genomewideline = F,
            p = "Prob.score",
            chrlabs=chr,ylim=c(0,max(-log10(df.chr$Prob.score))+2),
            ylab="-log10(PEP)")
}
dev.off()
