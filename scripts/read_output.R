results = as.matrix(read.table("glob_25.fasta.struct",fill=T,header=T))
this.name = "sigma2_5"; plot(results[which(results[,this.name]!=-1),this.name],type="l")
