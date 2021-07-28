matrix <- as.data.frame(lapply(Sys.glob("*seq2struct.dat"), read.csv,header=FALSE,sep=" "))
matrix32 <- matrix[,c(1,3,2,4)]
matrix32$V3 <- (matrix$V3)*-1
write.table(matrix32,"matrixFinalTransform.dat",sep=" ",na="",row.names=FALSE,quote=FALSE,col.names=FALSE)
