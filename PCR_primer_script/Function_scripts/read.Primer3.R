###name

###input

###output

###Description

###Read primer3 output
read.Primer3 <- function(file){
                        if(missing(file)){                                          stop("provide primer3 output file to read")}
                        input <- read.delim(file=file,header=F)
                        input <- as.character(unlist(input))

                        ###Identify separate entries "="
                        input[input =="="] <- "end=end"

                        input <- matrix(unlist(strsplit(input,split="=")),ncol=2,byrow=T)

                        primer.no <- as.numeric(as.character(input[grep("PRIMER_PAIR_NUM_RETURNED",input),2]))
                        seq.id <- as.character(input[grep("SEQUENCE_ID",input),2])

                        primer.stats <- cbind((grep("SEQUENCE_ID",input)+11),(grep("SEQUENCE_ID",input)+32)+(primer.no-1)*22)
                        primer.stats <- rbind(primer.stats,NA)

                        output <- NULL
                        for(i in 1:length(seq.id)){
                            tmp <- t(matrix(seq.id[i],ncol=23,nrow=primer.no[i]))
                            tmp[-1,] <- input[primer.stats[i,1]:primer.stats[i,2],2]
                            output <- rbind(output,t(tmp))
                        }
                        output <- output[,c(1:12,23)]
                        output <- as.data.frame(output)
                            colnames(output) <- c("Name","Penalty_pair","Penalty_left","Penalty_right","Sequence_left","Sequence_right",
                                                  "Annealing_left","Annealing_right","TM_left","TM_right","GC_percent_left","GC_percent_right",
                                                  "Product_size")

                            for(i in c(1,5:8)){
                                output[,i] <- as.character(unlist(output[,i]))
                            }
                            output <- cbind(output,matrix(unlist(strsplit(output[,7],split=",")),ncol=2,byrow=T),matrix(unlist(strsplit(output[,8],split=",")),ncol=2,byrow=T))
                            output <- output[,-7:-8]
                            colnames(output)[12:15] <- c("Annealing_left","Primer_size_left","Annealing_right","Primer_size_right")

                            for(i in c(2:4,7:15)){
                                output[,i] <- as.numeric(as.character(unlist(output[,i])))
                            }
                        return(output)
                       }