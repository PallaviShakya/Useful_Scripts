###name

###input

###output

###Description



###format Primer 3 can read
write.BoulderIO <- function(input,file,additives){
                            if(missing(input)){                                         stop("provide dataframe with columns: name, sequence, and comments (optional)")}
                            if(missing(file)){                                          stop("provide filename for output")}
                            if(missing(additives)){                                     additives <- c("PRIMER_PRODUCT_SIZE_RANGE=100-150",
                                                                                        "PRIMER_NUM_RETURN=3",
                                                                                        "SEQUENCE_TARGET=950,100",
                                                                                        "PRIMER_MIN_TM=58",
                                                                                        "PRIMER_MAX_TM=60","=")}
                            query_name <- as.character(unlist(input$name))
                            query_sequence <- as.character(unlist(input$sequence))

                            seq.id <- paste("SEQUENCE_ID=",query_name,sep="")
                            if("comments" %in% colnames(input)){
                              seq.id <- paste(seq.id,input$comments,sep="|")
                            }
                            seq.temp <- paste("SEQUENCE_TEMPLATE=",query_sequence,sep="")

                            output <- NULL
                            for(i in 1:nrow(input)){
                              output <- c(output,seq.id[i],seq.temp[i],additives)
                            }

                            write.table(output,file=file,quote=FALSE,sep="\r",col.names=FALSE,row.names=FALSE)
                           }