###name fa.return.selection; Mark Sterken 2021/06/16

###input
 # A file with a chromosome identifier, start and stop loction
 # An indexed fa file from read.fa

###output
 # a selection of regions, for example for designing primers

###Description


###Function to extract sequences from an fa database file
fa.return.selection <- function(input,fa.database){
                                if(missing(input)){                     stop("provide input, chromosome, start location, end location")}
                                if(missing(fa.database)){               stop("provide output from the read.fa function")}

                                query_chromosome <- as.character(unlist(input[1]))
                                query_start_bp <-   as.numeric(as.character(unlist(input[2])))
                                query_end_bp <-     as.numeric(as.character(unlist(input[3])))

                                sequence.nu <- NULL
                                for(i in 1:length(query_chromosome)){
                                    selection.start <- which(fa.database[,1] == query_chromosome[i] & fa.database[,2] <= query_start_bp[i] & fa.database[,3] >= query_start_bp[i])
                                    selection.end   <- which(fa.database[,1] == query_chromosome[i] & fa.database[,2] <= query_end_bp[i] & fa.database[,3] >= query_end_bp[i])
    
                                    if(length(selection.start:selection.end) > 0){
                                        output <- fa.database[selection.start:selection.end,]
                                        output <- unlist(strsplit(output[,4],split=""))[output[1,2]:output[nrow(output),3] %in% query_start_bp[i]:query_end_bp[i]]
                                    } else {
                                        warning(paste("selection",paste(input[i,],collapse=" "),"does not return results"))
                                        output <- "NA"
                                    }
                                    sequence.nu <- c(sequence.nu,paste(output,collapse=""))
                                }

                                output <- data.frame(cbind(chromosome=query_chromosome,
                                                           location_start=query_start_bp,
                                                           location_end=query_end_bp,
                                                           name=paste(">",query_chromosome,":",query_start_bp,":",query_end_bp,sep=""),
                                                           sequence=sequence.nu),
                                                     stringsAsFactors = FALSE
                                                     )

                                return(output)
                               }
