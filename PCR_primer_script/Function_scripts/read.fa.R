###name

###input

###output

###Description

###Returns an indexed .fa file for convenience


read.fa <- function(file,print.progress){
                    if(missing(file)){                                          stop("provide .fa file to read")}
                    if(missing(print.progress)){                                print.progress <- TRUE}
                    if(print.progress){message("loading .fa file")}
                    data.tmp <-read.table(file,sep="\r")
                    if(print.progress){message("indexing .fa file")}
                    data.tmp <- as.matrix(data.tmp)
                    names.nu <- data.tmp[grep(">",data.tmp)]
                    names.nu <- unlist(strsplit(names.nu,split=">"))[seq(2,length(names.nu)*2,by=2)]
                    names.nu <- names.nu[rep(c(1:length(names.nu)),times=diff(c(grep(">",data.tmp),length(data.tmp)+1)))]

                    names.nu <- names.nu[-grep(">",data.tmp)]
                    data.tmp <- data.tmp[-grep(">",data.tmp)]
                    counts <- nchar(data.tmp)
                    counts.last <- NULL
                    for(i in 1:length(unique(names.nu))){
                      counts.last <- c(counts.last,cumsum(counts[names.nu==unique(names.nu)[i]]))
                    }
                    counts.first <- counts.last-(counts.last[1]-1);

                    if(print.progress){message("compiling indexed database")}
                    database <- as.data.frame(cbind(names.nu,counts.first,counts.last,data.tmp))
                    rownames(database) <- 1:nrow(database)
                    colnames(database) <- c("Name","First","Last","Sequence")
                    database[,1] <- as.character(database[,1])
                    database[,2] <- as.numeric(as.character(database[,2]))
                    database[,3] <- as.numeric(as.character(database[,3]))
                    database[,4] <- as.character(database[,4])
                    return(database)
                   }
