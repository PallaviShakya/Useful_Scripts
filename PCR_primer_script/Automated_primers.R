################################################################################
###Make indel primers
################################################################################


setwd("E:/Nemawork/Projects/NIH_Ethanol/4pRIL/Indel_primers")
workwd <- getwd()

###Load pre-made functions (only for Mark; the required functions are added to your script)
    #uses eQTL pipeline functions https://git.wur.nl/mark_sterken/eQTL_pipeline
    #     transcriptomics functions https://git.wur.nl/mark_sterken/Transcriptomics.workbench
# git_dir <- "E:/Nemawork/Projects_R_zone/Git"
# source(paste(git_dir,"/NEMA_functions/Loader.R",sep=""))

    #Load custom functions
    for(i in 1:length(dir("./Function_scripts/"))){
        source(paste("./Function_scripts/", dir("./Function_scripts/")[i], sep = ""))
    }
    
################################################################################
###libraries etc
################################################################################
    install.packages("tidyverse")
    library(tidyverse)

################################################################################
###load required files
################################################################################
    
    ###the file should have at least 3 columns: chromosome, start (of indel) stop (of indel). There can be a 4th column with comments, etc. 
        gap.file <- read.delim("FH_Markerlist1.txt") %>%
                    ###we want to make sure to select sufficient region around the gap
                    dplyr::mutate(start_use=start-1000,end_use=start+1150)

    ###setup the fasta of the C. elegans genome (you can find it on ftp://ftp.wormbase.org/pub/wormbase/releases/WS276/species/c_elegans/PRJNA13758/)
    ###note that the version number matters (WSxxx), download the 'genomic.fa' file
    ###note this might take a while
        Ce_WS276 <- read.fa("c_elegans.PRJNA13758.WS276.genomic.fa")
        
################################################################################
###Get the input for primer 3 assembled
################################################################################
        
        
###this code fetches the correct sequence from the .fa file       
amplify.seq <- fa.return.selection(input=gap.file[,c(1,5,6)],fa.database=Ce_WS276)       
    amplify.seq <- data.frame(cbind(amplify.seq,comments=paste(apply(gap.file[,c(4,1,2,3)],1,paste,collapse="_"))))


    
additives <- c("PRIMER_PRODUCT_SIZE_RANGE=200-500",        ###what size you'd like your primer
               "PRIMER_NUM_RETURN=3",                      ###number of primers to return
               "SEQUENCE_TARGET=1000,150",                 ###where the gap starts and how many bases it runs ahead
               "PRIMER_MIN_TM=58",                         ###tm bottom/max range
               "PRIMER_MAX_TM=60","=")                     ##that last = is very important!!!


    amplicon <- write.BoulderIO(input=amplify.seq,file="Primer3_input_Fabries.txt",additives)
                            
###the file is also written to your work directory, there you need to find it
    
    
###at this point you need to install primer3; simply by unzipping the library. Then copy the file to the same directory on your hdd
###now you can open de command line (press start and type cmd and hit return) and navigate towards the directory, for example: 
    ###d:
    ###cd Primer3
    ###cd release-2.3.6
###now you can run primer 3, the command should look like this: primer3_core.exe < Primer3.input.txt > Primer3.output.txt
###this runs for a while and the output file you can copy again to your work directory
    
    
    selected.primers <- read.Primer3("Primer3_output_Fabries.txt"); head(selected.primers)
    
    write.table(selected.primers,file="primers_ordered.txt",sep="\t",quote=F)

        
        
        
        

    