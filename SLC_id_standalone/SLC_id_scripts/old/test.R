#!/usr/bin/Rscript
library(argparser)


p=arg_parser('SLC_pre_TMHMM')
p=add_argument(p, "--meta", help="path to metadat")
argv=parse_args(p)

 

getScriptPath <- function(){
    cmd.args <- commandArgs()
    m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
    script.dir <- dirname(regmatches(cmd.args, m))
    if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
    if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
    return(script.dir)
}


#scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
#sourcePath <- dirname(scriptPath)
scriptPath <- getScriptPath()


print('script path is')
print(scriptPath)

print('first argument is')
print(argv$meta)

#print('source path is')
#print(scriptPath)
