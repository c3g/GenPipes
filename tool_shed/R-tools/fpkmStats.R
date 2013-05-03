## genearte sats from cufflinks fpkm analysis
## 
## created by Maxim Caron
## adapted by Mathieu Bourgey
## 29/01/2013

covCor <- function(x, ...) UseMethod("pairsHalf")
## pairsHalf.R
## adaptation of the pairs method for fitting RNA-seq corelation graphics
## change made by mathieu Bourgey - mbourgey@genomequebec.com 
pairsHalf.formula <-
function(formula, data = NULL, ..., subset, na.action = stats::na.pass)
{
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m$... <- NULL
    m$na.action <- na.action # force in even if  default
    m[[1L]] <- as.name("model.frame")
    require(stats, quietly=TRUE)
    mf <- eval(m, parent.frame())
    pairs(mf, ...)
}

#################################################
## some of the changes are from code
## Copyright 1999 Dr. Jens Oehlschlaegel-Akiyoshi
## Others are by BDR and MM
## This version distributed under GPL (version 2 or later)
#################################################

pairsHalf.default <-
function (x, labels, panel = points, ...,
          lower.panel = panel, upper.panel = panel,
          diag.panel = NULL, text.panel = textPanel,
          label.pos = 0.5 + has.diag/3,
          cex.labels = NULL, font.labels = 1,
          row1attop = TRUE, gap = 1)
{
    textPanel <-
        function(x = 0.5, y = 0.5, txt, cex, font)
            text(x, y, txt, cex = cex, font = font)

    localAxis <- function(side, x, y, xpd, bg, col=NULL, main, oma, ...) {
      ## Explicitly ignore any color argument passed in as
      ## it was most likely meant for the data points and
      ## not for the axis.
        Axis(x, side=side, xpd=NA, ...)
	Axis(y, side=side, xpd=NA, ...)
        
    }

    localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
    localLowerPanel <- function(..., main, oma, font.main, cex.main)
        lower.panel(...)
    localUpperPanel <- function(..., main, oma, font.main, cex.main)
        upper.panel(...)

    localDiagPanel <- function(..., main, oma, font.main, cex.main)
        diag.panel(...)

    dots <- list(...); nmdots <- names(dots)
    if (!is.matrix(x)) {
        x <- as.data.frame(x)
        for(i in seq_along(names(x))) {
            if(is.factor(x[[i]]) || is.logical(x[[i]]))
               x[[i]] <- as.numeric(x[[i]])
            if(!is.numeric(unclass(x[[i]])))
                stop("non-numeric argument to 'pairs'")
        }
    } else if (!is.numeric(x)) stop("non-numeric argument to 'pairs'")
    panel <- match.fun(panel)
    if((has.lower <- !is.null(lower.panel)) && !missing(lower.panel))
        lower.panel <- match.fun(lower.panel)
    if((has.upper <- !is.null(upper.panel)) && !missing(upper.panel))
        upper.panel <- match.fun(upper.panel)
    if((has.diag  <- !is.null( diag.panel)) && !missing( diag.panel))
        diag.panel <- match.fun( diag.panel)

    if(row1attop) {
        tmp <- lower.panel; lower.panel <- upper.panel; upper.panel <- tmp
        tmp <- has.lower; has.lower <- has.upper; has.upper <- tmp
    }

    nc <- ncol(x)
    if (nc < 2) stop("only one column in the argument to 'pairs'")
    has.labs <- TRUE
    if (missing(labels)) {
        labels <- colnames(x)
        if (is.null(labels)) labels <- paste("var", 1L:nc)
    }
    else if(is.null(labels)) has.labs <- FALSE
    oma <- if("oma" %in% nmdots) dots$oma else NULL
    main <- if("main" %in% nmdots) dots$main else NULL
    if (is.null(oma)) {
        oma <- c(4, 4, 4, 4)
        if (!is.null(main)) oma[3L] <- 6
    }
    opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
    on.exit(par(opar))

    for (i in if(row1attop) 1L:nc else nc:1L)
        for (j in 1L:nc) {
            localPlot(x[, j], x[, i], xlab = "", ylab = "",
                      axes = FALSE, type = "n", ...)
            if(i == j || (i < j && has.lower) || (i > j && has.upper) ) {
                box()
			##1 row top axe
                if(i == 1  && j != 1)
                    localAxis(1 + 2*row1attop, x[, j], x[, i], ...)
			##4 row botom axe
#                 if(i == nc && (  j %% 2  || !has.upper || !has.lower ))
#                     localAxis(3 - 2*row1attop, x[, j], x[, i], ...)
# 			##1 column left axe
#                 if(j == 1  && (!(i %% 2) || !has.upper || !has.lower ))
#                     localAxis(2, x[, j], x[, i], ...)
			##last column rigth axe
                if(j == nc && i != nc)
                    localAxis(4, x[, j], x[, i], ...)
                mfg <- par("mfg")
                if(i == j) {
                    if (has.diag) localDiagPanel(as.vector(x[, i]), ...)
                    if (has.labs) {
                        par(usr = c(0, 1, 0, 1))
                        if(is.null(cex.labels)) {
                            l.wid <- strwidth(labels, "user")
                            cex.labels <- max(0.8, min(2, .9 / max(l.wid)))
                        }
                        text.panel(0.5, label.pos, labels[i],
                                   cex = cex.labels, font = font.labels)
                    }
                } else if(i < j)
                    localLowerPanel(as.vector(x[, j]), as.vector(x[, i]), ...)
                else
                    localUpperPanel(as.vector(x[, j]), as.vector(x[, i]), ...)
                if (any(par("mfg") != mfg))
                    stop("the 'panel' function made a new plot")
            } else par(new = FALSE)

        }
    if (!is.null(main)) {
        font.main <- if("font.main" %in% nmdots) dots$font.main else par("font.main")
        cex.main <- if("cex.main" %in% nmdots) dots$cex.main else par("cex.main")
        mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
    }
    invisible(NULL)
}


####################
### Main function
####################
args=commandArgs(TRUE)
paternFile=args[1]
resultDir=args[2]
outputBaseName=args[3]

master<-c()
listFile=list.files(resultDir,pattern=paternFile,recursive=T)
for(i in 1:length(listFile)) {
	
	
	nameS<-strsplit(listFile[i], "/")[1]
	d2<-read.table(paste(resultDir,listFile[i],sep="/"), header=T, sep="\t")
	colnames(d2)[9]<-paste("coverage_",nameS,sep="")
        colnames(d2)[11]<-paste("FPKM_", nameS,sep="")
	
	if(i == 1) {
		
		master<-d2[,c(1,7,8,9,11)]
		#master<-master[grep("ENS", master[,1]),]

	}
	else {
		subs<-d2[,c(1,9,11)]
		#subs<-subs[grep("ENS", subs[,1]),]
		master<-merge(master, subs, by.x=1, by.y=1)
		
	}
	cat(paste("processing",listFile[i],"\n",sep=" "))
}	
write.table(master, paste(outputBaseName,"Stat.txt",sep="."), quote=F, row.names=F, col.names=T,sep=",")



jpeg(paste(outputBaseName,"Correlation.jpeg",sep="."),1920,1080)


d1=master
if(ncol(d1)<=6) {
        stop("Not enough columns; have at least 2 samples?")
}
cd1<-d1[,c(seq(5,ncol(d1),2))]

upPan <- function(...) {
	points(..., col="black")
	abline(a=0,b=1,col="red")
}

lowPan<-function(x,y,...) {
	text(mean(par("usr")[1:2]), mean(par("usr")[3:4]), signif(cor(x,y),2),cex=2)
}

nameCol=strsplit(colnames(cd1), "_")

for(i in 1:length(nameCol)) {
        colnames(cd1)[i]<-paste(nameCol[[i]][-1],collapse="_")
}

covCor(cd1, xlim=c(0,100),ylim=c(0,100), pch=".", lower.panel=lowPan, upper.panel=upPan, main="Pearson correlation of transcript fpkm (tophat/cufflinks)")

dev.off()


