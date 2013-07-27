#' Read table
#' @param filename Filename of central database
#' @return Central database
#' @export
bring  <-  function(filename){
  data  <-  read.csv(paste0(dataDir,"/",filename,".csv"), h=TRUE, stringsAsFactors=FALSE, strip.white=TRUE, check.names=FALSE)
  data
}

#' Adds traits to the central dataset
#' @param data Central database
#' @param traits Traits database
#' @param matcH Name of column or vector from which to conduct the merge
#' @param targets Vector of columns names to be merged
#' @param verbose If TRUE, print messages to screen, good for isolating problems  
#' @return Central database with additional traits columns
#' @export
getTraits  <-  function(data, traits, matcH, targets, verbose=FALSE){
  if(missing(traits)) stop("Error: traits database not provided")
  if(missing(matcH)) stop("Error: matching element not provided")
  
  for(x in targets){if(verbose) cat("\nstart merging", x); data[[x]]<-traits[[x]][match(data[[matcH]], traits[[matcH]])]}
  data
}

#' Combines different traits into one
#' @param data Central database
#' @param targets Vector of columns names to be combined
#' @param unified Name of final unified trait
#' @return Central database with additional unified trait column
#' @export
combTraits  <-  function(data, unified, targets){
  un.vec  <-  vector(mode='character', length=nrow(data))
  for(k in targets) un.vec <- paste0(un.vec, "_", data[[k]])
  data[[unified]]  <-  sub("_", "", un.vec)
  
  data
}

#' Creates tables of traits (rows) and scale (columns)
#' @param data Central database
#' @param vals Values to be used in table
#' @param target Traits to be used
#' @param scale Scale to be used for columns
#' @param func Function to pe parsed to \code{tapply}
#' @param assign Should final table be assigned to an object?
#' @param tabname Name to be assigned to final table
#' @return assigns tabname a table of target (rows) x scale (columns) taking vals as input
#' @export
createTabs  <-  function(data, vals, target, scale, func, assign=FALSE, tabname){
  tab  <-  tapply(data[[vals]], list(data[[target]], data[[scale]]), func)
  tab[is.na(tab)]  <-  0
  if(assign) assign(tabname, tab, pos=1)
  tab
}

#' Plots rarefaction curves with 95 percent CI
#' @description Calculates sample-based rarefaction curves via function
#'              \code{specaccum} from package \code{vegan}.
#' @param data Central database
#' @param x an identity element
#' @param ID identity vector to be matched with x
#' @param permut Number of iterations
#' @return individual plots for each site via \code{quartz} 
#' @export
plotRarefac  <-  function(data, x, ID, permut){
  dat  <-  data[ID==x,]
  tab  <-  tapply(dat$abun, list(dat$transect_id, dat$code), sum); tab[is.na(tab)]<-0
  rar  <-  specaccum(tab, method="random", permutations=permut)
  if(Sys.info()[['sysname']]=="Darwin"){quartz()} else {windows()}
  plot(rar, ci.type="polygon", ci.col="red", las=1, main=x, xlab="Transects", ylab="Species", cex.lab=2, cex.axis=1.2, yaxs="i") 
}

#' Generate random null matrices
#' @description Uses algorithms from package \code{vegan} to generate
#'              random matrices
#' @param m A quantitative matrix
#' @param iter Number of iteractions, i.e. number of random matrices to be created 
#' @param model Function to be chosen. See \code{Details}.
#' @param ... Further arguments from \code{permatswap} or \code{permatfull}
#' @return a list with iter random matrices
#' @details If \code{model} = 1, algorithm from \code{permatswap} 
#' is called. If \code{model} = 2, algorithm from \code{permatfull} 
#' is called.
#' @export
nullMatsVegan <- function (m, iter, model, ...){
  if(model %in% 1:2 == FALSE)
    stop("Available models are numbered from 1 to 2")
  if(model == 1)
    return(permatswap(m, times=iter, ...))
  else
    return(permatfull(m, times=iter, ...))
}

#' Mgen algorithm (modified by Barneche on 2013)
#' @description This is a generic function to build null models. 
#' @param m quantitative matrix
#' @param n total number of occurrences in the matrix (if individuals, default is the total sum of the matrix)
#' @param zs logical. See \code{Details}
#' @param rep.cell logical. See \code{Details}
#' @return a randomized matrix
#' @details \code{Mgen} is general in the sense that it allows any type of probability matrix 
#'          to be used for constructing the simulated matrices. It does not, however, constrain rows 
#'          and column totals, nor does it constrain connectance. If rep.cell=TRUE, repeated interactions
#'          are added, thus generating a quantitative matrix with cell values as positive integers. 
#'          If rep.cell=FALSE, no repeated assignment of interactions is allowed, thus generating a binary 
#'          matrix of zeros and ones. Please note that when rep.cell=FALSE the number of interactions to be 
#'          assinged must be equal or lower than the number of cells in the matrix. if \code{zs} is FALSE 
#'          it does not allow any column or row to be empty, i.e. sum up zero (default);
#'          if TRUE, it does. If \code{rep.cell} is TRUE it returns a quantitative matrix (default);
#'          if FALSE, it returns a binary matrix.
#' @references "Vazquez DP, Chacoff N and Cagnolo L (2009) Evaluating multiple determinants of the structure of mutualistic networks. Ecology, 90:2039-2046"
#' @export
Mgen  <-  function (m, n=sum(m), zs = FALSE, rep.cell = TRUE) {
  if (rep.cell == FALSE & n > (nrow(m) * ncol(m))) {
    message("Argument n should be smaller than the number of cells in matrix")
  } else {
    mac  <- cumsum(m)
    mint <- matrix(0, nrow(m), ncol(m))
    if (zs == FALSE) {
      c1 <- sample(ncol(m), nrow(m), replace = TRUE, prob = colSums(m))
      for (i in 1:nrow(m)) {
        mint[i, c1[i]] <- 1
      }
      r1 <- sample(nrow(m), ncol(m), replace = TRUE, prob = rowSums(m))
      for (i in 1:ncol(m)) {
        mint[r1[i], i] <- 1
      }
    }
    rand  <-  runif((n-sum(mint)), min(mac), 1)
    for(i in rand){
      ri  <-  min(which(mac >= i))
      if (rep.cell == TRUE) 
        mint[ri] <- mint[ri] + 1
      if (rep.cell == FALSE) 
        mint[ri] <- 1
    }
    mint
  }
}

#' Iterates over Mgen algorithm
#' @description Iterates over \code{Mgen} for quantitative matrices. 
#' @param m A quantitative matrix
#' @param iter Number of iteractions, i.e. number of random matrices to be created
#' @param ... Further arguments from \code{Mgen}
#' @return a list with iter random matrices
#' @export
nullMatsMgen <- function(m, iter, ...) {
  # creating the reference probability matrix
  #on which the randomization process will be based
  pmat <- (rowSums(m)) %*% t(colSums(m))/(sum(m))^2 
  null <- list()
  for(i in 1:iter){
    null[[i]] <- Mgen(pmat, n=sum(m), ...) 
  }
  null  #all random matrices in a list
}

#' Generates null binary matrices following the null model 2 in Bascompte et al. 2003 PNAS. 
#' Written by M.M.Pires, optmized by Barneche (in progress)
#' @description This null created random matrices by randomly resorting the 1's 
#'              among the matrix cells according to marginal totals of rows and columns.
#' @param mat Binary matrix
#' @param iter Number of iterations (recommended 1000)
#' @return mat.t array of randomized matrices
#' @details Each cell has a probability of being filled that is proportional to the 
#'          number of occurrences of individuals in sites: cij = 1/2*(Pi/C + Pj/R) 
#'          where Pi= row sums; Pj = column sums; C = number of columns; and R = number of rows.
#' @references "Bascompte J, Jordano P, Melian CJ and Olesen JM (2003) The nested assembly of 
#'              plant-animal mutualistic networks. PNAS, 100(16):9383-9387"
#' @export
nullBinMats <- function(mat, iter=100){
  nR<-nrow(mat);    nC<-ncol(mat)
  mR<-rowSums(mat); mC<-colSums(mat)
  #generating a matrix with probablities for each cell
  tR    <-  rep(mR, each=nC); tC<-rep(mC, nR)
  prob  <-  matrix(((tR/nC)+(tC/nR))/2, nR, nC, byrow=TRUE)
  #filling theoretical matrices
  mat.t=array(0,c(m,n,iter)) #To store theoretical matrices
  s=1
  
  while (s<=iter){
    rand=matrix(runif(nR*nC),nR,nC)
    aux=mat.t[,,s] #avoid indexing problems
    aux[which(rand<prob)]=1 #filling empty matrix
    #avoid zeroed rows or columns 
    #rows
    rm.aux=rowSums(aux)
    cols=sample(1:nC,sum(rm.aux==0), replace=TRUE) #randomly selecting columns
    for (i in 1:sum(rm.aux==0)){
      aux[which(rm.aux==0)[i],cols[i]]=1
    }
    #columns
    cm.aux=colSums(aux)
    rows=sample(1:nR,sum(cm.aux==0), replace=TRUE) #randomly selecting rows
    for (i in 1:sum(cm.aux==0)){
      aux[rows[i],which(cm.aux==0)[i]]=1
    }
    #storing matrices
    if (sum(aux)==sum(mat)){ #whenever the resulting matrix has a different number of interactions the code runs again
      mat.t[,,s]=aux  #store the matrix within the array
      s=s+1
    }
  }  
  return(mat.t)
}

#' Performs a comparison between matrices via a given function
#' @description Are all outputs from a given function the same between two matrices?
#' @param x index of mat2
#' @param mat1 original matrix
#' @param mat2 matrix to be compared
#' @param func function used to compare 
#' @return lgical. If TRUE then matrices are identical according to the parameter requested
#' @export
compMats  <-  function(x, func, mat1, mat2){
  if(all(as.numeric(func(mat1)) == as.numeric(func(mat2[[x]]))))
    res  <-  TRUE
  else
    res  <-  FALSE
  res
}


#' Performs dignostics between original and null matrices
#' @description Performs dignostics for null models 1-4
#' @param original Original matrix
#' @param output List containing all random matrices generated via \code{nullMatsVegan}
#' @param iters Index of matrices to be compared (default uses the full length of output,
#'        i.e. original number of iterations)
#' @param func Function used to compare matrices via \code{compMats}
#' @param verbose If TRUE, print messages to screen, good for isolating problems   
#' @return Output messages from \code{compMats} via \code{sapply}
#' @export
nullDiagnostics <- function(original, output, iters, func, verbose=TRUE){
  if(class(original) != "matrix" & class(output) != "list") 
    stop("original must be a matrix and output must be a list")
  if(missing(iters)) 
    iters  <-  1:length(output)
  if(verbose)
    cat("\ntesting function", func, "\n")
  funct  <-  get(func)
  if(class(funct) != "function")
    stop(cat("Chosen function",func,"does not exist\n"))
  result  <-  invisible(unlist(lapply(iters, compMats, func=funct, mat1=original, mat2=output)))
  if(all(result)){
    if(length(iters)==length(output)){quant="All"}else{quant="Chosen"}
    cat(quant, "random matrices are identical to the original with respect to", func, "\n")
  }else{
    nums  <-  which(result != TRUE)
    cat("random matrices", nums, "differ from the original with respect to", func, "\n")
  }
}  

#' Counts number of Zeros in a matrix
#' @param matrix Input matrix
#' @return Number of zeros
#' @export
nZeros  <-  function(matrix){
  length(which(matrix !=0))
}

#' Plots differences between random and original matrix with respect to a function
#' @param x function used for comparison
#' @param original Original matrix
#' @param random One of the random matrices generated via \code{nullMatsVegan}
#' @param ... additional arguments to \code{plot}
#' @return A plot with the 1-1 line
#' @export
plotFuncDiff  <-  function(original, random, x, lim, ...){
  funct  <-  get(x)
  plot(funct(original), funct(random), ylim=c(0,lim), xlim=c(0,lim), main=x, xlab="Original matrix", ylab="1st Random matrix", ...)
  abline(0,1, lty=2, col="red")
}

#' Calculates summary of NODFs across all random matrices
#' @param matList list containing random matrices
#' @param ... optional arguments to \code{matNODF}
#' @return Summary statistics of NODF values across all random matrices in matList
#' @references ????
#' @export
summaryNODF  <-  function(matList, ...){
  sumnodf <- sapply(matList, matNODF, ...)
  quantile(sumnodf, probs=c(0.025, 0.975))
}

#' Calculates NODF of a matrix
#' @param x a matrix
#' @param ... optional arguments to \code{nestednodf}
#' @return NODF
#' @references ????
#' @export
matNODF  <-  function(x, summary=FALSE, ...){
  if(class(x)!="matrix" | length(x[x<0])>0)
    stop("x must be a matrix")
  nest <- nestednodf(x, ...)
  if(summary)
    nest  <-  nest$statistic[3]
  nest
}

#' Creates objects containing WNODF and NODF for a list of matrices
#' @param x A list of matrices
#' @param fgLevels levels of traits combinations
#' @param names Object tags to be pasted to each \code{fgLevels}
#'        which will results in the object names
#' @return assigns objects to global environment
#' @export
nestObjects  <-  function(x, fgLevels, names){
  for(i in seq_along(fgLevels)){
    obj  <-  x[[i]]
    assign(paste0("wnodf_",names[i]),        lapply(obj, matNODF, weighted=TRUE,  order=TRUE), pos=1)
    assign(paste0("wnodf_",names[i],"_sum"), lapply(obj, matNODF, weighted=TRUE,  order=TRUE, summary=TRUE), pos=1)
    assign(paste0("nodf_",names[i]),         lapply(obj, matNODF, weighted=FALSE, order=TRUE), pos=1)
    assign(paste0("nodf_",names[i],"_sum"),  lapply(obj, matNODF, weighted=FALSE, order=TRUE, summary=TRUE), pos=1)
  }
}

#' Creates objects containing WNODF and NODF for a list of matrices
#' @param matrices A list of matrices
#' @param names Tags to be pasted into each matrix object name
#' @param ... optional arguments to \code{nullMatsVegan}
#' @return assigns lists of null matrices to global environment
#' @export
createNullTran  <-  function(matrices, names, ...){
  assign(paste0("n",names,".quant"),lapply(matrices, function(mats)nullMatsVegan(m=mats, mtype="count", ...)), pos=1)
  assign(paste0("n",names,".bin"),  lapply(matrices, function(mats)nullMatsVegan(m=mats, mtype="prab", ...)), pos=1)
}

#' Reports observed WNODF and NODF against 95 percent CI from null models
#' @param bin A list of observed binary NODF summaries
#' @param nullbin A list of random-generated binary NODF summaries
#' @param quant A list of observed quantitative WNODF summaries
#' @param nullquant A list of observed quantitative WNODF summaries
#' @param output should results be written to output folder? (default=TRUE)
#' @param scales if \code{output} is TRUE, scales is incorporated into ouput name. 
#'        Must be either 'site' or 'samples'
#' @param addInfo if \code{output} is TRUE, addInfo is incorporated into ouput name.
#' @return a table of summaries comparing observed mean against 95 percent CI
#' @export
reportNestedness  <-  function(bin, nullbin, quant, nullquant, output=TRUE, scales, addInfo){
  if(!missing(scales) && scales %in% c("site","samples")==FALSE)
    stop("Error: argument scales must be either 'site' or 'samples'")
  if(!output && c(!missing(scales) | !missing(addInfo)))
    warning("Arguments scales and addInfo are only used if output=TRUE")
  obin  <-  unlist(bin)
  nbin  <-  sapply(nullbin, function(rands){summaryNODF(rands$perm, summary=TRUE)})
  oqun  <-  unlist(quant)
  nqun  <-  sapply(nullquant, function(rands){summaryNODF(rands$perm, summary=TRUE)})
  fin   <-  rbind(t(rbind(obin,nbin)), t(rbind(oqun,nqun)))
  names <-  as.character(sapply(rownames(fin), function(x)strsplit(x,"[.]")[[1]][1]))
  rownames(fin)   <- paste0(names, rep(c("_nodf","_wnodf"), each=length(obin)))
  colnames(fin)[1]<- "observed"
  if(output){
    if(missing(addInfo))
      finalName  <-  paste0(scales,"_scale_nestedness")
    else
      finalName  <-  paste0(scales,"_scale_nestedness", addInfo)
    write.csv(fin, paste0(OutsDir,"/", finalName, ".csv"))
  } else {
    fin
  }
}


#' Plot binary and quantitative matrices with a color code (in {phaget4})
#' @param x matrix
#' @param ... whatever parameters a plot function may receive
#' @return colored matrix plot
#' @details function could be retrieved by 'source("http://www.phaget4.org/R/myImagePlot.R")'
#' @references citation(package = "phaget4")
#' @export
myImagePlot <- function(x, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0,1,length=256),  # Red
                    seq(0,1,length=256),  # Green
                    seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7)
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
}