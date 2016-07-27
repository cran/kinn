


biasedCoinToss <-
  function(p)
  {
    toss <- as.numeric(runif(1))
    if (toss > p)
      return (0)
    return (1)
  }


buildGraph <-
  function(X,Y,method = "det",uue = "linear")
  {
    mat    <- buildSimilarityMatrix(X)
    nedges <- numEdgesBasedOnNoise(X,Y,uue)
    if (method == "det")
      g <- simMatToDetermAdjGraph(mat,nedges,directed = F)
    return(g)
  }

buildGraphModel <-
  function(X,Y,pre.process=F)
  {
    listModel <- list()
    simClust(X) -> mX
    XY = cbind(X,Y)
    if (pre.process == T)
       simClust(XY) -> mXY
    else
      simClust(XY,1) -> mXY
    c2d <- mXY$classification
    gx <- list()
    vx <- list()
    vy <- list()
    cxy <- list()

    for (i in 1:mX$G)
    {
      id = (mX$classification == i)
      gx[[i]] = buildGraph(X[id,],Y[id])
      vx[[i]] = X[id,]
      vy[[i]] = Y[id]
      cxy[[i]] = c2d[id]

    }

    listModel$mX = mX
    listModel$mXY = mXY
    listModel$gx <- gx
    listModel$vx <- vx
    listModel$vy <- vy
    listModel$X <- X
    listModel$Y <- Y
    listModel$cxy <- cxy

    return (listModel)

  }



copyLowerToUpperTriangle <-
  function(m)
  {
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    return (m)
  }


createGraphFileName <- function(filename,i)
{
  prefix = str_extract(filename, "([a-zA-Z0-9]+)")
  graphname = str_c(prefix,as.character(i),".","gml")
  return (graphname)
}

eliminateSelfSimilarity <- function(mat)
{
  mat <- mat - diag(mat) * diag(nrow(mat))
}

eliminateUpperTriangle <-
  function(m)
  {
    m[upper.tri(m)] <- 0
    return (m)
  }
entropy2d <-
  function(cxy)
  {
    if (length(cxy) == 0)
      return (0)

    table(cxy) -> t
    e <- 0
    pv <- 0
    pe <- list()
    for (i in 1:length(t))
    {
      p <- t[i] / sum(t)
      e <- e + log2(1 / p) * p
    }
    if (e != 0)
      e <- e / log2(length(t))
    return(round(e,3))
  }



calculateEstimator <-
  function (g,simv,Y,alpha)
  {
    neighboursInGraph(g,simv) -> ng
    apply(ng,1,function(x)
      (x %*% Y) / sum(x != 0)) -> ngy
    #correcting in case no neighbours are present.
    ngy[is.na(ngy)] <- (Y[simv])[is.na(ngy)]
    (1 - alpha) * ngy + alpha * Y[simv]
  }


buildSimilarityMatrix <-
  function(x,sigma = 1)
  {
    mat <- gausskernel(X = x,sigma = sigma)
    return (mat)
  }




mostSimilarIndices <-
  function(X,x)
  {
    X <- as.matrix(X)
    x <- as.matrix(x)
    s <- buildSimilarityMatrix(rbind(X,x))
    s <- s[(nrow(X) + 1):nrow(s),1:nrow(X),drop = F]
    apply(s,1,function(l) max(l)) -> maxv
    simv<-rep(0,nrow(s))
    for (i in 1:nrow(s))
    {
      simv[i]<-which(s[i,]==max(s[i,]))
    }

    #simvec<-simvec[2:length(simvec)]
    return (simv)
  }
estimate <-
  function(g,i,x)
  {
    simv <- mostSimilarIndices(g$vx[[i]],x)
    calculateEstimator(g$gx[[i]],simv,g$vy[[i]],0.8) -> est

    return (est)

  }

#' Builds kinn regression model
#' @import mclust
#' @importFrom stringr str_extract str_c
#' @import KRLS
#' @import caTools
#' @importFrom graphics plot
#' @importFrom stats as.formula  lm predict  runif smooth.spline var
#' @param  f formula in the form y~x1+x2+...
#' @param  data dataframe which contains training data
#' @param  pre.process boolean flag which cluster data by density distributions
#' @return kinn model object
#' @export
#' @examples
#' library(kinn)
#' x<-runif(100,min = 0,max=10)
#' e<-rnorm(100)
#' y<-2*x+3+e
#' df<-data.frame(x,y)
#' model=kinn.train("y~x",df)


kinn.train <-
  function (f,data,pre.process=T)
  {
    g <- as.formula(f)
    vars=all.vars(g)
    indep <- all.vars(g)[2]
    if (indep ==  ".")
      indep = setdiff(colnames(data),vars[-1])
    else
      indep = vars[-1]

    dep=vars[1]
    X <- as.matrix(data[,indep])
    Y <- as.matrix(data[,dep])
    graph <- buildGraphModel(X,Y,pre.process)
    graph$dep <- dep
    graph$indep <- indep
    return (graph)
  }

neighboursInGraph <- function(g,simv)
{
  #if (simv == NULL)
  #if ((length(simv)==0) || (simv == NULL))
  return(g[simv,])
}

normMatrixToExpEdges <-
  function (mat,nedges,directed)
  {
    s <- sum(mat)
    proportion <- nedges / s
    if (directed == FALSE)
      proportion <- (proportion * 2)
    mat <- proportion * mat
    print(sum(mat))
    print(mat)
    return(mat)
  }


unexpLinearNoise <- function(x,y)
{
  m <- lm(y ~ x)
  summary(m) -> s
  return (1 - s$r.squared)
}



unexpSplineNoise <-
  function(x,y)
  {
    smooth.spline(x, y) -> sp
    predict(sp,x) -> p

    return (var(y - p$y) / var(y))
  }


numEdgesBasedOnNoise <-
  function(x,y,uue = "linear")
  {
    if (uue == "linear")
      noise <- unexpLinearNoise(x,y)
    else
      noise <- unexpSplineNoise(x,y)

    n <- nrow(as.matrix(x))
    edges <- max(noise * n * (n - 1) / 2,n)
    return (edges)

  }


#' Plots model graphs and (optional) saves them as files in a gml format
#' @param  gmodel kinn model generated from kinn function
#' @param  filename prefix for the subgraphs files in gml format(optional)
#' @export
#' @importFrom igraph graph.adjacency write.graph layout_with_fr
#' @importFrom caret createDataPartition
#' @importFrom graphics plot
#' @importFrom stats as.formula lm predict  runif smooth.spline var
#' @examples
#' library(kinn)
#' library(caret)
#' x<-runif(100,min = 0,max=10)
#' e<-rnorm(100)
#' y<-2*x+3+e
#' df<-data.frame(x,y)
#' inTrain <- createDataPartition(y = df$y, p = 0.7, list = FALSE)
#' train <-df[inTrain, ]
#' test <- df[-inTrain, ]
#' model=kinn.train("y~x",train)
#' kinn.plot(model,'subgraphfile')

kinn.plot <-
  function(gmodel,filename=NULL)
  {
    for (i in 1:length(gmodel$gx))
    {
      m = gmodel$gx[[1]]
      g = graph.adjacency(m,mode = "undirected",weighted = NULL)
      if (! is.null(filename))
      {
        gfilename = createGraphFileName(filename,i)
        write.graph(g, gfilename, format = c("gml"))
      }
      plot(g, layout=layout_with_fr, vertex.size=3,vertex.label.cex=0.5,

           vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.5)
    }
  }


predict.graph <-
  function(g,X)
  {
    predict(g$mX,X)$classification -> zx
    Y <- rep(0,nrow(X))
    for (i in 1:g$mX$G)
    {
      if (i %in% zx)
      {
        e <- estimate(g,i,X[i == zx])
        Y[i == zx] <- e
      }

    }
    print(Y)
    return (Y)
  }


getEstimatorsVector <-
  function (g,simv,Y,C,alpha)
  {
    y2d <- list()
    entropyEst <- list()

    neighboursInGraph(g,simv) -> ng
    ent <- 0
    probv <- list()
    t(apply(ng,1,function(x)
      (x * Y))) -> ngvy
    t(apply(ng,1,function(x)
      (x * C))) -> ngvc
    for (i in 1:nrow(ngvc))
    {
      ngvc[i,] -> line
      line[line > 0] -> line
      ent[i] <- entropy2d(line)
      probv[[i]] <- prob2d(line)
    }


    for (i in 1:nrow(ngvc))
    {
      for (j  in 1:max(C))
      {
        Y[ngvc[i,] == j] -> y2d[[j]]
      }
      if (ent[i] > 0.5)
        sig <- "*"
      else
        if (ent[i] > 0.2)
          sig <- "."
        else
          sig <- " "

        ye <- vector()

        if (probv[[i]][1] > 0)
        {
          k = 1
          for (j  in 1:max(ngvc[i,]))
          {
            if (length(y2d[[j]]) > 0)
            {
              ye[k] <- mean(y2d[[j]])
              k = k + 1
            }
          }


          if (ent[i] == 0)
          {
            s <-
              paste0("p=",round(probv[[i]],3)," yhat=",round(as.numeric(ye),3))
            ests <- paste0(sig," ","entropy=",ent[i]," ",s)
          }
          else
          {
            s <-
              paste0("(p=",round(probv[[i]],3)," yhat=",round(ye,3),")",collapse = ",")
            ests <- paste0(sig," ","entropy=",ent[i]," ",s)
          }
        }
        else
        {
          sig = " "
          s <- paste0("p=",1)
          ests <-
            paste0(sig," ","entropy=",0," ",s,"-no neighbours in graph.")

        }

        entropyEst[[i]] <- ests
    }

    return (entropyEst)


  }





estimateWithEntropy <-
  function(g,i,x)
  {
    simv <- mostSimilarIndices(g$vx[[i]],x)
    ewe <-
      getEstimatorsVector(g$gx[[i]],simv,g$vy[[i]],g$cxy[[i]],0.8)
    return (ewe)

  }


#' Predicts response vector using kinn model
#' @param  g kinn model object
#' @param  data dataframe which holds preictors
#' @param  verbose printing detalied log
#' @return predicton response vector
#' @import mclust
#' @importFrom graphics plot
#' @importFrom stats  as.formula lm predict  runif smooth.spline var
#' @importFrom stringr str_extract str_c
#' @import KRLS
#' @import caTools
#' @export
#' @examples
#' library(kinn)
#' library(caret)
#' x<-runif(100,min = 0,max=10)
#' e<-rnorm(100)
#' y<-2*x+3+e
#' df<-data.frame(x,y)
#' inTrain <- createDataPartition(y = df$y, p = 0.7, list = FALSE)
#' train <-df[inTrain, ]
#' test <- df[-inTrain, ]
#' model=kinn.train("y~x",train)
#' yhat=kinn.predict(model,test)
kinn.predict <-
  function(g,data,verbose=F)
  {
    X = as.matrix(data[,g$indep])
    predict(g$mX,X)$classification -> zx
    Y <- rep(0,nrow(X))
    YE <- rep(0,nrow(X))

    for (i in 1:g$mX$G)
    {
      if (i %in% zx)
      {
        e <- estimate(g,i,X[i == zx,])
        ewe <- estimateWithEntropy(g,i,X[i == zx,])
        Y[i == zx] <- e
        YE[i == zx] <- ewe
      }

    }
    g$Yhat <- Y
    g$YhatP <- YE
    substring(YE,1,100)->ye
    if (verbose == T)
    { message("prediction")
      for (i in  Y)
        print (i)

      message("prediction entropy")
      for (i in  ye)
        print (i)

    }

    return (Y)
  }


prob2d <-
  function(cxy)
  {
    if (length(cxy) == 0)
      return (0)
    table(cxy) -> t
    pv <- c(0)
    for (i in 1:length(t))
    {
      pv[i] <- t[i] / sum(t)
    }

    return(pv)
  }


simClust <- function(data,ngraph = 10)
{
  mixclust = Mclust(data,G = 1:ngraph)
  return (mixclust)
}

topValuesDeterminsticEdges <-
  function(m,nedges,directed)
  {
    if (directed == FALSE)
      nedges <- nedges * 2

    r <- rank(m,ties.method = "first")
    topr <- (length(r) - nedges + 1)
    w <- which(r < topr)
    r[w] <- 0
    r[-w] <- 1
    g <- matrix(r,nrow = nrow(m))

    return (g)
  }

simMatToDetermAdjGraph <-
  function(mat,nedges,directed)
  {
    mat <- eliminateSelfSimilarity(mat)
    mat <- topValuesDeterminsticEdges(mat,nedges,directed)
    return(mat)
  }


simMatToRandomAdjGraph <- function(mat,nedges,directed)
{
  mat <- eliminateSelfSimilarity(mat)

  if (directed ==  FALSE)
    m <- eliminateUpperTriangle(mat);
  mat <- normMatrixToExpEdges(mat,nedges,directed)
  mat <- apply(mat, 1:2, function(x)
    biasedCoinToss(x))

  if (directed ==  FALSE)
    mat <- copyLowerToUpperTriangle(mat);

  return (mat)

}




topValuesDeterminsticEdges <-
  function(m,nedges,directed)
  {
    if (directed == FALSE)
      nedges <- nedges * 2

    r <- rank(m,ties.method = "first")
    topr <- (length(r) - nedges + 1)
    w <- which(r < topr)
    r[w] <- 0
    r[-w] <- 1
    g <- matrix(r,nrow = nrow(m))

    return (g)
  }


