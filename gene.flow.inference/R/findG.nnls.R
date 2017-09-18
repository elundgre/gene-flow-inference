findG.nnls <- function(H,G_adj,const_coal=TRUE){
  #H is a coalescence time matrix
  #G_struct is an adjacency matrix specifying the structure of the markov chain
  #require(Matrix)
  #require(MASS)
  #require(igraph)
  #require(linprog)
  #require(quadprog)
  #require(nnls)
  
  # make sure matrix is square and find dimentions
  if(NROW(H) == NCOL(H)) n <- NROW(H)
  else stop("Matrix is not square. Please enter a symmetric matrix")
  
  # make sure G and H have the same dimensions
  if((NROW(H) != NROW(G_adj)) || (NCOL(H) != NCOL(G_adj))) stop("H and G_adj must have the same dimensions")
  
  #make sure matrix is symmetric up to computational error
  if(max(H-t(H))/max(H) > 1e-10) stop("Matrix is not symmetric. Please enter a symmetric matrix")
  
  #make sure G_adj is in sparse matrix format and add diagonal for structure matrix
  G_adj <- as(G_adj,"dgCMatrix")
  G_struct <- G_adj
  Matrix::diag(G_struct) <- 1
  #future: add check to make sure all non-zero entries are 1
  
  #construct matrix A such that Ag = -1 where g is the unknown values of G and gamma
  #not quite as stupid as the first version, but still pretty awful
  #create empty vector for gamma
  if(const_coal == TRUE) gam <- 0 else gam <- rep(0,n)
  #create empty matrix for A
  A <- Matrix::Matrix(0,nrow=n^2,ncol=n^2+length(gam),sparse=TRUE)
  for(i in 1:n){
    for(j in 1:i){
      A[(i-1)*n+j,((i-1)*n+1):(i*n)] = A[(i-1)*n+j,((i-1)*n+1):(i*n)] + G_struct[i,]*H[,j]
      A[(i-1)*n+j,((j-1)*n+1):(j*n)] = A[(i-1)*n+j,((j-1)*n+1):(j*n)] + G_struct[j,]*H[,i]
      if(i == j){
        if(const_coal==TRUE){
          A[(i-1)*n+i,n^2+1] = A[(i-1)*n+i,n^2+1] - H[i,i]
        } else A[(i-1)*n+i,n^2+i] = A[(i-1)*n+i,n^2+i] - H[i,i]
      }
    }
  }
  
  #return(A)
  
  #subtract coefficients on diagonal elements from relevant columns #this is no longer wrong, but *really* slow
  for(i in 1:n){
    A[,((i-1)*n+1:n)] = A[,((i-1)*n+1:n)] - A[,(i-1)*n+i]*apply(A[,((i-1)*n+1:n)],c(1,2),
                                                                FUN=function(x){if(x>0) return(1) else return(0)})
  }
 # return(A)
  #remove cols for diagonal entries
  A <- A[,-((1:n-1)*n+1:n)]
  
  #get rid of all zero rows and columns
  A <- A[,!apply(A==0,2,all),drop=FALSE]
  A <- A[!apply(A==0,1,all),,drop=FALSE]
  
  #return(A)
  
  #create vector of right side of eq
  b <- rep(-1,NROW(A))
  
  #square or overdetermined system
  #if(NROW(A) >= NCOL(A)){
  if(FALSE){ #deprecated
    QR <- qr(A,LAPACK=F) 
    g <- as.vector(solve(QR,b))
    #return(list(A,g))
  }
  
  #if(TRUE){ #quadratic programming solution of least squares solution
  if(FALSE){ #deprecated
    Dmat <- as.matrix(t(A) %*% A)
    dvec <- as.vector(t(A) %*% b)
    qp <- solve.QP(Dmat,dvec,Amat=Matrix::diag(NCOL(Dmat)),bvec=rep(0,NCOL(Dmat)))
    g <- as.vector(qp$solution)
    return(qp)
    nulls <- Null(Dmat)
  }
  
  if(TRUE){ #solution with non-negative least squares package (nnls)
  #if(FALSE){ #seems to provide
    nnls <- nnls::nnls(A,b)
    g <- as.vector(nnls$x)
    nulls <- MASS::Null(Matrix::t(A))
  }
  
  #underdetermined system
  #if(NROW(A) < NCOL(A)){
  if(FALSE){ #also deprecated
    Amat <- as.matrix(A) #solveLP requires input to be a matrix; system must be small
    cvec <- rep(0,NCOL(Amat)) #all feasible solutions should be valid
    lp <- solveLP(cvec,b,Amat,const.dir = rep("=",NROW(Amat)),lpSolve=TRUE)
    g <- as.vector(lp$solution)
    nulls <- Null(t(Amat))
  }
  
  #return(list(A,g,G_adj))
  
  #extract gamma from g
  gam <- g[(length(G_adj@x)+1):length(g)]
  
  gs <- g
  
  #put g back into matrix form
  G <- G_adj
  G@x <- g[1:length(G_adj@x)]
  G <- Matrix::t(G)
  g <- G@x
  Matrix::diag(G) <- -Matrix::rowSums(G)
  
  #find ordering function 
  G_ord <- G_adj
  G_ord@x <- as.numeric(1:length(G_adj@x)) #because it doesn't like integers, i guess
  G_ord <- Matrix::t(G_ord)
  ord <- as.integer(G_ord@x)
  #print(g-gs[ord])
  
  #order null space and combined parameters
  gs[1:length(G_adj@x)] <- gs[ord]
  nulls[1:length(G_adj@x),] <- nulls[ord,]
  
  return(list(A=A,G=G,gam=gam,g=g,gs=gs,nulls=nulls))
  
}