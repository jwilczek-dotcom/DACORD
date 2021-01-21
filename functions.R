#..............................................................................................
# D A C O R D  -  Fuctions
#..............................................................................................

rm(list=ls(all=TRUE))

# (1) Functions from other or obsolete packages
# mesheR (obsolete package):
glVisible <- function (mesh, offset = 0.001) {
  mesh <- vcgUpdateNormals(mesh)
  mesh0 <- meshOffset(mesh, offset)
  viewpoint <- c(glView(), 0)
  normals <- viewpoint - mesh0$vb
  mesh0$normals <- normals
  tmp <- as.logical(vcgRaySearch(mesh0, mesh)$quality)
  out <- tmp
  out[tmp] <- FALSE
  out[!tmp] <- TRUE
  return(out)
}
meshOffset <- function (mesh, offset) {
  if (is.null(mesh$normals)) 
    mesh <- vcgUpdateNormals(mesh)
  mesh$vb[1:3, ] <- mesh$vb[1:3, ] + offset * mesh$normals[1:3, 
                                                           ]
  invisible(mesh)
}
glView <- function () {
  out <- c(0, 0, 0)
  info <- rgl.projection()
  mdl <- info$model
  out[1] <- -(mdl[1] * mdl[13] + mdl[2] * mdl[14] + mdl[3] * 
                mdl[15])
  out[2] <- -(mdl[5] * mdl[13] + mdl[6] * mdl[14] + mdl[7] * 
                mdl[15])
  out[3] <- -(mdl[9] * mdl[13] + mdl[10] * mdl[14] + mdl[11] * 
                mdl[15])
  return(out)
}

# conicfit package (for 'CircleFitByKasa' function):
CircleFitByKasa <- function (XY) {
  # dependencies: pracma
  P <- mldivide(cbind(XY, 1), matrix(XY[, 1]^2 + XY[, 2]^2, ncol = 1))
  Pout = cbind(P[1]/2, P[2]/2, sqrt((P[1]^2 + P[2]^2)/4 + P[3]))
  Pout
}

# cwhmisc package (for 'rotV' function):
rotV <- function (v, w = c(0, 0, 1)) {
  u <- vec.prod(v, w)
  if (lV(u) <= .Machine$double.eps * 16) {
    res <- diag(3)
  }
  else {
    phi <- angle(v, w)
    res <- rotA(phi, u)
  }
  return(res)
}
lV <- function (v) { as.vector(sqrt(v %*% v)) }
angle <- function (v, w) {
  cs <- v %*% w
  lp <- lV(v) * lV(w)
  pc <- cs/lp
  if (abs(pc - 1) >= .Machine$double.eps * 128) {
    res <- acos(pc)
  }
  else {
    vp <- lV(vecprod(v, w))
    res <- asin(vp/lp)
  }
  return(res)
}
rotA <- function (phi, P = c(0, 0, 1)) {
  r <- rotL(-acos(P[3]/sqrt(t(P) %*% P)), 1, 3) %*% rotL(atan2(P[2], P[1]), 1, 2)
  t(r) %*% rotL(phi, 1, 2) %*% r
}
rotL <- function (phi, k = 1, m = 2, N = 3) {
  res <- diag(N)
  if (k != m) {
    ss <- sin(phi)
    cc <- cos(phi)
    res[k, k] <- res[m, m] <- cc
    res[k, m] <- ss
    res[m, k] <- -ss
  }
  return(res)
}

# SDMTools package
wt.mean <- function (x, wt) {
  s = which(is.finite(x * wt))
  wt = wt[s]
  x = x[s]
  return(sum(wt * x)/sum(wt))
}
wt.var <- function (x, wt) {
  s = which(is.finite(x + wt))
  wt = wt[s]
  x = x[s]
  xbar = wt.mean(x, wt)
  return(sum(wt * (x - xbar)^2) * (sum(wt)/(sum(wt)^2 - sum(wt^2))))
}

# 'plotrix' package
rescale <- function (x, newrange) {
  if (missing(x) | missing(newrange)) {
    usage.string <- paste("Usage: rescale(x,newrange)\n", 
                          "\twhere x is a numerical object and newrange is the new min and max\n", 
                          sep = "", collapse = "")
    stop(usage.string)
  }
  if (is.numeric(x) && is.numeric(newrange)) {
    xna <- is.na(x)
    if (all(xna)) 
      return(x)
    if (any(xna)) 
      xrange <- range(x[!xna])
    else xrange <- range(x)
    if (xrange[1] == xrange[2]) 
      return(x)
    mfac <- (newrange[2] - newrange[1])/(xrange[2] - xrange[1])
    return(newrange[1] + (x - xrange[1]) * mfac)
  }
  else {
    warning("Only numerical objects can be rescaled")
    return(x)
  }
}



# (2) DACORD functions
# 0. Support functions
deg2rad <- function(deg) { radians <- deg*pi/180; return(radians) }

rad2deg <- function(rad) { degrees <- rad/pi*180; return(degrees) }

vec.prod <- function(a, b, plots=F) {
  # Vector/cross product of two 3D vectors (return vector which is perpendicular to both vectors)
  # Last Update: 2015; 2016/12/08
  # Dependencies:  plot3d {rgl}, lines3d {rgl}
  
  # Arguments:
  #   a:        3D vector ('vector', num, 3)
  #   b:        3D vector ('vector', num, 3)
  #   plots:    visualisation of the process ('logical')
  # Value:
  #   vec.prod: Vector/cross product
  
  # examples:
  # vec.prod(a=c(2,3,4),b=c(5,6,7),plots=T)
  
  cx = a[2]*b[3] - a[3]*b[2]
  cy = a[3]*b[1] - a[1]*b[3]
  cz = a[1]*b[2] - a[2]*b[1]
  vec.prod <- c(cx,cy,cz)
  if (plots==T) {
    plot3d(0,0,0)
    lines3d(rbind(c(0,0,0),a), col="blue")
    lines3d(rbind(c(0,0,0),b), col="blue")
    lines3d(rbind(c(0,0,0),vec.prod))
  }
  return(vec.prod)
}

mesh.transform <- function(mesh, par, method=("ptab"), plots=F) {
  # Transformation of the mesh by four parameters
  # Last update: 2015; 2016/12/08
  # Dependencies:  plot3d {rgl}; par3d {}; lines3d {rgl};
  #                rotate3d {rgl}; translate3d {rgl}; vcgUpdateNormals {Rvcg}
  
  # Arguments:
  #   mesh:     triangular mesh ('mesh3d')
  #   par:      ordered transformation parameters ('vector', num, 4)
  #             par=c('phi','theta','a','b')
  #             'phi' for rotation around x-axis (in radians)
  #             'theta' for rotation around y-axis (in radians)
  #             'a'for translation along x-axis (in mm)
  #             'b'for translation along y-axis (in mm)
  #   method:   order of transformations ('vector', char, 1)
  #             'ptab' for phi,theta,a,b
  #             'tpab' for theta,phi,a,b
  #   plots:    visualisation of the process/results ('logical')
  # Value:
  #   mesh:     transformed triangular mesh ('mesh3d')
  
  # examples:
  # mesh.transform(mesh, par=c(1,1,1,1), method=("ptab"), plots=T)
  
  if (plots==T) { par3d(windowRect=c(960,30,1920,1040)); plot3d(mesh, meshColor="legacy") }
  
  phi <- par[1]
  theta <- par[2]
  a <- par[3]
  b <- par[4]
  if (method=="ptab") {
    mesh <- rotate3d(mesh,phi,1,0,0)
    mesh <- rotate3d(mesh,theta,0,1,0)
  }
  if (method=="tpab") {
    mesh <- rotate3d(mesh,theta,0,1,0)
    mesh <- rotate3d(mesh,phi,1,0,0)
  }
  mesh <- translate3d(mesh,a,b,0)
  mesh <- vcgUpdateNormals(mesh)
  if (plots==T) { plot3d(mesh, meshColor="legacy", add=T) }
  return(mesh)
}

mesh.flip.vertically <- function(mesh, plots=F) {
  # Flip mesh by xy-plane
  # Last update: 2016/10/18
  # Dependencies:  par3d {rgl}, plot3d {rgl}, updateNormals {Morpho}

  # Arguments:
  #   mesh:     triangular mesh ('mesh3d')
  #   plots:    visualisation of the process/results ('logical')
  # Value:
  #   mesh:     flipped mesh ('mesh3d')
  
  # examples:
  # mesh.flip.vertically(mesh, plots=T)
  
  if (plots==T) { par3d(windowRect=c(960,30,1920,1040)); plot3d(mesh, asp="iso", color="grey") }
  mesh$vb[3,] <- mesh$vb[3,] * (-1)
  mesh <- updateNormals(mesh)
  if (plots==T) { plot3d(mesh, meshColor="legacy", add=T) }
  return(mesh)
}

mesh.align2xy <- function(mesh, plots=F) {
  # Align mesh with xy-plane passing through z=0
  # Last Update: 2016/10/28
  # Dependencies: par3d {rgl}, plot3d {rgl}, updateNormals {Morpho}
  
  # Arguments:
  #   mesh:       triangular mesh ('mesh3d')
  #   plots:      visualisation of the process/results ('logical')
  # Value:
  #   result:     aligned mesh ('mesh3d')
  
  # examples:
  #   mesh.align2xy(mesh, plots=T)
  
  if (plots==T) { par3d(windowRect=c(960,30,1920,1040)); plot3d(mesh, asp="iso", color="grey") }
  mesh$vb[3,] <- mesh$vb[3,]-max(mesh$vb[3,])
  mesh <- updateNormals(mesh)
  if (plots==T) { plot3d(mesh, meshColor="legacy", add=T) }
  return(mesh)
  
}

rechant <- function (M, pts=1000) {
  # Re-sampling points of the 2D outline
  # Last update: 2016/10/18
  # Dependencies: spsample {sp}
  
  # Arguments:
  #   M:        matrix of xy coordinates of points on the outline ('matrix', num, Nx3)
  #   pts:      number of sampled points
  # Value:
  #   M:        matrix of xy coordinates of resampled points ('mesh3d')
  
  # examples:
  #   rechant(matrix(1:100,50,2), pts=10)
  
  if (pts <= dim(M)[1]) {
    M <- M[seq(1,dim(M)[1],length=pts),]
  }
  else if (pts == dim(M)[1]) {
    M <- M
  }
  else if (pts > dim(M)[1]) {
    Ldig <- Line(M)
    M <- spsample(Ldig, pts, type="regular", offset=c(0,1))@coords
  }
  M <- as.matrix(M)
  return(M)
}

outline <- function (M) {
  Mm <- cbind(M,c(1:(dim(M)[1])))
  y.max <-  Mm[which(Mm[,2]==max(Mm[,2])),]
  y.max <- as.matrix(y.max)
  if (length(y.max)==3) {
    y.max.ind <- y.max[3]
  }
  else if (length(y.max)>3) {
    x.mean <- mean(y.max[,1])
    y.max.ind <- y.max[which.min(abs(y.max[,1]-x.mean)),3]
  }
  y.max.ind <- as.numeric(y.max.ind)
  OUT <- M[c(y.max.ind:dim(M)[1], 1:(y.max.ind-1)),]
  return(OUT)
}

sect.hor <- function (mesh, nb, col="blue", size=1, plots=F) {
  if (nb==0) { nb <- nb_planes(mesh) }
  else { nb <- nb }
  lim <- quantile(mesh$vb[3,],c(0.05,0.95))
  cz <- seq(from=lim[1],to=lim[2], length.out=nb+2)
  cz <- cz[-c(1,length(cz))]
  sect <- matrix(numeric(0), 0,3)
  for (k in 1:(length(cz))) {
    v1 <- c(0,0,cz[k])
    v2 <- c(0,1,cz[k])
    v3 <- c(1,0,cz[k])
    xy <- tryCatch(meshPlaneIntersect(mesh,v1,v2,v3), error = function(c) matrix(0,1,3))
    if (length(xy[,1])>5) { sect <- rbind(sect,xy) }
  }
  if (plots==T) { points3d(sect, col=col, size=size) }
  return (sect)
}

sect.ver <- function (mesh, nb, by.deg=1, col="blue", size=1, plots=F) {
  mima <- rho.range(mesh,plots=F)
  if (nb==0) { rho <- seq(from=mima[1], to=mima[2], by=deg2rad(by.deg)) }
  else { rho <- seq(from=mima[1], to=mima[2], length=nb+2) }
  rho <- -rho
  rho <- rho[-c(1, length(rho))]  # take out the first and the last value
  v1 <- c(0,0,min(mesh$vb[3,]))
  v2 <- c(0,0,max(mesh$vb[3,]))
  v3 <- c(1,0,0)
  sect <- matrix(numeric(0),0,3)
  for (k in 1:(length(rho))) {
    v3.rot <- rotate3d(c(1,0,0),rho[k],0,0,1)
    pla <- tryCatch(meshPlaneIntersect(mesh,v1,v2,v3.rot), error = function(c) 0)
    if (length(pla)!=1) { sect <- rbind(sect,pla) }
  }
  if (plots==T) { points3d(sect, col=col, size=size) }
  return (sect)
}


# 1. Acquisition, importation and preparation of the model
mesh.clean <- function(mesh, it=3) {
  # Mesh cleaning
  # Last update: 2016/12/10
  # Dependencies: 'vcgClean' {Rvcg}; 'vcgUpdateNormals {Rvcg}; 'rmVertex' {Morpho}
  
  # Arguments:
  #   mesh:       triangular mesh ('mesh3d')
  #   it:         number of iterations ('vector', num, 1)
  # Value:
  #   mesh:       cleaned triangular mesh ('mesh3d' class)
  
  for(i in 1:it) {
    mesh <- vcgClean(mesh, sel=c(1:7), silent=T)   # cleaning mesh (delete manifold, duplicated vertices, etc.)
    mesh <- vcgUpdateNormals(mesh, type = 1)       # updating normals
  }
  for(i in 1:it) {
    if (sum(apply(mesh$normals[1:3,],2,sum)==0)>0) {           # if erroneous normals (c(0,0,0)) are present, erase them
      erra <- which(apply(mesh$normals[1:3,],2,sum)==0)
      mesh <- rmVertex(mesh, index=erra, keep=F)
    }
    mesh <- vcgUpdateNormals(mesh, type=1)
  }
  return(mesh)
}

mesh.clean.borders <- function(mesh, it=2, plots=F) {
  # Erasing borders of the mesh
  # Last update: 2016/10/17
  # Dependencies: vcgBorder {Rvcg}; rmVertex {Morpho}; par3d {rgl}; plot3d {rgl}; points3d {rgl}
  
  # Arguments:
  #   mesh:       triangular mesh ('mesh3d')
  #   it:         number of iterations ('vector', num, 1)
  # Value:
  #   mesh:       triangular mesh with erased borders ('mesh3d')
  
  if (plots==T) { plot3d(mesh, meshColor="legacy", type="wire") ; par3d(windowRect=c(960,30,1920,1040)) }
  for(i in 1:it) {
    b.borders <- vcgBorder(mesh)
    if (plots==T) { points3d(t(mesh$vb[1:3,])[which(b.borders$bordervb==T),], col="red", size=5) }
    mesh <- rmVertex(mesh, index=which(b.borders$bordervb==F), keep=T)
  }
  if (plots==T) { plot3d(mesh, meshColor="legacy", col="grey", add=T) }
  return(mesh)
}

decim <- function(mesh, nb) {
  # Reduce the number of points on the mesh
  # Last Update: 2016/10/10
  # Dependencies: 'vcgQEdecim' {Rvcg}
  
  # Arguments:
  #   mesh:       triangular mesh ('mesh3d' class)
  #   nb:         number of desired points on the mesh ('vector', num, 1)
  # Value:
  #   mesh:       decimated triangular mesh ('mesh3d' class)
  
  dec <- nb/dim(mesh$vb)[2]
  mesh <- vcgQEdecim(mesh, scaleindi=T, percent=dec, silent=T)
  return(mesh)
}


# Extraction of surfaces relevant to axis estimation
distNormal2axis <- function (mesh, i) {
  # The function performs the calculation of the distances between the i-th normal of the (preoriented) model and z-axis
  # Last Update: 2016/10/16
  # Dependencies: none
  # More details: 
  
  # Arguments:
  #   mesh:       pre-alignd triangular mesh, decimated if possible ('mesh3d')
  #   i:          index of the point ('vector', num, 1)
  # Value:
  #   dist:       distance between normal and the z-axis ('mesh3d')
  
  # examples:
  #   distNormal2axis(mesh=mesh, i=10)
  
  cen <- c(0,0,0)
  cen.n <- c(0,0,1)
  p <- mesh$vb[1:3,i]
  p.n <- mesh$normals[1:3,i]
  dist <- dist3Dline(cen,cen.n,p,p.n)
  return(dist)
}

distNormals2axis.parallel <- function (mesh, clust=0, infos=T) {
  # The function calculates the distances between normals of the (preoriented) model
  # (i.e. model where the rotation axis coincides with the z-axis) and rotation axis (z-axis)
  # Last Update: 2016/12/10
  # Dependencies: foreach {foreach}; makeCluster {parallel}; registerDoParallel {doParallel}; stopImplicitCluster {doParallel}
  # More details: see 'dist3Dline' and 'distNormal2axis'
  
  # Arguments:
  #   mesh:       pre-aligned triangular mesh, decimated if possible ('mesh3d')
  #   clust:      number of cores for parallel ('vector', num, 1)
  #               '0' for automatic detection of cores
  #   infos:      print the information/result ('logical')
  # Value:
  #   dist:       distances between mesh normals and the z-axis ('vector', num, N)
  
  # examples:
  #   distNormals2axis.parallel(mesh=mesh, clust=8, infos=T)
  
  if (clust==0) { clust=detectCores() }
  
  cl <- makeCluster(clust)
  registerDoParallel(cl)
  dist3Dline <- function(x0, n0, x1, n1) {
    mm <- cbind(n0,n1)
    xy <- sapply(1:3,function(x)det(mm[-x,])*(2*(x%%2)-1))
    d <- as.vector(crossprod((x0-x1),xy)/sqrt(sum(xy^2)))
    if (is.nan(d)) { d=0 }
    dis <- abs(d)
    return(dis)
  }
  distNormal2axis <- function (mesh, i) {
    cen <- c(0,0,0)
    cen.n <- c(0,0,1)
    p <- mesh$vb[1:3,i]
    p.n <- mesh$normals[1:3,i]
    dist <- dist3Dline(cen,cen.n,p,p.n)
    return(dist)
  }
  if (infos==T) { print("Parallel started") ; st <- Sys.time()}
  dist <- foreach(i=1:dim(mesh$vb)[2],.combine=c,.packages=c('Morpho', 'Rvcg')) %dopar% distNormal2axis(mesh,i)
  stopImplicitCluster()
  stopCluster(cl)
  if (infos==T) { print(paste("Parallel finished in ",round(c(Sys.time()-st)),"secs.")) }
  return (dist)
}

distNpval <- function (mesh, i, plots=F) {
  # The function projects the i-th point of the mesh (point p) to the opposite side of the model (point x). Afterwards, it calculates:
  # (i) the distance between p and x
  # (ii) the p-value of the angle between normals of these two points (pn and xn)
  # Last Update: 2016/12/07
  # Dependencies: setRays {Rvcg}, angleTest {Morpho}, plot3d {rgl}, points3d {rgl}, lines3d {rgl}, text3d {rgl}

  # Arguments:
  #   mesh:       pre-alignd triangular mesh, decimated if possible ('mesh3d')
  #   i:          index of the point ('vector', num, 1)
  #   plots:      visualisation of the process ('logical')
  # Value:
  #   result:     vector containing 2 elements ('vector', num, 2)
  #               result=c(dist,p.val)
  #               'dist' is the distance between p and x
  #               'p.val' is the p-value of the angle between pn and xn
  
  # examples:
  #   distNpval(mesh=mesh, i=10, plots=T)
  
  p <- t(mesh$vb[1:3,i])
  pn <- mesh$normals[1:3,i]
  ray <- setRays(p-pn*0.001,-pn)
  x <- vcgRaySearch(ray, mesh)
  xn <- x$normals[1:3]
  dist <- sqrt(sum((x$vb[1:3]-p)^2))
  ptest <- tryCatch(angleTest(pn,xn), error = function(c) list(angle=0, p.value=0))
  p.val <- ptest$p.value
  result <- c(dist,p.val)
  if(plots==T) {
    plot3d(mesh, type="wire", col="lightgrey"); par3d(windowRect=c(960,30,1920,1040)) 
    points3d(rbind(p), size=10, col="red")
    text3d(rbind(p), texts="p", cex=2)
    lines3d(rbind(p,p+pn*10), lwd=2, col="red")
    text3d(apply(rbind(p,p+pn*10),2,mean), texts="pn", cex=2, col="red")
    lines3d(rbind(t(ray$vb),t(ray$vb+ray$normals*10)), lwd=2, col="green")
    text3d(apply(rbind(t(ray$vb),t(ray$vb+ray$normals*10)),2,mean), texts="ray", cex=2, col="green")
    points3d(t(x$vb[1:3]), size=10, col="blue")
    text3d(t(x$vb[1:3]), texts="x", cex=2)
    lines3d(rbind(t(x$vb),t(x$vb+x$normals*10)), lwd=2, col="blue")
    text3d(apply(rbind(t(x$vb[1:3]),t(x$vb[1:3]+x$normals[1:3]*10)),2,mean), texts="xn", cex=2, col="blue")
  }
  return(result)
}

distNpval.parallel <- function (mesh, clust=0, infos=T) {
  # The function projects points of the mesh (points p) to the opposite side of the model (points x). Afterwards, it calculates:
  # (i) the distances between p and x
  # (ii) the p-values of the angle between normals of these point pairs (pn and xn)
  # Last Update: 2016/12/10
  # Dependencies: foreach {foreach}; makeCluster {parallel}; registerDoParallel {doParallel}; stopImplicitCluster {doParallel}
  # More details: see 'distNpval'
  
  # Arguments:
  #   mesh:       non-aligned triangular mesh, decimated if possible ('mesh3d')
  #   clust:      number of cores for parallel ('vector', num, 1)
  #               '0' for automatical detection of cores
  #   infos:      print the information/result ('logical')
  # Value:
  #   result:     vector containing 6 elements ('vector', num, 6)
  #               result=c(phi,theta,a,b,F1,F2)
  #               'phi' and 'theta' are initial rotation parameters
  #               'a' and 'b' are x- and y- translation parameters
  #               'F1' and 'F2' are .....
  
  # examples:
  #   distNpval.parallel(mesh=mesh, clust=8, infos=T)
  
  if (clust==0) { clust=detectCores() }
  cl <- makeCluster(clust)
  registerDoParallel(cl)
  distNpval <- function (mesh, i) {
    p <- t(mesh$vb[1:3,i])
    pn <- mesh$normals[1:3,i]
    ray <- setRays(p-pn*0.001,-pn)
    x <- vcgRaySearch(ray, mesh)
    xn <- x$normals[1:3]
    dist <- sqrt(sum((x$vb[1:3]-p)^2))
    ptest <- tryCatch(angleTest(pn,xn), error = function(c) list(angle=0, p.value=0))
    p.val <- ptest$p.value
    result <- c(dist,p.val)
    return(result)
  }
  if (infos==T) { print("Parallel started") ; st <- Sys.time()}
  res <- foreach(i=1:dim(mesh$vb)[2],.combine=c,.packages=c('Morpho', 'Rvcg')) %dopar% distNpval(mesh,i)
  stopImplicitCluster()
  stopCluster(cl)
  if (infos==T) { print(paste("Parallel finished in ",round(c(Sys.time()-st)),"secs.")) }
  return (res)
}


# The separation of the inner and outer part of the fragment
separate.byAxis <- function (mesh, plots=T) {
  # The separation of the inner and outer part of the aligned model is performed by examining
  # whether a vertex shifted by its normal is closer to or more distant from the rotation axis
  # Last Update: 2016/10/16
  # Dependencies: rmVertex {Morpho}, plot3d {rgl}, lines3d {rgl}

  # Arguments:
  #   mesh:         pre-aligned triangular mesh, decimated and cleaned if possible ('mesh3d')
  #   plots:        visualisation of result ('logical')
  # Value:
  #   mesh.triple:    the 'list' containing 3 triangular meshes ('list')
  #                   'mesh$mesh, 'mesh$S1' and 'mesh$S2' - all ('mesh3d')
  
  # examples:
  #   str(separate.byAxis(mesh=mesh, plots=T))
  
  r <- sqrt(mesh$vb[1,]^2+mesh$vb[2,]^2)
  r1 <- mesh$vb[1:3,]-mesh$normals[1:3,]*0.001
  r1 <- sqrt(r1[1,]^2+r1[2,]^2)
  mesh1 <- rmVertex(mesh, index=which(r<r1), keep=T)
  mesh2 <- rmVertex(mesh, index=which(r>r1), keep=T)
  if (plots==T) {
    plot3d(mesh, type="wire", col="grey", asp="iso"); par3d(windowRect=c(960,30,1920,1040)) 
    plot3d(mesh1, col="red", add=T)
    plot3d(mesh2, col="blue", add=T)
    lines3d(rbind( c(0,0,min(mesh$vb[3,])), c(0,0,max(mesh$vb[3,])) ), lwd=2, col="red")
  }
  mesh.triple <- list(S1=mesh1, S2=mesh2, mesh=mesh)
  return(mesh.triple)
}





# 2.1. MANUAL MODEL PREORIENTATION
plane3p <- function (A, B, C, infos=F) {
  # The Cartesian equation of a plane using three non-collinear points
  # Last Update: 2016/10/10
  # Notes: function 'crossprod' {base} can be used instead of 'vec.prod'

  # Arguments:
  #   A:          3D coordinates of the first point ('vector')
  #   B:          3D coordinates of the second point ('vector')
  #   C:          3D coordinates of the third point ('vector')
  #   infos:      whatever the result is printed ('logical')
  # Value:
  #   plane3p:    the Cartesian equation of the plane
  
  # Examples:
  #   plane3p(A=c(1,1,1),B=c(-1,1,0),C=c(2,0,3),infos=T)
  
  Ax=A[1]; Ay=A[2]; Az=A[3]    
  Bx=B[1]; By=B[2]; Bz=B[3]
  Cx=C[1]; Cy=C[2]; Cz=C[3]
  i <- (Bx-Ax)
  j <- (By-Ay)
  k <- (Bz-Az)
  AB <- c(i,j,k)
  ii <- (Cx-Ax)
  jj <- (Cy-Ay)
  kk <- (Cz-Az)
  AC <- c(ii,jj,kk)
  norm.vect <- vec.prod(AB,AC)
  a <- norm.vect[1]
  b <- norm.vect[2]
  c <- norm.vect[3]
  d <- sum(norm.vect)
  if (infos==T) { print(paste(a,"x + ", b,"y + ", c,"z + ", d," = 0", sep="")) }
  plane3p <- c(a,b,c,d)
  return(plane3p)
}

read.pp <- function(FileName) {
  temp <- scan(FileName,what="list",sep="\n",quiet=TRUE)
  Row <- M <- Name <- tp <- Mtp <- c()
  cpt_pts <- 0
  entetes <- termin <- list()
  for (i in 1:length(temp)) {
    if (length(grep("<point",temp[[i]]))>0) {
      cpt_pts <- cpt_pts+1
      dx <- gregexpr("x=\"",temp[[i]])[[1]][1]+3
      fx <- gregexpr("\" y=",temp[[i]])[[1]][1]-1
      Row[1] <- as.numeric(substr(temp[[i]],dx,fx))
      tp[1] <- substr(temp[[i]],1,dx-1)
      dy <- gregexpr("y=\"",temp[[i]])[[1]][1]+3
      fy <- gregexpr("\" z=",temp[[i]])[[1]][1]-1
      Row[2] <- as.numeric(substr(temp[[i]],dy,fy))
      tp[2] <- substr(temp[[i]],fx+1,dy-1)
      dz <- gregexpr("z=\"",temp[[i]])[[1]][1]+3
      fz <- gregexpr("\" active=",temp[[i]])[[1]][1]-1
      Row[3] <- as.numeric(substr(temp[[i]],dz,fz))
      tp[3] <- substr(temp[[i]],fy+1,dz-1)
      dn <- gregexpr("name=\"",temp[[i]])[[1]][1]+6
      fn <- gregexpr("\"/>",temp[[i]])[[1]][1]-1
      tp[4] <- substr(temp[[i]],fz+1,dn-1)
      tp[5] <- substr(temp[[i]],fn+1,nchar(temp[[i]]))
      M <- rbind(M,Row)
      Name[cpt_pts] <- as.character(substr(temp[[i]],dn,fn))
      Mtp <- rbind(Mtp,tp)
      
    }else{
      if (cpt_pts==0){
        entetes[[i]] <- temp[[i]]
      } else {
        termin[[i-cpt_pts-length(entetes)]] <- temp[[i]]
      }
    }
  }
  rownames(M) <- Name
  return(list(entetes=entetes,M=M,termin=termin,Mtp=Mtp))
}

Clear3d <- function(type = c("shapes", "bboxdeco", "material"), defaults = getr3dDefaults(), subscene = 0) {
  d <- .check3d()
  rgl.clear(type, subscene = subscene)
  return(d)
}

vcgRaySearch2 <- function(x, mesh, mintol = 0, maxtol = 1e+15, mindist = FALSE) {
  if (!inherits(mesh, "mesh3d") || !inherits(x, "mesh3d"))
    stop("arguments 'x' and 'mesh' need to be objects of class 'mesh3d'")
  mesh <- meshintegrity(mesh, facecheck = TRUE)
  x <- meshintegrity(x, normcheck = TRUE)
  vb <- mesh$vb[1:3, , drop = FALSE]
  it <- mesh$it - 1
  dimit <- dim(it)[2]
  dimvb <- dim(vb)[2]
  storage.mode(it) <- "integer"
  clost <- x$vb[1:3, , drop = FALSE]
  normals <- x$normals[1:3, , drop = FALSE]
  clostDim <- ncol(clost)
  maxtol <- as.numeric(maxtol)
  mintol <- as.numeric(mintol)
  mindist <- as.logical(mindist)
  tmp <- .Call("Rintersect", vb, it, clost, normals, mintol,
               maxtol, mindist)
}

PlacePt <- function(x,y,verts,norms,mesh,Start) {
  temp<-rgl.user2window(x=verts[,1],y=verts[,2],z=verts[,3],projection=Start$projection)
  X<-temp[,1]*Start$viewport[3]
  Y<-(1-temp[,2])*Start$viewport[4]
  X1<-x/Start$viewport[3]
  Y1<-1-y/Start$viewport[4]
  mX<-matrix(X[mesh$it],nrow=3,byrow=FALSE)
  mY<-matrix(Y[mesh$it],nrow=3,byrow=FALSE)
  d<-mX
  d[1,]<-sqrt((mX[1,]-mX[2,])^2+(mY[1,]-mY[2,])^2)
  d[2,]<-sqrt((mX[1,]-mX[3,])^2+(mY[1,]-mY[3,])^2)
  d[3,]<-sqrt((mX[2,]-mX[3,])^2+(mY[2,]-mY[3,])^2)
  dm<-max(d)
  sqIdx<- X>=x-dm & X<=x+dm & Y>=y-dm & Y<=y+dm
  SelecMesh<-rmVertex(mesh,which(sqIdx),keep=TRUE)
  Q<-rgl.window2user(X1,Y1,0)
  normQ<-rgl.window2user(X1,Y1,1)-rgl.window2user(X1,Y1,0)
  normQ<-normQ/sqrt(sum(normQ^2))
  lQ<-list(vb=matrix(c(Q,1),4,1),normals=matrix(c(normQ,1),4,1))
  class(lQ)<-"mesh3d"
  t1<-Sys.time()
  int<-vcgRaySearch(lQ, mesh)
  t2<-Sys.time()
  print(t2-t1)
  cas<-1
  if (int$quality==0){
    t1<-Sys.time()
    temp2<-rgl.user2window(x=verts[,1]+norms[,1],y=verts[,2]+norms[,2],z=verts[,3]+norms[,3],projection=Start$projection)
    normals<-temp2-temp
    u<-par3d()$observer
    alpha<-acos((t(u)%*%t(normals))/(sqrt(rowSums(normals^2))*sqrt(sum(u^2))))
    Idx<-alpha>pi/2
    Xs<-X[Idx]
    Ys<-Y[Idx]
    Dist<-sqrt((Xs-x)^2+(Ys-y)^2)
    idx<-which.min(Dist)
    int<-rmVertex(mesh,which(Idx)[idx],keep=TRUE)
    t2<-Sys.time()
    print(t2-t1)
    cas<-2
  }
  if (length(SelecMesh$vb)==0 | !is.matrix(SelecMesh$it)){
    Visibles<-idx<-NULL
  } else {
    isol<-vcgIsolated(SelecMesh,split=TRUE)
    vd<-matrix(0,length(isol),2)
    for (i in 1:length(isol)){
      temp<-sqrt(colMeans((isol[[i]]$vb-matrix(int$vb,4,dim(isol[[i]]$vb)[2]))^2))
      vd[i,1]<-min(temp)
      vd[i,2]<-which.min(temp)
    }
    idx<-which(vd[,1]==min(vd[,1]))
    Visibles<-cbind(isol[[idx]]$vb[1,],isol[[idx]]$vb[2,],isol[[idx]]$vb[3,])
    idx<-vd[idx,2]
  }
  return(list(Visibles=Visibles,idx=idx))
}

rgl.select2<-function(button = c("left", "middle", "right"),verts,norms,mesh,modify,A=NULL,vSp=NULL,dev,IdxPts=NULL,vTx=NULL){
  Start <- list()
  Sp<-Tx<-idx<-NULL
  firstTime<-TRUE
  Begin <- function(x, y) {
    Start$viewport <<- par3d("viewport")
    Start$projection <<- rgl.projection()
    temp<-PlacePt(x,y,verts,norms,mesh,Start)
    Visibles<<-temp$Visibles
    idx<<-temp$idx
    if (modify){
      distPt<-sqrt(rowSums((A-matrix(Visibles[idx,],dim(A)[1],dim(A)[2],byrow=TRUE))^2))
      Idx<<-which(distPt==min(distPt))
      IdxPts<<-Idx
      if (firstTime){
        rgl.pop("shapes",vSp[Idx])
        rgl.pop("shapes",vTx[Idx])
        firstTime<<-FALSE
      }else{
        rgl.pop("shapes",Sp)
        rgl.pop("shapes",Tx)
      }
    }else{
      if (!is.null(Sp)){
        rgl.pop("shapes",Sp)
        rgl.pop("shapes",Tx)
      }
      Idx<<-NULL
    }
    Sp<<-spheres3d(x=Visibles[idx,1],y=Visibles[idx,2],z=Visibles[idx,3],alpha=0.5)
    Tx<<-text3d(x=Visibles[idx,1],y=Visibles[idx,2],z=Visibles[idx,3],texts=as.character(IdxPts))
  }
  Update <- function(x, y) {
    temp<-PlacePt(x,y,verts,norms,mesh,Start)
    Visibles<<-temp$Visibles
    idx<<-temp$idx
    if (!is.null(Sp)){
      rgl.pop("shapes",Sp)
      rgl.pop("shapes",Tx)
    }
    Sp<<-spheres3d(x=Visibles[idx,1],y=Visibles[idx,2],z=Visibles[idx,3],alpha=0.5)
    Tx<<-text3d(x=Visibles[idx,1],y=Visibles[idx,2],z=Visibles[idx,3],texts=as.character(IdxPts))
  }
  End <- function(x,y) {}
  button <- match.arg(button)
  newhandler <- par3d("mouseMode")
  newhandler[button] <- "selecting"
  oldhandler <- par3d(mouseMode = newhandler)
  on.exit(par3d(mouseMode = oldhandler))
  rMul<-rgl.setMouseCallbacks(2, Begin, Update, End)
  dev<-rgl.cur()
  while (dev==rgl.cur()){
    if (!is.null(idx)){
      result <- rgl:::rgl.selectstate()
      if (result$state >= rgl:::msDONE)
      {break}
    }
  }
  if (dev!=rgl.cur()) {
    if (modify) {
      isDone<-TRUE
      return(list(isDone=isDone,isClosed=TRUE))
    }
  }
  rgl.setMouseCallbacks(2)
  rgl:::rgl.setselectstate("none")
  if (result$state == rgl:::msDONE){
    isDone<-FALSE
    return(list(isDone=isDone,isClosed=FALSE))
  }else {
    isDone<-TRUE
    return(list(isDone=isDone,Pt=Visibles[idx,],sp=Sp,Idx=Idx,isClosed=FALSE,tx=Tx))
  }
}

Selectpoints3d <- function (objects, norms, mesh,modify,A=NULL,vSp=NULL,dev=NULL,IdxPts=NULL,vTx=NULL) {
  StopPts<-0
  while (StopPts==0) {
    verts <- rgl.attrib(objects, "vertices")
    temp<-rgl.select2(button="right",verts=verts,norms=norms,mesh=mesh,modify=modify,A=A,vSp=vSp,dev=dev,IdxPts=IdxPts,vTx=vTx)
    if (temp$isClosed | temp$isDone) {
      if (temp$isClosed) {
        rgl.close()
      }
      break
    }
  }
  res<-list(coords=temp$Pt,sp=temp$sp,Idx=temp$Idx,isClosed=temp$isClosed,tx=temp$tx)
  return(res)
}

submesh <- function(spec,keep) {
  spec2<-list()
  spec2$vb<-spec$vb[,keep]
  idxV<-is.element(spec$it,which(keep))
  idxV<-matrix(idxV,dim(spec$it)[1],dim(spec$it)[2],byrow=FALSE)
  M<-spec$it[,which(colSums(idxV)==3)]
  spec2$it<-matrix(match(M,which(keep)),dim(M)[1],dim(M)[2])
  spec2$normals<-spec$normals[,keep]
  spec2$material<-spec$material
  if (length(spec$material)>0){
    spec2$material$color<-spec$material$color[,which(colSums(idxV)==3)]
  }
  class(spec2) <- "mesh3d"
  return(spec2)
}

SetPtZoom <- function(dd,specFull,Trans,Pt,ptsize,modify=FALSE,A=NULL,vSp=NULL,dev=NULL,IdxPts=NULL,wR2=NULL) {
  keep<- dd<0.15*max(dd)
  specFull2<-submesh(specFull,keep)
  temp<-vcgIsolated(specFull2,split=TRUE)
  Min<-+Inf
  for (ii in 1:length(temp)){
    vb<-temp[[ii]]$vb[1:3,]
    cs<-sqrt(colSums((vb-(Pt+Trans))^2))
    MinT<-cs[which(cs==min(cs))]
    if (MinT<Min){
      specFull2<-temp[[ii]]
      Min<-MinT
    }
  }
  specimen2<-specFull2$vb[1:3,]
  Trans2<-rowMeans(specimen2)
  specimen2<-specimen2-Trans2
  specFull2$vb[1:3,]<-specimen2
  param3d<-par3d()
  d2<-open3d()
  par3d(windowRect=wR2)
  ids2<-plot3d(specimen2[1,],specimen2[2,],specimen2[3,],size = ptsize,aspect = FALSE); par3d(windowRect=c(960,30,1920,1040)) 
  mesh2<-specFull2
  if (is.null(mesh2$material)) {
    mesh2$material <- "gray"
  }
  shade3d(mesh2)
  rgl.viewpoint(userMatrix=param3d$userMatrix)
  res2<-Selectpoints3d(ids2["data"],norms=t(specFull2$normals[1:3,]),mesh=specFull2,modify,A,vSp,dev,IdxPts)
  rgl.close()
  return(list(coords=res2$coords,sp=res2$sp,Trans2=Trans2,tx=res2$tx))
}

Digit.Fixed.New <- function (spec, specFull, fixed, TemplateFile = NULL, IdxPtsTemplate = 1:3, index = 1:fixed, ptsize = 1, center = TRUE, wR1=c(0,50,830,904),wR2=c(840,50,1672,904)) {
  require(rgl)
  require(Rvcg)
  require(sp)
  require(Morpho)
  if (is.null(TemplateFile)){
    IdxPtsTemplate<-1:fixed
  }else{
    if (is.null(TemplateFile)){
      warning("Template data required. Execution stops...")
      return(0)
    }else{
      p1<-length(IdxPtsTemplate)
      Template<-list()
      Template$M<-TemplateFile
    }
  }
  spec.name <- deparse(substitute(spec))
  mesh <- NULL
  Trans <- 0
  if (center == TRUE) {
    specimen <- scale(as.matrix(t(spec$vb)[, -4]), scale = FALSE)
    spec$vb <- rbind(t(specimen), 1)
    Trans<-attr(specimen,"scaled:center")
  }
  if (center == FALSE) {
    specimen <- as.matrix(t(spec$vb)[, -4])
  }
  mesh <- spec
  if (is.null(mesh$material)) {
    mesh$material <- "gray"
  }
  d1<-Clear3d()
  par3d(windowRect=wR1)
  ids1<-plot3d(specimen[, 1], specimen[, 2], specimen[, 3], size = 1, col="white", alpha=0, aspect = FALSE); par3d(windowRect=c(960,30,1920,1040))     ################# size=ptsize
  if (!is.null(mesh)) {
    shade3d(mesh, add = TRUE)
  }
  A<-Adeci<-matrix(NA,fixed,3)
  rownames(A)<-1:fixed
  colnames(A)<-c("x","y","z")
  Idx<-setdiff(1:fixed,index[IdxPtsTemplate])
  vSp<-vTx<-Sp<-Tx<-rep(NA,fixed)
  for (i in 1:fixed){
    if (i<=length(IdxPtsTemplate)){
      idx_pts<-index[IdxPtsTemplate[i]]
      res<-Selectpoints3d(objects=ids1["data"],norms=t(spec$normals[1:3,]),mesh=spec,modify=FALSE,A=A,dev=rgl.cur(),IdxPts=idx_pts)
      Sp[idx_pts]<-res$sp
      Tx[idx_pts]<-res$tx
      Pt<-res$coords
      dd<-sqrt(colSums((specFull$vb[1:3,]-Trans-Pt)^2))
      res2<-SetPtZoom(dd,specFull,Trans,Pt,ptsize,dev=rgl.cur(),IdxPts=idx_pts,wR2=wR2)
      Trans2<-res2$Trans2
      A[idx_pts,]<-res2$coords+Trans2-Trans
      rgl.set(d1)
      rgl.pop("shapes",Sp[idx_pts])
      rgl.pop("shapes",Tx[idx_pts])
      pts<-projRead(t(A[idx_pts,]),spec,sign=FALSE)      # sign ??????????????????????????
      Adeci[idx_pts,]<-pts$vb[1:3]
      sp<-spheres3d(pts$vb[1:3,],color = "green",alpha=0.5)
      tx<-text3d(pts$vb[1:3,],texts=as.character(idx_pts),col="red",cex=2)
      vSp[idx_pts]<-sp
      vTx[idx_pts]<-tx
      if(is.null(TemplateFile)==FALSE & i==length(IdxPtsTemplate)){
        configA<-as.matrix(A[is.na(A[,1])==FALSE,])
        configB<-as.matrix(Template$M)
        p2<-dim(configB)[1]
        configC<-configB[IdxPtsTemplate,]
        transA<-colMeans(configA)
        AA<-configA-matrix(transA,p1,3,byrow=TRUE)
        scaleA<-1/sqrt(sum(AA^2))
        AA<-AA*scaleA
        transB<-colMeans(configC)
        BB<-configC-matrix(transB,p1,3,byrow=TRUE)
        scaleB<-1/sqrt(sum(BB^2))
        BB<-BB*scaleB
        sv<-svd(t(AA)%*%BB)
        U<-sv$v
        V<-sv$u
        sig<-sign(det(t(AA)%*%BB))
        V[,3]<-sig * V[,3]
        rot<- U%*%t(V)
        BB<-configB
        BB<-BB-matrix(transB,p2,3,byrow=TRUE)
        BB<-BB*scaleB
        BB<-BB%*%rot
        BB<-BB/scaleA
        BB<-BB+matrix(transA,p2,3,byrow=TRUE)
        B<-BB
        B[is.na(A[,1])==FALSE]<-A[is.na(A[,1])==FALSE]
        ptsB<-projRead(B,spec,sign=FALSE) # sign ??????????????????????????
        vv<-index[Idx]
        for (ii in 1:length(vv)){
          vSp[vv[ii]]<-spheres3d(t(ptsB$vb[1:3,vv[ii]]),color = "blue", radius = 1,alpha=0.5)
          vTx[vv[ii]]<-text3d(t(ptsB$vb[1:3,vv[ii]]),texts=as.character(vv[ii]),col="red",cex=2)
        }
      }
    }else{
      idx_pts<-index[Idx[i-length(IdxPtsTemplate)]]
      dd<-sqrt(colSums((specFull$vb[1:3,]-Trans-B[idx_pts,])^2))
      Pt<-B[idx_pts,]+Trans
      res2<-SetPtZoom(dd,specFull,Trans,Pt,ptsize,IdxPts=idx_pts,wR2=wR2)
      Trans2<-res2$Trans2
      A[idx_pts,]<-res2$coords+Trans2-Trans
      rgl.set(d1)
      rgl.pop("shapes",vSp[idx_pts])
      rgl.pop("shapes",vTx[idx_pts])
      pts<-projRead(t(A[idx_pts,]),spec,sign=FALSE)      # sign ??????????????????????????
      Adeci[idx_pts,]<-pts$vb[1:3]
      vSp[idx_pts]<-spheres3d(pts$vb[1:3,],color = "green", radius = 1,alpha=0.5)
      vTx[idx_pts]<-text3d(pts$vb[1:3,],texts=as.character(idx_pts),col="red",cex=2)
    }
  }
  Stop<-0
  while((Stop==0)){
    res<-Selectpoints3d(ids1["data"],norms=t(spec$normals[1:3,]),mesh=spec,modify=TRUE,Adeci,vSp,dev=d1,vTx=vTx)
    if (res$isClosed) {break}
    idx_pts<-res$Idx
    Sp[idx_pts]<-res$sp
    Tx[idx_pts]<-res$tx
    Pt<-res$coords
    dd<-sqrt(colSums((specFull$vb[1:3,]-Trans-Pt)^2))
    res2<-SetPtZoom(dd,specFull,Trans,Pt,ptsize,IdxPts=idx_pts,wR2=wR2)
    Trans2<-res2$Trans2
    A[idx_pts,]<-res2$coords+Trans2-Trans
    rgl.set(d1)
    rgl.pop("shapes",Sp[idx_pts])
    rgl.pop("shapes",Tx[idx_pts])
    pts<-projRead(t(A[idx_pts,]),spec,sign=FALSE)
    Adeci[idx_pts,]<-pts$vb[1:3]
    sp<-spheres3d(pts$vb[1:3,],color = "green",alpha=0.5)
    vSp[idx_pts]<-sp
    tx<-text3d(pts$vb[1:3,],texts=as.character(idx_pts),col="red",cex=2)
    vTx[idx_pts]<-tx
  }
  return(list(A=A))
}


digit.fixed.lite <- function (spec, fixed, ptsize = 1, center = FALSE) {
  # Slighly altered digit.fixed {geomorph} Digitize 3D landmarks on mesh3d object
  # Notes:  There was problem with 'Digit.Fixed.New' function. To solve this issue,
  #         the original 'digit.fixed' function from 'geomorph' package was used instead
  #         However this creates problems because clicked point is sometimes identified on another surface

  spec.name <- deparse(substitute(spec))
  mesh <- NULL
  if (inherits(spec, "shape3d") == TRUE || inherits(spec, "mesh3d") == TRUE) {
    if (center == TRUE) {
      specimen <- scale(as.matrix(t(spec$vb)[, -4]), scale = FALSE)
      spec$vb <- rbind(t(specimen), 1)
    }
    if (center == FALSE) {
      specimen <- as.matrix(t(spec$vb)[, -4])
    }
    mesh <- spec
    if (is.null(mesh$material)) {
      mesh$material <- "gray"
    }
  }
  else if (inherits(spec, "matrix") == FALSE) {
    stop("File is not a shape3d/mesh3d object or xyz matrix")
  }
  else if (inherits(spec, "matrix") == TRUE && dim(spec)[2] == 3) {
    if (center == TRUE) {
      specimen <- scale(spec, scale = FALSE)
    }
    if (center == FALSE) {
      specimen <- spec
    }
  }
  else {
    stop("File is not matrix in form: vertices by xyz")
  }
  clear3d()
  ids <- plot3d(specimen[, 1], specimen[, 2], specimen[, 3], size = ptsize, aspect = FALSE)
  if (!is.null(mesh)) {
    rgl.bringtotop(stay = TRUE)
    shade3d(mesh, meshColor = "legacy", add = TRUE)
  }
  selected <- matrix(NA, nrow = fixed, ncol = 3)
  fix <- NULL
  i=1
  for (i in 1:fixed) {
    keep <- ans <- NULL
    keep <- selectpoints3d(ids["data"], value = FALSE, button = "right")[2]
    points3d(specimen[keep, 1], specimen[keep, 2], specimen[keep, 3], size = 10, color = "red", add = TRUE)
    selected[i, ] <- as.numeric(specimen[keep, ])
    fix <- c(fix, keep)
    rgl.bringtotop(stay = FALSE)
  }
  rgl.close()
  rownames(selected) <- seq(from = 1, to = nrow(selected))
  colnames(selected) <- c("xpts", "ypts", "zpts")
  return(selected)
}


AOS_manual <- function(mesh, plots=F, plots.fin=F) {
  # Manual pre-orientation of the model
  # Last Update: 2019/09/27
  # Dependencies: plot3d {rgl}; lines3d {rgl}; points3d {rgl}; rotate3d {rgl}; translate3d {rgl};
  #               CircleFitByKasa {conicfit}; rotV {cwhmisc}
  
  # Notes:  There was problem with 'Digit.Fixed.New' function. To solve this issue,
  #         the original Digit.Fixed function from 'geomorph' package was used
  #         However this creates problems because clicked point is sometimes identified on another surface

  # Arguments:
  #   mesh:       triangular mesh ('mesh3d')
  #   plots:      visualisation of the process ('logical')
  #   plots.fin:  visualisation of the result ('logical')
  # Value:
  #   mesh:       pre-oriented triangular mesh ('mesh3d')
  
  # examples:
  #   AOS_manual(mesh=mesh, plots=T, plots.fin=F)

  par3d(windowRect=c(960,30,1920,1040)); plot3d(mesh, meshColor="legacy", asp="iso")
  # pts <- Digit.Fixed.New(mesh,mesh,3,center=FALSE)$A            # problems with Digit.Fixed.New
  pts <- digit.fixed.lite(mesh,fixed = 3,ptsize = 2,center=F)                      # newly added
  if (plots==T) { par3d(windowRect=c(960,30,1920,1040)); plot3d(mesh, meshColor="legacy", asp="iso"); points3d(rbind(pts), col="red", size=20) }
  plane <- plane3p(A=pts[1,],B=pts[2,],C=pts[3,])
  plane.norm <- plane[1:3]/sqrt(sum(plane[1:3]^2))
  mesh.M.rot <- rotV(plane.norm, c(0,0,1))
  mesh.rot <- rotate3d(mesh,matrix=mesh.M.rot)
  pts.rot <- rotate3d(pts,matrix=mesh.M.rot)
  if (plots==T) { par3d(windowRect=c(960,30,1920,1040)); plot3d(mesh.rot, meshColor="legacy"); points3d(pts.rot, col="red", size=20) }
  circle <- CircleFitByKasa(pts.rot[,1:2])
  cx <- circle[1]; cy <- circle[2]; r <- circle[3]
  mesh.preal <- translate3d(mesh.rot,-cx,-cy,0)
  pts.preal <- translate3d(pts.rot,-cx,-cy,0)
  if (plots==T) {
    plot3d(mesh.preal, meshColor="legacy"); par3d(windowRect=c(960,30,1920,1040)) 
    points3d(pts.preal, col="red", size=10)
    lines3d(rbind(c(0,0,min(mesh.preal$vb[3,])),c(0,0,max(mesh.preal$vb[3,]))), lwd=2, col="red")
    rr <- max(sqrt(mesh.preal$vb[1,]^2+mesh.preal$vb[2,]^2))
    plot(mesh.preal$vb[1,],mesh.preal$vb[2,], col="lightgrey", asp=1, xlim=c(-rr,rr), ylim=c(-rr,rr),
         main="xy projection", xlab="x", ylab="y")
    points(pts.preal[,1],pts.preal[,2], col="red", pch=19)  
    points(0,0, col="red", pch=19)
    draw.circle2d(0,0, r=r, col.border="red")
  }
  mid.pt <- apply(pts.preal[c(1,3),],2,mean); mid.pt[3] <- 0
  if (plots==T) {
    points3d(mid.pt[1],mid.pt[2],mid.pt[3], col="blue",size=10)
    points(mid.pt[1],mid.pt[2], col="blue", pch=19)
  }
  mid.pt.rot <- rotV(mid.pt,c(1,0,0))
  mesh.preal <- rotate3d(mesh.preal, matrix=mid.pt.rot)
  if (plots==T) {
    plot3d(mesh.preal, meshColor="legacy"); par3d(windowRect=c(960,30,1920,1040)) 
    lines3d(rbind(c(0,0,min(mesh.preal$vb[3,])),c(0,0,max(mesh.preal$vb[3,]))), lwd=2, col="red")
    
    rr <- max(sqrt(mesh.preal$vb[1,]^2+mesh.preal$vb[2,]^2))
    plot(mesh.preal$vb[1,],mesh.preal$vb[2,], col="lightgrey", asp=1, xlim=c(-rr,rr), ylim=c(-rr,rr),
         main="xy projection", xlab="x", ylab="y")
    points(pts.preal[,1],pts.preal[,2], col="red", pch=19)
    points(0,0, col="red", pch=19)
    draw.circle2d(0,0, r=r, col.border="red")
  }
  if (plots.fin==T) {
    plot3d(mesh.preal, meshColor="legacy"); par3d(windowRect=c(960,30,1920,1040)) 
    lines3d(rbind(c(0,0,min(mesh.preal$vb[3,])),c(0,0,max(mesh.preal$vb[3,]))), lwd=2, col="red")
  }
  mesh <- mesh.preal
  return(mesh)
}





# 2.2. AUTOMATIC MODEL PRE-ORIENTATION
AOS_normals.all <- function(mesh=mesh, par=par, plots=F, infos=F) {
  # The function calculates
  # (i) the optimal position of the rotation axis for a given axis direction
  # (ii) the sum of squared distances between normals and the optimal axis
  # Note that the model is fixed and only the z-axis is rotated
  # Last Update: 2016/10/10
  # Dependencies: 'rotate3d' {rgl}; 'prcomp' {stats}; lmRob {robust}; plot3d {rgl}; plotNormals {Morpho}; lines3d {rgl}
  # More details: Halir 1997, Halir 1999
  # Notes: function 'AOS_normals.halir.f1' returns only the fifth element of the result (i.e. 'SSD')
  
  # Arguments:
  #   mesh:       triangular mesh ('mesh3d')
  #   par:        2 rotation parameters (theta, phi) of the z-axis (here used as the rotation axis) ('vector')
  #               par=c(theta,phi)
  #               'phi' for the rotation around the x-axis (in radians)
  #               'theta' for the rotation around the y-axis (in radians)
  #   plots:      visualisation of the process/results ('logical')
  #   infos:      print the result ('logical')
  # Value:
  #   result:     vector containing 5 elements ('vector')
  #               result=c(phi,theta,a,b,SSD)
  #               'phi' and 'theta' are initial rotation parameters
  #               'a' and 'b' are x- and y- translation parameters
  #               'SSD' sum of squared distances between normals and the optimal axis

  # examples:
  #   AOS_normals.all(mesh=mesh, par=c(0,0), plots=T, infos=T)
  
  phi <- par[1]
  theta <- par[2]
  axe <- c(0,0,1)
  axe <- rotate3d(axe,phi,1,0,0)
  axe <- rotate3d(axe,theta,0,1,0)
  n0 <- axe
  ni <- t(mesh$normals[-4,])
  Ai <- matrix(NA,nrow(ni),3)
  for (i in 1:nrow(ni)) {
    mm <- cbind(n0,ni[i,])
    Ai.temp <- sapply(1:3,function(x) det(mm[-x,])*(2*(x%%2)-1))
    Ai[i,] <- Ai.temp/sqrt(sum(Ai.temp^2))
  }
  Xi <- t(mesh$vb[-4,])
  qi <- rep(NA,nrow(Ai))
  for (i in 1:nrow(Ai)) {
    qi[i] <- crossprod(Xi[i,],Ai[i,])
  }
  pcs <- prcomp(Ai)
  X0 <- lmRob(qi~pcs$x[,1:2])$coefficients[-1] %*% t(pcs$rotation[,1:2]) + pcs$center
  X0 <- as.vector(X0)
  sD <- rep(NA,nrow(Xi))
  for (i in 1:nrow(Xi)) { sD[i] <- dist3Dline(X0,n0,Xi[i,],ni[i,])^2 }
  F1 <- sum(sD)
  X0.t <- X0
  X0.t <- rotate3d(X0.t,-theta,0,1,0)
  X0.t <- rotate3d(X0.t,-phi,1,0,0)
  ab <- X0.t[1:2]
  a <- ab[1]
  b <- ab[2]
  if (plots==T) {
    meshx <- rotate3d(mesh,-theta,0,1,0)
    meshx <- rotate3d(meshx,-phi,1,0,0)
    
    plot3d(meshx, meshColor="legacy", asp="iso"); par3d(windowRect=c(960,30,1920,1040)) 
    plotNormals(meshx, long = 150, col="lightgrey", lwd=0.1)
    lines3d(rbind(c(ab[1],ab[2],min(meshx$vb[3,])),c(ab[1],ab[2],max(meshx$vb[3,]))), col="red", lwd=2)
  }
  if (infos==T) { print(paste("phi:", round(phi,4), "; theta:", round(theta,4),
                              "; a:", round(a,4), "; b:", round(b,4), "; F1:", round(F1,5)),sep="") }
  result <-c(phi,theta,a,b,F1)
  return(result)
}

AOS_normals.f1 <- function(mesh, par, plots, infos) {
  # function used in optimisation
  # returns only the fifth element of the 'AOS_normals.all' result
  # see 'AOS_normals.halir.all' for more details
  return(AOS_normals.all(mesh=mesh, par=par, plots=plots, infos=infos)[5])
}

dist3Dline <- function(x0, n0, x1, n1) {
  # Calculation the distance between two lines (X0 and X1) in 3D. Each line is defined by the point and the normal.
  # Last Update: 2016/10/10

  # Arguments:
  #   x0:         3D coordinates of the point on the first line X0 ('vector', num, 3)
  #   n0:         3D normal coordinates of the first line X0  ('vector', num, 3)
  #   x1:         3D coordinates of the point on the second line X1  ('vector', num, 3)
  #   n1:         3D normal coordinates of the second line X1  ('vector', num, 3)
  # Value:
  #   dis:        distance between two lines in 3D  ('vector', num, 1)
  
  # examples:
  # dist3Dline (x0=c(0,0,0), n0=c(0,1,0), x1=c(1,1,1), n1=c(0,1,0))
  
  mm <- cbind(n0,n1)
  xy <- sapply(1:3,function(x)det(mm[-x,])*(2*(x%%2)-1))
  d <- as.vector(crossprod((x0-x1),xy)/sqrt(sum(xy^2)))
  if (is.nan(d)) { d=0 }
  dis <- abs(d)
  return(dis)
}



# 3. OPTIMISATION OF THE ROTATION AXIS POSITION
# 3.1. Horizontal-section model orientation using Pareto front

nb_planes <- function(mesh) {
  # Calculates the maximum number of horizontal planes which can be used in horizontal-section methods
  # Last Update: 2016/10/10
  # Dependencies: 'meshres' {Morpho}
  # More details: Mara and Sablatnig 2006
  
  # Arguments:
  #   mesh:       triangular mesh ('mesh3d')
  # Value:
  #   nb:         the maximum number of horizontal planes ('vector', num, 1)
  
  # examples:
  #   nb_planes(mesh)
  
  mesh.res <- meshres(mesh)
  lz <- max(mesh$vb[3,])-min(mesh$vb[3,])
  nb <- 0; di1 <- mesh.res*2; di <- 0
  while (di1-di >= mesh.res*2) {
    nb <- nb+1
    di <- nb*lz/nb
    di1 <- (nb+1)*lz/nb
    di1-di >= mesh.res*2
  }
  return(nb)
}

circular_section <- function(centre,xy) {
  # Calculates the circular sector
  # Last Update: 2016/10/10
  # Dependencies: none
  
  # Arguments:
  #   centre:     xy coordinates of the circle centre ('vector', num, 2)
  #   xy:         xy coordinates of the points on the circular section ('matrix', num, Nx2)
  # Value:
  #   cs:         arc-length of the circular section in radians ('vector', num, 1)
  
  # examples:
  #   circular_section(c(0,0), xy=cbind(1:10,2:11))
  
  pl.cen <- xy[,1:2]-matrix(centre,dim(xy)[1],2,byrow=T)
  pol <- cart2pol(pl.cen)
  pol <- pol[order(pol[,1]),]
  pol <- pol[,1]
  pol <- c(pol,pol[1])
  d <- numeric()
  for (i in 1:(length(pol)-1)) {
    if (pol[i]>pol[i+1]) { d[i]=2*pi-(pol[i]-pol[i+1]) }
    else {
      d[i]=pol[i+1]-pol[i]
    }
  }
  cs <- 2*pi-max(d)
  return(cs)
}

AOS_sections.all <- function(mesh, par, nb=0, plots=F, infos=F) {
  # The function can be used only for one (i.e. inner OR outer) side of the fragment. It calculates:
  # (i) F1
  # (ii) F2
  # Notes: (i) function can be used only for one (i.e. inner OR outer) side of the fragment
  #        (ii) the model (not axis) is rotated
  #        (iii) at least 5 points in one intersection of the model with the horizontal plane are needed to perform the calculation
  #        (iv) at least 20% of the arc must be preserved
  #        (v) Discussion 2016/10/12: Instead of the model rotation, it may be faster to rotate the axis
  #        (vi) a corresponds to x bar, b to y bar in equation 3 in the text!!!
  # Last Update: 2016/10/10
  # Dependencies: 'rotate3d' {rgl}; 'quantile' {stats}; 'meshPlaneIntersect' {Morhpo}; 'CircleFitByKasa' {conicfit};
  #               'plot3d' {rgl}; 'points3d' {rgl}; 'wt.var' {SDMTools}; wt.mean {SDMTools}
  # More details:
  # Notes: function 'AOS_sections.f12' returns only 'F1' and 'F2'
  
  # Arguments:
  #   mesh:       triangular mesh ('mesh3d')
  #   par:        2 rotation parameters ('phi', 'theta') of the model transformation ('vector', num, 2)
  #               par=c(phi,theta)
  #               'phi' for the rotation around the x-axis (in radians)
  #               'theta' for the rotation around the y-axis (in radians)
  #   nb:         number of horizontal planes ('vector', num, 1)
  #               '0' for the maximum possible number of horizontal sections
  #   plots:      visualisation of the process/results ('logical')
  #   infos:      print the result ('logical')
  # Value:
  #   result:     vector containing 6 elements ('vector', num, 6)
  #               result=c(phi,theta,a,b,F1,F2)
  #               'phi' and 'theta' are initial rotation parameters
  #               'a' and 'b' are x- and y- translation parameters
  #               'F1' and 'F2' are .....
  
  # examples:
  #   AOS_sections.all(mesh=mesh, par=c(0,0), nb=0, plots=T, infos=T)
  
  phi <- par[1]
  theta <- par[2]
  mesh <- rotate3d(mesh,phi,1,0,0)
  mesh <- rotate3d(mesh,theta,0,1,0)
  if (nb==0) { nb <- nb_planes(mesh) }
  else { nb <- nb  }
  lim <- quantile(mesh$vb[3,],c(0.05,0.95))
  cz <- seq(from=lim[1],to=lim[2], length.out=nb+2)
  cz <- cz[-c(1,length(cz))]
  OUTPUTS <- matrix(numeric(0), ncol=4, nrow=0)
  colnames(OUTPUTS) <- c("cx","cy","wi=section_angle","res=sum((di-r)^2)/n")
  for (k in 1:(length(cz))) {
    v1 <- c(0,0,cz[k])
    v2 <- c(0,1,cz[k])
    v3 <- c(1,0,cz[k])
    xy <- tryCatch(meshPlaneIntersect(mesh,v1,v2,v3), error = function(c) matrix(0,1,3))
    if (length(xy[,1])>5) { 
      reg.circle <- as.vector(CircleFitByKasa(xy[,1:2]))
      cx <- reg.circle[1]
      cy <- reg.circle[2]
      r <- reg.circle[3]
      wi <- circular_section(centre=c(cx,cy), xy=xy[,1:2])
      di <- sqrt((xy[,1]-cx)^2+(xy[,2]-cy)^2)
      n <- length(di)
      res <- sum((di-r)^2)/n
      outputs <- c(cx,cy,wi,res)
      if (outputs[3]>20*pi/180) {
        OUTPUTS <- rbind(OUTPUTS, outputs)
        if (plots==T) {
          if (k==1) {
            plot3d(mesh, meshColor="legacy", asp="iso"); par3d(windowRect=c(960,30,1920,1040)) 
            limi <- r*2
            plot(0, type="n", xlab="x", ylab="y", xlim=c(min(xy[,1])-limi,max(xy[,1])+limi), ylim=c(min(xy[,2])-limi,max(xy[,2])+limi), asp=1)
          }
          points3d(xy, col="blue", size=8)
          points3d(cx,cy,cz[k], col="red", size=8)
          points(xy[,1:2], col="blue", cex=0.9)
          points(cx,cy, col="red", pch=19)
          draw.circle2d(cx,cy, r=r, col.border="red")
        }
      }
    }
  }
  F1 <- wt.var(OUTPUTS[,1],OUTPUTS[,3])+wt.var(OUTPUTS[,2],OUTPUTS[,3])
  F2 <- wt.var(OUTPUTS[,4],OUTPUTS[,3])
  a <- wt.mean(OUTPUTS[,1],OUTPUTS[,3])
  b <- wt.mean(OUTPUTS[,2],OUTPUTS[,3])
  if (plots==T) { lines3d(rbind(c(a,b,min(mesh$vb[3,])),c(a,b,max(mesh$vb[3,]))), lwd=2, col="red") }
  if (infos==T) { print(paste("phi:", round(phi,4), "; theta:", round(theta,4),
                              "; a:", round(a,4), "b:", round(b,4), "; F1:", round(F1,5), "; F2:", round(F2,5)),sep="") }
  result <- c(phi,theta,a,b,F1,F2)
  return(result)
}

AOS_sections.f12 <- function(mesh, par, nb=0, plots, infos) {
  # function used in optimisation
  # returns only the fifth and sixth elements of the 'AOS_sections.all' result
  # see 'AOS_sections.all' for more details
  return(AOS_sections.all(mesh=mesh, par=par, nb=nb, plots=plots, infos=infos)[5:6])
}

AOS_sections_pair.all <- function(mesh.pair, par, nb=0, plots=F, infos=F) {
  # The function is used for both (i.e. inner AND outer) sides of the fragment.
  # Notes: (i) the function is used for both (i.e. inner AND outer) sides of the fragment
  #        (ii) the model (not axis) is rotated
  #        (iii) at least 5 points in one intersection of the model with the horizontal plane are needed to perform the calculation
  #        (iv) at least 20% of the arc must be preservaed
  #        (v) Discussion 2016/10/12: Instead of the model rotation, it may be faster to rotate the axis
  # Last Update: 2016/10/12
  # Dependencies: 'rotate3d' {rgl}; 'quantile' {stats}; 'meshPlaneIntersect' {Morhpo}; 'CircleFitByKasa' {conicfit};
  #               'plot3d' {rgl}; 'points3d' {rgl}; 'wt.var' {SDMTools}; wt.mean {SDMTools}
  # Notes: function 'AOS_sections_pair.f12' returns only 'F1' and 'F2'
  
  # Arguments:
  #   mesh.pair:  the 'list' containing 2 triangular meshes - 'mesh$S1' and 'mesh$S2' - both ('mesh3d') 
  #   par:        2 rotation parameters ('phi', 'theta') of the model transformation ('vector', num, 1)
  #               par=c(phi,theta)
  #               'phi' for the rotation around the x-axis (in radians)
  #               'theta' for the rotation around the y-axis (in radians)
  #   nb:         number of horizontal planes ('vector', num, 1)
  #               '0' for the maximum possible number of horizontal sections
  #   plots:      visualisation of the process/results ('logical')
  #   infos:      print the result ('logical')
  # Value:
  #   result:     vector containing 6 elements ('vector', num, 6)
  #               result=c(phi,theta,a,b,F1,F2)
  #               'phi' and 'theta' are initial rotation parameters
  #               'a' and 'b' are x- and y- translation parameters
  #               'F1' and 'F2' are outputs to be minimised
  #               
  
  # examples:
  #   AOS_sections_pair.all(mesh.pair=mesh.pair, par=c(0,0), nb=0, plots=T, infos=T)
  

  phi <- par[1]
  theta <- par[2]
  OUTPUTS <- matrix(numeric(0), ncol=5, nrow=0)
  colnames(OUTPUTS) <- c("cx","cy","wi=section_angle","res=sum((di-r)^2)/n", "E/I")
  for (i in 1:2) {
    mesh <- mesh.pair[[i]]
    mesh <- rotate3d(mesh,phi,1,0,0)
    mesh <- rotate3d(mesh,theta,0,1,0)
    if (nb==0) { nb <- nb_planes(mesh) }
    else { nb <- nb }
    lim <- quantile(mesh$vb[3,],c(0.05,0.95))
    cz <- seq(from=lim[1],to=lim[2], length.out=nb+2)
    cz <- cz[-c(1,length(cz))]
    for (k in 1:(length(cz))) {
      v1 <- c(0,0,cz[k])
      v2 <- c(0,1,cz[k])
      v3 <- c(1,0,cz[k])
      xy <- tryCatch(meshPlaneIntersect(mesh,v1,v2,v3), error = function(c) matrix(0,1,3))
      if (length(xy[,1])>5) { 
        reg.circle <- as.vector(CircleFitByKasa(xy[,1:2]))
        cx <- reg.circle[1]
        cy <- reg.circle[2]
        r <- reg.circle[3]
        wi <- circular_section(centre=c(cx,cy), xy=xy[,1:2])
        di <- sqrt((xy[,1]-cx)^2+(xy[,2]-cy)^2)
        n <- length(di)
        res <- sum((di-r)^2)/n
        outputs <- c(cx,cy,wi,res,i)
        if (outputs[3]>20*pi/180) {
          OUTPUTS <- rbind(OUTPUTS, outputs)
          if (plots==T) {
            if (i==1 & k==1) {
              plot3d(mesh, meshColor="legacy", asp="iso"); par3d(windowRect=c(960,30,1920,1040)) 
              limi <- r*2
              plot(0, type="n", xlab="x", ylab="y", xlim=c(min(xy[,1])-limi,max(xy[,1])+limi), ylim=c(min(xy[,2])-limi,max(xy[,2])+limi), asp=T)
            }
            if (i==2 & k==1) { plot3d(mesh, meshColor="legacy", add=T) }
            
            points3d(xy, col="blue", size=8)
            points3d(cx,cy,cz[k], col="red", size=8)
            
            points(xy[,1:2], col="blue", cex=0.9)
            points(cx,cy, col="red", pch=19)
            draw.circle2d(cx,cy, r=r, col.border="red")
          }
        }
      }
    }
  }
  F1 <- wt.var(OUTPUTS[,1],OUTPUTS[,3])+wt.var(OUTPUTS[,2],OUTPUTS[,3])
  F2 <- wt.var(OUTPUTS[,4],OUTPUTS[,3])
  a <- wt.mean(OUTPUTS[,1],OUTPUTS[,3])
  b <- wt.mean(OUTPUTS[,2],OUTPUTS[,3])
  if (plots==T) { lines3d(rbind(c(a,b,min(c(mesh.pair[[1]]$vb[3,],mesh.pair[[2]]$vb[3,]))),c(a,b,max(c(mesh.pair[[1]]$vb[3,],mesh.pair[[2]]$vb[3,])))), lwd=2, col="red") }
  if (infos==T) { print(paste("phi:", round(phi,4), "; theta:", round(theta,4),
                              "; a:", round(a,4), "b:", round(b,4), "; F1:", round(F1,5), "; F2:", round(F2,5)),sep="") }
  result <- c(phi,theta,a,b,F1,F2)
  return(result)
}

AOS_sections_pair.f12 <- function(mesh.pair, par, nb=0, plots=T, infos) {
  # function used in optimisation
  # returns only the fifth and sixth elements of the 'AOS_sections_pair.all' result
  # see 'AOS_sections_pair.all' for more details
  return(AOS_sections_pair.all(mesh.pair=mesh.pair, par=par, nb=nb, plots=plots, infos=infos)[5:6])
}

select_best_pareto <- function(best.paretos, nb, plots=F) {
  # Calculation of the best Pareto front, i.e. the mean values of phi and theta for 'nb' points
  # Last Update: 2017/01/29
  # Dependencies: none
  
  # Arguments:
  #   best.paretos:   result of nsga2 optimisation function ('list')
  #   nb:             number of optimal values ('vector', num, 1)
  #   plots:          visualisation of the process/results ('logical')
  # Value:
  #   best.par:       vector containing 2 elements  ('vector', num, 2)
  #                   best.par=c(phi,theta)
  #                  'phi' and 'theta' are initial rotation parameters
  
  res <- best.paretos$value
  res.norm <- cbind((res[,1]-min(res[,1])) / (max(res[,1])-min(res[,1])),(res[,2]-min(res[,2])) / (max(res[,2])-min(res[,2])))
  ord <- order(res.norm[,1]+res.norm[,2])
  if (plots==T) {
    plot(best.paretos, xlab="F1", ylab="F2", main="Objective space")
    points(best.paretos$value[ord[1:nb],1], best.paretos$value[ord[1:nb],2], col="blue", pch=19)
    plot(best.paretos$par, xlab="phi", ylab="theta", main="Parameter space", xlim=c(-pi/2,pi/2), ylim=c(-pi/2,pi/2), asp=1)
    points(best.paretos$par[ord[1:nb],1], best.paretos$par[ord[1:nb],2], col="blue", pch=19)
    abline(v=0, lty=2); abline(h=0, lty=2)
  }
  final <- best.paretos$par[ord[1:nb],]
  if (length(final) == 2) { best.par <- final }
  else { best.par <- apply(final,2,mean) }
  return(best.par)
}



# 3.2.1 Horizontal-section model orientation using circle properties - optimisation of 4 parameters 
AOS_circle4DDL.all <- function(mesh, par, nb=0, plots=F, infos=F) {
  # The function can be used only for one (i.e. inner OR outer) side of the fragment.
  # Notes: (i) can be used only for one (i.e. inner OR outer) side of the fragment
  #        (ii) the model (not axis) is rotated and translated
  #        (iii) function 'AOS_circle4DDL.f12' returns only 'F1'
  #        (iv) Discussion 2016/10/12: Instead of the model rotation, it may be faster to rotate the axis
  # Last Update: 2016/10/10
  # Dependencies: rotate3d {rgl}, translate3d {rgl}, quantile {stats}, plot3d {rgl}, lines3d {rgl}, meshPlaneIntersect {Morpho}

  # Arguments:
  #   mesh:       triangular mesh ('mesh3d')
  #   par:        2 rotation ('phi', 'theta') and 2 translation ('a','b') parameters of the model transformation ('vector', num, 4)
  #               par=c(phi,theta,a,b)
  #               'phi' for the rotation around the x-axis (in radians)
  #               'theta' for the rotation around the y-axis (in radians)
  #               'a' for the translation along the x-axis (in mm)
  #               'b' for the translation along the y-axis (in mm)
  #   nb:         number of horizontal planes ('vector', num, 1)
  #               '0' for the maximum possible number of horizontal sections 
  #   plots:      visualisation of the process/results ('logical')
  #   infos:      print the result ('logical')
  # Value:
  #   result:     vector containing 5 elements ('vector', num, 5)
  #               result=c(phi,theta,a,b,F1)
  #               'phi' and 'theta' are initial rotation parameters
  #               'a' and 'b' are initial x- and y- translation parameters
  #               'F1'
  
  # examples:
  #   AOS_circle4DDL.all(mesh=mesh, par=c(0,0,0,0), nb=0, plots=T, infos=T)
  
  phi <- par[1]
  theta <- par[2]
  a <- par[3]
  b <- par[4]
  mesh <- rotate3d(mesh,phi,1,0,0)
  mesh <- rotate3d(mesh,theta,0,1,0)
  mesh <- translate3d(mesh,a,b,0)
  didivar <- numeric(0)
  if (nb==0) { nb <- nb_planes(mesh) }
  else { nb <- nb }
  lim <- quantile(mesh$vb[3,],c(0.05,0.95))
  cz <- seq(from=lim[1],to=lim[2], length.out=nb+2)
  cz <- cz[-c(1,length(cz))]
  if (plots==T) { plot3d(mesh, meshColor="legacy", asp="iso"); par3d(windowRect=c(960,30,1920,1040)) ; lines3d(rbind(c(0,0,min(mesh$vb[3,])),c(0,0,max(mesh$vb[3,]))), lwd=2, col="red") }
  for (k in 1:(length(cz))) {
    v1 <- c(0,0,cz[k])
    v2 <- c(0,1,cz[k])
    v3 <- c(1,0,cz[k])
    xy <- tryCatch(meshPlaneIntersect(mesh,v1,v2,v3), error = function(c) matrix(0,1,3))
    if (length(xy[,1])>5) { 
      di <- sqrt((xy[,1])^2+(xy[,2])^2)     # di is the distance of points from the circle centre (i.e. c(0,0))
      didivar <- c(didivar,var(di^2))
      if (plots==T) { points3d(xy, col="blue", size=8) }
    }
  }
  F1 <- sum(didivar)
  if (infos==T) { print(paste("phi:", round(phi,4), "; theta:", round(theta,4),
                              "; a:", round(a,4), "b:", round(b,4), "; F1:", round(F1,5)),sep="") }
  result <- c(phi,theta,a,b,F1)
  return(result)
}

AOS_circle4DDL.f1 <- function(mesh, par, nb, plots, infos) {
  # function used in optimisation
  # returns only the fifth element of the 'AOS_circle4DDL.all' result
  # see 'AOS_circle4DDL.all' for more details
  return(AOS_circle4DDL.all(mesh=mesh, par=par, nb=nb, plots=plots, infos=infos)[5])
}

AOS_circle4DDL_pair.all <- function(mesh.pair, par, nb=0, plots=F, infos=F) {
  # The function is used for both (i.e. inner AND outer) sides of the fragment.
  # Notes: (i) function is used for both (i.e. inner AND outer) sides of the fragment
  #        (ii) the model (not axis) is rotated
  #        (iii) function 'AOS_circle4DDL_pair.f1' returns only 'F1'
  #        (iv) Discussion 2016/10/12: Instead of the model rotation, it may be faster to rotate the axis
  # Last Update: 2016/10/13
  # Dependencies: rotate3d {rgl}, translate3d {rgl}, quantile {stats}, plot3d {rgl}, lines3d {rgl}, meshPlaneIntersect {Morpho}

  # Arguments:
  #   mesh.pair:  the 'list' containing 2 triangular meshes - 'mesh$S1' and 'mesh$S2' - both ('mesh3d') 
  #   par:        2 rotation ('phi', 'theta') and 2 translation ('a','b') parameters of the model transformation ('vector', num, 4)
  #               par=c(theta,phi,a,b)
  #               'phi' for the rotation around the x-axis (in radians)
  #               'theta' for the rotation around the y-axis (in radians)
  #               'a' for the translation along the x-axis (in mm)
  #               'b' for the translation along the y-axis (in mm)
  #   nb:         number of horizontal planes ('vector', num, 1)
  #               '0' for the maximum possible number of horizontal sections 
  #   plots:      visualisation of the process/results ('logical')
  #   infos:      print the result ('logical')
  # Value:
  #   result:     vector containing 5 elements ('vector', num, 5)
  #               result=c(phi,theta,a,b,F1)
  #               'phi' and 'theta' are initial rotation parameters
  #               'a' and 'b' are initial x- and y- translation parameters
  #               'F1'
  
  # examples:
  #   AOS_circle4DDL_pair.all(mesh.pair=mesh.pair, par=c(0,0,0,0), nb=0, plots=T, infos=T)

  phi <- par[1]
  theta <- par[2]
  a <- par[3]
  b <- par[4]
  didivar <- numeric(0)
  for (i in 1:2) {
    mesh <- mesh.pair[[i]]
    mesh <- rotate3d(mesh,phi,1,0,0)
    mesh <- rotate3d(mesh,theta,0,1,0)
    mesh <- translate3d(mesh,a,b,0)
    if (nb==0) { nb <- nb_planes(mesh) }
    else { nb <- nb }
    lim <- quantile(mesh$vb[3,],c(0.05,0.95))
    cz <- seq(from=lim[1],to=lim[2], length.out=nb+2)
    cz <- cz[-c(1,length(cz))]
    if(plots==T) {
      if(i==1) { plot3d(mesh, asp="iso", col="red") ; par3d(windowRect=c(960,30,1920,1040)) }
      if(i==2) { plot3d(mesh, col="blue", add=T) }
      lines3d(rbind(c(0,0,min(mesh$vb[3,])),c(0,0,max(mesh$vb[3,]))), lwd=2, col="red")
    }
    for (k in 1:(length(cz))) {
      v1 <- c(0,0,cz[k])
      v2 <- c(0,1,cz[k])
      v3 <- c(1,0,cz[k])
      xy <- tryCatch(meshPlaneIntersect(mesh,v1,v2,v3), error = function(c) matrix(0,1,3))
      if (length(xy[,1])>5) { 
        di <- sqrt((xy[,1])^2+(xy[,2])^2)
        didivar <- c(didivar,var(di^2))
        if (plots==T) { points3d(xy, col="blue", size=8) }
      }
    }
  }
  F1 <- sum(didivar)
  if (infos==T) { print(paste("phi:", round(phi,4), "; theta:", round(theta,4),
                              "; a:", round(a,4), "b:", round(b,4), "; F1:", round(F1,5)),sep="") }
  result <- c(phi,theta,a,b,F1)
  return(result)
}

AOS_circle4DDL_pair.f1 <- function(mesh.pair, par, nb, plots, infos) {
  # function used in optimisation
  # returns only the fifth element of the 'AOS_circle4DDL_pair.all' result
  # see 'AOS_circle4DDL_pair.all' for more details
  return(AOS_circle4DDL_pair.all(mesh.pair=mesh.pair, par=par, nb=nb, plots=plots, infos=infos)[5])
}



# 3.2.2 Horizontal-section model orientation using circle properties - optimisation of 3 parameters 
AOS_circle3DDL.all <- function(mesh, par, nb=0, plots=F, infos=F) {
  # The function can be used only for one (i.e. inner OR outer) side of the fragment. It calculates
  # the variation of the squared distances between model points of each intersection and the rotation axis
  # coincided (i.e. fixed) with the z-axis
  # Notes: (i) can be used only for one (i.e. inner OR outer) side of the fragment
  #        (ii) the model (not axis) is rotated and translated
  #        (iii) function 'AOS_circle3DDL.f12' returns only 'F1'
  #        (iv) Discussion 2016/10/12: Instead of the model rotation, it may be faster to rotate the axis
  # Last Update: 2016/10/13
  # Dependencies: rotate3d {rgl}, translate3d {rgl}, quantile {stats}, plot3d {rgl}, lines3d {rgl}, meshPlaneIntersect {Morpho}

  # Arguments:
  #   mesh:       triangular mesh ('mesh3d')
  #   par:        1 rotation ('theta') and 2 translation ('a','b') parameters of the model ('vector', num, 3)
  #               par=c(theta,a,b)
  #               'theta' for the rotation around the y-axis (in radians)
  #               'a' for the translation along the x-axis (in mm)
  #               'b' for the translation along the y-axis (in mm)
  #   nb:         number of horizontal planes ('vector', num, 1)
  #               '0' for the maximum possible number of horizontal sections 
  #   plots:      visualisation of the process/results ('logical')
  #   infos:      print the result ('logical')
  # Value:
  #   result:     vector containing 4 elements ('vector', num, 4)
  #               result=c(theta,a,b,F1)
  #               'theta' is the initial rotation parameter
  #               'a' and 'b' are the initial x- and y- translation parameters
  #               'F1'
  
  # examples:
  #   AOS_circle3DDL.all(mesh=mesh, par=c(0,0,0), nb=0, plots=T, infos=T)
  
  theta <- par[1]
  a <- par[2]
  b <- par[3]
  mesh <- rotate3d(mesh,theta,0,1,0)
  mesh <- translate3d(mesh,a,b,0)
  didivar <- numeric(0)
  if (nb==0) { nb <- nb_planes(mesh) }
  else { nb <- nb }
  lim <- quantile(mesh$vb[3,],c(0.05,0.95))
  cz <- seq(from=lim[1],to=lim[2], length.out=nb+2)
  cz <- cz[-c(1,length(cz))]
  if (plots==T) { plot3d(mesh, meshColor="legacy", asp="iso"); par3d(windowRect=c(960,30,1920,1040)) ; lines3d(rbind(c(0,0,min(mesh$vb[3,])),c(0,0,max(mesh$vb[3,]))), lwd=2, col="red") }
  for (k in 1:(length(cz))) {
    v1 <- c(0,0,cz[k])
    v2 <- c(0,1,cz[k])
    v3 <- c(1,0,cz[k])
    xy <- tryCatch(meshPlaneIntersect(mesh,v1,v2,v3), error = function(c) matrix(0,1,3))
    if (length(xy[,1])>5) { 
      di <- sqrt((xy[,1])^2+(xy[,2])^2)
      didivar <- c(didivar,var(di^2))
      if (plots==T) { points3d(xy, col="blue", size=8) }
    }
  }
  F1 <- sum(didivar)
  if (infos==T) { print(paste("theta:", round(theta,4), "; a:", roun(a,4), "b:", round(b,4), "; F1:", round(F1,5)),sep="") }
  result <- c(theta,a,b,F1)
  return(result)
}

AOS_circle3DDL.f1 <- function(mesh, par, nb, plots, infos) {
  # function used in optimisation
  # returns only the fourth element of the 'AOS_circle3DDL.all' result
  # see 'AOS_circle3DDL.all' for more details
  return(AOS_circle3DDL.all(mesh=mesh, par=par, nb=nb, plots=plots, infos=infos)[4])
}

AOS_circle3DDL_pair.all <- function(mesh.pair, par, nb=0, plots=F, infos=F) {
  # The function is used for both (i.e. inner AND outer) sides of the fragment. It calculates
  # the variation of the squared distances between model points of each intersection and the rotation axis
  # coincided (i.e. fixed) with the z-axis
  # Notes: (i) function is used for both (i.e. inner AND outer) sides of the fragment
  #        (ii) the model (not axis) is rotated
  #        (iii) function 'AOS_circle3DDL_pair.f1' returns only the fourth element (i.e. 'F1')
  #        (iv) Instead of the model rotation, it may be faster to rotate the axis
  # Last Update: 2016/10/13
  # Dependencies: rotate3d {rgl}, translate3d {rgl}, quantile {stats}, plot3d {rgl}, lines3d {rgl}, meshPlaneIntersect {Morpho}

  # Arguments:
  #   mesh.pair:  the 'list' containing 2 triangular meshes - 'mesh$S1' and 'mesh$S2' - both ('mesh3d') 
  #   par:        1 rotation ('theta') and 2 translation ('a','b') parameters of the model transformation ('vector', num, 1)
  #               par=c(theta,a,b)
  #               'theta' for the rotation around the y-axis (in radians)
  #               'a' for the translation along the x-axis (in mm)
  #               'b' for the translation along the y-axis (in mm)
  #   nb:         number of horizontal planes; nb=0 for the maximum possible number of horizontal sections 
  #   plots:      visualisation of the process/results ('logical')
  #   infos:      print the result ('logical')
  # Value:
  #   result:     vector containing 4 elements ('vector', num, 4)
  #               result=c(theta,a,b,F1)
  #               'theta' is the initial rotation parameter
  #               'a' and 'b' are the initial x- and y- translation parameters
  #               'F1'
  
  # examples:
  #   AOS_circle3DDL_pair.all(mesh.pair=mesh.pair, par=c(0,0,0), nb=0, plots=T, infos=T)

  theta <- par[1]
  a <- par[2]
  b <- par[3]
  didivar <- numeric(0)
  for (i in 1:2) {
    mesh <- mesh.pair[[i]]
    mesh <- rotate3d(mesh,theta,0,1,0)
    mesh <- translate3d(mesh,a,b,0)
    if (nb==0){ nb <- nb_planes(mesh) }
    else { nb <- nb }
    lim <- quantile(mesh$vb[3,],c(0.05,0.95))
    cz <- seq(from=lim[1],to=lim[2], length.out=nb+2)
    cz <- cz[-c(1,length(cz))]
    if (plots==T) {
      if(i==1) { plot3d(mesh, asp="iso", col="red"); par3d(windowRect=c(960,30,1920,1040))  }
      if(i==2) { plot3d(mesh, col="blue", add=T) }
      lines3d(rbind(c(0,0,min(mesh$vb[3,])),c(0,0,max(mesh$vb[3,]))), lwd=2, col="red")
    }
    for (k in 1:(length(cz))) {
      v1 <- c(0,0,cz[k])
      v2 <- c(0,1,cz[k])
      v3 <- c(1,0,cz[k])
      xy <- tryCatch(meshPlaneIntersect(mesh,v1,v2,v3), error = function(c) matrix(0,1,3))
      if (length(xy[,1])>5) { 
        di <- sqrt((xy[,1])^2+(xy[,2])^2)
        didivar <- c(didivar,var(di^2))
        if (plots==T) { points3d(xy, col="blue", size=8) }
      }
    }
  }
  F1 <- sum(didivar)
  if (infos==T) { print(paste("theta:", round(theta,4),"; a:", round(a,4), "b:", round(b,4), "; F1:", round(F1,5)),sep="") }
  result <- c(theta,a,b,F1)
  return(result)
}

AOS_circle3DDL_pair.f1 <- function(mesh.pair, par, nb, plots, infos) {
  # function used in optimisation
  # returns only the fourth element of the 'AOS_circle3DDL_pair.all' result
  # see 'AOS_circle3DDL_pair.all' for more details
  return(AOS_circle3DDL_pair.all(mesh.pair=mesh.pair, par=par, nb=nb, plots=plots, infos=infos)[4])
}



# 3.3. Vertical-section model orientation using profile superimposition
rho.range <- function (mesh, plots=F) {
  # Calculates range values of the circular section of the model seen in xy-projection
  # Last Update: 2016/12/08
  # Dependencies: none
  # Notes: angle is calculated from the origin of the coordinate system (x=0, y=0)

  # Arguments:
  #   mesh:       triangular mesh ('mesh3d')
  #   plots:      visualisation of the process/results ('logical')
  # Value:
  #   result:     vector containing 4 elements ('vector', num, 4)
  #               result=c(theta,a,b,F1)
  #               'theta' is the initial rotation parameter
  #               'a' and 'b' are the initial x- and y- translation parameters
  #               'F1'
  
  # examples:
  #   AOS_circle3DDL_pair.all(mesh.pair=mesh.pair, par=c(0,0,0), nb=0, plots=T, infos=T)
  
  
  centre <- c(0,0)
  xy <- t(mesh$vb[1:2,])
  pl.cen <- xy[,1:2]-matrix(centre,dim(xy)[1],2,byrow=T)
  pol <- cart2pol(pl.cen)
  pol.temp <- cbind(pol[,1],pol)
  #pol <- pol[order(pol[,1]),]
  pol <- pol[,1]
  if ( sum(pol>pi/2)>0 && sum(pol<(-pi/2))>0 ) {
    pol[pol<0] <- 2*pi + pol[pol<0]
    pol.temp[,1] <- pol
  }
  mi <- min(pol)
  ma <- max(pol)
  if (plots==T) {
    p.mi <- xy[which(pol.temp[,1]==mi),]
    p.ma <- xy[which(pol.temp[,1]==ma),]
    
    plot(xy, main="xy-plane", xlab="x", ylab="y", asp=1, cex=0.3, col="grey",
         xlim=c(-max(xy[,1])*(-1),max(xy[,1])), ylim=c(-max(xy[,2]),max(xy[,2])))
    points(centre[1], centre[2], pch=1, col="red")
    points(centre[1], centre[2], pch=3, col="red")
    points(p.mi[1],p.mi[2], col="red", pch=19)
    points(p.ma[1],p.ma[2], col="red", pch=19)
    lines(rbind(centre,p.mi), col="red")
    lines(rbind(centre,p.ma), col="red")
  }
  info <- c(mi,ma)
  return(info)
}

AOS_profile.all <- function(mesh, par, nb=6, by.deg=1, ref="longest", plots=F, infos=F) {
  # The function can be used for the entire model. It calculates:
  # the sum of squared distances between the referential intersection (profile) points and their closest correspondences within the rz plane
  # Notes: (i) the model (not axis) is rotated and translated
  #        (ii) instead of the model rotation, it may be faster to rotate the axis
  #        (iii) function 'AOS_profile.f1' returns only 'F1'
  # Last Update: 2017/03/21
  # Dependencies: rotate3d {rgl}, translate3d {rgl}, plot3d {rgl}, lines3d {rgl}, meshPlaneIntersect {Morpho}; planes3d {rgl}

  # Arguments:
  #   mesh:       triangular mesh ('mesh3d')
  #   par:        2 rotation ('phi', 'theta') and 2 translation ('a','b') parameters of the model transformation ('vector', num, 4)
  #               par=c(phi,theta,a,b)
  #               'phi' for the rotation around the x-axis (in radians)
  #               'theta' for the rotation around the y-axis (in radians)
  #               'a' for the translation along the x-axis (in mm)
  #               'b' for the translation along the y-axis (in mm)
  #   nb:         number of vertical planes ('vector', num, 1)
  #               nb=0 if the vertical planes are defined by equally spaced angles (see 'by.deg')
  #   by.deg:     number of degrees based on which the model will be sectioned ('vector', num, 1)
  #   ref:        definition of the referential profile ('vector', char, 1)
  #               'longest' for the section with the highest range of z-coordinates
  #               'max.pts' for the section possessing the maximum nombre of points
  #               'middle' for the section which is situated "in the middle" of the fragment
  #   plots:      visualisation of the process/results ('logical')
  #   infos:      print the result ('logical')
  # Value:
  #   result:     vector containing 5 elements ('vector', num, 5)
  #               result=c(phi,theta,a,b,F1)
  #               'phi' and 'theta' are the initial rotation parameters
  #               'a' and 'b' are the initial x- and y- translation parameters
  #               'F1'
  
  # examples:
  #   AOS_profile.all(mesh=mesh, par=c(0,0,0,0), nb=60, plots=T, infos=T)
  
  phi <- par[1]
  theta <- par[2]
  a <- par[3]
  b <- par[4]
  mesh <- rotate3d(mesh,phi,1,0,0)
  mesh <- rotate3d(mesh,theta,0,1,0)
  mesh <- translate3d(mesh,a,b,0)
  mima <- rho.range(mesh,plots=F)
  if (nb==0) { rho <- seq(from=mima[1], to=mima[2], by=deg2rad(by.deg)) }
  else { rho <- seq(from=mima[1], to=mima[2], length=nb+2) }
  rho <- -rho
  rho <- rho[-c(1, length(rho))]
  v1 <- c(0,0,min(mesh$vb[3,]))
  v2 <- c(0,0,max(mesh$vb[3,]))
  v3 <- c(1,0,0)
  if (plots==T) { par3d(windowRect=c(960,30,1920,1040)); plot3d(mesh, meshColor="legacy", asp="iso"); lines3d(rbind(c(0,0,min(mesh$vb[3,])),c(0,0,max(mesh$vb[3,]))), lwd=2, col="red") }
  RZ <- matrix(numeric(0),0,2)
  ind <- numeric(0)
  pr.range <- numeric(0)
  for (k in 1:(length(rho))) {
    v3.rot <- rotate3d(c(1,0,0),rho[k],0,0,1)
    pla <- tryCatch(meshPlaneIntersect(mesh,v1,v2,v3.rot),error = function(c) 0)
    if (length(pla)!=1) {
      rz <- cbind(sqrt(pla[,1]^2+pla[,2]^2),pla[,3])
      RZ <- rbind(RZ,rz)
      ind <- c(ind,dim(rz)[1])
      pr.range <- c(pr.range,max(rz[,2])-min(rz[,2]))
      if (plots==T) {
        v3.rot.p <- rotate3d(v3.rot,pi/2,0,0,1)
        planes3d(v3.rot.p[1],v3.rot.p[2],v3.rot.p[3], alpha=0.1, col="red")
      }
    }
  }
  if (ref=="longest") { r.ind <- which.max(pr.range) }
  if (ref=="max.pts") { r.ind <- which.max(ind) }
  if (ref=="middle")  { r.ind <- round(nb/2) }
  ind <- c(0,cumsum(ind))
  TARGET <- RZ[(ind[r.ind]+1):(ind[r.ind+1]),]
  SOURCE <- RZ[-((ind[r.ind]+1):(ind[r.ind+1])),]
  if (plots==T) {
    plot(SOURCE, asp=1, col="blue", cex=0.5, xlab="r", ylab="z", main="Profil plane")
    points(TARGET, col="red", cex=0.5)
  }
  di = numeric(dim(SOURCE)[1])
  for (p in 1:dim(SOURCE)[1]) {
    SOURCE.M <- matrix(SOURCE[p,],dim(TARGET)[1],2,byrow=T)         
    di[p] <- min(sqrt(apply((SOURCE.M-TARGET)^2,1,sum)))
  }
  F1 = sum(di^2)
  if (infos==T) { print(paste("phi:", round(phi,4), "; theta:", round(theta,4),
                              "; a:", round(a,4), "b:", round(b,4), "; F1:", round(F1,5)),sep="") }
  result <- c(phi,theta,a,b,F1)
  return(result)
}

AOS_profile.f1 <- function(mesh, par, nb, by.deg, ref, plots, infos) {
  # function used in optimisation
  # returns only the fifth element of the 'AOS_profile.all' result
  # see 'AOS_profile.all' for more details
  return(AOS_profile.all(mesh=mesh, par=par, nb=nb, by.deg=by.deg, ref=ref, plots=plots, infos=infos)[5])
}



# 3.4. Vertical-section model orientation using profile curve fitting
AOS_polynomial.all <- function(mesh, par, po=10, plots=F, infos=F) {
  # The function can be used only for one (i.e. inner OR outer) side of the fragment. It calculates the sum of squared residuals of profile points from the polynomial
  # Notes: (i) the model (not axis) is rotated and translated
  #        (ii) function 'AOS_polynomial.f1' returns only 'F1'
  #        (iii) Discussion 2016/10/12: Instead of the model rotation, it may be faster to rotate the axis
  # Last Update: 2017/03/21
  # Dependencies: rotate3d {rgl}, translate3d {rgl}, lm {stats}, predict {stats}, plot3d {rgl},
  #               lines3d {rgl}, turn3d {rgl}
  # More details: Willis et al. 2000
  
  # Arguments:
  #   mesh:       triangular mesh ('mesh3d')
  #   par:        2 rotation ('phi', 'theta') and 2 translation ('a','b') parameters of the model transformation ('vector', num, 4)
  #               par=c(phi,theta,a,b)
  #               'phi' for the rotation around the x-axis (in radians)
  #               'theta' for the rotation around the y-axis (in radians)
  #               'a' for the translation along the x-axis (in mm)
  #               'b' for the translation along the y-axis (in mm)
  #   po:         polynomial order
  #   plots:      visualisation of the process/results ('logical')
  #   infos:      print the result ('logical')
  # Value:
  #   result:     vector containing 5 elements ('vector', num, 5)
  #               result=c(phi,theta,a,b,F1)
  #               'phi' and 'theta' are the initial rotation parameters
  #               'a' and 'b' are the initial x- and y- translation parameters
  #               'F1'
  
  # examples:
  #   AOS_polynomial.all(mesh=mesh, par=c(0,0,0,0), po=12, plots=T, infos=T)
  
  phi <- par[1]
  theta <- par[2]
  a <- par[3]
  b <- par[4]
  mesh <- rotate3d(mesh,phi,1,0,0)
  mesh <- rotate3d(mesh,theta,0,1,0)
  mesh <- translate3d(mesh,a,b,0)
  xyz <- t(mesh$vb[1:3,])
  xyz <- xyz[order(xyz[,3]),]
  x <- xyz[,1]
  y <- xyz[,2]
  z <- xyz[,3]
  r <- sqrt(x^2+y^2)
  rz <- cbind(r,z)
  RZ <- as.data.frame(rz)
  colnames(RZ) <- c("r","z")
  polymod <- lm(r ~ poly(z, po, raw=T), data = RZ)
  F1 <- sum(polymod$residuals^2)
  if (plots==T) {
    d <- seq(min(z), max(z), length=200)
    pred <- predict(polymod, data.frame(z=d))
    d.pred <- cbind(pred,d)
    plot(RZ, asp=1, xlim=c(0,max(pred)), xlab="r", ylab="z", cex=0.5, main="Profile plane")
    lines(d.pred, col="red", lwd=2)
    lines(rbind(c(0,min(z)),c(0,max(z))), col="red", lwd=2)
    profil <- d.pred
    profil <- profil[,2:1]
    reconst <- turn3d(profil, n=64, smooth=T)
    reconst <- rotate3d(reconst,0.5*pi,0,1,0)
    plot3d(reconst, type="wire", col="lightgrey", asp="iso"); par3d(windowRect=c(960,30,1920,1040)) 
    plot3d(mesh, meshColor="legacy", add=T)
    lines3d(rbind(c(0,0,min(mesh$vb[3,])),c(0,0,max(mesh$vb[3,]))), lwd=2, col="red")
  }
  if (infos==T) { print(paste("phi:", round(phi,4), "; theta:", round(theta,4),
                              "; a:", round(a,4), "b:", round(b,4), "; F1:", round(F1,5)),sep="") }
  result <- c(phi,theta,a,b,F1)
  return(result)
}

AOS_polynomial.f1 <- function(mesh, par, po, plots, infos) {
  # function used in optimisation
  # returns only the fifth element of the 'AOS_polynomial.all' result
  # see 'AOS_polynomial.all' for more details
  return(AOS_polynomial.all(mesh=mesh, par=par, po=po, plots=plots, infos=infos)[5])
}

AOS_polynomial_pair.all <- function(mesh.pair, par, po=10, plots=F, infos=F) {
  # The function is used for both (i.e. inner AND outer) sides of the fragment. It calculates
  # the sum of squared residuals of profile points from the polynomials
  # Notes: (i) the model (not axis) is rotated and translated
  #        (ii) function 'AOS_polynomial_pair.f1' returns only 'F1'
  #        (iii) Discussion 2016/10/12: Instead of the model rotation, it may be faster to rotate the axis
  # Last Update: 2016/10/13
  # Dependencies: rotate3d {rgl}, translate3d {rgl}, lm {stats}, predict {stats}, plot3d {rgl},
  #               lines3d {rgl}, turn3d {rgl}
  # More details: Willis et al. 2000
  
  # Arguments:
  #   mesh:       triangular mesh ('mesh3d')
  #   mesh.pair:  the 'list' containing 2 triangular meshes - 'mesh$S1' and 'mesh$S2' - both ('mesh3d') 
  #               par=c(phi,theta,a,b)
  #               'phi' for the rotation around the x-axis (in radians)
  #               'theta' for the rotation around the y-axis (in radians)
  #               'a' for the translation along the x-axis
  #               'b' for the translation along the y-axis
  #   po:         polynomial order
  #   plots:      visualisation of the process/results ('logical')
  #   infos:      print the result ('logical')
  # Value:
  #   result:     vector containing 5 elements ('vector')
  #               result=c(phi,theta,a,b,F1)
  #               'phi' and 'theta' are the initial rotation parameters
  #               'a' and 'b' are the initial x- and y- translation parameters
  #               'F1'
  
  # examples:
  #   AOS_polynomial_pair.all(mesh.pair=mesh.pair, par=c(0,0,0,0), po=12, plots=T, infos=T)
  
  phi <- par[1]
  theta <- par[2]
  a <- par[3]
  b <- par[4]
  sumdidi <- numeric(0)
  for (i in 1:2) {
    mesh <- mesh.pair[[i]]
    mesh <- rotate3d(mesh,phi,1,0,0)
    mesh <- rotate3d(mesh,theta,0,1,0)
    mesh <- translate3d(mesh,a,b,0)
    xyz <- t(mesh$vb[1:3,])
    xyz <- xyz[order(xyz[,3]),]
    x <- xyz[,1]
    y <- xyz[,2]
    z <- xyz[,3]
    r <- sqrt(x^2+y^2)
    rz <- cbind(r,z)
    RZ <- as.data.frame(rz)
    colnames(RZ) <- c("r","z")
    polymod <- lm(r ~ poly(z, po, raw=T), data = RZ)
    sumdidi <- c(sumdidi,sum(polymod$residuals^2))
    if (plots==T) {
      d <- seq(min(z), max(z), length=200)
      pred <- predict(polymod, data.frame(z=d))
      d.pred <- cbind(pred,d)
      if (i==1) {
        mesh.xyz <- rbind(t(mesh.pair[[1]]$vb[1:3,]),t(mesh.pair[[2]]$vb[1:3,]))
        mesh.rz <- cbind(sqrt(mesh.xyz[,1]^2+mesh.xyz[,2]^2),mesh.xyz[,3])
        plot(mesh.rz, asp=1, xlim=c(0,max(pred)), xlab="r", ylab="z", cex=0.5, main="Profile plane")
        lines(rbind(c(0,min(mesh.xyz[,3])),c(0,max(mesh.xyz[,3]))), col="red", lwd=2)
        plot3d(mesh.pair[[1]], col="red", asp="iso"); par3d(windowRect=c(960,30,1920,1040)) 
      }
      if (i==2) { plot3d(mesh.pair[[2]], col="blue", add=T) }
      lines(d.pred, col="red", lwd=2)
      profil <- d.pred
      profil <- profil[,2:1]
      reconst <- turn3d(profil, n=64, smooth=T)
      reconst <- rotate3d(reconst,0.5*pi,0,1,0)
      plot3d(reconst, type="wire", col="lightgrey", add=T)
    }
  }
  F1 <- sum(sumdidi)
  if (infos==T) { print(paste("phi:", round(phi,4), "; theta:", round(theta,4),
                              "; a:", round(a,4), "b:", round(b,4), "; F1:", round(F1,5)),sep="") }
  result <- c(phi,theta,a,b,F1)
  return(result)
}

AOS_polynomial_pair.f1 <- function(mesh.pair, par, po, plots, infos) {
  # function used in optimisation
  # returns only the fifth element of the 'AOS_profile.all' result
  # see 'AOS_polynomial_pair.all' for more details
  return(AOS_polynomial_pair.all(mesh.pair=mesh.pair, par=par, po=po, plots=plots, infos=infos)[5])
}


# 3.5. Rim/base-tangent model orientation
AOS_rimbase.all <- function(mesh, par, tresh=1, nb=20, by.deg=1, part="rim", three.p=F, plots=F, infos=F) {
  # The function can be used for the entire model containing RIM or BASE. It calculates
  # the sum of squared distances between z-coordinates of the highest/lowest points (defined by threshold)
  # and the maximum/minimum z-coordinate
  # Notes: (i) the model (not axis) is rotated and translated
  #        (ii) Discussion 2016/10/12: Instead of the model rotation, it may be faster to rotate the axis
  # Last Update: 2017/03/21
  # Dependencies: rotate3d {rgl}, translate3d {rgl}, plot3d {rgl}, lines3d {rgl}, meshPlaneIntersect {Morpho}; planes3d {rgl}
  # Notes: function 'AOS_rimbase.f1' returns only 'F1'
  
  # Arguments:
  #   mesh:       triangular mesh ('mesh3d')
  #   par:        2 rotation ('phi', 'theta') and 2 translation ('a','b') parameters of the model transformation ('vector', num, 4)
  #               par=c(phi,theta,a,b)
  #               'phi' for the rotation around the x-axis (in radians)
  #               'theta' for the rotation around the y-axis (in radians)
  #               'a' for the translation along the x-axis (in mm)
  #               'b' for the translation along the y-axis (in mm)
  #   tresh:      the distance (in mm) defining the "points belonging to the rim" (in mm) ('vector', num, 1)
  #   nb:         number of vertical planes ('vector', num, 1)
  #               nb=0 if the vertical planes are defined by equally spaced angles (see 'by.deg')
  #   by.deg:     number of degrees based on which the model will be sectioned ('vector', num, 1)
  #   part:        definition of the part on which the axis will be searched ('vector', char, 1)
  #               'rim' for rim
  #               'base' for base
  #   # three.p   if only 3 points defining the horizontal plane will be used in calculation ('logical')
  #   plots:      visualisation of the process/results ('logical')
  #   infos:      print the result ('logical')
  # Value:
  #   result:     vector containing 5 elements ('vector', num, 5)
  #               result=c(phi,theta,a,b,F1)
  #               'phi' and 'theta' are the initial rotation parameters
  #               'a' and 'b' are the initial x- and y- translation parameters
  #               'F1'
  
  # examples:
  #   AOS_rimbase.all(mesh=mesh, par=c(0,0,0,0), tresh=0.5, nb=0, by.deg=2, part="rim", three.p=T, plots=T, infos=T)

  phi <- par[1]
  theta <- par[2]
  a <- par[3]
  b <- par[4]
  mesh <- rotate3d(mesh,phi,1,0,0)
  mesh <- rotate3d(mesh,theta,0,1,0)
  mesh <- translate3d(mesh,a,b,0)
  mima <- rho.range(mesh,plots=F)
  if (nb==0) { rho <- seq(from=mima[1], to=mima[2], by=deg2rad(by.deg)) }
  else {
    rho <- seq(from=mima[1], to=mima[2], length=nb+2)
  }
  rho <- -rho
  rho <- rho[-c(1, length(rho))]
  v1 <- c(0,0,min(mesh$vb[3,]))
  v2 <- c(0,0,max(mesh$vb[3,]))
  v3 <- c(1,0,0)
  if (plots==T) { plot3d(mesh, meshColor="legacy", asp="iso"); par3d(windowRect=c(960,30,1920,1040)) ; lines3d(rbind(c(0,0,min(mesh$vb[3,])),c(0,0,max(mesh$vb[3,]))), lwd=2, col="red") }
  pts <- matrix(numeric(0),0,3)
  for (k in 1:(length(rho))) {
    v3.rot <- rotate3d(c(1,0,0),rho[k],0,0,1)
    pla <- tryCatch(meshPlaneIntersect(mesh,v1,v2,v3.rot), error = function(c) matrix(0,1,3))
    if (dim(pla)[1]>1) {
      if (part=="rim") { p <- pla[which.max(pla[,3]),] }
      if (part=="base") { p <- pla[which.min(pla[,3]),] }
      pts <- rbind(pts,p)
      if (plots==T) {
        v3.rot.p <- rotate3d(v3.rot,pi/2,0,0,1)
        planes3d(v3.rot.p[1],v3.rot.p[2],v3.rot.p[3], alpha=0.1, col="red")
        points3d(p[1],p[2],p[3], col="blue",size=5)
      }
    }
  }
  if (part=="rim") {
    while ( sum(abs(pts[,3]-max(pts[,3]))<tresh) < 3 ) { tresh <- tresh + 0.1 }
    pts.tresh <- pts[abs(pts[,3]-max(pts[,3]))<tresh,]
    if (three.p==T) { pts.tresh <- pts.tresh[order(pts.tresh[,3])[1:3],] }
    F1 <- mean( (pts.tresh[,3]-max(pts.tresh[,3]))^2 )
  }
  if (part=="base") {
    while ( sum(abs(pts[,3]-min(pts[,3]))<tresh) < 3 ) { tresh <- tresh + 0.1 }
    pts.tresh <- pts[abs(pts[,3]-min(pts[,3]))<tresh,]
    if (three.p==T) { pts.tresh <- pts.tresh[order(pts.tresh[,3],decreasing=T)[1:3],] }
    F1 <- mean( (pts.tresh[,3]-min(pts.tresh[,3]))^2 )
  }
  if (plots==T) { points3d(pts.tresh, col="red",size=10) }
  if (infos==T) { print(paste("phi:", round(phi,4), "; theta:", round(theta,4),
                              "; a:", round(a,4), "b:", round(b,4), "; F1:", round(F1,5)),sep="") }
  result <- c(phi,theta,a,b,F1)
  return(result)
}

AOS_rimbase.f1 <- function(mesh, par, tresh, nb, by.deg, part, three.p, plots, infos) {
  # function used in optimisation
  # returns only the fifth element of the 'AOS_rimbasel.all' result
  # see 'AOS_rimbase.all' for more details
  return(AOS_rimbase.all(mesh=mesh, par=par, tresh=tresh, nb=nb, by.deg=by.deg, part=part, three.p=three.p, plots=plots, infos=infos)[5])
}



# 4.2. Profile extraction
profile2d <- function(mesh, method="envelop", longest.by.deg=0.2, brush.by.deg=0.6, brush.tresh=5, envelop.alpha=0.1,
                      envelop.alpha.add.by=0.1, plots=F) {
  # Calculation of the profile
  # Last Update: 2017/01/22
  # Dependencies: identify3d {rgl}, plot3d {rgl}, points3d {rgl}, lines3d {rgl}, rgl.close {rgl},
  #               rotate3d {rgl}, meshPlaneIntersect {Morpho}
  #               ashape {alphahull}, graph.edgelist {igraph}, degree {igraph}, clusters {igraph},
  #               E {igraph}, get.shortest.paths {igraph}, V {igraph}

  # Arguments:
  #   mesh:                 triangular mesh ('mesh3d')
  #   method:               method used for profile extraction ('vector', char, 1)
  #                         'arbitrary' for the profile selected by the operator
  #                         'middle' for the profile located in the middle of the fragment's horizontal arc-section
  #                         'longest' for the longest preserved profile
  #                         'brush' for the cleaning profile. Cleaning is by left mouse click, to exit
  #                         click inside the blue circle. Note that when this method is selected, alpha shape
  #                         is calculated so 'envelop.alpha' and 'envelp.alpha.add.by' parameters can be adjusted
  #                         'envelop' for the whole envelope profile
  #   longest.by.deg:       used for 'longest' method, indicates by how many degrees the fragment is vertically sectioned (in degrees) ('vector', num, 1)
  #   envelop.alpha:        the initial alpha value for ashape calculation; used in 'envelop' method ('vector', num, 1)
  #   envelop.alpha.add.by: used in 'envelop' method. If the problem with the edge connection is found during
  #                         the calculation of the ashape with the initial alpha set to 'envelop.alpha',
  #                         the 'envelop.alpha' is iteratively augmented by 'envelop.alpha.add.by' until the
  #                         ashape is properly calculated ('vector', num, 1)
  #   plots:                visualisation of the process/results ('logical')
  # Value:
  #   profil:               list containing 3 elements  ('list')
  #                         '$profil' contains all profile coordinates ('matrix', num, Nx2)
  #                         '$extern' contains outer profile coordinates ('matrix', num, Nx2)
  #                         '$intern' contains inner profile coordinates ('matrix', num, Nx2)
  
  # examples:
  # profile2d(mesh, method="arbitrary", plots=T)
  # profile2d(mesh, method="middle", plots=T)
  # profile2d(mesh, method="longest", plots=T)
  # profile2d(mesh, method="envelop", plots=T)
  # profile2d(mesh, method="brush", plots=T)
  
  
  if (method=="arbitrary") {
    plotOrientedMesh(mesh, method="3d")
    idx <- identify3d(t(mesh$vb[1:3,]), n=1)
    rgl.close()
    pt <- mesh$vb[1:3,idx]
  }
  if (method=="middle") {
    mid <- sum((rho.range(mesh, plots=F)))/2
    pt <- c(pol2cart(c(mid,1)),0)
  }
  if (method=="longest") {
    mima <- rho.range(mesh,plots=F)
    rho <- seq(from=mima[1], to=mima[2], by=deg2rad(longest.by.deg))
    v1 <- c(0,0,min(mesh$vb[3,]))
    v2 <- c(0,0,max(mesh$vb[3,]))
    v3 <- c(1,0,0)
    RES <- matrix(numeric(0),0,2); colnames(RES) <- c("abs(max(z)-min(z)),", "rho")
    for (k in 1:(length(rho))) {
      v3.rot <- rotate3d(c(1,0,0),rho[k],0,0,1)
      pla <- tryCatch(meshPlaneIntersect(mesh,v1,v2,v3.rot),
                      error = function(c) NULL)
      if (length(pla)!=0) {
        RES <- rbind(RES,c(abs(max(pla[,3])-min(pla[,3])), rho[k]))
      }
    }
    longest <- RES[which.max(RES[,1]),2]
    pt <- c(pol2cart(c(-longest,1)),0)
  }
  if (method=="arbitrary" || method=="middle" || method=="longest") {
    sect <- meshPlaneIntersect(mesh,c(0,0,0),c(0,0,1),pt)
    rz <- cbind(sqrt(sect[,1]^2+sect[,2]^2),sect[,3])
    if (plots==T) { plotOrientedMesh(mesh, method="3d"); points3d(sect, col="red", size=10) }
  }
  if (method=="envelop") {
    r <- sqrt(mesh$vb[1,]^2+mesh$vb[2,]^2)
    z <- mesh$vb[3,]
  }
  if (method=="brush") {
    r <- sqrt(mesh$vb[1,]^2+mesh$vb[2,]^2)
    z <- mesh$vb[3,]
    pts <- unique(matrix(c(r,z),length(r),2))
    plots=F
    longest.by.deg=4
    mima <- rho.range(mesh,plots=F)
    rho <- seq(from=mima[1], to=mima[2], by=deg2rad(brush.by.deg))
    v1 <- c(0,0,min(mesh$vb[3,]))
    v2 <- c(0,0,max(mesh$vb[3,]))
    v3 <- c(1,0,0)
    PLA <- matrix(numeric(0),0,2); colnames(PLA) <- c("r", "z")
    for (k in 1:(length(rho))) {
      v3.rot <- rotate3d(c(1,0,0),rho[k],0,0,1)
      pla <- tryCatch(meshPlaneIntersect(mesh,v1,v2,v3.rot),
                      error = function(c) NULL)
      if (length(pla)!=0) {
        pla.rz <- cbind(sqrt(pla[,1]^2+pla[,2]^2),pla[,3])
        PLA <- rbind(PLA,pla.rz)
        if (plots==T) {
          v3.rot.p <- rotate3d(v3.rot,pi/2,0,0,1)
          planes3d(v3.rot.p[1],v3.rot.p[2],v3.rot.p[3], alpha=0.1, col="red")
        }
      }
    }
    # errase
    p.exit <- c(min(pts[,1])-10,max(pts[,2]))
    p.back <- c(min(pts[,1])-40,max(pts[,2]))
    p <- c(10000,100000)
    pl.lim <- c(0,max(pts[,1]))
    pts.old <- pts
    while (sum((p-p.exit)^2)>brush.tresh) {
      plot(pts, type="n", asp=1, xlim=pl.lim, xlab=c("r"), ylab=c("z"))
      points(PLA[,1], PLA[,2], col="red", pch=1, cex=0.5)
      points(pts, pch=19, cex=0.5)
      points(p.exit[1],p.exit[2], col="blue", cex=brush.tresh)
      points(p.exit[1],p.exit[2], col="blue", cex=brush.tresh, pch=3)
      text(p.exit[1],p.exit[2], labels="done", col="blue")
      points(p.back[1],p.back[2], col="red", cex=brush.tresh)
      points(p.back[1],p.back[2], col="red", cex=brush.tresh, pch=3)
      text(p.back[1],p.back[2], labels="back", col="red")
      
      p <- locator(n=1)
      p <- c(p$x,p$y)
      points(p[1],p[2], col="blue", cex=brush.tresh)
      
      if (sum((p-p.back)^2)<brush.tresh) {
        pts <- pts.old
      }
      
      # try
      else {
        pts.old <- pts
        dis <- sqrt(apply((pts-matrix(p,dim(pts)[1],dim(pts)[2],byrow=T))^2,1,sum))
        err <- which(dis<brush.tresh)
        if (length(err)!=0) { pts <- pts[-err,] }
      }
    }
    
    r <- pts[,1]
    z <- pts[,2]
    print("Brushing finished")
  }
  if (method=="envelop" || method=="brush") {
    pts <- unique(matrix(c(r,z),length(r),2))
    err <- 1
    while (err>0) {
      a <- ashape(pts, alpha = envelop.alpha)
      a.g <- graph.edgelist(cbind(as.character(a$edges[, "ind1"]), as.character(a$edges[,"ind2"])), directed = FALSE)
      err <- 0
      if (!is.connected(a.g)) { err <- err+1 }
      if (any(degree(a.g) != 2)) { err <- err+1 }
      if (clusters(a.g)$no > 1) { err <- err+1 }
      if (err>0) { envelop.alpha <- envelop.alpha + envelop.alpha.add.by }
    }
    cutg <- a.g - E(a.g)[1]
    ends <- names(which(degree(cutg) == 1))
    path <- get.shortest.paths(cutg, ends[1], ends[2])[[1]]  
    pathX <- as.numeric(V(a.g)[path[[1]]]$name)
    pathX <- c(pathX, pathX[1])
    rz <- a$x[pathX, ]
  }
  rz <- sort.coords.profile(rz) # if (rz[1,1]==rz[2,1]) { rz <- rz[-1,] } # if the coordinates of the first and the second point are the same
  if (rz[1,1]<rz[2,1]) {
    extern <- rz[1:which.min(rz[,2]),]
    intern <- rz[(which.min(rz[,2])+1):dim(rz)[1],]
  }
  if (rz[1,1]>rz[2,1]) {
    extern <- rz[c(1,dim(rz)[1]:which.min(rz[,2])),]
    intern <- rz[(which.min(rz[,2])-1):2,]
  }
  profil <- rbind(extern,intern)
  if (intern[which.min(abs((intern[,2]-mean(profil[,2])))),1] > extern[which.min(abs((extern[,2]-mean(profil[,2])))),1]) {
    intern.temp <- intern; extern.temp <- extern
    intern <- extern.temp; extern <- intern.temp
  }
  if (intern[1,2]<intern[dim(intern)[1],2]) { intern <- intern[(dim(intern)[1]):1,] }
  if (extern[1,2]<extern[dim(extern)[1],2]) { extern <- extern[(dim(extern)[1]):1,] }
  if (extern[1,2]<intern[1,2]) { extern <- rbind(intern[1,],extern) }
  if (plots==T) {
    plot(profil, type="n", xlim=c(min(profil[,1]),0), asp=1, axes=T, xlab="", ylab="")
    lines(rbind(c(0,0),c(0,min(profil[,2]))), col="red", lwd=2)
    polygon(profil, col="grey", border="darkgrey")
    lines(extern, col="blue")
    points(rbind(extern[1,]), col="blue", pch=3)
    lines(intern, col="red")
    points(rbind(intern[1,]), col="red")
  }
  profil <- list(profil=profil, extern=extern, intern=intern)
  return (profil)
}

profile2d.repair <- function(profil, pt, tolerance=10, plots=F) {
  # Repair the profile
  # Last Update: 2018/01/08

  # Arguments:
  #   profil:     profile object ('list')
  #   pt:         xy-coordinates of the point ('vector', num, 2)
  #   tolerance:  maximum distance used for point identification (in mm) ('vector', num, 1)
  #   plots:      visualisation of the process/results ('logical')
  # Value:
  #   profil:     list containing 3 elements  ('list')
  #               '$profil' contains all profile coordinates ('matrix', num, Nx2)
  #               '$extern' contains outer profile coordinates ('matrix', num, Nx2)
  #               '$intern' contains inner profile coordinates ('matrix', num, Nx2)
  
  # examples:
  # profile2d.repair(profil, c(-48,-2), tolerance=5, plots=T)

  pro <- profil$profil
  extern <- profil$extern
  intern <- profil$intern
  distP2O <- function(pt,M){ sqrt(apply((M-matrix(pt,dim(M)[1],2,byrow=T))^2,1,sum))  }
  if (min(distP2O(pt,pro))>tolerance) { warning("No point identified within tolerated distance"); return(profil);  }
  p <- which.min(distP2O(pt,pro))
  if (sum(sqrt(apply((extern-matrix(pro[p,],dim(extern)[1],2,byrow=T))^2,1,sum))==0)>0) {
    pp <- which(sqrt(apply((extern-matrix(pro[p,],dim(extern)[1],2,byrow=T))^2,1,sum))==0)
    if (abs(max(extern[,2])-extern[pp,2]) < abs(min(extern[,2])-extern[pp,2])) { extern <- extern[pp:(dim(extern)[1]),] }
    else { extern <- extern[1:pp,] }
  }
  if (sum(sqrt(apply((intern-matrix(pro[p,],dim(intern)[1],2,byrow=T))^2,1,sum))==0)>0) {
    pp <- which(sqrt(apply((intern-matrix(pro[p,],dim(intern)[1],2,byrow=T))^2,1,sum))==0)
    if (abs(max(intern[,2])-intern[pp,2]) < abs(min(intern[,2])-intern[pp,2])) { intern <- intern[pp:(dim(intern)[1]),] }
    else { intern <- intern[1:pp,] }
  }
  else {
    va <- numeric(4)
    va[1] <- sqrt(sum((extern[1,]-pt)^2))
    va[2] <- sqrt(sum((extern[dim(extern)[1],]-pt)^2))
    va[3] <- sqrt(sum((intern[1,]-pt)^2))
    va[4] <- sqrt(sum((intern[dim(intern)[1],]-pt)^2))
    va <- which.min(va)
    if (va==1) {
      p2 <- which(distP2O(extern[1,],pro)==0)-1
      extern <- rbind(pro[p:p2,], extern)
    }
    if (va==2) {
      p2 <- which(distP2O(extern[dim(extern)[1],],pro)==0)+1
      extern <- rbind(extern, pro[p2:p,])
    }
    if (va==3) {
      p2 <- which(distP2O(intern[1,],pro)==0)+1
      intern <- rbind(pro[p:p2,], intern)
    }
    if (va==4) {
      p2 <- which(distP2O(intern[dim(intern)[1],],pro)==0)-1
      intern <- rbind(intern, pro[p2:p,])
    }
  }
  if (plots==T) {
    plot(profil$profil, type="n", xlim=c(min(profil$profil[,1]),0), ylim=c(min(profil$profil[,2])-5,max(profil$profil[,2])+5), asp=1, axes=T, xlab="", ylab="")
    lines(rbind(c(0,0),c(0,min(pro[,2]))), col="red", lwd=2)
    polygon(pro, col="grey", border="darkgrey")
    lines(extern, col="blue", lwd=2)
    lines(intern, col="red", lwd=2)
  }
  profil <- list(profil=pro, extern=extern, intern=intern)
  return(profil)
}

cumchord <- function(M) {
  cumsum(sqrt(apply((M-rbind(M[1,], M[-(dim(M)[1]),]))^2,1,sum)))
}

profile2d.rec <- function (profil, ptt, rec.len=6, rec.off=2, plots=F) {
  # Calculate "reconstructed part of vessel" by using polynomials
  # Last Update: 2018/01/22
  
  # Arguments:
  #   profil:     profile object ('list')
  #   pt:         xy-coordinates of the point ('vector', num, 2)
  #   plots:      visualisation of the process/results ('logical')
  # Value:
  #   lin:        coordinates of the line ('matrix')
  
  if (plots==T) {
    plot(profil$profil, type="n", xlim=c(min(profil$profil[,1]),0), ylim=c(min(profil$profil[,2])-5,max(profil$profil[,2])+5), asp=1, axes=T, xlab="", ylab="")
    lines(rbind(c(0,0),c(0,min(profil$profil[,2]))), col="red", lwd=2)
    polygon(profil$profil, col="grey", border="darkgrey")
    lines(profil$extern, col="blue")
    lines(profil$intern, col="red")
  }
  
  pro <- profil$profil
  extern <- profil$extern
  intern <- profil$intern
  extern1 <- Line(extern)
  extern1 <- spsample(extern1, 10, type="regular", offset=c(0,0))@coords
  extern <- rbind(extern1,extern[dim(extern)[1],])
  intern1 <- Line(intern)
  intern1 <- spsample(intern1, 10, type="regular", offset=c(0,0))@coords
  intern <- rbind(intern1,intern[dim(intern)[1],])
  pt <- ptt[1,]
  
  va <- numeric(4)
  va[1] <- sqrt(sum((extern[1,]-pt)^2))
  va[2] <- sqrt(sum((extern[dim(extern)[1],]-pt)^2))
  va[3] <- sqrt(sum((intern[1,]-pt)^2))
  va[4] <- sqrt(sum((intern[dim(intern)[1],]-pt)^2))
  va <- which.min(va)
  
  pt <- ptt[2,]
  if (va==1) { rec <- rbind(pt,extern) }
  if (va==2) { rec <- rbind(extern,pt) }
  if (va==3) { rec <- rbind(pt,intern) }
  if (va==4) { rec <- rbind(intern,pt) }
  M <- rec
  row.names(M) <- NULL
  M1 <- as.data.frame(M)
  z <- cumchord(M1)
  fo1 <- splinefun(z,M1[,1],method="natural")
  fo2 <- splinefun(z,M1[,2], method="natural")
  x <- spline(z,M1[,1],method="natural", n=100)$y
  y <- spline(z,M1[,2],method="natural", n=100)$y
  rec <- cbind(x,y)
  if (va==1) {
    ind <- which.min((sqrt(apply((rec-matrix(extern[1,],dim(rec)[1],2,byrow=T))^2,1,sum))))
    rec <- rec[ind:1,]
  }
  if (va==2) {
    ind <- which.min((sqrt(apply((rec-matrix(extern[dim(extern)[1],],dim(rec)[1],2,byrow=T))^2,1,sum))))
    rec <- rec[ind:dim(rec)[1],]
  }
  if (va==3) {
    ind <- which.min((sqrt(apply((rec-matrix(intern[1,],dim(rec)[1],2,byrow=T))^2,1,sum))))
    rec <- rec[ind:1,]
  }
  if (va==4) {
    ind <- which.min((sqrt(apply((rec-matrix(intern[dim(intern)[1],],dim(rec)[1],2,byrow=T))^2,1,sum))))
    rec <- rec[ind:dim(rec)[1],]
  }
  ind <- which((sqrt(apply((rec-matrix(rec[1,],dim(rec)[1],2,byrow=T))^2,1,sum)))<(rec.len+rec.off))
  rec <- rec[ind,]
  ind <- which((sqrt(apply((rec-matrix(rec[1,],dim(rec)[1],2,byrow=T))^2,1,sum)))>rec.off)
  rec <- rec[ind,]
  lines(rec, col="grey", lty=2, lwd=1)
  
  rec <- Line(rec)
  rec <- spsample(rec, 40, type="regular", offset=c(0,1))@coords

  return(rec)
}

sort.coords.profile<- function(M, plots=F) {
  # Ordonation of the unorganised points in the (profile) outline
  # Last Update: 2016/10/30
  # Dependencies: none
  # Notes:  (i) points should be as close as possible
  #         (ii) function cannot handle incurved shapes
  
  # Arguments:
  #   M:        xy-coordinates of the unorganised points of the (profile) outline ('matrix', num, Nx2)
  #   plots:    visualisation of the process/results ('logical')
  # Value:
  #   res:      xy-coordinates of the organised points in the outline  ('matrix', num, Nx2)
  
  # examples:
  # sort.coords(M, plots=T)
  
  if (plots==T) { plot(M, asp=1) }
  ind <- which.max(M[,2])
  val <- (matrix(M[ind,],1,2,byrow=T))
  M <- M[-ind,]
  res <- val
  if (plots==T) { points(val, col="red", pch=19) }
  for (i in 1:dim(M)[1]) {
    if (length(M)==2) { M <- matrix(M,1,2,byrow=T) }
    VAL <- matrix(val,dim(M)[1],2,byrow=T)
    ind <- which.min(sqrt(apply((M-VAL)^2,1,sum)))
    val <- M[ind,]
    M <- M[-ind,]
    res <- rbind(res,val)
    if (plots==T) { points(val[1],val[2], col="red") }
  }
  if (plots==T) { plot(res, type="n", asp=1); lines(res, col="red") }
  return(res)
}

sort.coords.profile_old <- function(M, plots=F) {
  # Ordonation of the unorganised points in the (profile) outline
  # Last Update: 2016/10/30
  # Dependencies: none
  # Notes:  (i) points should be as close as possible
  #         (ii) function cannot handle incurved shapes
  
  # Arguments:
  #   M:        xy-coordinates of the unorganised points of the (profile) outline ('matrix', num, Nx2)
  #   plots:    visualisation of the process/results ('logical')
  # Value:
  #   res:      xy-coordinates of the organised points in the outline  ('matrix', num, Nx2)
  
  # examples:
  # sort.coords(M, plots=T)

  if (plots==T) { plot(M, asp=1) }
  ind <- which.max(M[,2])
  val <- (matrix(M[ind,],1,2,byrow=T))
  M <- M[-ind,]
  res <- val
  if (plots==T) { points(val, col="red", pch=19) }
  mea <- numeric(0)
  for (i in 1:dim(M)[1]) {
    if (length(M)==2) { M <- matrix(M,1,2,byrow=T) }
    VAL <- matrix(val,dim(M)[1],2,byrow=T)
    dis <- sqrt(apply((M-VAL)^2,1,sum))
    ind <- which.min(dis)
    mea <- c(mea,dis[ind])
    if (dis[ind]>mean(mea)*10) {
      print(i)
      M <- M[-ind,]
    }
    if (dis[ind]<mean(mea)*10) {
      val <- M[ind,]
      M <- M[-ind,]
      res <- rbind(res,val)
      if (plots==T) { points(val[1],val[2], col="red") }
    }
  }
  if (plots==T) { plot(res, type="n", asp=1); lines(res, col="red") }
  return(res)
}
switch.profile <- function(profil, plots=F) {
  # Switch profile side
  # Last Update: 2016/10/30

  # Arguments:
  #   profil:   list containing 3 elements ('list')
  #             'profil$profil', 'profil$extern' and 'profil$intern'
  #   plots:    visualisation of the process/results ('logical')
  # Value:
  #   profil:   list containing 3 elements  ('list')
  #             '$profil' contains all profile coordinates ('matrix', num, Nx2)
  #             '$extern' contains outer profile coordinates ('matrix', num, Nx2)
  #             '$intern' contains inner profile coordinates ('matrix', num, Nx2)
  
  # examples:
  # plot(profil$profil, asp=1)
  # plot(switch.profile(profil)$profil, asp=1)
  
  profil$profil[,1] <- -profil$profil[,1]
  profil$extern[,1] <- -profil$extern[,1]
  profil$intern[,1] <- -profil$intern[,1]
  profil <- list(profil=profil$profil, extern=profil$extern, intern=profil$intern)
  return (profil)
}




# Illustrations 

plot3d.bl <- function() {
  # Blank 3D plot used to plot objects
  # Dependencies: par3d {rgl}; view3d {rgl}; clear3d {rgl}; light3d {rgl}
  
  par3d(windowRect=c(960-3000,30-3000,1920-3000,1040-3000))
  view3d(theta=-10, phi=-60, zoom=0.66)
  # view3d(zoom=0.66)
  clear3d(type="light")
  light3d(ambient="white", diffuse="white", specular="black")
}

plot.blank <- function() {
  # Blank plot
  plot(0, type="n", axes=F, xlab="", ylab="")
}

plotOrientedMesh <- function(mesh, method="both", col="black", asp=1, cex=0.5, pch=19, add=F, alpha=1) {
  # The function that draws the oriented model in: (i) plot (rz-plane) and (ii) rgl (3D)
  # Last Update: 2016/10/28
  # Dependencies: plot3d {rgl}, lines3d {rgl}
  
  # Arguments:
  #   mesh:       triangular oriented mesh ('mesh3d')
  #   method:     method to draw ('vector', char, 1)
  #               'both' for RZ-plane and RGL
  #               'rz' for rz-plane
  #               '3d' for 3d in rgl
  #   asp, cex, pch, ...  graphical (see 'plot' or 'plot3d' for more info)
  #   add:        whatever the plot should be added ('logical')
  #   alpha:      model transparency ('vector', num, 1)
  # Value:
  #   graphic
  
  # examples:
  #   plotOrientedMesh(mesh, method="both", add=F)
  
  if (method=="3d" || method=="both") {
    plot3d(mesh, asp="iso", col=col, alpha=alpha, add=add)
    lines3d(rbind(c(0,0,min(mesh$vb[3,])),c(0,0,max(mesh$vb[3,]))), lwd=2, col="red")  
  }
  if (method=="rz" || method=="both") {
    r <- sqrt(mesh$vb[1,]^2+mesh$vb[2,]^2)
    if (add==F) { plot(cbind(r,mesh$vb[3,]), type="n", asp=asp, main="Profile plane", xlab="r", ylab="z", xlim=c(0,max(r))) }
    points(cbind(r,mesh$vb[3,]), col=col, cex=cex, pch=pch)
    lines( rbind(c(0,min(mesh$vb[3,])),c(0,max(mesh$vb[3,]))), col="red", lwd=2)
  }
}

reconst.2d <- function(pottery, add=F, ae.scale, CSPI, LIN, REC, main="",
                       profil.fill.col="darkgrey", lines.col="black", lines.lwd=1, 
                       add.rim=F, add.base=F,
                       add.lines=F, l.lines.col="black", l.lines.lwd=1,
                       add.diam=F,
                       add.scale=F, scale.fill.col="black",
                       add.cspi=F, cspi.fill.col="black", cspi.bord.col="black",
                       add.rec=F) {
  
  # 2D archaeological drawing
  # Last Update: 2016/10/30

  # Arguments:
  #   pottery:            pottery object ('list') - see 'profil2pottery' function for more details
  #   add:                add this graphic to layout ('logical')
  #   ae.scale:           'scale' object
  #   CSPI:               'SCPI' object
  #   LIN:                array with the lines coordinates ('array', num, 2x2xN)
  #   REC:                array with the reconstruction lines coordinates ('array', num, 40x2xN)
  #   main:               plot title ('vector', char, 1)
  #   profil.fill.col:    color of profile ('vector', char, 1)
  #   lines.col:          color of profile lines ('vector', char, 1)
  #   lines.lwd:          thickness of profile lines ('vector', num, 1)
  #   add.rim:            add rim line ('logical')
  #   add.base:           add base line ('logical')
  #   add.diam:           add diameter on the top of the drawing ('logical')
  #   add.lines:          add.lines on the reconstruction part ('logical')
  #   l.lines.col:        color of lines ('vector', char, 1)
  #   l.lines.lwd:        thickness of lines ('vector', num, 1)
  #   add.scale           add.scale ('logical')
  #   scale.fill.col:     scale filling color ('vector', char, 1)
  #   add.cspi:           add circle sector preservation indicator (CSPI) ('logical')
  #   cspi.fill.col:      CSPI filling color ('vector', char, 1)
  #   cspi.bord.col:      CSPI border color ('vector', char, 1)
  #   add.rec:            add reconstruction lines ('logical')
  # Value:
  #   graphic
  
  if (add==F) { plot(0, type="n", xlim=pottery$x.lim, ylim=pottery$y.lim, main=main, asp=1, axes=F, xlab="", ylab="") }
  if (add.lines==T) { if(dim(LIN)[3]>0) { for (i in 1:dim(LIN)[3]) { lines(LIN[,,i], col=l.lines.col, lwd=l.lines.lwd) } } }
  polygon(pottery$profil, col=profil.fill.col, border=lines.col, lwd=lines.lwd)
  lines(pottery$axi, col=lines.col, lwd=lines.lwd)
  lines(pottery$right, col=lines.col, lwd=lines.lwd)
  if (add.rim==T) { lines(pottery$rim, col=lines.col, lwd=lines.lwd) }
  if (add.base==T) { lines(pottery$base, col=lines.col, lwd=lines.lwd) }
  if (add.diam==T) { text(pottery$rupt[1],pottery$rupt[2], label=paste("d.", round((pottery$rupt[1]*2)/10,1), "cm"), pos=3, cex=1) }
  if (add.scale==T) {
    for (i in 1:dim(ae.scale$coords)[3]) {
      if (i %% 2 != 0) { polygon(ae.scale$coords[,,i], col="black") }
      else { polygon(ae.scale$coords[,,i]) }
    }
    text(ae.scale$spt[1],ae.scale$spt[2], label="0", pos=3)
    text(ae.scale$ept[1],ae.scale$ept[2], label=paste(ae.scale$value/10, "cm"), pos=3)
  }
  if (add.cspi==T) {
    polygon(CSPI$c.coords, border = cspi.bord.col, lwd=1)
    polygon(CSPI$cs.coords, col=cspi.fill.col, border=cspi.bord.col, lwd=1)
    text(CSPI$cs.val.pos[1], CSPI$cs.val.pos[2], labels=paste(round(CSPI$cs.val.prc),"%"), cex=0.7)
  }
  if (add.rec==T) { if (dim(REC)[3]>1) { for(i in 1:dim(REC)[3]) { lines(REC[,,i], lty=5, col=lines.col) } } }
}

lines.2d <- function (pottery, click.coo, x.off=5) {
  # Add a line to the 2D drawing
  # Last Update: 2016/11/22

  # Arguments:
  #   pottery:    pottery object ('list') - see 'profil2pottery' function for more details
  #   click.coo:  coordinates of the clicked point ('2d vector')
  # Value:
  #   lin:        coordinates of the line ('matrix')
  
  if (click.coo[1]>0) {
    pt <- pottery$right[which.min(abs(pottery$right[,2]-click.coo[2])),]
    lin <- rbind(c(0,pt[2]),pt)
  }
  if (click.coo[1]<0) {
    pt1 <- pottery$extern[which.min(abs(pottery$extern[,2]-click.coo[2])),]
    pt2 <- pottery$intern[which.min(abs(pottery$intern[,2]-click.coo[2])),]
    if (abs(pt1[1]-click.coo[1]) > abs(pt2[1]-click.coo[1])) {
      pt <- pt2
      pt[1] <- pt[1]+x.off
      lin <- rbind(pt,c(0,pt[2]))
    }
    else {
      pt <- pt1
      ext <- unique(pottery$extern)
      pts <- matrix(numeric(0),0,2)
      A1 <- pt
      A2 <- c(pt[1]+3, pt[2])
      for (j in 1:(dim(ext)[1]-1)) {
        B1 <- ext[j,]
        B2 <- ext[j+1,]
        x1 <- A1[1]; y1 <- A1[2]
        x2 <- A2[1]; y2 <- A2[2]
        x3 <- B1[1]; y3 <- B1[2]
        x4 <- B2[1]; y4 <- B2[2]
        if (  (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)==0  ) { }
        if (  (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)!=0  ) {
          Px <- ( (x1*y2-y1*x2)*(x3-x4) - (x1-x2)*(x3*y4-y3*x4) ) /
            ( (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4) )
          Py <- ( (x1*y2-y1*x2)*(y3-y4) - (y1-y2)*(x3*y4-y3*x4) ) /
            ( (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4) )
          P <- c(Px,Py)
          B1c <- B1-P
          B2c <- B2-P
          xmin <- min(c(B1c[1],B2c[1])); xmax <- max(c(B1c[1],B2c[1]))
          ymin <- min(c(B1c[2],B2c[2])); ymax <- max(c(B1c[2],B2c[2]))
          if (xmin<=0 & xmax>=0 & ymin<=0 & ymax>=0) {    #} & i!=j) {
            pts <- rbind(pts,P)
          }
        }
      }
      if (dim(pts)[1]<=3) {
        pt <- pt2
        pt[1] <- pt[1]+x.off
        lin <- rbind(pt,c(0,pt[2]))
      }
      else {
        lin <- pts[2:1,]
        lin[1,1] <- lin[1,1]+x.off
      }
    }
  }
  return(lin)
}




# Objects used for illustrations

profil2pottery <- function(profil, type="2d", x.off=5, off.side=10) {
  # Transform 'profil' object to 'pottery' object
  # The goal is to have one single object containing all the information needed for plotting the pottery
  # and not to calculate information within the 'reconstruct2d' and 'reconstruct2d.rgl' plot functions
  # Last Update: 2016/11/30
  # Notes: Profile needs to be oriented on left side. If is not # profil <- switch.profile(profil)
  
  # Arguments:
  #   profil:     list containing 3 elements  ('list')
  #               'profil$profil' contains all profile coordinates ('matrix', num, Nx2)
  #               'profil$extern' contains outer profile coordinates ('matrix', num, Nx2)
  #               'profil$intern' contains inner profile coordinates ('matrix', num, Nx2)
  #   type:       'NULL' for transformation to 2D/3D pottery object ('vector', char, 1)
  #               '2d' for transformation only to 2D pottery object
  #   x.off:      offset of the rim and base lines from the profile (in mm) ('vector', num, 1)
  #   off.side:   offset around the drawing (in mm) ('vector', num, 1)
  #   po.extern:  polynomial order of extern outline ('vector', num, 1)
  #   po.intern:  polynomial order of intern outline ('vector', num, 1)
  #   rec.len:    length of the reconstructed dotted part ('vector', num, 1)
  # Value:
  #   pottery:    list containing 12 elements  ('list')
  #               '$profil' contains all profile coordinates ('matrix', num, Nx2)
  #               '$extern' contains outer profile coordinates ('matrix', num, Nx2)
  #               '$intern' contains inner profile coordinates ('matrix', num, Nx2)
  #               '$axi' contains coordinates of the rotation axis ('matrix', num, Nx2)
  #               '$right' contains outer profile coordinates plotted on the right side ('matrix', num, Nx2)
  #               '$left' contains outer profile coordinates plotted on the left (==$extern) ('matrix', num, Nx2)
  #               '$rim' contains rim coordinates ('matrix', num, Nx2)
  #               '$base' contains base coordinates ('matrix', num, Nx2)
  #               '$x.lim' contains x limit of the plot ('vector', num, 2)
  #               '$y.lim' contains y limit of the plot ('vector', num, 2)
  #               '$rupt' contains xy coordinates of ru point ('vector', num, 2)
  #               '$rpt' contains xy coordinates of y value of rim ('vector', num, 2)
  #               '$bpt' contains xy coordinates of y value of base ('vector', num, 2)
  #               '$mask' contains coordinates of the rectangle mask ('matrix', num, Nx2)

  pro <- profil$profil
  out <- profil$extern
  int <- profil$intern

  left <- out[1:which.min(out[,2]),]
  right <- cbind(-out[,1],out[,2])
  rownames(right) <- NULL
  axi <- rbind(c(0,min(right[,2])),c(0,max(right[,2])))
  nb <- dim(left)[1]
  lim <- c(max(right[,2]),min(right[,2]))
  z <- seq(from=lim[1],to=lim[2], length.out=nb)
  right2 <- matrix(numeric(0),0,2)
  for (i in 1:(length(z)-1)) {
    sel <- right[right[,2]<z[i] & right[,2]>z[i+1],]
    if (length(sel)==2) { p <- sel }
    else { p <- sel[which.max(sel[,1]),] }
    right2 <- rbind(right2,p)
  }
  right <- rbind(right[which.max(right[,2]),],right2,right[which.min(right[,2]),])
  rim <- rbind(right[which.max(right[,2]),],right[which.max(right[,2]),]); rim[1,1] <- -rim[1,1]+x.off
  base <- rbind(right[which.min(right[,2]),],right[which.min(right[,2]),]); base[1,1] <- -base[1,1]
  if (base[1,1]+x.off < 0) { base[1,1] <- base[1,1]+x.off }
  x.lim <- c(min(pro[,1])-off.side, -min(pro[,1])+off.side)
  y.lim <- c(min(pro[,2])-off.side, max(pro[,2])+off.side)
  mask <- rbind(c(x.lim[1],y.lim[1]), c(0,y.lim[1]), c(0,y.lim[2]), c(x.lim[1],y.lim[2]))
  rupt <- rim[2,]
  rpt <- rim[1,2]
  bpt <- base[1,2]
  profil3d <- cbind(pro[,1],0,pro[,2])
  profil3d <- rbind(profil3d,profil3d[1,])
  extern3d <- cbind(out[,1],0,out[,2])
  intern3d <- cbind(int[,1],0,int[,2])
  axi3d <- cbind(axi[,1],0,axi[,2])
  left3d <- cbind(left[,1],0,left[,2])
  right3d <- cbind(right[,1],0,right[,2])
  rim3d <- cbind(rim[,1],0,rim[,2])
  base3d <- cbind(base[,1],0,base[,2])
  mask3d <- cbind(mask[,1],0,mask[,2])
  rupt3d <- c(rupt[1],0,rupt[2])
  
  if (type=="all") {
    pottery <- list(profil=profil$profil, extern=profil$extern, intern=profil$intern,
                    axi=axi,
                    left=left, right=right, rim=rim, base=base,
                    x.lim=x.lim, y.lim=y.lim,
                    mask=mask,
                    rupt=rupt,
                    rpt=rpt,
                    bpt=bpt,
                    
                    profil3d=profil3d, extern3d=extern3d, intern3d=intern3d,
                    axi3d=axi3d,
                    left3d=left3d, right3d=right3d, rim3d=rim3d, base3d=base3d,
                    mask3d=mask3d,
                    rupt3d=rupt3d,
                    rpt3d=rpt,
                    bpt3d=bpt
    )
  }
  if (type=="2d") {
    pottery <- list(profil=profil$profil, extern=profil$extern, intern=profil$intern,
                    axi=axi,
                    left=left, right=right, rim=rim, base=base,
                    x.lim=x.lim, y.lim=y.lim,
                    mask=mask,
                    rupt=rupt,
                    rpt=rpt,
                    bpt=bpt
    )
  }
  return(pottery)
}

archeo.scale <- function(scale.val=50, scale.seg=5, scale.height=10/4, scale.y=0, plots=F) {
  # Create 'scale' object for pottery drawing
  # Last Update: 2016/11/13
  
  # Arguments:
  #   scale.val:        maximum value shown on the scale. Corresponds to scale length (in mm) ('vector', num, 1)
  #   scale.height:     ratio between the length and height of one scale segment ('vector', num, 1)
  #                     '10/10' for square (i.e. 1/1 or 10 mm)
  #                     '10/2' for rectangle (i.e. 1/2 or 5 mm)
  #                     '10/3' for longer rectangle (i.e. 1/3 or 3.33 mm)
  #                     '10/4' for longer rectangle (i.e. 1/4 or 2.5 mm)
  #   scale.y:          position of the scale in the drawing (in mm) ('vector', num, 1)
  #                     calculated as 'pottery$bpt'-'x' mm
  #   plots:            visualisation of the process/results ('logical')
  # Value:
  #   archeo.scale:     list containing 3 elements  ('list')
  #                     '$scale.coords' array whose matrices contain coordinates of scale rectangles ('array', num, 4x2xN)
  #                     '$spt' contains xy coordinates of the position of the 0 value ('vector', num, 2)
  #                     '$ept' contains xy coordinates of the position of the scale.val ('vector', num, 2)
  
  # example: archeo.scale()
  
  d <- seq(from=-scale.val/2, to=scale.val/2, length=scale.seg+1)
  A <- array(numeric(0), dim=c(4,2,scale.seg))
  A3d <- array(numeric(0), dim=c(4,3,scale.seg))
  for (i in 1:scale.seg) {
    A[,,i] <- rbind(  c(d[i],scale.y),
                      c(d[i+1],scale.y),
                      c(d[i+1],scale.y-scale.height),
                      c(d[i],scale.y-scale.height)  )
    A3d[,,i] <- rbind(c(d[i],0,scale.y),
                      c(d[i+1],0,scale.y),
                      c(d[i+1],0,scale.y-scale.height),
                      c(d[i],0,scale.y-scale.height)  )
  }
  spt <- c(d[1],scale.y)
  ept <- c(d[scale.seg+1],scale.y)
  value <- scale.val/10
  if (plots==T) {
    for (i in 1:dim(A)[3]) {
      if (i %% 2 != 0) { polygon(A[,,i], col="black") }
      else { polygon(A[,,i]) }
    }
    text(spt[1],spt[2], label="0", pos=3)
    text(ept[1],ept[2], label=paste(value, "cm"), pos=3)
  }
  spt3d <- c(spt[1],0,spt[2])
  ept3d <- c(ept[1],0,ept[2])
  archeo.scale <- list(value=scale.val, coords=A, spt=spt, ept=ept,
                       coords3d=A3d, spt3d=spt3d, ept3d=ept3d)
  return(archeo.scale)
}

archeo.colo.scale <- function(x, colo, scale.length=100, scale.height=10/4, scale.y=0, method="segments", nb.seg=5, nb.by=0.5, plots=F) {
  # Create 'archeo.colo.scale' object for pottery drawing
  # Last Update: 2016/11/28
  # More details: created replicating 'meshDist' scale of Morpho package
  # Notes: not yet implemented for rgl.
  
  # Arguments:
  #   x:                values used in scales ('vector', num, N)
  #   scale.length:     scale length (in mm) ('vector', num, 1)
  #   scale.height:     ratio between the length and height of one scale segment ('vector', num, 1)
  #                     '10/10' for square (i.e. 1/1 or 10 mm)
  #                     '10/2' for rectangle (i.e. 1/2 or 5 mm)
  #                     '10/3' for longer rectangle (i.e. 1/3 or 3.33 mm)
  #                     '10/4' for longer rectangle (i.e. 1/4 or 2.5 mm)
  #   scale.y:          position of the scale in the drawing (in mm) ('vector', num, 1)
  #   method:           method by which numbers are drawn ('vector', char, 1)
  #                     'segments' scale will be segmented into 'nb.seg' values
  #                     'by' scale is segmented by 'nb.by' intervals
  #   nb.seg:           number of segments used for plot (see 'method') ('vector', num, 1)
  #   nb.by:            number by which the values are segmented (see 'method') ('vector', num, 1)
  #   plots:            visualisation of the process/results ('logical')
  # Value:
  #   archeo.colo.scale:     list containing 5 elements  ('list')
  #                         '$coords' array whose matrices contain coordinates of scale rectangles ('array', num, 4x2xN)
  #                         '$colo' vector containing scale colors ('vector', char, N)
  #                         '$values' matrix containing xy coordinates and values of the scale (in mm, rounded to 2 decimal) ('matrix', num, Nx2)
  #                         '$values.lin' array containing xy coordinates of the line ('array', num, 2x2xN)
  #                         '$lin' matrix containing xy coordinates of the line ('matrix', num, Nx2)
  
  # examples:
  # archeo.colo.scale(1:1000, colo=1:8)
  # archeo.colo.scale(1:1000, colo=gray.colors(10))
  
  scale.seg <- length(colo)
  d <- seq(from=-scale.length/2, to=scale.length/2, length=scale.seg+1)
  A <- array(numeric(0), dim=c(4,2,scale.seg))
  for (i in 1:scale.seg) {
    A[,,i] <- rbind(  c(d[i],scale.y),
                      c(d[i+1],scale.y),
                      c(d[i+1],scale.y-scale.height),
                      c(d[i],scale.y-scale.height)  )
  }
  if (method=="segments") {
    values <- cbind(d,scale.y-scale.height,round(x))
    values <- values[seq(1,length(x), length.out = nb.seg),]
  }
  if (method=="by") {
    nb <- seq(min(x),max(x),by=nb.by)
    uni <- (max(d)-min(d))/max(x)
    values <- cbind(min(d)+nb*uni,scale.y-scale.height,nb)
  }
  B <- array(numeric(0), dim=c(2,2,dim(values)[1]))
  for (i in 1:dim(values)[1]) {
    B[,,i] <- rbind(c(values[i,1],scale.y-scale.height),
                    c(values[i,1],scale.y-scale.height-1))
  }
  lin <- rbind(c(d[1],scale.y-scale.height),c(d[length(d)],scale.y-scale.height))
  if (plots==T) {
    plot(0, type="n", xlim=range(d), asp=1, axes=F, xlab="", ylab="")
    for (i in 1:dim(A)[3]) {
      polygon(A[,,i], col=colo[i], border=colo[i])
    }
    for (i in 1:dim(B)[3]) {
      lines(B[,,i], col="black")
    }
    text(values[,1:2], label=values[,3], pos=1)
    lines(lin, col="black")
  }
  archeo.colo.scale <- list(coords=A, colo=colo, values=values, values.lin=B, lin=lin)
  return(archeo.colo.scale)
}

add.archeo.colo.scale <- function(ae.colo.scale, text.cex=1) {
  # Add 'archeo.colo.scale' object to plot
  # Last Update: 2016/12/09
  # Dependencies: none
  # More details: see 'archeo.colo.scale' for more details
  # Notes: not yet implemented for rgl.
  
  # Arguments:
  #   ae.colo.scale:    'archeo.colo.scale' object
  #   text.cex:         size of label in the plot ('vector', num, 1)
  # Value:
  #   Graphical
  
  # examples:
  # plot(0, type="n", xlim=c(-100,100), asp=1, axes=F, xlab="", ylab="")
  # add.archeo.colo.scale(archeo.colo.scale(1:1000, colo=gray.colors(10)), text.cex=1)
  
  # plot(0, type="n", xlim=range(d), asp=1, axes=F, xlab="", ylab="")
  for (i in 1:dim(ae.colo.scale$coords)[3]) {
    polygon(ae.colo.scale$coords[,,i], col=ae.colo.scale$colo[i], border=ae.colo.scale$colo[i])
  }
  for (i in 1:dim(ae.colo.scale$values.lin)[3]) {
    lines(ae.colo.scale$values.lin[,,i], col="black")
  }
  text(ae.colo.scale$values[1:(dim(ae.colo.scale$values)[1]-1),1:2], label=ae.colo.scale$values[1:(dim(ae.colo.scale$values)[1]-1),3], pos=1, cex=text.cex)
  text(rbind(ae.colo.scale$values[dim(ae.colo.scale$values)[1],1:2]), label=paste(ae.colo.scale$values[dim(ae.colo.scale$values)[1],3], " mm"), pos=1, cex=text.cex)
  lines(ae.colo.scale$lin, col="black")
}

cspi <- function(cs.val, cs.val.x=0, cs.val.y=0, cs.radius=30, cs.x=0, cs.y=0, n=100, plots=F) {
  # Create 'CSPI' (circle sector preservation indicator)
  # Last Update: 2016/11/13
  # Dependencies: extrude3d {rgl}; rotate3d {rgl}; plot3d {rgl}; lines3d {rgl}
  
  # Arguments:
  #   cs.val:     CSPI value (in degrees) ('vector', num, 1)
  #               circular_section(centre = c(0,0), xy = t(mesh.drawing.decim$vb[1:3,]))*180/pi
  #   cs.radius:  CSPI radius (in mm) ('vector', num, 1)
  #   cs.x:       x-coordinate of the SCPI centre in the drawing (in mm) ('vector', num, 1)
  #   cs.y:       y-coordinate of the SCPI centre in the drawing (in mm) ('vector', num, 1)
  #               calculated as 'pottery$bpt'-'x' mm
  #   n:          nb of points defying the circle ('vector', num, 1)
  #   plots:      visualisation of the process/results ('logical')
  # Value:
  #   CSPI:       list containing 8 elements  ('list')
  #               '$cs.val' CSPI value (in degrees) ('vector', num, 1)
  #               '$cs.val.pos' xy-coordinate of the text of the SCPI value in the drawing (in mm) ('vector', num, 2)
  #               '$c.coords' matrix of xy-coordinates of the circle outline ('matrix', num, Nx2)
  #               '$cs.coords' matrix of xy-coordiantes of the circle sector outline ('matrix', num, Nx2)
  #               '$c.fill3d' mesh of the circle outline in 3D ('mesh')
  #               '$c.bord3d' matrix of xyz-coordinates of the circle outline in 3D ('matrix', num, Nx3)
  #               '$cc.fill3d' mesh of the circle sector in 3D ('mesh')
  #               '$cc.bord3d' matrix of xyz-coordinates of the circle sector outline in 3D ('matrix', num, Nx3)
  
  # examples:
  # cspi(20, plots=T)
  # cspi(90, plots=T)
  # cspi(210, plots=T)
  
  cs.val.pos <- c(cs.val.x, cs.val.y)
  cs.val.prc <- cs.val*100/360
  r <- cs.radius/10
  c.coords <- circle2d(n=n, r=r, cx=cs.x, cy=cs.y, plots=F)
  st <- 90+cs.val/2
  en <- 90-cs.val/2
  cs.coords <- circle.sect.2d(n=n, r=r, cx=cs.x, cy=cs.y, st=st, en=en)
  c.fill3d <- extrude3d(c.coords, thickness=0.001)
  c.fill3d <- rotate3d(c.fill3d, -pi/2, 1,0,0)
  c.bord3d <- circle3d.xz(n=n, r=r, cx=cs.x, cz=cs.y)
  c.bord3d <- rbind(c.bord3d,c.bord3d[1,])
  c.bord3d[,2] <- c.bord3d[,2]-0.004
  cs.fill3d <- extrude3d(cs.coords, thickness=0.001)
  cs.fill3d <- rotate3d(cs.fill3d, -pi/2, 1,0,0)
  cs.fill3d$vb[2,] <- cs.fill3d$vb[2,]-0.002
  cs.bord3d <- cbind(cs.coords[,1],0,cs.coords[,2])
  cs.bord3d <- rbind(cs.bord3d,cs.bord3d[1,])
  if (plots==T) {
    plot(c.coords, asp=1, type="n")
    polygon(c.coords, col="red", border = "black", lwd=2)
    polygon(cs.coords, col="green", border="violet", lwd=4)
    plot3d(c.fill3d, col="red", add=T)
    lines3d(c.bord3d, col="black", lwd=2)
    plot3d(cs.fill3d, col="green", add=T)
    lines3d(cs.bord3d, col="violet", lwd=4, add=T)
  }
  cspi <- list(c.coords=c.coords, cs.coords=cs.coords, cs.val=cs.val, cs.val.prc=cs.val.prc, cs.val.pos=cs.val.pos,
               c.fill3d=c.fill3d, c.bord3d=c.bord3d, cs.fill3d=cs.fill3d, cs.bord3d=cs.bord3d)
  return(cspi)
}

pottery.outline <- function(mesh, method="as", axa=1, axb=3, as.alpha=0.1, as.alpha.add.by=0.1, infos=F, plots=F) {
  # Calculation of the fragment outline, i.e. outline of the mesh which is projected to a single plane
  # Last Update: 2016/11/04
  # Dependencies: ashape {alphahull}, graph.edgelist {igraph}, degree {igraph}, clusters {igraph},
  #               E {igraph}, get.shortest.paths {igraph}, V {igraph}; vcgBorder {Rvcg}
  # Notes:  (i)   mesh should be (at least for 'vcgBorder') cutted by 'glVisible' function
  #         (ii)  'vcgBorder' is faster, but (i) gives lower nb of points which (ii) are not connected
  #         (iii) 'as' is slower, but (i) gives lots of points which (ii) are connected
  
  # Arguments:
  #   mesh:               triangular mesh ('mesh3d')
  #                       mesh should be cut in the front view
  #   method:             method used for outline extraction ("vector', char, 1)
  #                       'as' calculate border as the alpha shape of the vertices projected onto 'axa' and 'axb' plane
  #                       'vcgBorder' calculate border by 'vcgBorder' function
  #   as.axa:             the first axis in which coordinates are projected ('vector', num, 1)
  #   as.axb:             the second axis in which coordinates are projected  ('vector', num, 1)
  #   as.alpha:           the initial alpha value for ashape calculation  ('vector', num, 1)
  #   as.alpha.add.by:    If the problem with the edge connection is found during  ('vector', num, 1)
  #                       the calculation of the ashape with the initial alpha set to 'as.alpha',
  #                       the 'as.alpha' is iteratively augmented by 'as.alpha.add.by' until the
  #                       ashape is properly calculated ('vector')
  #   plots:              visualisation of the process/results ('logical')
  # Value:
  #   bord:               xy-coordinates of the outline  ('matrix', num, Nx2)
  
  # examples:
  # pottery.outline(mesh, method="as", plots=T)
  # pottery.outline(mesh, method="vcgBorder", plots=T)
  
  if (method=="as") {
    pla <- t(mesh$vb[c(axa,axb),])
    pts <- unique(pla)
    err <- 1
    while (err>0) {
      if (infos==T) { print(paste("alpha=",as.alpha)) }
      a <- ashape(pts, alpha = as.alpha)
      a.g <- graph.edgelist(cbind(as.character(a$edges[, "ind1"]), as.character(a$edges[,"ind2"])), directed = FALSE)
      err <- 0
      if (!is.connected(a.g)) { err <- err+1 }
      if (any(degree(a.g) != 2)) { err <- err+1 }
      if (clusters(a.g)$no > 1) { err <- err+1 }
      if (err>0) { as.alpha <- as.alpha + as.alpha.add.by }
    }
    cutg <- a.g - E(a.g)[1]
    ends <- names(which(degree(cutg) == 1))
    path <- get.shortest.paths(cutg, ends[1], ends[2])[[1]]  
    pathX <- as.numeric(V(a.g)[path[[1]]]$name)
    pathX <- c(pathX, pathX[1])
    bord <- a$x[pathX, ]
  }
  if (method=="vcgBorder") {
    pts <- t(mesh$vb[c(axa,axb),])
    bord.mesh <- vcgBorder(mesh)
    bord <- pts[which(bord.mesh$bordervb == 1),]
  }
  if (plots==T) {
    plot(pts, lwd = 1, col = "gray", asp=1)
    if (method=="as") {         lines(bord, lwd = 2, col="red") }
    if (method=="vcgBorder") {  lines(bord, lwd = 1, col="red"); points(bord, col="red", pch=19) }
  }
  return(bord)
}


# Support functions for illustration
mesh.align2front <- function(mesh, plots=F) {
  # Align the centre of the model with the y axis (i.e. with c(0,1,0))
  # Created: 2016/11/01
  # Last Update: 2016/11/01
  # Dependencies: plot3d {rgl}, lines3d {rgl}, rotate3d {rgl},
  #               rotV {cwhmisc}, updateNormals {Morpho}
  
  # Arguments:
  #   mesh:     triangular mesh ('mesh3d')
  #   plots:    visualisation of the process/results ('logical')
  # Value:
  #   mesh:     triangular mesh whose centre is aligned with the y-axis ('mesh3d')
  
  # examples:
  # mesh.align2front(mesh, plots=T)
  
  if (plots==T) { plotOrientedMesh(mesh) }
  mid <- sum((rho.range(mesh, plots=F)))/2
  mid.pt <- c(pol2cart(c(mid,1)),0)
  mid.pt.rot <- rotV(mid.pt,c(0,-1,0))
  mesh <- rotate3d(mesh, matrix=mid.pt.rot)
  mesh <- updateNormals(mesh)
  if (plots==T) { plot3d(mesh, meshColor="legacy", add=T) }
  return(mesh)
}

mesh.transform4drawing <- function(mesh, rho, x, plots=F) {
  # Transformation used in drawing
  # Dependencies: par3d {rgl}, plot3d {rgl}, lines3d {rgl}, rotate3d {rgl}, updateNormals {Morpho}
  
  # Arguments:
  #   mesh:     triangular mesh ('mesh3d')
  #   rho:      rotation parameter (in radians) ('vector', num, 1)
  #   x:        translation parameter (in mm) ('vector, num, 1)
  #   plots:    visualisation of the process/results ('logical')
  # Value:
  #   mesh:     transformed triangular mesh ('mesh3d')
  
  # examples:
  # mesh.transform4drawing(mesh, rho=0.3, x=0.2, plots=T)
  
  if (plots==T) { plot3d(mesh, meshColor="legacy") ; par3d(windowRect=c(960,30,1920,1040)) }
  mesh <- rotate3d(mesh,rho,0,0,1)
  mesh <- translate3d(mesh,x,0,0)
  mesh <- vcgUpdateNormals(mesh)
  if (plots==T) { plot3d(mesh, meshColor="legacy", add=T) }
  return(mesh)
}

mesh.centroid <- function(mesh) {
  # Calculate coordinates of the mesh centre
  centroid <- apply(mesh$vb[1:3,],1,mean)
  return(centroid)
}

cut.img <- function(img, plots=F) {
  # Clip png image (erase white parts)
  # Last Update: 2016/11/24
  # Dependencies: as.raster {grDevices}
  
  # Arguments:
  #   img:        image matrix ('raster array'), see 'readPNG' {png} for more information
  # Value:
  #   img.cutted: clipped image ('raster array') 
  #   plots:      visualisation of the process/results ('logical')
  
  # examples:
  # cut.img(img, plots=T)
  
  img.trip <- img[,,1]+img[,,2]+img[,,3]
  xxs <- which(apply(img.trip,2,sum)!=dim(img.trip)[1]*3)
  yys <- which(apply(img.trip,1,sum)!=dim(img.trip)[2]*3)
  img.cutted <- img[yys[1]:yys[length(yys)],xxs[1]:xxs[length(xxs)],]
  if (plots==T) {
    x <- img.trip[,1]
    y <- img.trip[,2]
    xx <- dim(img.trip)[1]
    yy <- dim(img.trip)[2]
    plot(0, xlim=c(0,xx),ylim=c(0,yy), type="n", asp=1)
    for (i in 1:dim(img.trip)[1]) {
      for (j in 1:dim(img.trip)[2]) {
        if (img.trip[i,j]==3) { points(i,j, col="black") }
        if (img.trip[i,j]!=3) { points(i,j, col="black", pch=19) }
      }
    }
    plot(as.raster(img.cutted))
  }
  return(img.cutted)
}

circle2d <- function(n=300, cx=0, cy=0, r=1, plots=F) {
  # 2D circle
  # Last Update: 2016/11/12

  # Arguments:
  #   n:          number of points on the circle ('vector', num, 1)
  #   cx:         x-coordinate of the circle centre ('vector', num, 1)
  #   cy:         y-coordinate of the circle centre ('vector', num, 1)
  #   r:          radius ('vector')
  #   plots:      visualisation of the process/results ('logical')
  # Value:
  #   circle2d:   coordinates of the circle ('matrix', num, Nx2)
  
  # examples:
  # circle2d(n=20, cx=0, cy=0, r=1, plots=T)
  
  theta <- seq(0, 2*pi, len=n) 
  x <- cos(theta)*r+cx
  y <- sin(theta)*r+cy
  if (plots==T) { plot(x,y, asp=1) ; lines(x,y) }
  circle2d <- cbind(x,y)
  return(circle2d)
}

draw.circle2d <- function(n=300, cx=0, cy=0, r=1, col=NA, col.border="black", lty=1, lwd=1) {
  # 2D circle
  # Last Update: 2016/11/12
  
  # Arguments:
  #   n:          number of points on the circle ('vector', num, 1)
  #   cx:         x-coordinate of the circle centre ('vector', num, 1)
  #   cy:         y-coordinate of the circle centre ('vector', num, 1)
  #   r:          radius ('vector')
  #   plots:      visualisation of the process/results ('logical')
  # Value:
  #   graphic
  
  # examples:
  # plot(0,0, asp=1); draw.circle2d(col="green", col.border="red", lty=3, lwd=3)
  
  xy <- circle2d(n=n, cx=cx, cy=cy, r=r, plots=F)
  x <- xy[,1]; y <- xy[,2]
  
  polygon(x,y, border=col.border, col=col, lty=lty, lwd=lwd)
}

circle3d.xz <- function(n=300, cx=0, cz=0, r=1, plots=F) {
  # 3D circle
  # Last Update: 2016/11/12
  # Dependencies: plot3d {rgl}, lines3d {rgl}
  
  # Arguments:
  #   n:          number of points on the circle ('vector', num, 1)
  #   cx:         x-coordinate of the circle centre ('vector', num, 1)
  #   cz:         y-coordinate of the circle centre ('vector', num, 1)
  #   r:          radius ('vector', num, 1)
  #   plots:      visualisation of the process/results ('logical')
  # Value:
  #   circle3d:   xyz-coordinates of the circle ('matrix', num, Nx3)
  
  # examples:
  # circle3d.xz(n=20, cx=0, cz=0, r=1, plots=T)
  
  theta <- seq(0, 2*pi, len=n) 
  x <- cos(theta)*r+cx
  z <- sin(theta)*r+cz
  y <- rep(0, n) 
  if (plots==T) { plot3d(x,y,z, asp="iso") ; lines3d(x,y,z) }
  circle3d <- cbind(x,y,z)
  return(circle3d)
}

circle3d <- function(n=1000, axi="x", cen=c(0,0,0), r=1, plots=F, add=T, cir.col="black", cir.lwd=1) {
  # 3D circle which is perpendicular to one of the axes (x/y/z)
  # Last Update: 2016/12/11
  # Dependencies: plot3d {rgl}, lines3d {rgl}
  
  # Arguments:
  #   n:          number of points on the circle ('vector', num, 1)
  #   axi:        axis perpendicular to the the circle ('vector', char, 1)
  #               'x' for x-axis, 'y' for y-axis, 'z' for z-axis
  #   cen:        xyz-coordinates of the circle centre ('vector', num, 3)
  #   r:          radius ('vector', num, 1)
  #   plots:      visualisation of the process/results ('logical')
  #   add:        if the circle is drawn ('logical')
  #   cir.col:    color of the circle ('vector', num/char, 1)
  #   cir.lwd:    thickness of the circle ('vector', num, 1)
  # Value:
  #   circle3d:   xyz-coordinates of the circle ('matrix', num, Nx3)
  
  # examples:
  # circle3d(n=30, axi="x", cen=c(10,20,30), r=30, plots=T, add=F, cir.col="red", cir.lwd=1)
  # circle3d(n=60, axi="y", cen=c(10,20,30), r=60, plots=T, add=T, cir.col="green", cir.lwd=3)
  # circle3d(n=90, axi="z", cen=c(10,20,30), r=90, plots=T, add=T, cir.col="blue", cir.lwd=6)
  
  if (add==T) { plots=T }
  theta <- seq(0, 2*pi, len=n)
  if (axi=="x") {
    x <- rep(0,n)
    y <- cos(theta)*r
    z <- sin(theta)*r
  }
  if (axi=="y") {
    x <- cos(theta)*r
    y <- rep(0,n)
    z <- sin(theta)*r
  }
  if (axi=="z") {
    x <- cos(theta)*r
    y <- sin(theta)*r
    z <- rep(0,n)
  }
  circle3d <- cbind(x,y,z)
  circle3d <- circle3d+matrix(cen,dim(circle3d)[1],dim(circle3d)[2],byrow=T)
  if (plots==T) {
    if (add==F) { plot3d(rbind(cen), asp="iso") }
    # points3d(rbind(cen), col=cir.col, size=cir.lwd)
    lines3d(circle3d, col=cir.col, size=cir.lwd)
  }
  return(circle3d)
}

circle.sect.2d <- function(n=300, r=1, cx=0, cy=0, st=0, en=90, plots=F) {
  # 2D circle sector
  # Last Update: 2016/11/13

  # Arguments:
  #   n:            number of points on the circle ('vector', num, 1)
  #   cx:           x-coordinate of the circle centre ('vector', num, 1)
  #   cy:           y-coordinate of the circle centre ('vector', num, 1)
  #   st:           starting angle of the sector (in degrees) ('vector', num, 1)
  #   en:           ending angle of the sector (in degrees) ('vector', num, 1)
  #   r:            radius ('vector', num, 1)
  # Value:
  #   circle.sect:  xy-coordinates of the circle sector ('matrix', num, Nx2)
  
  theta <- seq(st/180*pi, en/180*pi, len=n) 
  x <- cos(theta)*r+cx
  y <- sin(theta)*r+cy
  circle.sect <- rbind(c(cx,cy),cbind(x,y))
  if (plots==T) {
    plot(circle.sect, asp=1)
    polygon(circle.sect, col="red")
  }
  return(circle.sect)
}

A2M <- function (A) {
  # Array to matrix
  M <- matrix(numeric(), 0, dim(A)[2])
  for (i in 1:dim(A)[3]){ M <- rbind(M,A[,,i]) }
  return (M)
}

directional.lighting <- function(mesh, light, method="vertices") {
  # Calculation of the directional lighting
  # Last Update: 2016/11/05
  # Dependencies: facenormals {Morpho}
  # More details:
  # Notes: (i) it is faster if the mesh is preliminary cut in the front view (by 'glVisible' function)
  
  # Arguments:
  #   mesh:         triangular mesh ('mesh3d')
  #   light:        xyz-coordinates of the light source ('vector', num, 3)
  #   method:       feature on which the calculation is performed ('vector', char, 1)
  #                 'vertices'
  #                 'faces' for face normals
  # Value:
  #   DL:           directional lighting values ('vector', num, N)
  
  # examples:
  # directional.lighting(mesh, c(0,0,0), method="vertices")
  # directional.lighting(mesh, c(0,0,0), method="faces")
  
  if (method=="vertices") {
    pts <- t(mesh$vb[1:3,])
    normals <- t(mesh$normals[1:3,])
  }
  if (method=="faces") {
    pts <- t(facenormals(mesh)$vb)[,1:3]
    normals <- t(facenormals(mesh)$normals)[,1:3]
  }
  theta <- rep(NA, nrow(normals))
  light.pts <- matrix(NA,nrow=nrow(normals),ncol=3)
  for (i in 1:nrow(normals)) {
    light.pts[i,] <- pts[i,]-light
    theta[i] <- acos( sum(light.pts[i,]*normals[i,]) / ( sqrt(sum(light.pts[i,] * light.pts[i,])) *
                                                           sqrt(sum(normals[i,] * normals[i,])) ) )
  }
  DL <- cos(theta)
  return(DL)
}

val.filter <- function(val, plots=F) {
  # Filter used in Ambiant occlusion visualisation
  # return logical indicating if the value is found within the uniform distribution
  # Last Update: 2016/11/21
  # Dependencies: rescale {plotrix}
  # Notes: function needs to be checked
  
  # Arguments:
  #   val:        values ('vector', num, N)
  #   plots:      visualisation of the process/results ('logical')
  # Value:
  #   val.filter: whatever the value is within the given filter ('vector', logical, N)
  
  # examples:
  # val.filter(val=1:1000, plots=T)
  
  diff <- rescale(val, newrange=c(0, 1))
  prob <- runif(length(diff))
  val.filter <- prob<diff
  if (plots==T) { hist(prob); hist(diff, col="red", add=T); hist(diff[val.filter], col="orange", add=T) }
  return(val.filter)
}

mesh.error <- function(mesh, method="profile", profil=NULL, clust=0, plots=F, infos=F) {
  # Calculate the rotation axis quality/vase regularity
  # Last Update: 2016/11/11
  # Dependencies: foreach {foreach}; makeCluster {parallel}; registerDoParallel {doParallel}; stopImplicitCluster {doParallel}
  # Notes: (i) axis2 works only for one side
  
  # Arguments:
  #   mesh:       triangular mesh ('mesh3d')
  #   method:     method of calculation ('vector', char, 1)
  #               'axis' for calculation the radius of each point
  #               'axis2' for calculation of radius of points regrouped by z-values
  #               'profile' for calculation of the closest distances between each point of the mesh
  #               from the form given by the profile. Whole envelope profile is priviledged
  #   profil:     list containing 3 elements required for 'profile' method ('list')
  #               '$profil' contains all profile coordinates ('matrix', num, Nx2)
  #               '$extern' contains outer profile coordinates ('matrix', num, Nx2)
  #               '$intern' contains inner profile coordinates ('matrix', num, Nx2)
  #   clust:      number of cores for parallel ('vector', num, 1)
  #               '0' for automatical detection of cores
  #   plots:      visualisation of the process/results ('logical')
  #   infos:      print the information/result ('logical')
  # Value:
  #   err:        quality of alignment values for each point ('vector', num, N)
  
  # examples:
  # mesh.error(mesh, method="axis")
  # mesh.error(mesh, method="axis2")
  # mesh.error(mesh, method="profile", profil=profil)
  
  if (method=="axis") {
    xyz <- t(mesh$vb[1:3,])
    R <- sqrt(xyz[,1]^2+xyz[,2]^2)
    err <- R
  }
  if (method=="axis2") {
    xyz <- t(mesh$vb[1:3,])
    z <- xyz[,3]
    err <- rep(NA,dim(xyz)[1])
    brks <- seq(min(z),max(z),length.out = round(max(z)-min(z))/meshres(mesh))
    for (i in 1:(length(brks)-1)) {
      sel <- which((z>=brks[i])&(z<brks[i+1]))
      d <- sqrt(xyz[sel,1]^2+xyz[sel,2]^2)
      mean.d <- mean(d)
      r <- d-mean.d
      err[sel] <- r
    }
    err[which(z==max(z))] <- 0
  }
  if (method=="profile") {
    xyz <- t(mesh$vb[1:3,])
    R <- sqrt(xyz[,1]^2+xyz[,2]^2)
    RZ <- cbind(R,xyz[,3])
    pr <- profil$profil
    pr[,1] <- abs(pr[,1])
    if (plots==T) { plot(RZ, asp=1); lines(pr, col="red", lwd=3) }
    
    # Can be calculated also as...
    # err <- rep(NA,dim(xyz)[1])
    # for (i in 1:length(err)){ err[i] <- min(sqrt(apply((matrix(RZ[i,],dim(pr)[1],dim(pr)[2],byrow=T)-pr)^2,1,sum))) }
    
    if (clust==0) { clust=detectCores() }
    cl <- makeCluster(clust)
    registerDoParallel(cl)
    dista <- function(i,RZ,pr) {
      return(min(sqrt(apply((matrix(RZ[i,],dim(pr)[1],dim(pr)[2],byrow=T)-pr)^2,1,sum))))
    }
    if (infos==T) { print("Parallel started") ; st <- Sys.time()}
    err <- foreach(i=1:dim(xyz)[1]) %dopar% dista(i,RZ,pr)
    stopImplicitCluster()
    stopCluster(cl)
    if (infos==T) { print(paste("Parallel finished in ",round(c(Sys.time()-st)),"secs.")) }
    err <- unlist(err)
  }
  
  return(err)  
}

volume <- function(vol.perc) {
  inter <- pottery$intern
  vol.perc[1] <- vol.perc[1]+1
  zz <- seq(min(inter[,2]),max(inter[,2]),length=100)
  ind <- which(inter[,2]>zz[vol.perc[1]] & inter[,2]<zz[vol.perc[2]])
  M <- inter[ind,]
  if (length(ind)<3) { return(NULL) }
  MM <- rbind(M[dim(M)[1]:1,],cbind(-M[,1],M[,2]))
  polygon(MM, col="lightblue", border = "black")
  V <- numeric(dim(inter)[1]-1)
  for (i in 1:dim(M)[1]-1) {
    v <- abs(M[i+1,2]-M[i,2])
    r1 <- abs(M[i,1])
    r2 <- abs(M[i+1,1])
    VV <- pi*((r1+r2)/2)^2 * v       # VV <- 1/3*pi*v*(r1^2+r1*r2+r2^2)
    V[i] <- VV
  }
  V <- sum(V)/10^6
  pos <- apply(MM, 2, mean)
  text(pos[1], pos[2], labels=paste(round(V,3),"L"))
}

rulerF <- function() {
  points(ruler, pch=3, col="red")
  if (dim(ruler)[1]>1) {
    i=1
    while(i<dim(ruler)[1]){
      dista <- sqrt(sum((ruler[i,]-ruler[i+1,])^2))
      ce <- apply(ruler[i:(i+1),],2,mean) 
      lines(ruler[i:(i+1),], col="red")
      text(ce[1],ce[2],labels=paste(round(dista),"mm"), cex=0.8, col="red", pos=3)
      i=i+2
    }
  }
}

updateLIN <- function() {
  LIN <- array(numeric(0), dim=c(2,2,dim(pts)[1]))
  for (i in 1:dim(pts)[1]) { LIN[,,i] <- lines.2d(pottery, click.coo=pts[i,])  }
  LIN <<- LIN
}

updateREC <- function() {
  REC <<- array(numeric(0), dim=c(40,2,round(dim(ptt)[1]/2)))
  j = seq(from=1, to=dim(ptt)[1], by=2)
  k <- 1:length(j)
  for (i in 1:length(k)) {
    REC[,,k[i]] <- profile2d.rec(profil, ptt=ptt[j[i]:(j[i]+1),])
  }
  REC <<- REC
}


# 3D Reconstruction
reconst.3d.slices <- function(profil, plots=F, seg=360) {
  # Create the array containing the 3D pointcloud of the rotationally symmetrical object defined by its profile (here 'profile' object)
  # Last Update: 2016/12/09
  # Dependencies: rotate3d {rgl}, lines3d {rgl}
  
  # Arguments:
  #   profil:     profil object, i.e. list containing 3 elements ('list')
  #               '$profil' contains all profile coordinates ('matrix', num, Nx2)
  #               '$extern' contains outer profile coordinates ('matrix', num, Nx2)
  #               '$intern' contains inner profile coordinates ('matrix', num, Nx2)
  #   plots:      visualisation of the process/results ('logical')
  #   seg:        number of segments ('vector', num, 1)
  # Value:
  #   reconst:    array containing xyz-coordinates of the reconstructed model ('array', num, Nx3xseg)
  
  # examples:
  # reconst.3d.slices(profil, plots=T, seg=20)
  
  pro <- profil$profil
  pro <- cbind(pro[,1],0,pro[,2])
  reconst <- matrix(numeric(0),0,3)
  sig <- seq(from=-1, to=1, length=seg)
  reconst <- array(numeric(0),dim=c(dim(pro)[1], 3, seg))
  for(i in 1:seg) { 
    one.seg <- rotate3d(pro,sig[i]*pi,0,0,1)
    reconst[,,i] <- one.seg
    if (plots==T) { lines3d(one.seg, col="brown",size=3) }
  }
  return(reconst)
}

reconst.3d.slices.plot <- function(reconst, type="l", col="brown", size=3) {
  # Plot the 3D reconstruction of the symmetrical model
  # Last Update: 2016/12/09
  # Dependencies: points3d {rgl}, lines3d {rgl}
  # More details: see 'reconst.3d.slices' for more information
  
  # Arguments:
  #   reconst:    array containing xyz-coordinates of the reconstructed model ('array', num, Nx3xseg)
  #   type:       type of the reconstruction ('vector', char, 1)
  #               'p' for pointcloud
  #               'l' for lines of each slice
  #   col:        color of the reconstruction ('vector', char, 1)
  #   size:       size of the point/line in the plot ('vector', num, 1)
  # Value:
  #   Graphical
  
  # examples:
  # reconst.3d.slices.plot(reconst.3d.slices(profil, plots=F, seg=20), type="l")
  # reconst.3d.slices.plot(reconst.3d.slices(profil, plots=F, seg=20), type="p")
  
  if (type=="p") {
    M <- matrix(numeric(0),0,3)
    for (i in 1:dim(reconst)[3]) { M <- rbind(M,reconst[,,i]) }  
    points3d(M, col=col, size=size)
  }
  if (type=="l") {
    for (i in 1:dim(reconst)[3]) { lines3d(reconst[,,i], col=col, size=size) }  
  }
}

reconst.3d.mesh <- function (profil, plots=T, seg=128, col="grey") {
  # Create the mesh reconstruction of the rotationally symmetrical object defined by its profile (here 'profile' object)
  # Last Update: 2017/01/21
  # Dependencies: turn3d {rgl}; rotate3d {rgl}; shade3d {rgl}; updateNormals {Morpo}; 
  
  # Arguments:
  #   profil:     profil object, i.e. list containing 3 elements ('list')
  #               '$profil' contains all profile coordinates ('matrix', num, Nx2)
  #               '$extern' contains outer profile coordinates ('matrix', num, Nx2)
  #               '$intern' contains inner profile coordinates ('matrix', num, Nx2)
  #   plots:      visualisation of the process/results ('logical')
  #   seg:        number of segments ('vector', num, 1)
  #   col:        color of the mesh ('vector', char, 1)
  # Value:
  #   reconst:    triangular mesh of the reconstructed model ('mesh')
  
  # examples:
  # reconst.3d.mesh(profil, plots=T, seg=20)
  
  pro <- profil$profil
  pro <- pro[,2:1]
  if (sum(pro[,2]<0) >0) { pro[,2] <- -pro[,2] }
  reconst <- turn3d(pro, n=seg, smooth=T)
  reconst <- rotate3d(reconst,0.5*pi,0,1,0)
  reconst <- updateNormals(reconst)
  reconst$material$color <- matrix(col,3,dim(reconst$it)[2])
  if (plots==T) { shade3d(reconst, col=col, override=FALSE) }
  return(reconst)
}


