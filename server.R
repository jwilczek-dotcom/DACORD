#..............................................................................................
# D A C O R D  (Dessin Assiste de la CERamique par ORDinateur) -  Server
#..............................................................................................
# Last update: 2019/09/27

# SET PATH TO DACORD FOLDER:
WD <- c("D:/OneDrive/DACORD")

#..............................................................................................
# To increase the calculation time, some calculations are performed externally via MeshLab:
# '_01_clean_mesh.mlx' for mesh cleaning
# '_02_reorient_normals.mlx' for updating normals
# '_03_ambient_occlusion.mlx' for ambient occlusion
# '_04_reorient_normals_mesh_export.mlx' for updating normals
# The vessel regularity is calculated externally via Rscript ('_support_ME.R')

#..............................................................................................
# Update 2019/09/27:
# 1) Placing 3 points (fixed)
#     'Digit.Fixed.New' function stopped to work properly. Functon was replaced by altered 'digit.fixed'
#     function {geomorph}
#     However this may create issues because closest point on the mesh can sometimes be identified incorrectly
#     on the other side of the mesh. To minimise occurence of this situation, 3 points should be clicked precisely
#     on 3 points of the mesh 
# 2) Problem with texture (fixed)
#     Since {rgl} version later then 0.100.1, a 'meshColor' argument must be specified in the call of 'plot3d' function


#..............................................................................................
# LIBRARIES
# 'server.R' file:
library(shiny)
library(png)
library(mco)
library(optimx)
# 'DACORD_Functions.R' file:
library(rgl)
library(Rvcg)
library(Morpho)
library(sp)
library(parallel)
library(foreach)
library(iterators); library(doParallel)
library(fit.models); library(MASS); library(robustbase) ;library(rrcov); library(robust)
library(pracma)
library(alphahull)
library(igraph)
# ..............................................................................................

#..............................................................................................
# Do not alter
setwd(WD)
source("functions.R")
set.seed(12345)
options(shiny.maxRequestSize=200*1024^2)
PlotPath <<- paste(WD,"/temp/plot.png", sep="")
PlotQualPath <<- paste(WD,"/temp/plotQual.png", sep="")
cs <<- NULL
Rot.count <<- Rot.count2 <<- Rot.count3 <<- Rot.count4 <<- Rot.count.AODL <<- 1000
Trans.count <<- Trans.count2 <<- Trans.count3 <<- Trans.count4 <<- Trans.count.AODL <<- 1000
Lat.count <<- Lat.count2 <<- Lat.count3 <<- Lat.count.AODL <<- 1000
Lon.count <<- Lon.count2 <<- Lon.count3 <<- Lon.count.AODL <<- 1000
Dis.count <<- Dis.count2 <<- Dis.count3 <<- Dis.count.AODL <<- 1000
Lims.count <<- c(-1,1000)
Lims.count.DL <<- Lims.count.AO <<- c(-1,1000)
Add.light.count3 <<- T
Col.count3 <<- ""
Meth.count4 <<- ""
pts <<- matrix(numeric(0),0,2)
ptt <<- matrix(numeric(0),0,2)
ruler <<- matrix(numeric(0),0,2)
REC <<- array(numeric(0), dim=c(40,2,0))
LIN <<- array(numeric(0), dim=c(2,2,0))
Add.model.count <<- F
Slice.count <<- 0
#..............................................................................................

#..............................................................................................
# Parameters (Can be altered)
Infos <<- F       # whether information is printed in R Console
dur <<- 5         # duration of the message in Shiny
Itnmax <<- 100    # maximum iteration number in 'optimix'
Ps.PtS <<- 1      # postscript point size
#..............................................................................................

#..............................................................................................
# SHINY SERVER
shinyServer(function(input, output, session) {
  # # DACORD - ORIENTATION
  dataInputMesh <- reactive({
    inFile <- input$fileMesh
    if (is.null(inFile)) return(NULL)
    FileName <- inFile$name
    FileExt <- substr(FileName, nchar(FileName)-2, nchar(FileName))
    if (FileExt!="ply") { showNotification(paste("Please load model in PLY format"), duration = dur); return(NULL) }
    file <- inFile$datapath
    filen <- gsub("[\\]", "/", file)
    filen <- paste(substr(filen, 1, nchar(filen)-1), inFile$name, sep="")
    file.rename(file, filen)
    ### DEBUG ###
    # filen <<- "D:/OneDrive/DACORD/examples/A011.ply"
    ### END DEBUG ###
    mesh <- vcgImport(filen, clean=F, readcolor=T)
    if (is.null(mesh$material$color)) { mesh$material$color <- matrix(grey.colors(3)[2],3,dim(mesh$it)[2]) }
    mesh <<- mesh
    mesh.decim <<- decim(mesh,10000)
    return(mesh)
  })
  
  observeEvent(input$fileMesh, {
    withProgress(message = 'Preparing model', value = 0.1, {
      setProgress(0.5)
      mesh <- dataInputMesh()
      plot3d.bl()
      plot3d(mesh.decim, col="grey", asp="iso")
      widget <- rglwidget()
      output$plotOrientation3d <- renderRglwidget(widget)
      rgl.close()
      output$plotOrientation <- renderPlot({
        # plot(mesh.decim$vb[1,],mesh.decim$vb[2,], asp=1)
      })
      setProgress(1)
    })
  })
  
  observeEvent(input$btnPreorientation, {
    if (input$selPreorientationMethod=="") { showNotification(paste("No pre-orientation method selected. Please select pre-orientation method."), duration = dur);  return(NULL) }
    if (!exists("mesh")) { showNotification(paste("No file loaded. Please load file."), duration = dur); return(NULL) }
    if (input$selPreorientationMethod=="Automatic") { dataInputPreorientAutomatic() }
    if (input$selPreorientationMethod=="Manual") { dataInputPreorientManual() }
    mesh.preorient.decim <- decim(mesh.preorient,1000)
    plot3d.bl()
    plot3d(mesh.preorient.decim, col="lightgrey", asp="iso", alpha=0.6)
    plot3d(mesh.trip$S1, col="red", add=T)
    plot3d(mesh.trip$S2, col="blue", add=T)
    lines3d(rbind( c(0,0,min(mesh.preorient$vb[3,])), c(0,0,max(mesh.preorient$vb[3,])) ), lwd=2, col="red")
    widget <- rglwidget()
    output$plotOrientation3d <- renderRglwidget(widget)
    rgl.close()
    output$plotOrientation <- renderPlot({
      layout(matrix(1:3,1,3))
      r <- sqrt(mesh.preorient$vb[1,]^2+mesh.preorient$vb[2,]^2)
      plot(cbind(r,mesh.preorient$vb[3,]), asp=1, cex=0.5, main="Profile plane", xlab="r", ylab="z", xlim=c(0,max(r)))
      points(cbind(sqrt(mesh.trip$S1$vb[1,]^2+mesh.trip$S1$vb[2,]^2),mesh.trip$S1$vb[3,]), col="red", pch=19, cex=0.5)
      points(cbind(sqrt(mesh.trip$S2$vb[1,]^2+mesh.trip$S2$vb[2,]^2),mesh.trip$S2$vb[3,]), col="blue", pch=19, cex=0.5)
      lines(rbind(c(0,min(mesh.preorient$vb[3,])),c(0,max(mesh.preorient$vb[3,]))), col="red", lwd=2)
      plot.blank()
      plot.blank()
    })
    mesh.preorient.decim <<- mesh.preorient.decim
    showNotification(paste("Pre-orientation finished"), duration = dur)
  })
  
  observeEvent(input$btnInverse, {
    if (!exists("mesh")) { showNotification(paste("No file loaded. Please load file."), duration = dur); return(NULL) }
    if (!exists("mesh.preorient")) { showNotification(paste("No pre-oriented model. Please pre-orient model."), duration = dur); return(NULL) }
    if (!exists("mesh.trip")) { showNotification(paste("No pre-oriented model. Please pre-orient model."), duration = dur);  return(NULL) }
    mesh.preorient <- mesh.flip.vertically(mesh.preorient, plots=F)
    mesh.trip <- list(S1=mesh.flip.vertically(mesh.trip$S1, plots=F), S2=mesh.flip.vertically(mesh.trip$S2, plots=F), mesh=mesh.flip.vertically(mesh.trip$mesh, plots=F))
    mesh.preorient.decim <- decim(mesh.preorient,1000)
    plot3d.bl()
    plot3d(mesh.preorient.decim, col="lightgrey", asp="iso", alpha=0.6)
    plot3d(mesh.trip$S1, col="red", add=T)
    plot3d(mesh.trip$S2, col="blue", add=T)
    lines3d(rbind( c(0,0,min(mesh.preorient$vb[3,])), c(0,0,max(mesh.preorient$vb[3,])) ), lwd=2, col="red")
    widget <- rglwidget()
    output$plotOrientation3d <- renderRglwidget(widget)
    rgl.close()
    output$plotOrientation <- renderPlot({
      layout(matrix(1:3,1,3))
      r <- sqrt(mesh.preorient$vb[1,]^2+mesh.preorient$vb[2,]^2)
      plot(cbind(r,mesh.preorient$vb[3,]), asp=1, cex=0.5, main="Profile plane", xlab="r", ylab="z", xlim=c(0,max(r)))
      points(cbind(sqrt(mesh.trip$S1$vb[1,]^2+mesh.trip$S1$vb[2,]^2),mesh.trip$S1$vb[3,]), col="red", pch=19, cex=0.5)
      points(cbind(sqrt(mesh.trip$S2$vb[1,]^2+mesh.trip$S2$vb[2,]^2),mesh.trip$S2$vb[3,]), col="blue", pch=19, cex=0.5)
      lines(rbind(c(0,min(mesh.preorient$vb[3,])),c(0,max(mesh.preorient$vb[3,]))), col="red", lwd=2)
      plot.blank()
      plot.blank()
    })
    mesh.preorient <<- mesh.preorient
    mesh.trip <<- mesh.trip
    mesh.preorient.decim <<- mesh.preorient.decim
    showNotification(paste("Model inversed"), duration = dur)
  })
  
  observeEvent(input$btnOrientation, {
    if (input$selOrientationMethod=="") { showNotification(paste("No Orientation method selected. Please select Orientation method"), duration = dur);  return(NULL) }
    if (!exists("mesh")) { showNotification(paste("No file loaded. Please load file."), duration = dur); return(NULL) }
    if (!exists("mesh.preorient")) { showNotification(paste("No pre-oriented model. Please pre-orient model."), duration = dur); return(NULL) }
    if (!exists("mesh.trip")) { showNotification(paste("No pre-oriented model. Please pre-orient model."), duration = dur);  return(NULL) }
    showNotification(paste("Orientation running"), duration = dur)
    if (input$selOrientationMethod=="AOS_none") { dataInput_AOS_none() }
    if (input$selOrientationMethod=="AOS_sections") { dataInput_AOS_sections() }
    if (input$selOrientationMethod=="AOS_circle3DDL") { dataInput_AOS_circle3DDL() }
    if (input$selOrientationMethod=="AOS_circle4DDL") { dataInput_AOS_circle4DDL() }
    if (input$selOrientationMethod=="AOS_profile") { dataInput_AOS_profile() }
    if (input$selOrientationMethod=="AOS_polynomial") { dataInput_AOS_polynomial() }
    if (input$selOrientationMethod=="AOS_rimbase_rim") { dataInput_AOS_rimbase_rim() }
    if (input$selOrientationMethod=="AOS_rimbase_base") { dataInput_AOS_rimbase_base() }
    mesh.orient <- mesh.orient
    mesh.orient.decim <- decim(mesh.orient, 10000)
    plot3d.bl()
    plotOrientedMesh(mesh.orient.decim, method="3d", col="lightgrey")
    widget <- rglwidget()
    output$plotOrientation3d <- renderRglwidget(widget)
    rgl.close()
    output$plotOrientation <- renderPlot({
      layout(matrix(1:3,1,3))
      plotOrientedMesh(mesh, method="rz")
      plot.blank()
      plot.blank()
    })
    showNotification(paste("Orientation finished"), duration = dur)
  })
  
  observeEvent(input$btnSave, {
    if (!exists("mesh")) { showNotification(paste("No file loaded. Please load file."), duration = dur); return(NULL) }
    if (!exists("mesh.preorient")) { showNotification(paste("No pre-oriented model. Please pre-orient model."), duration = dur); return(NULL) }
    if (input$selOrientationMethod!="none") {
      if (!exists("mesh.trip")) { showNotification(paste("No pre-oriented model. Please pre-orient model."), duration = dur);  return(NULL) }
      if (!exists("mesh.orient")) { showNotification(paste("No oriented model. Please orient model."), duration = dur);  return(NULL) }
    }
    if (input$selOrientationMethod=="none") {
      mesh.orient <- mesh.preorient
    }
    SavePath <- input$savePath
    if (!file.exists(SavePath)) { dir.create(file.path(SavePath),recursive=T) }
    FileName <- substr(input$fileMesh$name, 1, nchar(input$fileMesh$name)-4)
    FullFilePath <- paste(SavePath,FileName,"_oriented",sep="")
    vcgPlyWrite(mesh.orient, FullFilePath, binary = T)
    meshi <- paste(FullFilePath,".ply",sep="")
    mesho <- paste(FullFilePath,".ply",sep="")
    meshs <- paste(WD,"/support/_02_reorient_normals.mlx",sep="")
    meshm <- "vn vc"
    MeshFullCmd <- paste("meshlabserver -i", meshi, "-o", mesho, "-m", meshm, "-s", meshs); setProgress(0.3)
    shell(MeshFullCmd, wait=T)
    showNotification(paste("Oriented model saved"), duration = dur)
  })
  
  observeEvent(input$selOrientationMethod, {
    if (input$selOrientationMethod=="AOS_none")         { updateSelectInput(session, 'selSide', label=NULL, c("None"="")) }
    if (input$selOrientationMethod=="AOS_circle4DDL")   { updateSelectInput(session, 'selSide', label=NULL, c("Both" = "both", "Outer (Blue)" = "blue", "Inner (Red)" = "red")) }
    if (input$selOrientationMethod=="AOS_circle3DDL")   { updateSelectInput(session, 'selSide', label=NULL, c("Both" = "both", "Outer (Blue)" = "blue", "Inner (Red)" = "red")) }
    if (input$selOrientationMethod=="AOS_sections")     { updateSelectInput(session, 'selSide', label=NULL, c("Both" = "both", "Outer (Blue)" = "blue", "Inner (Red)" = "red")) }
    if (input$selOrientationMethod=="AOS_profile")      { updateSelectInput(session, 'selSide', label=NULL, c("Whole fragment" = "whole", "Outer (Blue)" = "blue", "Inner (Red)" = "red")) }
    if (input$selOrientationMethod=="AOS_polynomial")   { updateSelectInput(session, 'selSide', label=NULL, c("Both" = "both", "Outer (Blue)" = "blue", "Inner (Red)" = "red")) }
    if (input$selOrientationMethod=="AOS_rimbase_rim")  { updateSelectInput(session, 'selSide', label=NULL, c("Rim"="")) }
    if (input$selOrientationMethod=="AOS_rimbase_base") { updateSelectInput(session, 'selSide', label=NULL, c("Base"="")) }
  })
  
  observeEvent(c(input$sldNb_AOS_circle4DDL, input$sldNb_AOS_circle3DDL), {
    if (!exists("mesh.preorient.decim")) { return(NULL) }
    if (input$selOrientationMethod=="AOS_circle4DDL") { Nb <- input$sldNb_AOS_circle4DDL }
    if (input$selOrientationMethod=="AOS_circle3DDL") { Nb <- input$sldNb_AOS_circle3DDL }
    mesh <- mesh.trip$mesh
    S1 <- mesh.trip$S1
    S2 <- mesh.trip$S2
    S1 <- vcgIsolated(S1, silent=T)
    S2 <- vcgIsolated(S2, silent=T)
    plot3d.bl()
    plot3d(mesh.preorient.decim, col="lightgrey", asp="iso", alpha=0.6)
    lines3d(rbind( c(0,0,min(mesh.preorient$vb[3,])), c(0,0,max(mesh.preorient$vb[3,])) ), lwd=2, col="red")
    if (input$selSide=="red" || input$selSide=="both") {
      plot3d(S1, col="red", add=T)
      sect.hor(S1, nb=Nb, col="blue", size=4, plots=T)
    }
    if (input$selSide=="blue" || input$selSide=="both") {
      plot3d(S2, col="blue", add=T)
      sect.hor(S2, nb=Nb, col="red", size=4, plots=T)
    }
    widget <- rglwidget()
    output$plotOrientation3d <- renderRglwidget(widget)
    rgl.close()
  })
  
  observeEvent(input$sldParetoNb_AOS_sections, {
    if (!exists("mesh")) { showNotification(paste("No file loaded. Please load file."), duration = dur); return(NULL) }
    if (!exists("mesh.preorient")) { showNotification(paste("No pre-oriented model. Please pre-orient model."), duration = dur); return(NULL) }
    if (!exists("mesh.trip")) { showNotification(paste("No pre-oriented model. Please pre-orient model."), duration = dur);  return(NULL) }
    if (!exists("best.paretos")) { showNotification(paste("No Pareto result. Please orient model by Pareto."), duration = dur);  return(NULL) }
    Pareto.Nb <- input$sldParetoNb_AOS_sections
    Plo <- input$cbxPlots_AOS_sections
    best.paretos <- best.paretos
    mesh <- mesh.trip$mesh
    S1 <- mesh.trip$S1
    S2 <- mesh.trip$S2
    S1 <- vcgIsolated(S1, silent=T)
    S2 <- vcgIsolated(S2, silent=T)
    res <- best.paretos$value
    res.norm <- cbind((res[,1]-min(res[,1])) / (max(res[,1])-min(res[,1])),(res[,2]-min(res[,2])) / (max(res[,2])-min(res[,2])))
    ord <- order(res.norm[,1]+res.norm[,2])
    total <- order(res.norm[,1]+res.norm[,2])
    final <- best.paretos$par[total[1:Pareto.Nb],]
    if (length(final) == 2) { chosen.solution <- final }
    else { chosen.solution <- apply(final,2,mean) }
    if (input$selSide=="both") {
      mesh.pair <- list(S1=S1, S2=S2)
      best.par <- AOS_sections_pair.all(mesh.pair=mesh.pair, par=chosen.solution, nb=0, plots=Plo)
    }
    if (input$selSide=="red") {
      best.par <- AOS_sections.all(mesh=S1, par=chosen.solution, nb=0, plots=Plo)
    }
    if (input$selSide=="blue") {
      best.par <- AOS_sections.all(mesh=S2, par=chosen.solution, nb=0, plots=Plo)
    }
    best.par <- c(best.par[1:2],-best.par[3:4])
    mesh <- mesh.transform(mesh, best.par)
    mesh.orient <- mesh.transform(mesh.preorient, best.par)
    mesh.orient.decim <- decim(mesh.orient, 10000)
    plot3d.bl()
    plotOrientedMesh(mesh.orient.decim, method="3d", col="lightgrey")
    widget <- rglwidget()
    output$plotOrientation3d <- renderRglwidget(widget)
    rgl.close()
    output$plotOrientation <- renderPlot({
      layout(matrix(1:3,1,3))
      plotOrientedMesh(mesh.orient, method="rz")
      plot(best.paretos, xlab="F1", ylab="F2", main="Objective space")
      points(best.paretos$value[ord[1:Pareto.Nb],1],best.paretos$value[ord[1:Pareto.Nb],2], col="blue", pch=19)
      plot(best.paretos$par, xlab="phi", ylab="theta", main="Parameter space", xlim=c(-pi/2,pi/2), ylim=c(-pi/2,pi/2), asp=1)
      points(best.paretos$par[ord[1:Pareto.Nb],1],best.paretos$par[ord[1:Pareto.Nb],2], col="blue", pch=19)
      abline(v=0, lty=2); abline(h=0, lty=2)
    })
    mesh <<- mesh
    mesh.orient <<- mesh.orient
  })
  
  observeEvent(input$sldNb_AOS_profile, {
    if (!exists("mesh.preorient.decim")) { return(NULL) }
    Nb <- input$sldNb_AOS_profile
    mesh <- mesh.trip$mesh
    S1 <- mesh.trip$S1
    S2 <- mesh.trip$S2
    S1 <- vcgIsolated(S1, silent=T)
    S2 <- vcgIsolated(S2, silent=T)
    plot3d.bl()
    plot3d(mesh.preorient.decim, col="lightgrey", asp="iso", alpha=0.6)
    lines3d(rbind( c(0,0,min(mesh.preorient$vb[3,])), c(0,0,max(mesh.preorient$vb[3,])) ), lwd=2, col="red")
    if (input$selSide=="whole") {
      sect.ver(mesh.preorient.decim, nb=Nb, col="blue", size=4, plots=T)
    }
    if (input$selSide=="red") {
      plot3d(S1, col="red", add=T)
      sect.ver(S1, nb=Nb, col="blue", size=4, plots=T)    
    }
    if (input$selSide=="blue") {
      plot3d(S2, col="blue", add=T)
      sect.ver(S2, nb=Nb, col="red", size=4, plots=T)
    }
    widget <- rglwidget()
    output$plotOrientation3d <- renderRglwidget(widget)
    rgl.close()
  })
  
  observeEvent(c(input$sldNb_AOS_rimbase_rim, input$sldTresh_AOS_rimbase_rim, input$sldNb_AOS_rimbase_base, input$sldTresh_AOS_rimbase_base), {
    if (!exists("mesh.preorient.decim")) { return(NULL) }
    if (input$selOrientationMethod=="AOS_rimbase_rim") {
      Tresh <- input$sldTresh_AOS_rimbase_rim
      Nb <- input$sldNb_AOS_rimbase_rim
      v1 <- c(0,0,(max(mesh.preorient.decim$vb[3,])-Tresh))
      v2 <- c(0,1,(max(mesh.preorient.decim$vb[3,])-Tresh))
      v3 <- c(1,0,(max(mesh.preorient.decim$vb[3,])-Tresh))
      mesh.preorient.decim.cut <- cutMeshPlane(mesh.preorient.decim, v1, v2, v3, keep.upper=F)
    }
    if (input$selOrientationMethod=="AOS_rimbase_base") {
      Tresh <- input$sldTresh_AOS_rimbase_base
      Nb <- input$sldNb_AOS_rimbase_base
      v1 <- c(0,0,(min(mesh.preorient.decim$vb[3,])+Tresh))
      v2 <- c(0,1,(min(mesh.preorient.decim$vb[3,])+Tresh))
      v3 <- c(1,0,(min(mesh.preorient.decim$vb[3,])+Tresh))
      mesh.preorient.decim.cut <- cutMeshPlane(mesh.preorient.decim, v1, v2, v3, keep.upper=T)
    }
    plot3d.bl()
    plotOrientedMesh(mesh.preorient.decim, method="3d", col="lightgrey", alpha=0.6)
    lines3d(rbind( c(0,0,min(mesh.preorient.decim$vb[3,])), c(0,0,max(mesh.preorient.decim$vb[3,])) ), lwd=2, col="red")
    plot3d(mesh.preorient.decim.cut, color="red", add=T)
    sect.ver(mesh.preorient.decim, nb=Nb, col="blue", size=4, plots=T)
    widget <- rglwidget()
    output$plotOrientation3d <- renderRglwidget(widget)
    rgl.close()
  })
  
  dataInputPreorientAutomatic <- reactive({
    withProgress(message = 'Automatic pre-orientation', value = 0.1, {
      # ### DEBUG ###
      # Decim = 5000
      # SmIter <- 5
      # PVal <- 0.95
      # Plo <- F
      # ### DEBUG ###
      Decim <- input$sldDecim_PA
      SmIter <- input$sldSmIter_PA
      PVal <- input$sldPVal_PA
      Plo <- input$cbxPlots_PA
      mesh.orig <- mesh
      mesh <- decim(mesh, Decim)
      mesh <- mesh.clean(mesh, it=3)
      mesh <- vcgSmooth(mesh, type="laplace", iteration=SmIter)
      while (sum(apply(mesh$normals[1:3,],2,sum)==0)>1) { mesh <- mesh.clean(mesh, it=1) }
      setProgress(0.2, "Extraction of unwanted features")
      res <- distNpval.parallel(mesh)
      dist <- res[seq(1,length(res),by=2)]        #  not included in the calculation
      p.val <- res[seq(2,length(res),by=2)]
      mesh <- rmVertex(mesh, index=which(p.val>PVal), keep=T)
      while (sum(apply(mesh$normals[1:3,],2,sum)==0)>1) { mesh <- mesh.clean(mesh, it=1) }
      mesh <- mesh.clean.borders(mesh, it=2)
      mesh <- vcgIsolated(mesh, facenum=5, silent=T)
      setProgress(0.5, "Rotational axis estimation (this may take a while)")
      par <- c(0.05,0.05)
      set.seed(1234)
      best.par <- optimx(par=par, mesh=mesh, plots=Plo, infos=Infos, fn=AOS_normals.f1, method="Nelder-Mead")
      best.par <- AOS_normals.all(mesh=mesh, par=as.numeric(best.par[1:4]), plots=Plo)
      setProgress(0.9, "Rotational axis found")
      mesh <- mesh.transform(mesh, par=-best.par, method="tpab", plots=Plo)
      mesh.preorient <- mesh.transform(mesh.orig, par=-best.par, method="tpab", plots=Plo)
      mesh.trip <- separate.byAxis(mesh, plots=Plo)
      mesh.preorient <<- mesh.preorient
      mesh.trip <<- mesh.trip
    })
  })
  
  dataInputPreorientManual <- reactive({
    withProgress(message = 'Manual pre-orientation', value = 0.1, {
      Decim <- input$sldDecim_PM
      SmIter <- input$sldSmIter_PM
      Perc <- input$sldPerc_PM
      Plo <- input$cbxPlots_PM
      setProgress(0.1, "Please define three points on the model lying in the same horizontal plane.")
      mesh <- mesh.preorient <- AOS_manual(mesh)
      setProgress(0.3, "Rotational axis found")
      mesh <- decim(mesh, Decim)
      mesh <- mesh.clean(mesh, it=1)
      mesh <- vcgSmooth(mesh, type="laplace", iteration=SmIter)
      while (sum(apply(mesh$normals[1:3,],2,sum)==0)>1) { mesh <- mesh.clean(mesh, it=1) }
      setProgress(0.5, "Extraction of unwanted features")
      d <- distNormals2axis.parallel(mesh)
      mesh <- rmVertex(mesh, index=which(d<quantile(d,Perc)), keep=T)
      mesh.trip <- separate.byAxis(mesh, plots=Plo)
      mesh.preorient <<- mesh.preorient
      mesh.trip <<- mesh.trip
    })
  })
  
  dataInput_AOS_none <- reactive({
    mesh <<- mesh.orient <<- mesh.preorient
  })
  
  dataInput_AOS_circle4DDL <- reactive({
    withProgress(message = '4DDL orientation', value = 0.1, {
      Nb <- input$sldNb_AOS_circle4DDL
      Plo <- input$cbxPlots_AOS_circle4DDL
      mesh <- mesh.trip$mesh
      S1 <- mesh.trip$S1
      S2 <- mesh.trip$S2
      S1 <- vcgIsolated(S1, silent=T)
      S2 <- vcgIsolated(S2, silent=T)
      setProgress(0.4, "Rotational axis estimation (this may take a while)")
      if (input$selSide=="both") {
        mesh.pair <- list(S1=S1, S2=S2)
        par <- c(0,0,0,0)
        set.seed(1234)
        best.par <- optimx(par=par, mesh.pair=mesh.pair, nb=Nb, plots=Plo, infos=Infos, fn=AOS_circle4DDL_pair.f1, method="Nelder-Mead", itnmax=Itnmax)
        best.par <- AOS_circle4DDL_pair.all(mesh.pair=mesh.pair, par=as.numeric(best.par[1:4]), plots=Plo)
      }
      if (input$selSide=="red") {
        par <- c(0,0,0,0)
        set.seed(1234)
        best.par <- optimx(par=par, mesh=S1, nb=Nb, plots=Plo, infos=Infos, fn=AOS_circle4DDL.f1, method="Nelder-Mead", itnmax=Itnmax)
        best.par <- AOS_circle4DDL.all(mesh=S1, par=as.numeric(best.par[1:4]), plots=Plo)
      }
      if (input$selSide=="blue") {
        par <- c(0,0,0,0)
        set.seed(1234)
        best.par <- optimx(par=par, mesh=S2, nb=Nb, plots=Plo, infos=Infos, fn=AOS_circle4DDL.f1, method="Nelder-Mead", itnmax=Itnmax)
        best.par <- AOS_circle4DDL.all(mesh=S2, par=as.numeric(best.par[1:4]), plots=Plo)
      }
      setProgress(0.8, "Rotational axis found")
      mesh <- mesh.transform(mesh, best.par)
      mesh.orient <- mesh.transform(mesh.preorient, best.par)
      mesh <<- mesh
      mesh.orient <<- mesh.orient
      setProgress(1, "Rotational axis found")
    })
  })
  
  dataInput_AOS_circle3DDL <- reactive({
    withProgress(message = '3DDL orientation', value = 0.1, {
      Nb <- input$sldNb_AOS_circle3DDL
      Plo <- input$cbxPlots_AOS_circle3DDL
      mesh <- mesh.trip$mesh
      S1 <- mesh.trip$S1
      S2 <- mesh.trip$S2
      S1 <- vcgIsolated(S1, silent=T)
      S2 <- vcgIsolated(S2, silent=T)
      setProgress(0.4, "Rotational axis estimation (this may take a while)")
      if (input$selSide=="both") {
        mesh.pair <- list(S1=S1, S2=S2)
        par <- c(0.1,0,0)
        set.seed(1234)
        best.par=optimx(par=par, mesh.pair=mesh.pair, nb=Nb, plots=Plo, infos=Infos, fn=AOS_circle3DDL_pair.f1, method="Nelder-Mead", itnmax=Itnmax)
        best.par=AOS_circle3DDL_pair.all(mesh.pair=mesh.pair, par=as.numeric(best.par[1:3]), plots=Plo)
      }
      if (input$selSide=="red") {
        par <- c(0.1,0,0)
        set.seed(1234)
        best.par <- optimx(par=par, mesh=S1, nb=Nb, plots=Plo, infos=Infos, fn=AOS_circle3DDL.f1, method="Nelder-Mead", itnmax=Itnmax)
        best.par <- AOS_circle3DDL.all(mesh=S1, par=as.numeric(best.par[1:3]), plots=Plo)
      }
      if (input$selSide=="blue") {
        par <- c(0.1,0,0)
        set.seed(1234)
        best.par <- optimx(par=par, mesh=S2, nb=Nb, plots=Plo, infos=Infos, fn=AOS_circle3DDL.f1, method="Nelder-Mead", itnmax=Itnmax)
        best.par <- AOS_circle3DDL.all(mesh=S2, par=as.numeric(best.par[1:3]), plots=Plo)
      }
      setProgress(0.8, "Rotational axis found")
      best.par <- c(0,best.par[1:3])
      mesh <- mesh.transform(mesh, best.par)
      mesh.orient <- mesh.transform(mesh.preorient, best.par)
      mesh <<- mesh
      mesh.orient <<- mesh.orient
      setProgress(1, "Rotational axis found")
    })
  })
  
  dataInput_AOS_sections <- reactive({
    withProgress(message = 'Section orientation', value = 0.1, {
      Nb <- input$sldNb_AOS_sections
      Plo <- input$cbxPlots_AOS_sections
      Low.lim <- rep(input$sldLimits_AOS_sections[1]*pi/180,2)
      Up.lim <- rep(input$sldLimits_AOS_sections[2]*pi/180,2)
      Gener <- input$sldGener_AOS_sections
      Popsize <- input$sldPopsize_AOS_sections
      Pareto.Nb <- input$sldParetoNb_AOS_sections
      mesh <- mesh.trip$mesh
      S1 <- mesh.trip$S1
      S2 <- mesh.trip$S2
      S1 <- vcgIsolated(S1, silent=T)
      S2 <- vcgIsolated(S2, silent=T)
      setProgress(0.4, "Rotational axis estimation (this may take a while)")
      if (input$selSide=="both") {
        mesh.pair <- list(S1=S1, S2=S2)
        par <- c(0,0)
        set.seed(1234)
        best.paretos <- nsga2(AOS_sections_pair.f12, 2, 2, mesh.pair=mesh.pair, nb=Nb, plots=Plo, infos=Infos,
                              lower.bounds=Low.lim, upper.bounds=Up.lim, generations=Gener, popsize=Popsize)
        best.par <- AOS_sections_pair.all(mesh.pair=mesh.pair, par=select_best_pareto(best.paretos, nb=Pareto.Nb, plots=Plo), nb=0, plots=Plo)
      }
      if (input$selSide=="red") {
        par <- c(0,0)
        set.seed(1234)
        best.paretos <- nsga2(AOS_sections.f12, 2, 2, mesh=S1, nb=Nb, plots=Plo, infos=Infos,
                              lower.bounds=Low.lim, upper.bounds=Up.lim, generations=Gener, popsize=Popsize)
        best.par <- AOS_sections.all(mesh=S1, par=select_best_pareto(best.paretos, nb=Pareto.Nb, plots=Plo), nb=0, plots=Plo)
      }
      if (input$selSide=="blue") {
        par <- c(0,0)
        set.seed(1234)
        best.paretos <- nsga2(AOS_sections.f12, 2, 2, mesh=S2, nb=Nb, plots=Plo, infos=Infos,
                              lower.bounds=Low.lim, upper.bounds=Up.lim, generations=Gener, popsize=Popsize)
        best.par <- AOS_sections.all(mesh=S2, par=select_best_pareto(best.paretos, nb=Pareto.Nb, plots=Plo), nb=0, plots=Plo)
      }
      setProgress(0.8, "Rotational axis found")
      res <- best.paretos$value
      res.norm <- cbind((res[,1]-min(res[,1])) / (max(res[,1])-min(res[,1])),(res[,2]-min(res[,2])) / (max(res[,2])-min(res[,2])))
      ord <- order(res.norm[,1]+res.norm[,2])
      total <- order(res.norm[,1]+res.norm[,2])
      best.par <- c(best.par[1:2],-best.par[3:4])
      mesh <- mesh.transform(mesh, best.par)
      mesh.orient <- mesh.transform(mesh.preorient, best.par)
      mesh.orient.decim <- decim(mesh.orient, 10000)
      plot3d.bl()
      plotOrientedMesh(mesh.orient.decim, method="3d", col="lightgrey")
      widget <- rglwidget()
      output$plotOrientation3d <- renderRglwidget(widget)
      rgl.close()
      output$plotOrientation <- renderPlot({
        layout(matrix(1:3,1,3))
        plotOrientedMesh(mesh.orient, method="rz")
        plot(best.paretos, xlab="F1", ylab="F2", main="Objective space")
        points(best.paretos$value[ord[1:Pareto.Nb],1],best.paretos$value[ord[1:Pareto.Nb],2], col="blue", pch=19)
        plot(best.paretos$par, xlab="phi", ylab="theta", main="Parameter space", xlim=c(-pi/2,pi/2), ylim=c(-pi/2,pi/2), asp=1)
        points(best.paretos$par[ord[1:Pareto.Nb],1],best.paretos$par[ord[1:Pareto.Nb],2], col="blue", pch=19)
        abline(v=0, lty=2); abline(h=0, lty=2)
      })
      mesh <<- mesh
      mesh.orient <<- mesh.orient
      best.paretos <<- best.paretos
      setProgress(1, "Rotational axis found")
    })
  })
  
  dataInput_AOS_profile <- reactive({
    withProgress(message = 'Profile orientation', value = 0.1, {
      Nb <- input$sldNb_AOS_profile
      Plo <- input$cbxPlots_AOS_profile
      mesh <- mesh.trip$mesh
      S1 <- mesh.trip$S1
      S2 <- mesh.trip$S2
      S1 <- vcgIsolated(S1, silent=T)
      S2 <- vcgIsolated(S2, silent=T)
      setProgress(0.4, "Rotational axis estimation (this may take a while)")
      if (input$selSide=="whole") {
        par <- c(0,0,0,0)
        set.seed(1234)
        best.par <- optimx(par=par, mesh=mesh, nb=Nb, by.deg=1, ref="longest", plots=Plo, infos=Infos, fn=AOS_profile.f1, method="Nelder-Mead", itnmax=Itnmax)
        best.par <- AOS_profile.all(mesh=mesh, nb=0, by.deg=1, par=as.numeric(best.par[1:4]), plots=Plo)
      }
      if (input$selSide=="red") {
        par <- c(0,0,0,0)
        set.seed(1234)
        best.par <- optimx(par=par, mesh=S1, nb=Nb, by.deg=1, ref="longest", plots=Plo, infos=Infos, fn=AOS_profile.f1, method="Nelder-Mead", itnmax=Itnmax)
        best.par <- AOS_profile.all(mesh=S1, nb=0, par=as.numeric(best.par[1:4]), plots=Plo)
      }
      if (input$selSide=="blue") {
        par <- c(0,0,0,0)
        set.seed(1234)
        best.par <- optimx(par=par, mesh=S2, nb=Nb, by.deg=1, ref="longest", plots=Plo, infos=Infos, fn=AOS_profile.f1, method="Nelder-Mead", itnmax=Itnmax)
        best.par <- AOS_profile.all(mesh=S2, nb=0, par=as.numeric(best.par[1:4]), plots=Plo)
      }
      setProgress(0.8, "Rotational axis found")
      mesh <- mesh.transform(mesh, best.par)
      mesh.orient <- mesh.transform(mesh.preorient, best.par)
      mesh <<- mesh
      mesh.orient <<- mesh.orient
      setProgress(1, "Rotational axis found")
    })
  })
  
  dataInput_AOS_polynomial <- reactive({
    withProgress(message = 'Polynomial orientation', value = 0.1, {
      PO <- input$sldPO_AOS_polynomial
      Plo <- input$cbxPlots_AOS_polynomial
      mesh <- mesh.trip$mesh
      S1 <- mesh.trip$S1
      S2 <- mesh.trip$S2
      setProgress(0.4, "Rotational axis estimation (this may take a while)")
      if (input$selSide=="both") {
        mesh.pair <- list(S1=S1, S2=S2)
        par <- c(0,0,0,0)
        set.seed(1234)
        best.par <- optimx(par=par, mesh.pair=mesh.pair, po=PO, plots=Plo, infos=Infos, fn=AOS_polynomial_pair.f1, method="Nelder-Mead", itnmax=Itnmax)
        best.par <- AOS_polynomial_pair.all(mesh.pair=mesh.pair, par=as.numeric(best.par[1:4]), plots=Plo)
      }
      if (input$selSide=="red") {
        par <- c(0,0,0,0)
        set.seed(1234)
        best.par <- optimx(par=par, mesh=S1, po=PO, plots=Plo, infos=Infos, fn=AOS_polynomial.f1, method="Nelder-Mead", itnmax=Itnmax)
        best.par <- AOS_polynomial.all(mesh=S1, par=as.numeric(best.par[1:4]), plots=Plo)
      }
      if (input$selSide=="blue") {
        par <- c(0,0,0,0)
        set.seed(1234)
        best.par <- optimx(par=par, mesh=S2, po=PO, plots=Plo, infos=Infos, fn=AOS_polynomial.f1, method="Nelder-Mead", itnmax=Itnmax)
        best.par <- AOS_polynomial.all(mesh=S2, par=as.numeric(best.par[1:4]), plots=Plo)
      }
      setProgress(0.8, "Rotational axis found")
      mesh <- mesh.transform(mesh, best.par)
      mesh.orient <- mesh.transform(mesh.preorient, best.par)
      mesh <<- mesh
      mesh.orient <<- mesh.orient
      setProgress(1, "Rotational axis found")
    })
  })
  
  dataInput_AOS_rimbase_rim <- reactive({
    withProgress(message = 'Rim orientation', value = 0.1, {
      Nb <- input$sldNb_AOS_rimbase_rim
      Tresh <- input$sldTresh_AOS_rimbase_rim
      Plo <- input$cbxPlots_AOS_rimbase_rim
      mesh <- mesh.preorient
      setProgress(0.4, "Rotational axis estimation (this may take a while)")
      par <- c(0,0,0,0)
      set.seed(1234)
      best.par <- optimx(par=par, mesh=mesh, tresh=Tresh, nb=Nb, by.deg=1, part="rim", three.p=F, plots=Plo, infos=Infos, fn=AOS_rimbase.f1, method="Nelder-Mead", itnmax=Itnmax)
      best.par <- AOS_rimbase.all(mesh=mesh, tresh=Tresh, nb=Nb, by.deg=1, part="rim", three.p=F, plots=Plo, infos=F, par=as.numeric(best.par[1:4]))
      setProgress(0.8, "Rotational axis found")
      mesh <- mesh.transform(mesh, best.par)
      mesh.orient <- mesh.transform(mesh.preorient, best.par)
      mesh <<- mesh
      mesh.orient <<- mesh.orient
      setProgress(1, "Rotational axis found")
    })
  })
  
  dataInput_AOS_rimbase_base <- reactive({
    withProgress(message = 'Base orientation', value = 0.1, {
      Nb <- input$sldNb_AOS_rimbase_base
      Tresh <- input$sldTresh_AOS_rimbase_base
      Plo <- input$cbxPlots_AOS_rimbase_base
      mesh <- mesh.preorient
      setProgress(0.4, "Rotational axis estimation (this may take a while)")
      par <- c(0,0,0,0)
      set.seed(1234)
      best.par <- optimx(par=par, mesh=mesh, tresh=1, nb=Nb, by.deg=1, part="base", three.p=F, plots=Plo, infos=Infos, fn=AOS_rimbase.f1, method="Nelder-Mead", itnmax=Itnmax)
      best.par <- AOS_rimbase.all(mesh=mesh, tresh=1, nb=Nb, by.deg=1, part="base", three.p=F, plots=Plo, infos=F, par=as.numeric(best.par[1:4]))
      setProgress(0.8, "Rotational axis found")
      mesh <- mesh.transform(mesh, best.par)
      mesh.orient <- mesh.transform(mesh.preorient, best.par)
      mesh <<- mesh
      mesh.orient <<- mesh.orient
      setProgress(1, "Rotational axis found")
    })
  })
  
  
  # # DACORD - ILLUSTRATION

  dataInputOrientedMesh <- reactive({
    inFile <- input$fileOrientedMesh
    if (is.null(inFile)) return(NULL)
    FileName <<- inFile$name
    FileExt <- substr(FileName, nchar(FileName)-2, nchar(FileName))
    if (FileExt!="ply") { showNotification(paste("Please load model in PLY format"), duration = dur); return(NULL) }
    file <- inFile$datapath
    filen <- gsub("[\\]", "/", file)
    filen <- paste(substr(filen, 1, nchar(filen)-1), inFile$name, sep="")
    file.rename(file, filen)
    mesh <- vcgImport(filen, clean=F, readcolor=T)
    if (is.null(mesh$material$color)) { mesh$material$color <- matrix(grey.colors(3)[2],3,dim(mesh$it)[2]) }
    mesh.oriented <<- mesh
    return(mesh.oriented)
  })
  
  observeEvent(input$fileOrientedMesh, {
    mesh.drawing <<- dataInputOrientedMesh()
    dataInputPreparation()
    showNotification(paste("Oriented model loaded"), duration = dur)
  })
  
  dataInputPreparation <- reactive({
    withProgress(message = 'Data preparation', value=0.1, {
      FullFilePath <- tempfile("")
      vcgPlyWrite(mesh.drawing, FullFilePath, binary = T)
      FullFilePath <- paste(gsub("[\\]", "/", FullFilePath),".ply",sep="")
      meshi <- FullFilePath
      mesho <- FullFilePath
      meshs <- paste(WD,"/support/_02_reorient_normals.mlx",sep="")
      meshm <- "vn vc"
      MeshFullCmd <- paste("meshlabserver -i", meshi, "-o", mesho, "-m", meshm, "-s", meshs)
      shell(MeshFullCmd, wait=T)
      mesh.drawing <- vcgImport(FullFilePath, clean=F, readcolor=T)
      setProgress(0.2)
      mesh.drawing <- mesh.align2front(mesh.drawing)
      mesh.drawing <- mesh.align2xy(mesh.drawing)
      mesh.drawing.grey <- mesh2grey(mesh.drawing)
      mesh.drawing.decim.mini <- decim(mesh.drawing, 2000)
      mesh.drawing.decim.mini <- vcgUpdateNormals(mesh.drawing.decim.mini)
      mesh.drawing.centroid <- mesh.centroid(mesh.drawing)
      cs.val <- circular_section(centre = c(0,0), xy = t(mesh.drawing.decim.mini$vb[1:3,]))*180/pi
      output$plotDrawing <- renderPlot({
        plotOrientedMesh(mesh.drawing.decim.mini, method="rz")
      })
      mesh.drawing.decim.maxi <<- vcgQEdecim(mesh.drawing, tarface = 50000, normcheck=F, silent=T)
      mesh.drawing.decim.maxi <<- vcgUpdateNormals(mesh.drawing.decim.maxi)
      setProgress(0.4)
      FullFilePath <- tempfile("")
      vcgPlyWrite(mesh.drawing, FullFilePath, binary = T)
      FullFilePath <- paste(gsub("[\\]", "/", FullFilePath),".ply",sep="")
      meshi <- FullFilePath
      mesho <- FullFilePath
      meshs <- paste(WD,"/support/_03_ambient_occlusion.mlx",sep="")
      meshm <- "vn vc vq"
      MeshFullCmd <- paste("meshlabserver -i", meshi, "-o", mesho, "-m", meshm, "-s", meshs)
      shell(MeshFullCmd, wait=T)
      meshao <- vcgImport(FullFilePath, clean=F, readcolor=T)
      mesh.drawing$ao <- meshao$quality
      setProgress(0.4)
      profil <- profile2d(mesh.drawing.decim.mini, method = "envelop", plots=F)
      vcgPlyWrite(mesh.drawing.decim.maxi, filename=paste(WD,"/temp/temp", sep=""), binary=F)
      write.table(profil$profil, paste(WD,"/temp/temp_profil.txt",sep=""), sep=";", row.names=F, col.names=F)
      shell(paste("Rscript ",WD,"/_support_ME.R",sep=""), wait=T)
      mesh.drawing <<- mesh.drawing
      mesh.drawing.grey <<- mesh.drawing.grey
      mesh.drawing.decim.mini <<- mesh.drawing.decim.mini
      mesh.drawing.centroid <<- mesh.drawing.centroid
      cs.val <<- cs.val
      setProgress(1)
    })
  })

  observeEvent(input$btnGetProfile, {
    if (!exists("mesh.drawing")) { showNotification(paste("No oriented model loaded. Please load the model."), duration = dur);  return(NULL) }
    if (input$selProfileMethod=="") {  showNotification(paste("No profile method selected. Please select the method."), duration = dur);  return(NULL)  }
    withProgress(message = 'Profile extraction', value=0.1, {
      if (input$selProfileMethod=="Profile_envelop")    { profil <<- profile2d(mesh.drawing.decim.mini, method = "envelop", plots=F) ; mai <- c("Whole envelop profile") }
      if (input$selProfileMethod=="Profile_middle")     { profil <<- profile2d(mesh.drawing, method = "middle", plots=F) ; mai <- c("In the middle of the fragment profile") }
      if (input$selProfileMethod=="Profile_longest")    { profil <<- profile2d(mesh.drawing, method = "longest", longest.by.deg=2, plots=F) ; mai <- c("Longest-preserved profile") }
      if (input$selProfileMethod=="Profile_arbitrary")  { profil <<- profile2d(mesh.drawing, method = "arbitrary", plots=F) ; mai <- c("Arbitrary selected profile") }
      setProgress(0.8)
      profil <<- switch.profile(profil)
      pottery <<- profil2pottery(profil, x.off = 5, off.side=10)
      CSPI <<- cspi(cs.val=cs.val, cs.val.x=7, cs.val.y=pottery$bpt-5, cs.radius=30, cs.y=pottery$bpt-5)
      ae.scale <<- archeo.scale(scale.y=pottery$bpt-10)
      RepairProfile_F()
      setProgress(1)
    })
    ptt <<- matrix(numeric(0),0,2)
    REC <<- array(numeric(0), dim=c(40,2,0))
  })
  
  observeEvent(input$btnRepairProfile, {
    if (!exists("mesh.drawing")) { showNotification(paste("No oriented model loaded. Please load the model."), duration = dur);  return(NULL) }
    Errase_REC()
    RepairProfile_F()
  })
  
  RepairProfile_F <- function() {
    output$plotDrawing <- renderPlot({
      plot(profil$profil, type="n", xlim=c(min(profil$profil[,1]),0), ylim=c(min(profil$profil[,2])-5,max(profil$profil[,2])+5), asp=1, axes=T, xlab="", ylab="")
      lines(rbind(c(0,0),c(0,min(profil$profil[,2]))), col="red", lwd=2)
      polygon(profil$profil, col="grey", border="darkgrey")
      lines(profil$extern, col="blue")
      lines(profil$intern, col="red")
      pt <- c(input$plotDrawing_dblclick$x, input$plotDrawing_dblclick$y)
      if (!is.null(pt)) { profil <<- profile2d.repair(profil, pt=pt, plots=T) }
      ptt <<- rbind(ptt,c(input$plotDrawing_click$x, input$plotDrawing_click$y))
      if (dim(ptt)[1]>0 & (dim(ptt)[1])%%2==0) { updateREC() }
      pottery <<- profil2pottery(profil, x.off = 5, off.side=10)
    })
    updateSliderInput(session,"sldTrans_Drawing_linear", value=0, min=-round(pottery$rupt[1]), max=round(pottery$rupt[1]), step=1)
    updateSliderInput(session,"sldTrans_Drawing_photo", value=0, min=-round(pottery$rupt[1]), max=round(pottery$rupt[1]), step=1)
    updateSliderInput(session,"sldTrans_Drawing_DL", value=0, min=-round(pottery$rupt[1]), max=round(pottery$rupt[1]), step=1)
    updateSliderInput(session,"sldTrans_Drawing_AO", value=0, min=-round(pottery$rupt[1]), max=round(pottery$rupt[1]), step=1)
    updateSliderInput(session,"sldTrans_Drawing_symmetry", value=0, min=-round(pottery$rupt[1]), max=round(pottery$rupt[1]), step=1)
  }
  
  Errase_REC <- function(){
    if (length(ptt)==4) {
      ptt <<- matrix(numeric(0),0,2)
      REC <<- array(numeric(0), dim=c(40,2,0))
    }
    if (length(ptt)>4) {
      if (dim(ptt)[1]%%2==0) { ptt <<- ptt[-c((dim(ptt)[1]-1):dim(ptt)[1]),] }
      else { ptt <<- ptt[-c((dim(ptt)[1]-2):dim(ptt)[1]),] }
    }
  }
  
  Drawing_linear <- function(){
    Add.rim <- input$cbxAddRim_Drawing_linear
    Add.base <- input$cbxAddBase_Drawing_linear
    Add.rec <- input$cbxAddRec_Drawing_linear
    Add.lines <- input$cbxAddLines_Drawing_linear
    Add.diam <- input$cbxAddDiameter_Drawing_linear
    Add.scale <- input$cbxAddScale_Drawing_linear
    Add.cspi <- input$cbxAddCSPI_Drawing_linear
    Vol.perc <- input$sldVolume_Drawing_photo
    Cs.val <- cs.val
    pts <<- rbind(pts,c(input$plotDrawing_click$x, input$plotDrawing_click$y))
    ruler <<- rbind(ruler, c(input$plotDrawing_dblclick$x, input$plotDrawing_dblclick$y))
    if (dim(pts)[1]>0) { updateLIN() }
    # Add.rim=T;Add.base=F; Add.rec=T; Add.lines=T; Add.diam=T; Add.scale=T; Add.cspi=T; Vol.perc=c(0,100)
    reconst.2d(pottery=pottery, ae.scale=ae.scale, LIN=LIN, REC=REC, CSPI=CSPI, main="",
               profil.fill.col="darkgrey", lines.col="black", lines.lwd=1,
               add.rim=Add.rim, add.base=Add.base,
               add.lines=Add.lines, l.lines.col="grey", l.lines.lwd=1,
               add.diam=Add.diam,
               add.scale=Add.scale, scale.fill.col="black",
               add.cspi=Add.cspi, cspi.fill.col="black", cspi.bord.col="black",
               add.rec=Add.rec)
    if (input$cbxAddVolume_Drawing_linear==T) { if(Vol.perc[1]!=Vol.perc[2]) { volume(Vol.perc) } }
    if (input$cbxAddRuler_Drawing_linear==T) { if (dim(ruler)[1]>0) { rulerF() } }
  }
  
  Drawing_linear_F <- function(){
    output$plotDrawing <- renderPlot({
      Drawing_linear()
    })
    output$plotDrawingHist<- renderPlot({})
  }
  
  Drawing_photo <- function(){
    Add.rim <- input$cbxAddRim_Drawing_photo
    Add.base <- input$cbxAddBase_Drawing_photo
    Add.rec <- input$cbxAddRec_Drawing_photo
    Add.lines <- input$cbxAddLines_Drawing_photo
    Add.diam <- input$cbxAddDiameter_Drawing_photo
    Add.scale <- input$cbxAddScale_Drawing_photo
    Add.cspi <- input$cbxAddCSPI_Drawing_photo
    Cs.val <- cs.val
    Col <- input$selColour_Drawing_photo
    Add.mask <- input$cbxAddMask_Drawing_photo
    Add.light <- input$cbxAddLight_Drawing_photo
    Lat <<- input$sldLightLat_Drawing_photo
    Lon <<- input$sldLightLon_Drawing_photo
    Dis <<- input$sldLightDis_Drawing_photo*10
    Rot <- -deg2rad(input$sldRot_Drawing_photo)
    Trans <- input$sldTrans_Drawing_photo
    
    if (Rot==Rot.count3 & Trans==Trans.count3 & Lat==Lat.count3 & Lon==Lon.count3 & Dis==Dis.count3 &
        Add.light==Add.light.count3 & Col==Col.count3) { }
    else {
      if (Col=="selColour_Drawing_photo_color") {
        mesh.drawing.positioned <<- mesh.transform4drawing(mesh.drawing, rho=Rot, x=Trans)
      }
      if (Col=="selColour_Drawing_photo_greyscale") {
        mesh.drawing.positioned <<- mesh.transform4drawing(mesh.drawing.grey, rho=Rot, x=Trans)
      }
      while (length(rgl.dev.list()) > 0) { rgl.close() }
      save <- par3d(skipRedraw=T)
      par3d(windowRect=c(960-3000,30-3000,1920-3000,1040-3000))
      view3d(theta = 0, phi = -90, fov=0, interactive=F)
      clear3d(type="light")
      light3d(ambient="white", diffuse="white", specular="black")
      plot3d(mesh.drawing.positioned, meshColor="legacy", add=T)
      ### DEBUG
      if (Add.light==T) {
        light <- sph2cart(c(deg2rad(-90+Lon),deg2rad(Lat),Dis))+mesh.drawing.centroid
        # points3d(rbind(light), col="orange", size=10); lines3d(rbind(mesh.drawing.centroid,light))
        light3d(x = rbind(light), ambient="black", diffuse = "white",specular = "black", viewpoint.rel = F) 
      }
      par3d(save)
      rgl.snapshot(PlotPath)
      rgl.close()
      xyz <- t(mesh.drawing.positioned$vb[1:3,])
      x <- xyz[,1]
      y <- xyz[,3]
      recta <- c(min(x),min(y),max(x),max(y))
      Img <- readPNG(PlotPath)
      Img <- cut.img(Img, plots=F)
      Img <- as.raster(Img)
      Rot.count3 <<- Rot
      Trans.count3 <<- Trans
      Lat.count3 <<- Lat
      Lon.count3 <<- Lon
      Dis.count3 <<- Dis
      Add.light.count3 <<- Add.light
      Col.count3 <<- Col
      Img <<- Img
      recta <<- recta
    }
    plot(0, type="n", xlim=pottery$x.lim, ylim=pottery$y.lim, asp=1, axes=F, xlab="", ylab="")
    rasterImage(Img, xleft=recta[1], ybottom=recta[2], xright=recta[3], ytop=recta[4])
    pts <<- rbind(pts,c(input$plotDrawing_click$x, input$plotDrawing_click$y))
    if (dim(pts)[1]>0) { updateLIN() }
    if (Add.mask==T) { polygon(pottery$mask, col="white", border="white") }
    reconst.2d(pottery=pottery, ae.scale=ae.scale, LIN=LIN, REC=REC, CSPI=CSPI, main="", add=T,
               profil.fill.col="darkgrey", lines.col="black", lines.lwd=1,
               add.rim=Add.rim, add.base=Add.base,
               add.lines=Add.lines, l.lines.col="grey", l.lines.lwd=1,
               add.diam=Add.diam,
               add.scale=Add.scale, scale.fill.col="black",
               add.cspi=Add.cspi, cspi.fill.col="black", cspi.bord.col="black",
               add.rec=Add.rec)
  }
  
  Drawing_photo_F <- function(){
    output$plotDrawing <- renderPlot({
      Drawing_photo()
    })
  }
  
  Drawing_DL <- function(){
    Add.rim <- input$cbxAddRim_Drawing_DL
    Add.base <- input$cbxAddBase_Drawing_DL
    Add.rec <- input$cbxAddRec_Drawing_DL
    Add.diam <- input$cbxAddDiameter_Drawing_DL
    Add.lines <- input$cbxAddLines_Drawing_DL
    Add.scale <- input$cbxAddScale_Drawing_DL
    Add.cspi <- input$cbxAddCSPI_Drawing_DL
    Add.mask <- input$cbxAddMask_Drawing_DL
    Cs.val <- cs.val
    Rot <- -deg2rad(input$sldRot_Drawing_DL)
    Trans <- input$sldTrans_Drawing_DL
    Lat <<- input$sldLightLat_Drawing_DL
    Lon <<- input$sldLightLon_Drawing_DL
    Dis <<- input$sldLightDis_Drawing_DL*10
    Filt <<- input$cbxFilt_Drawing_DL
    Lims <<- input$sldLims_Drawing_DL; Lims <<- sort(Lims)
    Tresh1 <<- input$sldTresh1_Drawing_DL; Tresh1 <<- sort(Tresh1)
    Tresh2 <<- input$sldTresh2_Drawing_DL; Tresh2 <<- sort(Tresh2)
    Tresh3 <<- input$sldTresh3_Drawing_DL; Tresh3 <<- sort(Tresh3)
    Proc <<- input$sldProc_Drawing_DL
    if (Rot==Rot.count & Trans==Trans.count & Lat==Lat.count & Lon==Lon.count & Dis==Dis.count) { }
    else {
      mesh <- mesh.transform4drawing(mesh.drawing, rho=Rot, x=Trans)
      mesh.decim.maxi <- mesh.transform4drawing(mesh.drawing.decim.maxi, rho=Rot, x=Trans)
      while (length(rgl.dev.list()) > 0) { rgl.close() }
      par3d(windowRect=c(960-3000,30-3000,1920-3000,1040-3000))
      plot3d(mesh, meshColor="legacy")
      rgl.viewpoint(0,-90,fov=0,interactive = T)
      visi <- glVisible(mesh)
      rgl.close()
      mesh.visi <<- rmVertex(mesh, index=which(visi), keep=T)
      while (length(rgl.dev.list()) > 0) { rgl.close() }
      par3d(windowRect=c(960-3000,30-3000,1920-3000,1040-3000))
      rgl.viewpoint(0,-90,fov=0,interactive = F)
      plot3d(mesh.decim.maxi, meshColor="legacy")
      visi2 <- glVisible(mesh.decim.maxi)
      rgl.close()
      mesh.visi.decim.maxi <<- rmVertex(mesh.decim.maxi, index=which(visi2), keep=T)
      if (Rot==Rot.count2 & Trans==Trans.count2) { }
      else {
        mesh.bd.pts <<- pottery.outline(mesh.visi.decim.maxi, method="as", as.alpha=2, plots=F, infos=F)
      }
      mesh.drawing.centroid <<- mesh.centroid(mesh.visi)
      light <- sph2cart(c(deg2rad(-90+Lon),deg2rad(Lat),Dis))+mesh.drawing.centroid
      DL <- directional.lighting(mesh.visi,light)
      if ((sum(is.nan(DL))>0)) { DL[is.nan(DL)] <- 1 }
      mesh.visi.pts <- t(mesh.visi$vb[c(1,3),])
      Rot.count <<- Rot
      Trans.count <<- Trans
      Lat.count <<- Lat
      Lon.count <<- Lon
      Dis.count <<- Dis
      mesh.visi.pts <<- mesh.visi.pts
      DL.visi <<- DL
    }
    if (Filt==T) { DL.visi <- DL.visi[val.filter(DL.visi)] }
    if (Lims[1]==Lims.count[1] & Lims[2]==Lims.count[2]) { }
    else {
      low.lim <- round(as.numeric(c(quantile(DL.visi,probs=Lims[1]/100))),3)
      up.lim <- round(as.numeric(c(quantile(DL.visi,probs=Lims[2]/100))),3)
      updateSliderInput(session, "sldTresh1_Drawing_DL", value=c(0.3,1), min=low.lim, max=up.lim, step=0.05)
      updateSliderInput(session, "sldTresh2_Drawing_DL", value=c(low.lim,low.lim), min=low.lim, max=up.lim, step=0.05)
      updateSliderInput(session, "sldTresh3_Drawing_DL", value=c(low.lim,low.lim), min=low.lim, max=up.lim, step=0.05)
      Lims.count <<- Lims
      low.lim <<- low.lim
      up.lim <<- up.lim
    }
    output$plotDrawingHist<- renderPlot({
      h <- hist(DL.visi, breaks=201, main="Select treshold value", xlim=c(low.lim,up.lim))
      colo <- rep("white", length(h$counts))
      colo[ h$breaks>Tresh1[1] & h$breaks<Tresh1[2] ] <- "grey"
      colo[ h$breaks>Tresh2[1] & h$breaks<Tresh2[2] ] <- "grey"
      colo[ h$breaks>Tresh3[1] & h$breaks<Tresh3[2] ] <- "grey"
      plot(h, col=colo, add=T)
    })
    mesh.visi.pts.drawn1 <<- mesh.visi.pts[DL.visi>Tresh1[1] & DL.visi<Tresh1[2],]
    mesh.visi.pts.drawn2 <<- mesh.visi.pts[DL.visi>Tresh2[1] & DL.visi<Tresh2[2],]
    mesh.visi.pts.drawn3 <<- mesh.visi.pts[DL.visi>Tresh3[1] & DL.visi<Tresh3[2],]
    if (Proc!=100) {
      sam <- sample(1:dim(mesh.visi.pts.drawn1)[1],round(dim(mesh.visi.pts.drawn1)[1]*Proc/100))
      mesh.visi.pts.drawn1 <<- mesh.visi.pts.drawn1[sam,]
      sam <- sample(1:dim(mesh.visi.pts.drawn2)[1],round(dim(mesh.visi.pts.drawn2)[1]*Proc/100))
      mesh.visi.pts.drawn2 <<- mesh.visi.pts.drawn2[sam,]
      sam <- sample(1:dim(mesh.visi.pts.drawn3)[1],round(dim(mesh.visi.pts.drawn3)[1]*Proc/100))
      mesh.visi.pts.drawn3 <<- mesh.visi.pts.drawn3[sam,]
    }
    plot(0, type="n", xlim=pottery$x.lim, ylim=pottery$y.lim, asp=1, axes=F, xlab="", ylab="")
    points(mesh.visi.pts.drawn1, col="grey", cex=0.1)
    points(mesh.visi.pts.drawn2, col="grey", cex=0.1)
    points(mesh.visi.pts.drawn3, col="grey", cex=0.1)
    lines(mesh.bd.pts)
    pts <<- rbind(pts,c(input$plotDrawing_click$x, input$plotDrawing_click$y))
    if (dim(pts)[1]>0) { updateLIN() }
    if (Add.mask==T) { polygon(pottery$mask, col="white", border="white")  }
    reconst.2d(pottery=pottery, ae.scale=ae.scale, LIN=LIN, REC=REC, CSPI=CSPI, main="", add=T,
               profil.fill.col="darkgrey", lines.col="black", lines.lwd=1,
               add.rim=Add.rim, add.base=Add.base,
               add.lines=Add.lines, l.lines.col="grey", l.lines.lwd=1,
               add.diam=Add.diam,
               add.scale=Add.scale, scale.fill.col="black",
               add.cspi=Add.cspi, cspi.fill.col="black", cspi.bord.col="black",
               add.rec=Add.rec)
  }
  
  Drawing_DL_F <- function(){
    output$plotDrawing <- renderPlot({ Drawing_DL() })
    output$plotDrawingHist<- renderPlot({})
  }
  
  Drawing_AO <- function(){
    
    Add.rim <- input$cbxAddRim_Drawing_AO
    Add.base <- input$cbxAddBase_Drawing_AO
    Add.rec <- input$cbxAddRec_Drawing_AO
    Add.lines <- input$cbxAddLines_Drawing_AO
    Add.diam <- input$cbxAddDiameter_Drawing_AO
    Add.scale <- input$cbxAddScale_Drawing_AO
    Add.cspi <- input$cbxAddCSPI_Drawing_AO
    Add.mask <- input$cbxAddMask_Drawing_AO
    Cs.val <- cs.val
    Rot <- -deg2rad(input$sldRot_Drawing_AO)
    Trans <- input$sldTrans_Drawing_AO
    Filt <<- input$cbxFilt_Drawing_AO
    Lims <- input$sldLims_Drawing_AO; Lims <<- sort(Lims)
    Tresh1 <<- input$sldTresh1_Drawing_AO; Tresh1 <<- sort(Tresh1)
    Tresh2 <<- input$sldTresh2_Drawing_AO; Tresh2 <<- sort(Tresh2)
    Tresh3 <<- input$sldTresh3_Drawing_AO; Tresh3 <<- sort(Tresh3)
    Proc <<- input$sldProc_Drawing_AO
    if (Rot==Rot.count2 & Trans==Trans.count2) { }
    else {
      mesh <- mesh.transform4drawing(mesh.drawing, rho=Rot, x=Trans)
      mesh.decim.maxi <- mesh.transform4drawing(mesh.drawing.decim.maxi, rho=Rot, x=Trans)
      while (length(rgl.dev.list()) > 0) { rgl.close() }
      par3d(windowRect=c(960-3000,30-3000,1920-3000,1040-3000))
      plot3d(mesh, meshColor="legacy")
      rgl.viewpoint(0,-90,fov=0,interactive = F)
      visi <- glVisible(mesh)
      rgl.close()
      mesh.visi <<- rmVertex(mesh, index=which(visi), keep=T)
      while (length(rgl.dev.list()) > 0) { rgl.close() }
      par3d(windowRect=c(960-3000,30-3000,1920-3000,1040-3000))
      rgl.viewpoint(0,-90,fov=0,interactive = F)
      plot3d(mesh.decim.maxi, meshColor="legacy")
      visi2 <- glVisible(mesh.decim.maxi)
      rgl.close()
      mesh.visi.decim.maxi <<- rmVertex(mesh.decim.maxi, index=which(visi2), keep=T)
      mesh.bd.pts <<- pottery.outline(mesh.visi.decim.maxi, method="as", as.alpha=2, plots=F, infos=F)
      mesh.visi.pts <- t(mesh.visi$vb[c(1,3),])
      AO <- mesh.drawing$ao
      AO.visi <<- AO[visi]
      Rot.count2 <<- Rot
      Trans.count2 <<- Trans
      mesh.visi.pts <<- mesh.visi.pts
      AO.visi <<- AO.visi
    }
    if (Filt==T) { AO.visi <- AO.visi[val.filter(AO.visi)] }
    if (Lims[1]==Lims.count[1] & Lims[2]==Lims.count[2]) { }
    else {
      low.lim <- round(as.numeric(c(quantile(AO.visi,probs=Lims[1]/100))),3)
      up.lim <- round(as.numeric(c(quantile(AO.visi,probs=Lims[2]/100))),3)
      updateSliderInput(session, "sldTresh1_Drawing_AO", value=c(low.lim,low.lim), min=low.lim, max=up.lim, step=0.0005)
      updateSliderInput(session, "sldTresh2_Drawing_AO", value=c(low.lim,low.lim), min=low.lim, max=up.lim, step=0.0005)
      updateSliderInput(session, "sldTresh3_Drawing_AO", value=c(low.lim,low.lim), min=low.lim, max=up.lim, step=0.0005)
      Lims.count <<- Lims
      low.lim <<- low.lim
      up.lim <<- up.lim
    }
    output$plotDrawingHist<- renderPlot({
      h <- hist(AO.visi, breaks=401, main="Select treshold value", xlim=c(low.lim,up.lim))
      colo <- rep("white", length(h$counts))
      colo[ h$breaks>Tresh1[1] & h$breaks<Tresh1[2] ] <- "grey"
      colo[ h$breaks>Tresh2[1] & h$breaks<Tresh2[2] ] <- "grey"
      colo[ h$breaks>Tresh3[1] & h$breaks<Tresh3[2] ] <- "grey"
      plot(h, col=colo, add=T)
    })
    mesh.visi.pts.drawn1 <<- mesh.visi.pts[AO.visi>Tresh1[1] & AO.visi<Tresh1[2],]
    mesh.visi.pts.drawn2 <<- mesh.visi.pts[AO.visi>Tresh2[1] & AO.visi<Tresh2[2],]
    mesh.visi.pts.drawn3 <<- mesh.visi.pts[AO.visi>Tresh3[1] & AO.visi<Tresh3[2],]
    if (Proc!=100) {
      sam <- sample(1:dim(mesh.visi.pts.drawn1)[1],round(dim(mesh.visi.pts.drawn1)[1]*Proc/100))
      mesh.visi.pts.drawn1 <<- mesh.visi.pts.drawn1[sam,]
      sam <- sample(1:dim(mesh.visi.pts.drawn2)[1],round(dim(mesh.visi.pts.drawn2)[1]*Proc/100))
      mesh.visi.pts.drawn2 <<- mesh.visi.pts.drawn2[sam,]
      sam <- sample(1:dim(mesh.visi.pts.drawn3)[1],round(dim(mesh.visi.pts.drawn3)[1]*Proc/100))
      mesh.visi.pts.drawn3 <<- mesh.visi.pts.drawn3[sam,]
    }
    plot(0, type="n", xlim=pottery$x.lim, ylim=pottery$y.lim, asp=1, axes=F, xlab="", ylab="")
    points(mesh.visi.pts.drawn1, col="grey", cex=0.1)
    points(mesh.visi.pts.drawn2, col="grey", cex=0.1)
    points(mesh.visi.pts.drawn3, col="grey", cex=0.1)
    lines(mesh.bd.pts)
    pts <<- rbind(pts,c(input$plotDrawing_click$x, input$plotDrawing_click$y))
    if (dim(pts)[1]>0) {
      updateLIN()
    }
    if (Add.mask==T) { polygon(pottery$mask, col="white", border="white")  }
    reconst.2d(pottery=pottery, ae.scale=ae.scale, LIN=LIN, REC=REC, CSPI=CSPI, main="", add=T,
               profil.fill.col="darkgrey", lines.col="black", lines.lwd=1,
               add.rim=Add.rim, add.base=Add.base,
               add.lines=Add.lines, l.lines.col="grey", l.lines.lwd=1,
               add.diam=Add.diam,
               add.scale=Add.scale, scale.fill.col="black",
               add.cspi=Add.cspi, cspi.fill.col="black", cspi.bord.col="black",
               add.rec=Add.rec)
  }
  
  Drawing_AO_F <- function(){
    output$plotDrawing <- renderPlot({ Drawing_AO() })
  }
  
  Drawing_AODL <- function(){
    Add.rim <- input$cbxAddRim_Drawing_AODL
    Add.base <- input$cbxAddBase_Drawing_AODL
    Add.rec <- input$cbxAddRec_Drawing_AODL
    Add.diam <- input$cbxAddDiameter_Drawing_AODL
    Add.lines <- input$cbxAddLines_Drawing_AODL
    Add.scale <- input$cbxAddScale_Drawing_AODL
    Add.cspi <- input$cbxAddCSPI_Drawing_AODL
    Add.mask <- input$cbxAddMask_Drawing_AODL
    Cs.val <- cs.val
    Rot <- -deg2rad(input$sldRot_Drawing_AODL)
    Trans <- input$sldTrans_Drawing_AODL
    Lat <<- input$sldLightLat_Drawing_AODL_DL
    Lon <<- input$sldLightLon_Drawing_AODL_DL
    Dis <<- input$sldLightDis_Drawing_AODL_DL*10
    Lims.DL <<- input$sldLims_Drawing_AODL_DL; Lims.DL <<- sort(Lims.DL)
    Tresh1.DL <<- input$sldTresh1_Drawing_AODL_DL; Tresh1.DL <<- sort(Tresh1.DL)
    Tresh2.DL <<- input$sldTresh2_Drawing_AODL_DL; Tresh2.DL <<- sort(Tresh2.DL)
    Tresh3.DL <<- input$sldTresh3_Drawing_AODL_DL; Tresh3.DL <<- sort(Tresh3.DL)
    Proc.DL <<- input$sldProc_Drawing_AODL_DL
    Lims.AO <<- input$sldLims_Drawing_AODL_AO; Lims.AO <<- sort(Lims.AO)
    Tresh1.AO <<- input$sldTresh1_Drawing_AODL_AO; Tresh1.AO <<- sort(Tresh1.AO)
    Tresh2.AO <<- input$sldTresh2_Drawing_AODL_AO; Tresh2.AO <<- sort(Tresh2.AO)
    Tresh3.AO <<- input$sldTresh3_Drawing_AODL_AO; Tresh3.AO <<- sort(Tresh3.AO)
    Proc.AO <<- input$sldProc_Drawing_AODL_AO
    if (Rot==Rot.count.AODL & Trans==Trans.count.AODL & Lat==Lat.count.AODL & Lon==Lon.count.AODL & Dis==Dis.count.AODL) { }
    else {
      mesh <- mesh.transform4drawing(mesh.drawing, rho=Rot, x=Trans)
      mesh.decim.maxi <- mesh.transform4drawing(mesh.drawing.decim.maxi, rho=Rot, x=Trans)
      while (length(rgl.dev.list()) > 0) { rgl.close() }
      par3d(windowRect=c(960-3000,30-3000,1920-3000,1040-3000))
      plot3d(mesh, meshColor="legacy")
      rgl.viewpoint(0,-90,fov=0,interactive = T)
      visi <- glVisible(mesh)
      rgl.close()
      mesh.visi <<- rmVertex(mesh, index=which(visi), keep=T)
      while (length(rgl.dev.list()) > 0) { rgl.close() }
      par3d(windowRect=c(960-3000,30-3000,1920-3000,1040-3000))
      rgl.viewpoint(0,-90,fov=0,interactive = F)
      plot3d(mesh.decim.maxi, meshColor="legacy")
      visi2 <- glVisible(mesh.decim.maxi)
      rgl.close()
      mesh.visi.decim.maxi <<- rmVertex(mesh.decim.maxi, index=which(visi2), keep=T)
      if (Rot==Rot.count2 & Trans==Trans.count2) { }
      else {mesh.bd.pts <<- pottery.outline(mesh.visi.decim.maxi, method="as", as.alpha=2, plots=F, infos=F) }
      mesh.drawing.centroid <<- mesh.centroid(mesh.visi)
      light <- sph2cart(c(deg2rad(-90+Lon),deg2rad(Lat),Dis))+mesh.drawing.centroid
      DL <- directional.lighting(mesh.visi,light)
      if ((sum(is.nan(DL))>0)) { DL[is.nan(DL)] <- 1 }
      mesh.visi.pts <- t(mesh.visi$vb[c(1,3),])
      AO <- mesh.drawing$ao
      AO.visi <<- AO[visi]
      Rot.count.AODL <<- Rot
      Trans.count.AODL <<- Trans
      Lat.count.AODL <<- Lat
      Lon.count.AODL <<- Lon
      Dis.count.AODL <<- Dis
      mesh.visi.pts <<- mesh.visi.pts
      DL.visi <<- DL
      AO.visi <<- AO.visi
    }
    if (Lims.DL[1]==Lims.count.DL[1] & Lims.DL[2]==Lims.count.DL[2]) { }
    else {
      low.lim.DL <- round(as.numeric(c(quantile(DL.visi,probs=Lims.DL[1]/100))),3)
      up.lim.DL <- round(as.numeric(c(quantile(DL.visi,probs=Lims.DL[2]/100))),3)
      updateSliderInput(session, "sldTresh1_Drawing_AODL_DL", value=c(0.3,1), min=low.lim.DL, max=up.lim.DL, step=0.05)
      updateSliderInput(session, "sldTresh2_Drawing_AODL_DL", value=c(low.lim.DL,low.lim.DL), min=low.lim.DL, max=up.lim.DL, step=0.05)
      updateSliderInput(session, "sldTresh3_Drawing_AODL_DL", value=c(low.lim.DL,low.lim.DL), min=low.lim.DL, max=up.lim.DL, step=0.05)
      Lims.count.DL <<- Lims.DL
      low.lim.DL <<- low.lim.DL
      up.lim.DL <<- up.lim.DL
    }
    if (Lims.AO[1]==Lims.count.AO[1] & Lims.AO[2]==Lims.count.AO[2]) { }
    else {
      low.lim.AO <- round(as.numeric(c(quantile(AO.visi,probs=Lims.AO[1]/100))),3)
      up.lim.AO <- round(as.numeric(c(quantile(AO.visi,probs=Lims.AO[2]/100))),3)
      updateSliderInput(session, "sldTresh1_Drawing_AODL_AO", value=c(low.lim.AO,low.lim.AO), min=low.lim.AO, max=up.lim.AO, step=0.0005)
      updateSliderInput(session, "sldTresh2_Drawing_AODL_AO", value=c(low.lim.AO,low.lim.AO), min=low.lim.AO, max=up.lim.AO, step=0.0005)
      updateSliderInput(session, "sldTresh3_Drawing_AODL_AO", value=c(low.lim.AO,low.lim.AO), min=low.lim.AO, max=up.lim.AO, step=0.0005)
      Lims.count.AO <<- Lims.AO
      low.lim.AO <<- low.lim.AO
      up.lim.AO <<- up.lim.AO
    }
    output$plotDrawingHist<- renderPlot({
      layout(matrix(1:2,1,2))
      h <- hist(DL.visi, breaks=201, main="Select treshold value", xlim=c(low.lim.DL,up.lim.DL))
      colo <- rep("white", length(h$counts))
      colo[ h$breaks>Tresh1.DL[1] & h$breaks<Tresh1.DL[2] ] <- "grey"
      colo[ h$breaks>Tresh2.DL[1] & h$breaks<Tresh2.DL[2] ] <- "grey"
      colo[ h$breaks>Tresh3.DL[1] & h$breaks<Tresh3.DL[2] ] <- "grey"
      plot(h, col=colo, add=T)
      h <- hist(AO.visi, breaks=401, main="Select treshold value", xlim=c(low.lim.AO,up.lim.AO))
      colo <- rep("white", length(h$counts))
      colo[ h$breaks>Tresh1.AO[1] & h$breaks<Tresh1.AO[2] ] <- "grey"
      colo[ h$breaks>Tresh2.AO[1] & h$breaks<Tresh2.AO[2] ] <- "grey"
      colo[ h$breaks>Tresh3.AO[1] & h$breaks<Tresh3.AO[2] ] <- "grey"
      plot(h, col=colo, add=T)
    })
    mesh.visi.pts.drawn1.DL <<- mesh.visi.pts[DL.visi>Tresh1.DL[1] & DL.visi<Tresh1.DL[2],]
    mesh.visi.pts.drawn2.DL <<- mesh.visi.pts[DL.visi>Tresh2.DL[1] & DL.visi<Tresh2.DL[2],]
    mesh.visi.pts.drawn3.DL <<- mesh.visi.pts[DL.visi>Tresh3.DL[1] & DL.visi<Tresh3.DL[2],]
    mesh.visi.pts.drawn1.AO <<- mesh.visi.pts[AO.visi>Tresh1.AO[1] & AO.visi<Tresh1.AO[2],]
    mesh.visi.pts.drawn2.AO <<- mesh.visi.pts[AO.visi>Tresh2.AO[1] & AO.visi<Tresh2.AO[2],]
    mesh.visi.pts.drawn3.AO <<- mesh.visi.pts[AO.visi>Tresh3.AO[1] & AO.visi<Tresh3.AO[2],]
    if (Proc.DL!=100) {
      sam <- sample(1:dim(mesh.visi.pts.drawn1.DL)[1],round(dim(mesh.visi.pts.drawn1.DL)[1]*Proc.DL/100))
      mesh.visi.pts.drawn1.DL <<- mesh.visi.pts.drawn1.DL[sam,]
      sam <- sample(1:dim(mesh.visi.pts.drawn2.DL)[1],round(dim(mesh.visi.pts.drawn2.DL)[1]*Proc.DL/100))
      mesh.visi.pts.drawn2.DL <<- mesh.visi.pts.drawn2.DL[sam,]
      sam <- sample(1:dim(mesh.visi.pts.drawn3.DL)[1],round(dim(mesh.visi.pts.drawn3.DL)[1]*Proc.DL/100))
      mesh.visi.pts.drawn3.DL <<- mesh.visi.pts.drawn3.DL[sam,]
    }
    if (Proc.AO!=100) {
      sam <- sample(1:dim(mesh.visi.pts.drawn1.AO)[1],round(dim(mesh.visi.pts.drawn1.AO)[1]*Proc.AO/100))
      mesh.visi.pts.drawn1.AO <<- mesh.visi.pts.drawn1.AO[sam,]
      sam <- sample(1:dim(mesh.visi.pts.drawn2.AO)[1],round(dim(mesh.visi.pts.drawn2.AO)[1]*Proc.AO/100))
      mesh.visi.pts.drawn2.AO <<- mesh.visi.pts.drawn2.AO[sam,]
      sam <- sample(1:dim(mesh.visi.pts.drawn3.AO)[1],round(dim(mesh.visi.pts.drawn3.AO)[1]*Proc.AO/100))
      mesh.visi.pts.drawn3.AO <<- mesh.visi.pts.drawn3.AO[sam,]
    }
    plot(0, type="n", xlim=pottery$x.lim, ylim=pottery$y.lim, asp=1, axes=F, xlab="", ylab="")
    points(mesh.visi.pts.drawn1.DL, col="black", cex=0.1)
    points(mesh.visi.pts.drawn2.DL, col="black", cex=0.1)
    points(mesh.visi.pts.drawn3.DL, col="black", cex=0.1)
    points(mesh.visi.pts.drawn1.AO, col="black", cex=0.1)
    points(mesh.visi.pts.drawn2.AO, col="black", cex=0.1)
    points(mesh.visi.pts.drawn3.AO, col="black", cex=0.1)
    lines(mesh.bd.pts)
    pts <<- rbind(pts,c(input$plotDrawing_click$x, input$plotDrawing_click$y))
    if (dim(pts)[1]>0) { updateLIN() }
    if (Add.mask==T) { polygon(pottery$mask, col="white", border="white")  }
    reconst.2d(pottery=pottery, ae.scale=ae.scale, LIN=LIN, REC=REC, CSPI=CSPI, main="", add=T,
               profil.fill.col="darkgrey", lines.col="black", lines.lwd=1,
               add.rim=Add.rim, #add.base=Add.base,
               add.lines=Add.lines, l.lines.col="grey", l.lines.lwd=1,
               add.diam=Add.diam,
               add.scale=Add.scale, scale.fill.col="black",
               add.cspi=Add.cspi, cspi.fill.col="black", cspi.bord.col="black",
               add.rec=Add.rec)
  }
  
  Drawing_AODL_F <- function(){
    output$plotDrawing <- renderPlot({ Drawing_AODL() })
    output$plotDrawingHist<- renderPlot({})
  }
  
  Drawing_symmetry <- function(){
    Add.rim <- input$cbxAddRim_Drawing_symmetry
    Add.base <- input$cbxAddBase_Drawing_symmetry
    Add.rec <- input$cbxAddRec_Drawing_symmetry
    Add.lines <- input$cbxAddLines_Drawing_symmetry
    Add.diam <- input$cbxAddDiameter_Drawing_symmetry
    Add.scale <- input$cbxAddScale_Drawing_symmetry
    Add.cspi <- input$cbxAddCSPI_Drawing_symmetry
    Add.mask <- input$cbxAddMask_Drawing_symmetry
    Cs.val <- cs.val
    Rot <- -deg2rad(input$sldRot_Drawing_symmetry)
    Trans <- input$sldTrans_Drawing_symmetry
    Meth <- input$selMethod_Drawing_symmetry
    if (Rot==Rot.count4 & Trans==Trans.count4 & Meth==Meth.count4) { }
    else {
      mesh <- mesh.drawing.decim.maxi
      if (Meth=="selMethod_Drawing_symmetry_axis") {
        err <- meshDist(mesh, distvec=mesh.error(mesh,method="axis"), sign=T, steps=100, plot=F)
        ae.colo.scale <- archeo.colo.scale(x=err$colramp[[2]], colo=err$colramp$col, scale.length=50, method="segments", nb.seg=5, scale.y=10)
      }
      if (Meth=="selMethod_Drawing_symmetry_profile") {
        err <- meshDist(mesh, distvec=ME, sign=T, steps=100, plot=F)
        ae.colo.scale <- archeo.colo.scale(x=err$colramp[[2]], colo=err$colramp$col, scale.length=50, method="by", nb.by=0.5, scale.y=10)
      }
      mesh.drawing.decim.maxi.quality.positioned <<- mesh.transform4drawing(err$colMesh, rho=Rot, x=Trans)
      while (length(rgl.dev.list()) > 0) { rgl.close() }
      save <- par3d(skipRedraw=T)
      par3d(windowRect=c(960-3000,30-3000,1920-3000,1040-3000))
      view3d(theta = 0, phi = -90, fov=0, interactive=F)
      clear3d(type="light")
      light3d(ambient="white", diffuse="white", specular="black")
      plot3d(mesh.drawing.decim.maxi.quality.positioned, meshColor="legacy", add=T)
      par3d(save)
      rgl.snapshot(PlotQualPath)
      rgl.close()
      xyz <- t(mesh.drawing.decim.maxi.quality.positioned$vb[1:3,])
      x <- xyz[,1]
      y <- xyz[,3]
      recta <- c(min(x),min(y),max(x),max(y))
      Img <- readPNG(PlotQualPath)
      Img <- cut.img(Img, plots=F)
      Img <- as.raster(Img)
      Rot.count4 <<- Rot
      Trans.count4 <<- Trans
      Meth.count4 <<- Meth
      Img <<- Img
      recta <<- recta
      mesh.drawing.decim.maxi.quality.positioned <<- mesh.drawing.decim.maxi.quality.positioned
      ae.colo.scale <<- ae.colo.scale
    }
    plot(0, type="n", xlim=pottery$x.lim, ylim=pottery$y.lim, asp=1, axes=F, xlab="", ylab="")
    rasterImage(Img, xleft=recta[1], ybottom=recta[2], xright=recta[3], ytop=recta[4])
    pts <<- rbind(pts,c(input$plotDrawing_click$x, input$plotDrawing_click$y))
    if (dim(pts)[1]>0) { updateLIN() }
    if (Add.mask==T) { polygon(pottery$mask, col="white", border="white") }
    reconst.2d(pottery=pottery, ae.scale=ae.scale, LIN=LIN, REC=REC, CSPI=CSPI, main="", add=T,
               profil.fill.col="darkgrey", lines.col="black", lines.lwd=1,
               add.rim=Add.rim, add.base=Add.base,
               add.lines=Add.lines, l.lines.col="grey", l.lines.lwd=1,
               add.diam=Add.diam,
               add.scale=Add.scale, scale.fill.col="black",
               add.cspi=Add.cspi, cspi.fill.col="black", cspi.bord.col="black",
               add.rec=Add.rec)
    add.archeo.colo.scale(ae.colo.scale, text.cex=0.8)
  }
  
  Drawing_symmetry_F <- function(){
    MEFile <- paste(WD,"/temp/temp_ME.txt", sep="")
    if (!file.exists(MEFile)) { showNotification(paste("Pottery regularity is not calculated yet. Please wait a while."), duration = dur);  return(NULL) }
    ME <<- as.numeric(unlist(read.table(MEFile, sep=";")))
    if (length(ME)!=dim(mesh.drawing.decim.maxi$vb)[2]) { showNotification(paste("ME for diferent model stored in ME file. ME was probably not calculated yet."), duration = dur);  return(NULL) }
    if (input$selMethod_Drawing_symmetry=="") { showNotification(paste("No method selected. Please choose the method."), duration = dur); return(NULL) }
    output$plotDrawing <- renderPlot({ Drawing_symmetry() })
    output$plotDrawingHist<- renderPlot({})
  }
  
  Errase_line <- function(){
    if (length(pts)==2) {
      pts <<- matrix(numeric(0),0,2)
      LIN <<- array(numeric(0), dim=c(2,2,0))
    }
    if (length(pts)>2) { pts <<- pts[-dim(pts)[1],] }
  }
  
  observeEvent(input$btnDraw, {
    if (!exists("mesh.drawing")) { showNotification(paste("No oriented model loaded. Please load the model."), duration = dur);  return(NULL) }
    if (!exists("profil")) { showNotification(paste("No profile found. Please get the profile first."), duration = dur);  return(NULL) }
    if (input$selDrawingMethod=="") {  showNotification(paste("No drawing method selected. Please select the method."), duration = dur);  return(NULL)  }
    while (length(rgl.dev.list()) > 0) { rgl.close() }
    
    if (input$selDrawingMethod=="Drawing_linear") {
      
      Drawing_linear_F()
      
      showNotification(paste("Linear drawing performed"), duration = dur)
      
      }
    if (input$selDrawingMethod=="Drawing_photo") { Drawing_photo_F() ; showNotification(paste("Photographic drawing performed"), duration = dur)}
    if (input$selDrawingMethod=="Drawing_DL") { Drawing_DL_F() ; showNotification(paste("Shaded drawing (type 1) performed"), duration = dur)}
    if (input$selDrawingMethod=="Drawing_AO") { Drawing_AO_F() ; showNotification(paste("Shaded drawing (type 2) performed"), duration = dur)}
    if (input$selDrawingMethod=="Drawing_AODL") { Drawing_AODL_F() ; showNotification(paste("Shaded drawing (combined) performed"), duration = dur) }
    if (input$selDrawingMethod=="Drawing_symmetry") { Drawing_symmetry_F() ; showNotification(paste("Pottery regularity drawing performed"), duration = dur) }
  })
  
  observeEvent(input$btnSave_Drawing, {
    if (!exists("mesh.drawing")) { showNotification(paste("No oriented model loaded. Please load the model."), duration = dur);  return(NULL) }
    if (!exists("profil")) { showNotification(paste("No profile found. Please get the profile first."), duration = dur);  return(NULL) }
    if (input$selDrawingMethod=="") {  showNotification(paste("No drawing method selected. Please select the method."), duration = dur);  return(NULL)  }
    
    SaveDrawingPath <- input$saveDrawingPath
    FileName <- substr(FileName, 1, nchar(FileName)-4)
    if (!file.exists(SaveDrawingPath)) { dir.create(file.path(SaveDrawingPath),recursive=T) }

    setEPS()
    
    if (input$selDrawingMethod=="Drawing_linear") {
      postscript(paste(SaveDrawingPath,FileName,"_drawing_linear.eps", sep=""), pointsize=Ps.PtS)
      Drawing_linear()
    }
    if (input$selDrawingMethod=="Drawing_photo") {
      if (input$selColour_Drawing_photo=="") { showNotification(paste("No color mode selected. Please select color mode."), duration = dur);  return(NULL) }
      if (input$selColour_Drawing_photo=="selColour_Drawing_photo_color") {
        postscript(paste(SaveDrawingPath,FileName,"_drawing_photo_colo.eps", sep=""), pointsize=Ps.PtS)
      }
      if (input$selColour_Drawing_photo=="selColour_Drawing_photo_greyscale") {
        postscript(paste(SaveDrawingPath,FileName,"_drawing_photo_grey.eps", sep=""), pointsize=Ps.PtS)
      }
      Drawing_photo()
    }
    if (input$selDrawingMethod=="Drawing_DL") {
      postscript(paste(SaveDrawingPath,FileName,"_drawing_shaded_DL.eps", sep=""), pointsize=Ps.PtS)
      Drawing_DL()
    }
    if (input$selDrawingMethod=="Drawing_AO") {
      postscript(paste(SaveDrawingPath,FileName,"_drawing_shaded_AO.eps", sep=""), pointsize=Ps.PtS)
      Drawing_AO()
    }
    if (input$selDrawingMethod=="Drawing_AODL") {
      postscript(paste(SaveDrawingPath,FileName,"_drawing_shaded_AODL.eps", sep=""), pointsize=Ps.PtS)
      Drawing_AODL()
    }
    if (input$selDrawingMethod=="Drawing_symmetry") {
      postscript(paste(SaveDrawingPath,FileName,"_drawing_regularity.eps", sep=""), pointsize=Ps.PtS)
      Drawing_symmetry()
    }
  
    dev.off()
    showNotification(paste("Illustration saved"), duration = dur)
  })
  
  observeEvent(input$btnErraseLine_Drawing_linear, {
    Errase_line()
    Drawing_linear_F()
  })
  
  observeEvent(input$btnErraseLine_Drawing_photo, {
    Errase_line()
    Drawing_photo_F()
  })
  
  observeEvent(input$btnErraseLine_Drawing_DL, {
    Errase_line()
    Drawing_DL_F()
  })
  
  observeEvent(input$btnErraseLine_Drawing_AO, {
    Errase_line()
    Drawing_AO_F()
  })
  
  observeEvent(input$btnErraseLine_Drawing_AODL, {
    Errase_line()
    Drawing_AODL_F()
  })
    
  observeEvent(input$btnErraseLine_Drawing_symmetry, {
    Errase_line()
    Drawing_symmetry_F()
  })
  
  observeEvent(input$btnErraseRuler_Drawing_linear, {
    if (length(ruler)==2) { ruler <<- matrix(numeric(0),0,2) }
    if (length(ruler)>2) { ruler <<- ruler[-dim(ruler)[1],] }
    Drawing_linear_F()
  })
  
  
  # # 3D RECONSTRUCTION
  
  dataInputOrientedMesh2 <- reactive({
    inFile <- input$fileOrientedMesh2
    if (is.null(inFile)) return(NULL)
    FileName <<- inFile$name
    FileExt <- substr(FileName, nchar(FileName)-2, nchar(FileName))
    if (FileExt!="ply") { showNotification(paste("Please load model in PLY format"), duration = dur); return(NULL) }
    file <- inFile$datapath
    filen <- gsub("[\\]", "/", file)
    filen <- paste(substr(filen, 1, nchar(filen)-1), inFile$name, sep="")
    file.rename(file, filen)
    mesh <- vcgImport(filen, clean=F, readcolor=T)
    mesh.oriented <<- mesh
    return(mesh.oriented)
  })
  
  observeEvent(input$fileOrientedMesh2, {
    mesh.drawing <<- dataInputOrientedMesh2()
    dataInputPreparation2()
    showNotification(paste("Oriented model loaded"), duration = dur)
  })
  
  dataInputPreparation2 <- reactive({
    mesh.drawing <<- mesh.align2front(mesh.drawing)
    mesh.drawing <<- mesh.align2xy(mesh.drawing)
    mesh.drawing.decim.mini <<- decim(mesh.drawing, 2000)
    mesh.drawing.decim.mini <<- vcgUpdateNormals(mesh.drawing.decim.mini)
    output$plot_Reconst <- renderPlot({
      plotOrientedMesh(mesh.drawing.decim.mini, method="rz")
    })
  })
  
  observeEvent(input$btnGetProfile2, {
    if (!exists("mesh.drawing")) { showNotification(paste("No oriented model loaded. Please load the model."), duration = dur);  return(NULL) }
    if (input$selProfileMethod2=="") {  showNotification(paste("No profile method selected. Please select the method."), duration = dur);  return(NULL)  }
    withProgress(message = 'Profile extraction', value=0.1, {
      if (input$selProfileMethod2=="Profile_envelop")    { profil <<- profile2d(mesh.drawing.decim.mini, method = "envelop", plots=T) ; mai <- c("Whole envelop profile") }
      if (input$selProfileMethod2=="Profile_middle")     { profil <<- profile2d(mesh.drawing, method = "middle", plots=T) ; mai <- c("In the middle of the fragment profile") }
      if (input$selProfileMethod2=="Profile_longest")    { profil <<- profile2d(mesh.drawing, method = "longest", longest.by.deg=2, plots=T) ; mai <- c("Longest-preserved profile") }
      if (input$selProfileMethod2=="Profile_arbitrary")  { profil <<- profile2d(mesh.drawing, method = "arbitrary", plots=T) ; mai <- c("Arbitrary selected profile") }
      setProgress(0.8)
      profil <<- switch.profile(profil)
      pottery <<- profil2pottery(profil, x.off = 5, off.side=10)
      output$plot_Reconst <- renderPlot({
        plot(pottery$profil, main=mai, xlab="", ylab="", type="n", asp=1)
        lines(pottery$profil)
        lines(pottery$extern, col="blue")
        lines(pottery$intern, col="red")
      })
      setProgress(1)
    })
  })
  
  observeEvent(input$btnReconst3d, {
    if (!exists("mesh.drawing")) { showNotification(paste("No oriented model loaded. Please load the model."), duration = dur);  return(NULL) }
    if (!exists("profil")) { showNotification(paste("No profile found. Please get the profile first."), duration = dur);  return(NULL) }
    if (input$selReconst3dMethod=="") { showNotification(paste("No 3D reconstruction method selected. Please select method."), duration = dur);  return(NULL) }
    if (!exists("reconst.sect")) { reconst.sect <<- reconst.3d.slices(profil, plots=F, seg=720) }
    if (!exists("reconst.plas")) { reconst.plas <<- reconst.3d.mesh(profil, seg=128, col="grey") }
    Add.model <<- input$cbxAddModel_Reconst3d
    while (length(rgl.dev.list()) > 0) { rgl.close() }
    view3d(theta = 15, phi = -60, interactive=T)
    par3d(windowRect=c(960,30,1920,1040))
    rgl.bringtotop()
    if (input$selReconst3dMethod=="Model_pointcloud") {
      save <- par3d(skipRedraw=T)
      reconst.3d.slices.plot(reconst.sect, type="p")
      par3d(save)
    }
    if (input$selReconst3dMethod=="Model_sliced") {
      if (Slice.count==0) { reconst.3d.slices.plot(reconst.sect, type="l") ; Slice.count <<- 1}
      else {
        save <- par3d(skipRedraw=T)
        reconst.3d.slices.plot(reconst.sect, type="l")
        par3d(save)
        rgl.bringtotop()
      }
    }
    if (input$selReconst3dMethod=="Model_wireframe") { plot3d(reconst.plas, type="wire", col="lightgrey", add=T) }
    if (input$selReconst3dMethod=="Model_plastic") { plot3d(reconst.plas, col="lightgrey", add=T) }
    if (Add.model==T) { plot3d(mesh.drawing, meshColor="legacy", add=T) }
    rgl.bringtotop()
    Add.model.count <<- Add.model
  })
  
  observeEvent(input$cbxAddModel_Reconst3d, {
    if (!exists("mesh.drawing")) { showNotification(paste("No oriented model loaded. Please load the model."), duration = dur);  return(NULL) }
    if (!exists("profil")) { showNotification(paste("No profile found. Please get the profile first."), duration = dur);  return(NULL) }
    if (input$selReconst3dMethod=="") { showNotification(paste("No 3D reconstruction method selected. Please select method."), duration = dur);  return(NULL) }
    Add.model <<- input$cbxAddModel_Reconst3d
    if (Add.model==Add.model.count) {}
    else {
      if (!exists("reconst.sect")) { reconst.sect <<- reconst.3d.slices(profil, plots=F, seg=720) }
      if (!exists("reconst.plas")) { reconst.plas <<- reconst.3d.mesh(profil, seg=128, col="grey") }
      while (length(rgl.dev.list()) > 0) { rgl.close() }
      par3d(windowRect=c(960,30,1920,1040))
      view3d(theta = 15, phi = -60, interactive=T)
      rgl.bringtotop()
      if (input$selReconst3dMethod=="Model_pointcloud") {
        save <- par3d(skipRedraw=T)
        reconst.3d.slices.plot(reconst.sect, type="p")
        par3d(save)
      }
      if (input$selReconst3dMethod=="Model_sliced") {
        if (Slice.count==0) { rgl.bringtotop(); reconst.3d.slices.plot(reconst.sect, type="l") ; Slice.count <<- 1}
        else {
          save <- par3d(skipRedraw=T)
          reconst.3d.slices.plot(reconst.sect, type="l")
          par3d(save)
        }
      }
      if (input$selReconst3dMethod=="Model_wireframe") { plot3d(reconst.plas, type="wire", col="lightgrey", add=T) }
      if (input$selReconst3dMethod=="Model_plastic") { plot3d(reconst.plas, col="lightgrey", add=T) }
      if (input$cbxAddModel_Reconst3d==T) { plot3d(mesh.drawing, meshColor="legacy", add=T) }
      rgl.bringtotop()
      Add.model.count <<- Add.model
    }
  })
  
  observeEvent(input$btnSave_Reconstruction, {
    if (!exists("mesh.drawing")) { showNotification(paste("No oriented model loaded. Please load the model."), duration = dur);  return(NULL) }
    if (!exists("profil")) { showNotification(paste("No profile found. Please get the profile first."), duration = dur);  return(NULL) }
    if (!exists("reconst.sect")) { showNotification(paste("No 3D reconstruction found. Please perform 3D reconstruction first."), duration = dur);  return(NULL) }
    if (!exists("reconst.plas")) { showNotification(paste("No 3D reconstruction found. Please perform 3D reconstruction first."), duration = dur);  return(NULL) }
    SavePath <- input$saveReconstructionPath
    if (!file.exists(SavePath)) { dir.create(file.path(SavePath),recursive=T) }
    FileName <- substr(FileName, 1, nchar(FileName)-4)
    if (input$selReconst3dMethod=="Model_pointcloud" || input$selReconst3dMethod=="Model_sliced") {
      FullFilePath <- paste(SavePath,FileName,"_pointCloud.xyz", sep="")
      M <- A2M(reconst.sect)
      write.table(M, FullFilePath, row.names=F, col.names=F)
      showNotification(paste("Pointcloud of 3D reconstruction saved."), duration = dur)
    }
    if (input$selReconst3dMethod=="Model_wireframe" || input$selReconst3dMethod=="Model_plastic") {
      FullFilePath <- paste(SavePath,FileName,"_reconst",sep="")
      vcgPlyWrite(reconst.plas, FullFilePath, binary = T)
      meshi <- paste(FullFilePath,".ply",sep="")
      mesho <- paste(FullFilePath,".ply",sep="")
      meshs <- paste(WD,"/support/_04_reorient_normals_mesh_export.mlx",sep="")
      meshm <- "vn vc"
      MeshFullCmd <- paste("meshlabserver -i", meshi, "-o", mesho, "-m", meshm, "-s", meshs)
      shell(MeshFullCmd, wait=T)
      showNotification(paste("Mesh of 3D reconstruction saved."), duration = dur)
    }
  })
  
  
  # Updaters
  observeEvent(input$cbxAddRim_Drawing_linear, {
    val <- input$cbxAddRim_Drawing_linear
    updateCheckboxInput(session, "cbxAddRim_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddRim_Drawing_photo, {
    val <- input$cbxAddRim_Drawing_photo
    updateCheckboxInput(session, "cbxAddRim_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddRim_Drawing_DL, {
    val <- input$cbxAddRim_Drawing_DL
    updateCheckboxInput(session, "cbxAddRim_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddRim_Drawing_AO, {
    val <- input$cbxAddRim_Drawing_AO
    updateCheckboxInput(session, "cbxAddRim_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddRim_Drawing_AODL, {
    val <- input$cbxAddRim_Drawing_AODL
    updateCheckboxInput(session, "cbxAddRim_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddRim_Drawing_symmetry, {
    val <- input$cbxAddRim_Drawing_symmetry
    updateCheckboxInput(session, "cbxAddRim_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddRim_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddBase_Drawing_linear, {
    val <- input$cbxAddBase_Drawing_linear
    updateCheckboxInput(session, "cbxAddBase_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddBase_Drawing_photo, {
    val <- input$cbxAddBase_Drawing_photo
    updateCheckboxInput(session, "cbxAddBase_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddBase_Drawing_DL, {
    val <- input$cbxAddBase_Drawing_DL
    updateCheckboxInput(session, "cbxAddBase_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddBase_Drawing_AO, {
    val <- input$cbxAddBase_Drawing_AO
    updateCheckboxInput(session, "cbxAddBase_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddBase_Drawing_AODL, {
    val <- input$cbxAddBase_Drawing_AODL
    updateCheckboxInput(session, "cbxAddBase_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddBase_Drawing_symmetry, {
    val <- input$cbxAddBase_Drawing_symmetry
    updateCheckboxInput(session, "cbxAddBase_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddBase_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddRec_Drawing_linear, {
    val <- input$cbxAddRec_Drawing_linear
    updateCheckboxInput(session, "cbxAddRec_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddRec_Drawing_photo, {
    val <- input$cbxAddRec_Drawing_photo
    updateCheckboxInput(session, "cbxAddRec_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddRec_Drawing_DL, {
    val <- input$cbxAddRec_Drawing_DL
    updateCheckboxInput(session, "cbxAddRec_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddRec_Drawing_AO, {
    val <- input$cbxAddRec_Drawing_AO
    updateCheckboxInput(session, "cbxAddRec_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddRec_Drawing_AODL, {
    val <- input$cbxAddRec_Drawing_AODL
    updateCheckboxInput(session, "cbxAddRec_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddRec_Drawing_symmetry, {
    val <- input$cbxAddRec_Drawing_symmetry
    updateCheckboxInput(session, "cbxAddRec_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddRec_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddScale_Drawing_linear, {
    val <- input$cbxAddScale_Drawing_linear
    updateCheckboxInput(session, "cbxAddScale_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddScale_Drawing_photo, {
    val <- input$cbxAddScale_Drawing_photo
    updateCheckboxInput(session, "cbxAddScale_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddScale_Drawing_DL, {
    val <- input$cbxAddScale_Drawing_DL
    updateCheckboxInput(session, "cbxAddScale_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddScale_Drawing_AO, {
    val <- input$cbxAddScale_Drawing_AO
    updateCheckboxInput(session, "cbxAddScale_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddScale_Drawing_AODL, {
    val <- input$cbxAddScale_Drawing_AODL
    updateCheckboxInput(session, "cbxAddScale_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddScale_Drawing_symmetry, {
    val <- input$cbxAddScale_Drawing_symmetry
    updateCheckboxInput(session, "cbxAddScale_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddScale_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddMask_Drawing_linear, {
    val <- input$cbxAddMask_Drawing_linear
    updateCheckboxInput(session, "cbxAddMask_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddMask_Drawing_photo, {
    val <- input$cbxAddMask_Drawing_photo
    updateCheckboxInput(session, "cbxAddMask_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddMask_Drawing_DL, {
    val <- input$cbxAddMask_Drawing_DL
    updateCheckboxInput(session, "cbxAddMask_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddMask_Drawing_AO, {
    val <- input$cbxAddMask_Drawing_AO
    updateCheckboxInput(session, "cbxAddMask_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddMask_Drawing_AODL, {
    val <- input$cbxAddMask_Drawing_AODL
    updateCheckboxInput(session, "cbxAddMask_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddMask_Drawing_symmetry, {
    val <- input$cbxAddMask_Drawing_symmetry
    updateCheckboxInput(session, "cbxAddMask_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddMask_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddCSPI_Drawing_linear, {
    val <- input$cbxAddCSPI_Drawing_linear
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddCSPI_Drawing_photo, {
    val <- input$cbxAddCSPI_Drawing_photo
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddCSPI_Drawing_DL, {
    val <- input$cbxAddCSPI_Drawing_DL
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddCSPI_Drawing_AO, {
    val <- input$cbxAddCSPI_Drawing_AO
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddCSPI_Drawing_AODL, {
    val <- input$cbxAddCSPI_Drawing_AODL
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddCSPI_Drawing_symmetry, {
    val <- input$cbxAddCSPI_Drawing_symmetry
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddCSPI_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddDiameter_Drawing_linear, {
    val <- input$cbxAddDiameter_Drawing_linear
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddDiameter_Drawing_photo, {
    val <- input$cbxAddDiameter_Drawing_photo
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddDiameter_Drawing_DL, {
    val <- input$cbxAddDiameter_Drawing_DL
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddDiameter_Drawing_AO, {
    val <- input$cbxAddDiameter_Drawing_AO
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddDiameter_Drawing_AODL, {
    val <- input$cbxAddDiameter_Drawing_AODL
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_symmetry", value=val)
  })
  observeEvent(input$cbxAddDiameter_Drawing_symmetry, {
    val <- input$cbxAddDiameter_Drawing_symmetry
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_linear", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_photo", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_DL", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_AO", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_AODL", value=val)
    updateCheckboxInput(session, "cbxAddDiameter_Drawing_symmetry", value=val)
  })
  observeEvent(input$sldTrans_Drawing_photo, {
    val <- input$sldTrans_Drawing_photo
    # updateSliderInput(session, "sldTrans_Drawing_photo", value=val)
    updateSliderInput(session, "sldTrans_Drawing_DL", value=val)
    updateSliderInput(session, "sldTrans_Drawing_AO", value=val)
    updateSliderInput(session, "sldTrans_Drawing_AODL", value=val)
    updateSliderInput(session, "sldTrans_Drawing_symmetry", value=val)
  })
  observeEvent(input$sldTrans_Drawing_DL, {
    val <- input$sldTrans_Drawing_DL
    updateSliderInput(session, "sldTrans_Drawing_photo", value=val)
    # updateSliderInput(session, "sldTrans_Drawing_DL", value=val)
    updateSliderInput(session, "sldTrans_Drawing_AO", value=val)
    updateSliderInput(session, "sldTrans_Drawing_AODL", value=val)
    updateSliderInput(session, "sldTrans_Drawing_symmetry", value=val)
  })
  observeEvent(input$sldTrans_Drawing_AO, {
    val <- input$sldTrans_Drawing_AO
    updateSliderInput(session, "sldTrans_Drawing_photo", value=val)
    updateSliderInput(session, "sldTrans_Drawing_DL", value=val)
    # updateSliderInput(session, "sldTrans_Drawing_AO", value=val)
    updateSliderInput(session, "sldTrans_Drawing_AODL", value=val)
    updateSliderInput(session, "sldTrans_Drawing_symmetry", value=val)
  })
  observeEvent(input$sldTrans_Drawing_AODL, {
    val <- input$sldTrans_Drawing_AODL
    updateSliderInput(session, "sldTrans_Drawing_photo", value=val)
    updateSliderInput(session, "sldTrans_Drawing_DL", value=val)
    updateSliderInput(session, "sldTrans_Drawing_AO", value=val)
    # updateSliderInput(session, "sldTrans_Drawing_AODL", value=val)
    updateSliderInput(session, "sldTrans_Drawing_symmetry", value=val)
  })
  observeEvent(input$sldTrans_Drawing_symmetry, {
    val <- input$sldTrans_Drawing_symmetry
    updateSliderInput(session, "sldTrans_Drawing_photo", value=val)
    updateSliderInput(session, "sldTrans_Drawing_DL", value=val)
    updateSliderInput(session, "sldTrans_Drawing_AO", value=val)
    updateSliderInput(session, "sldTrans_Drawing_AODL", value=val)
    # updateSliderInput(session, "sldTrans_Drawing_symmetry", value=val)
  })
  observeEvent(input$sldRot_Drawing_photo, {
    val <- input$sldRot_Drawing_photo
    # updateSliderInput(session, "sldRot_Drawing_photo", value=val)
    updateSliderInput(session, "sldRot_Drawing_DL", value=val)
    updateSliderInput(session, "sldRot_Drawing_AO", value=val)
    updateSliderInput(session, "sldRot_Drawing_AODL", value=val)
    updateSliderInput(session, "sldRot_Drawing_symmetry", value=val)
  })
  observeEvent(input$sldRot_Drawing_DL, {
    val <- input$sldRot_Drawing_DL
    updateSliderInput(session, "sldRot_Drawing_photo", value=val)
    # updateSliderInput(session, "sldRot_Drawing_DL", value=val)
    updateSliderInput(session, "sldRot_Drawing_AO", value=val)
    updateSliderInput(session, "sldRot_Drawing_AODL", value=val)
    updateSliderInput(session, "sldRot_Drawing_symmetry", value=val)
  })
  observeEvent(input$sldRot_Drawing_AO, {
    val <- input$sldRot_Drawing_AO
    updateSliderInput(session, "sldRot_Drawing_photo", value=val)
    updateSliderInput(session, "sldRot_Drawing_DL", value=val)
    # updateSliderInput(session, "sldRot_Drawing_AO", value=val)
    updateSliderInput(session, "sldRot_Drawing_AODL", value=val)
    updateSliderInput(session, "sldRot_Drawing_symmetry", value=val)
  })
  observeEvent(input$sldRot_Drawing_AODL, {
    val <- input$sldRot_Drawing_AODL
    updateSliderInput(session, "sldRot_Drawing_photo", value=val)
    updateSliderInput(session, "sldRot_Drawing_DL", value=val)
    updateSliderInput(session, "sldRot_Drawing_AO", value=val)
    # updateSliderInput(session, "sldRot_Drawing_AODL", value=val)
    updateSliderInput(session, "sldRot_Drawing_symmetry", value=val)
  })
  observeEvent(input$sldRot_Drawing_symmetry, {
    val <- input$sldRot_Drawing_symmetry
    updateSliderInput(session, "sldRot_Drawing_photo", value=val)
    updateSliderInput(session, "sldRot_Drawing_DL", value=val)
    updateSliderInput(session, "sldRot_Drawing_AO", value=val)
    updateSliderInput(session, "sldRot_Drawing_AODL", value=val)
    # updateSliderInput(session, "sldRot_Drawing_symmetry", value=val)
  })
  
  
})

