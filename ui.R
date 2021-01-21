#..............................................................................................
# D A C O R D  -  UI
#..............................................................................................
# Last update: 2018/01/22

library(shiny)
library(rgl)

width.adj <- c(2, 10)                       # width panels adjustement (side panel, main panel)
mainPanel.adj <- c("800px", "520px")        # main panel adjustement (width, height)


shinyUI(
  navbarPage("DACORD",
    tabPanel("Orientation",
       fluidPage(
         sidebarLayout(
           sidebarPanel(
             strong('1. Load file'),
             br(),
             fileInput('fileMesh', label=NULL, accept=c('.ply', '.PLY')),
             strong('2. Model pre-orientation'),
             selectInput('selPreorientationMethod', label=NULL, c("Select pre-orientation method" = "", "Manual" = "Manual", "Automatic" = "Automatic")),
             checkboxInput('cbxAdvancedPreorientation', label="Advanced options", value=FALSE),
             conditionalPanel(
               condition = "input.cbxAdvancedPreorientation",
               conditionalPanel(
                 condition = "input.selPreorientationMethod == 'Manual'",
                 sliderInput('sldDecim_PM', label="Initial decimation", value=5000, min=1000, max=20000, step=1000),
                 sliderInput('sldSmIter_PM', label="Smoothing iterations", value=5, min=1, max=40, step=1),
                 sliderInput('sldPerc_PM', label="Percentile for extraction", value=0.50, min=0.01, max=1.0, step=0.10),
                 checkboxInput('cbxPlots_PM', label="Visualisations", value=FALSE)
               ),
               conditionalPanel(
                 condition = "input.selPreorientationMethod == 'Automatic'",
                 sliderInput('sldDecim_PA', label="Initial decimation", value=5000, min=1000, max=20000, step=1000),
                 sliderInput('sldSmIter_PA', label="Smoothing iterations", value=5, min=1, max=40, step=1),
                 sliderInput('sldPVal_PA', label="p-value for extraction", value=0.95, min=0, max=1, step=0.01),
                 checkboxInput('cbxPlots_PA', label="Visualisations", value=FALSE)
               )
             ),
             actionButton('btnPreorientation', 'Run Pre-orientation'),
             br(),
             br(),
             br(),
             actionButton('btnInverse', 'Inverse'),
             br(),
             br(),
             br(),
             strong('3. Optimisation of the rotation axis position'),
             selectInput('selOrientationMethod', label=NULL, c("Select orientation method" = "",
                                                               "0. None" = "AOS_none",
                                                               "1a. Horizontal circle adjustment using radius (4DDL)" = "AOS_circle4DDL",
                                                               "1b. Horizontal circle adjustment using radius (3DDL)" = "AOS_circle3DDL",
                                                               "2.	Horizontal circle adjustment using multi-criteria approach" = "AOS_sections",
                                                               "3a. Vertical profile superimposition" = "AOS_profile",
                                                               "3b. Vertical profile curve fitting" = "AOS_polynomial",
                                                               "4a. Tangent plane to rim" = "AOS_rimbase_rim",
                                                               "4b. Tangent plane to base" = "AOS_rimbase_base")
             ),
             selectInput('selSide', label=NULL, c("Both" = "both", "Outer (Blue)" = "blue", "Inner (Red)" = "red")),

             conditionalPanel(
               condition = "input.selOrientationMethod == 'AOS_sections'",
               sliderInput('sldParetoNb_AOS_sections', label="Number of best solutions", value=2, min=1, max=20, step=1)
             ),
             checkboxInput('cbxAdvancedOrientation', label="Advanced options", value=FALSE),
             conditionalPanel(
               condition = "input.cbxAdvancedOrientation",
               conditionalPanel(
                 condition = "input.selOrientationMethod == 'AOS_circle4DDL'",
                 sliderInput('sldNb_AOS_circle4DDL', label="Number of sections", value=12, min=0, max=30, step=1),
                 checkboxInput('cbxPlots_AOS_circle4DDL', label="Visualisations", value=FALSE)
               ),
               conditionalPanel(
                 condition = "input.selOrientationMethod == 'AOS_circle3DDL'",
                 sliderInput('sldNb_AOS_circle3DDL', label="Number of sections", value=12, min=0, max=30, step=1),
                 checkboxInput('cbxPlots_AOS_circle3DDL', label="Visualisations", value=FALSE)
               ),
               conditionalPanel(
                 condition = "input.selOrientationMethod == 'AOS_sections'",
                 sliderInput('sldNb_AOS_sections', label="Number of sections", value=12, min=0, max=30, step=1),
                 sliderInput('sldLimits_AOS_sections', label="Search parameter limits", value=c(-10,10), min=-90, max=90, step=1),
                 sliderInput('sldGener_AOS_sections', label="Search generations", value=20, min=12, max=100, step=4),
                 sliderInput('sldPopsize_AOS_sections', label="Search popsize", value=20, min=12, max=100, step=4),
                 checkboxInput('cbxPlots_AOS_sections', label="Visualisations", value=FALSE)
               ),
               conditionalPanel(
                 condition = "input.selOrientationMethod == 'AOS_profile'",
                 sliderInput('sldNb_AOS_profile', label="Number of sections", value=12, min=0, max=30, step=1),
                 checkboxInput('cbxPlots_AOS_profile', label="Visualisations", value=FALSE)
               ),
               conditionalPanel(
                 condition = "input.selOrientationMethod == 'AOS_polynomial'",
                 sliderInput('sldPO_AOS_polynomial', label="Polynomial order", value=12, min=2, max=20, step=1),
                 checkboxInput('cbxPlots_AOS_polynomial', label="Visualisations", value=FALSE)
               ),
               conditionalPanel(
                 condition = "input.selOrientationMethod == 'AOS_rimbase_rim'",
                 sliderInput('sldNb_AOS_rimbase_rim', label="Number of sections", value=12, min=0, max=30, step=1),
                 sliderInput('sldTresh_AOS_rimbase_rim', label="Z-threshold", value=1, min=0.2, max=10, step=0.1),
                 checkboxInput('cbxPlots_AOS_rimbase_rim', label="Visualisations", value=FALSE)
               ),
               conditionalPanel(
                 condition = "input.selOrientationMethod == 'AOS_rimbase_base'",
                 sliderInput('sldNb_AOS_rimbase_base', label="Number of sections", value=12, min=0, max=30, step=1),
                 sliderInput('sldTresh_AOS_rimbase_base', label="Z-threshold", value=1, min=0.2, max=10, step=0.1),
                 checkboxInput('cbxPlots_AOS_rimbase_base', label="Visualisations", value=FALSE)
               )
             ),
             actionButton('btnOrientation', 'Run Orientation'),
             br(),
             br(),
             br(),
             strong('4. Save oriented model'),
             textInput('savePath', label=NULL, value="C:/DACORD/oriented/"),
             actionButton('btnSave', 'Save'),
             width=width.adj[1]
           ),
           mainPanel(
             rglwidgetOutput('plotOrientation3d', width = mainPanel.adj[1], height = mainPanel.adj[2]),
             plotOutput('plotOrientation'),
             width=width.adj[2]
           )
         )
       )
    ),

    tabPanel("Drawing",
      fluidPage(
        sidebarLayout(
          sidebarPanel(
            strong('1. Load file'),
            fileInput('fileOrientedMesh', label=NULL, accept=c('.ply')),
            strong('2. Profile type selection'),
            selectInput('selProfileMethod', label=NULL,
                        c("Select profile method" = "",
                          "Whole envelope profile" = "Profile_envelop",
                          "In the middle of the fragment" = "Profile_middle",
                          "Longest-preserved profile" = "Profile_longest",
                          "Arbitrarily selected profile" = "Profile_arbitrary")
            ),
            actionButton('btnGetProfile', 'Get profile'),
            actionButton('btnRepairProfile', 'Erase last line'),
            br(),
            br(),
            br(),
            strong('3. Drawing type selection'),
            selectInput('selDrawingMethod', label=NULL,
                        c("Select drawing method" = "",
                          "Linear drawing" = "Drawing_linear",
                          "Photographic drawing" = "Drawing_photo",
                          "Shaded drawing (DL)" = "Drawing_DL",
                          "Shaded drawing (AO)" = "Drawing_AO",
                          "Shaded drawing (combined)" = "Drawing_AODL",
                          "Pottery regularity" = "Drawing_symmetry")
            ),
            conditionalPanel(
              condition = "input.selDrawingMethod == 'Drawing_linear'",
              checkboxInput('cbxFeatures_Drawing_linear', label=strong('Additional features'), value=FALSE),
              conditionalPanel(
                condition = "input.cbxFeatures_Drawing_linear",
                checkboxInput('cbxAddRim_Drawing_linear', label="Rim"),
                checkboxInput('cbxAddBase_Drawing_linear', label="Base"),
                checkboxInput('cbxAddRec_Drawing_linear', label="Missing profile parts"),
                checkboxInput('cbxAddLines_Drawing_linear', label="Lines"),
                actionButton('btnErraseLine_Drawing_linear', label="Back (Erase last line)"),
                checkboxInput('cbxAddScale_Drawing_linear', label="Scale"),
                checkboxInput('cbxAddCSPI_Drawing_linear', label="Circle sector preservation indicator"),
                checkboxInput('cbxAddDiameter_Drawing_linear', label="Rim diameter"),
                checkboxInput('cbxAddRuler_Drawing_linear', label="Measurements"),
                actionButton('btnErraseRuler_Drawing_linear', label="Back (Erase last point)"),
                checkboxInput('cbxAddVolume_Drawing_linear', label="Volume"),
                sliderInput('sldVolume_Drawing_photo', label="Boundary", value=c(0,100), min=0, max=100, step=1)
              )
            ),
            conditionalPanel(
              condition = "input.selDrawingMethod == 'Drawing_photo'",
              selectInput('selColour_Drawing_photo', label=NULL, c("Select colour mode" = "", "Colour mode" = "selColour_Drawing_photo_color", "Greyscale mode" = "selColour_Drawing_photo_greyscale"), selected = "selColour_Drawing_photo_color"),
              checkboxInput('cbxAddLight_Drawing_photo', label=strong('Additional light'), value=FALSE),
              checkboxInput('cbxFeatures_Drawing_photo', label=strong('Additional features'), value=FALSE),
              conditionalPanel(
                condition = "input.cbxFeatures_Drawing_photo",
                checkboxInput('cbxAddMask_Drawing_photo', label="Mask"),
                checkboxInput('cbxAddRim_Drawing_photo', label="Rim"),
                checkboxInput('cbxAddBase_Drawing_photo', label="Base"),
                checkboxInput('cbxAddRec_Drawing_photo', label="Missing profile parts"),
                checkboxInput('cbxAddLines_Drawing_photo', label="Lines"),
                actionButton('btnErraseLine_Drawing_photo', label="Back (Erase last line)"),
                checkboxInput('cbxAddScale_Drawing_photo', label="Scale"),
                checkboxInput('cbxAddCSPI_Drawing_photo', label="Circle sector preservation indicator"),
                checkboxInput('cbxAddDiameter_Drawing_photo', label="Diameter")
              ),
              checkboxInput('cbxFragmentRotTrans_Drawing_photo', label=strong("Fragment position"), value=FALSE),
              conditionalPanel(
                condition = "input.cbxFragmentRotTrans_Drawing_photo",
                sliderInput('sldRot_Drawing_photo', label="Rotation", value=0, min=-90, max=90, step=1),
                sliderInput('sldTrans_Drawing_photo', label="Translation", value=0, min=-90, max=90, step=1)
              ),
              checkboxInput('cbxAdvancedLight_Draw_photographic', label=strong("Advanced light"), value=FALSE),
              conditionalPanel(
                condition = "input.cbxAdvancedLight_Draw_photographic",
                sliderInput('sldLightLat_Drawing_photo', label="Latitude", value=45, min=-90, max=90, step=1),
                sliderInput('sldLightLon_Drawing_photo', label="Longitude", value=-45, min=-90, max=90, step=1),
                sliderInput('sldLightDis_Drawing_photo', label="Distance", value=30, min=5, max=300, step=1)
              )
            ),
            conditionalPanel(
              condition = "input.selDrawingMethod == 'Drawing_DL'",
              checkboxInput('cbxFeatures_Drawing_DL', label=strong('Additional features'), value=FALSE),
              conditionalPanel(
                condition = "input.cbxFeatures_Drawing_DL",
                checkboxInput('cbxAddMask_Drawing_DL', label="Mask"),
                checkboxInput('cbxAddRim_Drawing_DL', label="Rim"),
                checkboxInput('cbxAddBase_Drawing_DL', label="Base"),
                checkboxInput('cbxAddRec_Drawing_DL', label="Missing profile parts"),
                checkboxInput('cbxAddLines_Drawing_DL', label="Lines"),
                actionButton('btnErraseLine_Drawing_DL', label="Back (Erase last line)"),
                checkboxInput('cbxAddScale_Drawing_DL', label="Scale"),
                checkboxInput('cbxAddCSPI_Drawing_DL', label="Circle sector preservation indicator"),
                checkboxInput('cbxAddDiameter_Drawing_DL', label="Diameter")
              ),
              checkboxInput('cbxFragmentRotTrans_Drawing_DL', label=strong("Fragment position"), value=FALSE),
              conditionalPanel(
                condition = "input.cbxFragmentRotTrans_Drawing_DL",
                sliderInput('sldRot_Drawing_DL', label="Rotation", value=0, min=-90, max=90, step=1),
                sliderInput('sldTrans_Drawing_DL', label="Translation", value=0, min=-90, max=90, step=1)
              ),
              checkboxInput('cbxAdvancedLight_Drawing_DL', label=strong("Advanced light"), value=FALSE),
              conditionalPanel(
                condition = "input.cbxAdvancedLight_Drawing_DL",
                sliderInput('sldLightLat_Drawing_DL', label="Latitude", value=45, min=-90, max=90, step=1),
                sliderInput('sldLightLon_Drawing_DL', label="Longitude", value=-45, min=-90, max=90, step=1),
                sliderInput('sldLightDis_Drawing_DL', label="Distance", value=30, min=5, max=300, step=1)
              ),
              checkboxInput('cbxAdvancedShading_Drawing_DL', label=strong("Advanced shading"), value=FALSE),
              conditionalPanel(
                condition = "input.cbxAdvancedShading_Drawing_DL",
                checkboxInput('cbxFilt_Drawing_DL', label="Filter", value=F),
                sliderInput('sldLims_Drawing_DL', label="Percentiles", value=c(0,100), min=0, max=100, step=1),
                sliderInput('sldTresh1_Drawing_DL', label="Boundaries 1", value=c(-1,-1), min=-1, max=1, step=0.1),
                sliderInput('sldTresh2_Drawing_DL', label="Boundaries 2", value=c(-1,-1), min=-1, max=1, step=0.1),
                sliderInput('sldTresh3_Drawing_DL', label="Boundaries 3", value=c(-1,-1), min=-1, max=1, step=0.1),
                sliderInput('sldProc_Drawing_DL', label="Percentage", value=100, min=0, max=100, step=5)
              )
            ),
            conditionalPanel(
              condition = "input.selDrawingMethod == 'Drawing_AO'",
              checkboxInput('cbxFeatures_Drawing_AO', label=strong('Additional features'), value=FALSE),
              conditionalPanel(
                condition = "input.cbxFeatures_Drawing_AO",
                checkboxInput('cbxAddMask_Drawing_AO', label="Mask"),
                checkboxInput('cbxAddRim_Drawing_AO', label="Rim"),
                checkboxInput('cbxAddBase_Drawing_AO', label="Base"),
                checkboxInput('cbxAddRec_Drawing_AO', label="Missing profile parts"),
                checkboxInput('cbxAddLines_Drawing_AO', label="Lines"),
                actionButton('btnErraseLine_Drawing_AO', label="Back (Erase last line)"),
                checkboxInput('cbxAddScale_Drawing_AO', label="Scale"),
                checkboxInput('cbxAddCSPI_Drawing_AO', label="Circle sector preservation indicator"),
                checkboxInput('cbxAddDiameter_Drawing_AO', label="Diameter")
              ),
              checkboxInput('cbxFragmentRotTrans_Drawing_AO', label=strong("Fragment position"), value=FALSE),
              conditionalPanel(
                condition = "input.cbxFragmentRotTrans_Drawing_AO",
                sliderInput('sldRot_Drawing_AO', label="Rotation", value=0, min=-90, max=90, step=1),
                sliderInput('sldTrans_Drawing_AO', label="Translation", value=0, min=-90, max=90, step=1)
              ),
              checkboxInput('cbxAdvancedShading_Drawing_AO', label=strong("Advanced shading"), value=FALSE),
              conditionalPanel(
                condition = "input.cbxAdvancedShading_Drawing_AO",
                checkboxInput('cbxFilt_Drawing_AO', label="Filter", value=F),
                sliderInput('sldLims_Drawing_AO', label="Percentiles", value=c(5,95), min=0, max=100, step=1),
                sliderInput('sldTresh1_Drawing_AO', label="Boundaries 1", value=c(0,1), min=0, max=1, step=0.005),
                sliderInput('sldTresh2_Drawing_AO', label="Boundaries 2", value=c(0,1), min=0, max=1, step=0.005),
                sliderInput('sldTresh3_Drawing_AO', label="Boundaries 3", value=c(0,1), min=0, max=1, step=0.005),
                sliderInput('sldProc_Drawing_AO', label="Percentage", value=100, min=0, max=100, step=5)
              )
            ),
            conditionalPanel(
              condition = "input.selDrawingMethod == 'Drawing_AODL'",
              checkboxInput('cbxFeatures_Drawing_AODL', label=strong('Additional features'), value=FALSE),
              conditionalPanel(
                condition = "input.cbxFeatures_Drawing_AODL",
                checkboxInput('cbxAddMask_Drawing_AODL', label="Mask"),
                checkboxInput('cbxAddRim_Drawing_AODL', label="Rim"),
                checkboxInput('cbxAddRec_Drawing_AODL', label="Missing profile parts"),
                checkboxInput('cbxAddRecDown_Drawing_AODL', label="Reconstruction down"),
                checkboxInput('cbxAddLines_Drawing_AODL', label="Lines"),
                actionButton('btnErraseLine_Drawing_AODL', label="Back (Erase last line)"),
                checkboxInput('cbxAddScale_Drawing_AODL', label="Scale"),
                checkboxInput('cbxAddCSPI_Drawing_AODL', label="Circle sector preservation indicator"),
                checkboxInput('cbxAddDiameter_Drawing_AODL', label="Diameter")
              ),
              checkboxInput('cbxFragmentRotTrans_Drawing_AODL', label=strong("Fragment position"), value=FALSE),
              conditionalPanel(
                condition = "input.cbxFragmentRotTrans_Drawing_AODL",
                sliderInput('sldRot_Drawing_AODL', label="Rotation", value=0, min=-90, max=90, step=1),
                sliderInput('sldTrans_Drawing_AODL', label="Translation", value=0, min=-90, max=90, step=1)
              ),
              checkboxInput('cbxAdvancedLight_Drawing_AODL_DL', label=strong("Advanced light (DL)"), value=FALSE),
              conditionalPanel(
                condition = "input.cbxAdvancedLight_Drawing_AODL_DL",
                sliderInput('sldLightLat_Drawing_AODL_DL', label="Latitude", value=45, min=-90, max=90, step=1),
                sliderInput('sldLightLon_Drawing_AODL_DL', label="Longitude", value=-45, min=-90, max=90, step=1),
                sliderInput('sldLightDis_Drawing_AODL_DL', label="Distance", value=30, min=5, max=300, step=1)
              ),
              checkboxInput('cbxAdvancedShading_Drawing_AODL_DL', label=strong("Advanced shading (DL)"), value=FALSE),
              conditionalPanel(
                condition = "input.cbxAdvancedShading_Drawing_AODL_DL",
                sliderInput('sldLims_Drawing_AODL_DL', label="Percentiles", value=c(0,100), min=0, max=100, step=1),
                sliderInput('sldTresh1_Drawing_AODL_DL', label="Boundaries 1", value=c(-1,-1), min=-1, max=1, step=0.1),
                sliderInput('sldTresh2_Drawing_AODL_DL', label="Boundaries 2", value=c(-1,-1), min=-1, max=1, step=0.1),
                sliderInput('sldTresh3_Drawing_AODL_DL', label="Boundaries 3", value=c(-1,-1), min=-1, max=1, step=0.1),
                sliderInput('sldProc_Drawing_AODL_DL', label="Percentage", value=100, min=0, max=100, step=5)
              ),
              checkboxInput('cbxAdvancedShading_Drawing_AODL_AO', label=strong("Advanced shading (AO)"), value=FALSE),
              conditionalPanel(
                condition = "input.cbxAdvancedShading_Drawing_AODL_AO",
                sliderInput('sldLims_Drawing_AODL_AO', label="Percentiles", value=c(5,95), min=0, max=100, step=1),
                sliderInput('sldTresh1_Drawing_AODL_AO', label="Boundaries 1", value=c(0,1), min=0, max=1, step=0.005),
                sliderInput('sldTresh2_Drawing_AODL_AO', label="Boundaries 2", value=c(0,1), min=0, max=1, step=0.005),
                sliderInput('sldTresh3_Drawing_AODL_AO', label="Boundaries 3", value=c(0,1), min=0, max=1, step=0.005),
                sliderInput('sldProc_Drawing_AODL_AO', label="Percentage", value=100, min=0, max=100, step=5)
              )
            ),
            conditionalPanel(
              condition = "input.selDrawingMethod == 'Drawing_symmetry'",
              selectInput('selMethod_Drawing_symmetry', label=NULL, c("Select drawing quality method" = "", "According to rotation axis" = "selMethod_Drawing_symmetry_axis", "According to ideal vessel" = "selMethod_Drawing_symmetry_profile"), selected="selMethod_Drawing_symmetry_profile"),
              checkboxInput('cbxFeatures_Drawing_symmetry', label=strong('Additional features'), value=FALSE),
              conditionalPanel(
                condition = "input.cbxFeatures_Drawing_symmetry",
                checkboxInput('cbxAddMask_Drawing_symmetry', label="Mask"),
                checkboxInput('cbxAddRim_Drawing_symmetry', label="Rim"),
                checkboxInput('cbxAddBase_Drawing_symmetry', label="Base"),
                checkboxInput('cbxAddRec_Drawing_symmetry', label="Missing profile parts"),
                checkboxInput('cbxAddLines_Drawing_symmetry', label="Lines"),
                actionButton('btnErraseLine_Drawing_symmetry', label="Back (Erase last line)"),
                checkboxInput('cbxAddScale_Drawing_symmetry', label="Scale"),
                checkboxInput('cbxAddCSPI_Drawing_symmetry', label="Circle sector preservation indicator"),
                checkboxInput('cbxAddDiameter_Drawing_symmetry', label="Diameter")
              ),
              checkboxInput('cbxFragmentRotTrans_Drawing_symmetry', label=strong("Fragment position"), value=FALSE),
              conditionalPanel(
                condition = "input.cbxFragmentRotTrans_Drawing_symmetry",
                sliderInput('sldRot_Drawing_symmetry', label="Rotation", value=0, min=-90, max=90, step=1),
                sliderInput('sldTrans_Drawing_symmetry', label="Translation", value=0, min=-90, max=90, step=1)
              )
            ),
            actionButton('btnDraw', 'Draw'),
            br(),
            br(),
            br(),
            strong('4. Save current drawing'),
            textInput('saveDrawingPath', label=NULL, value="C:/DACORD/images/"),
            actionButton('btnSave_Drawing', 'Save'),
            width=width.adj[1]
          ),
          mainPanel(
            plotOutput('plotDrawing',click="plotDrawing_click",dblclick="plotDrawing_dblclick"),
            plotOutput('plotDrawingHist'),
            width=width.adj[2]
          )
        )
      )
    )
    ,

    tabPanel("3D reconstruction",
       fluidPage(
         sidebarLayout(
           sidebarPanel(
             strong('1. Load file'),
             fileInput('fileOrientedMesh2', label=NULL, accept=c('.ply')),

             strong('2. Profile type selection'),
             selectInput('selProfileMethod2', label=NULL,
                         c("Select profile method" = "",
                           "Whole envelope profile" = "Profile_envelop",
                           "In the middle of the fragment" = "Profile_middle",
                           "Longest-preserved profile" = "Profile_longest",
                           "Arbitrarily selected profile" = "Profile_arbitrary")
             ),
             actionButton('btnGetProfile2', 'Get profile'),
             br(),
             br(),
             br(),
             strong('3. 3D reconstruction selection'),
             selectInput('selReconst3dMethod', label=NULL,
                         c("Select 3D reconstruction method" = "",
                           "Plastic model" = "Model_plastic",
                           "Sliced model" = "Model_sliced",
                           "Wireframe model" = "Model_wireframe",
                           "Point cloud" = "Model_pointcloud"
                         )
             ),
             checkboxInput('cbxAddModel_Reconst3d', 'Add model'),
             actionButton('btnReconst3d', 'Reconstruct'),
             br(),
             br(),
             br(),
             strong('4. Save 3D reconstruction'),
             textInput('saveReconstructionPath', label=NULL, value="C:/DACORD/images/"),
             actionButton('btnSave_Reconstruction', 'Save'),
             width=width.adj[1]),
           mainPanel(
             rglwidgetOutput('plotReconst3d', width = mainPanel.adj[1], height = mainPanel.adj[2]),
             plotOutput('plot_Reconst'),
             width=width.adj[2]
           )
         )
       )
    )
  )
)
