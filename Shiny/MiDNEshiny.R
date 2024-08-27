#' MiDNE Shiny app
#'
#' @import shiny
#' @import shinydashboard
#' @import shinyalert
#' @import shinyFiles
#' @import shinyMatrix
#' @import shinyjs
#' @import shinycssloaders
#' @import readr
#' @import glue
#' @import conflicted
#' @import fpc
#' @import parallel
#' @import plotly
#' @import ggplot2
#' @import gprofiler2
#' @import visNetwork
#' @import tidyverse
#' @import dendextend
#' @import crosstalk
#' @import magrittr
#' @import dplyr
#' @import igraph
#' @import vroom
#' @import zip
#' @import MoNETA
#' @importFrom utils data read.csv
#' @importFrom stats as.dendrogram dist kmeans hclust
#' @importFrom graphics abline par
#' @param MAXreq shiny max request size
#' @return Shiny app
#' @export


source('MiDNE/R/network_inference.R')
source('MiDNE/R/gen_sim_mat_MH.R')
source('MiDNE/R/myfunctions.R')


MiDNEshiny = function(MAXreq = 10000) {
  options(shiny.maxRequestSize = MAXreq * 1024^2)
  shiny::shinyApp(ui, server)
  
}


netSummary <- function(network){
  
  nodes <- c(network[[1]], network[[2]])
  degree <- table(nodes) %>%
    dplyr::as_tibble(.) %>%
    dplyr::arrange(dplyr::desc(n))
  top10 <- c(degree[c(1:10), 'nodes'])
  netSummary_list <- list( 'dim'= dim(network),
                           'edges'= dim(network)[1],
                           'nodes'= length(unique(nodes)),
                           'maxConnections'= degree[[1, 'n']]#,
                           #'maxDegreeNodes'= top10$nodes
  )
  return(netSummary_list)
}


ui <- shinydashboard::dashboardPage(
  
  shinydashboard::dashboardHeader(title = shiny::span("MiDNE ", style = "color: white; font-size: 28px")),
  
  shinydashboard::dashboardSidebar(
    shinyjs::useShinyjs(),
    shinydashboard::sidebarMenu(id = 'tabs',
                                shiny::tags$head(shiny::tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")),
                                shiny::tags$head(shiny::tags$style(".inactiveLink {
                                                                         pointer-events: none;
                                                                        cursor: default;
                                                                        }")),
                                style = "position: fixed; height: 90vh; overflow-y: auto;",
                                shinydashboard::menuItem("Data", tabName = "mat_tab", icon = shiny::icon("table"),
                                                         shinydashboard::menuItem("Loading", tabName = "mat_sub_1"),
                                                         shinydashboard::menuItem("Pre-processing", tabName = "mat_sub_2")),
                                shinydashboard::menuItem("Network", tabName = "net_tab", icon = shiny::icon("circle-nodes"),
                                                         shinydashboard::menuItem("Network inference", tabName = "net_sub_1"),
                                                         shinydashboard::menuItem("Network Filtering", tabName =  "net_sub_2")),
                                shinydashboard::menuItem("RWR", tabName = "rwr_tab", icon = shiny::icon("person-running"),
                                                         shinydashboard::menuItem("Loading", tabName = "rwr_tab_1"),
                                                         shinydashboard::menuItem("Multiplex network", tabName = "rwr_tab_2"),
                                                         shinydashboard::menuItem("RWR Parameters", tabName = "rwr_tab_3")),
                                shinydashboard::menuItem("Dimensionality Reduction", tabName = "dr_tab", icon = shiny::icon("filter"),
                                                         shinydashboard::menuItem("Loading", tabName = "dr_tab_1"),
                                                         shinydashboard::menuItem("DR", tabName = "dr_tab_2")
                                ),
                                shinydashboard::menuItem("Clustering", tabName = "cl_tab", icon = shiny::icon("layer-group")),
                                shinydashboard::menuItem("Enrichment Analysis", icon = shiny::icon("fingerprint"),
                                                         shinydashboard::menuItem("Cluster2Pathway", tabName = "clupath_tab"),
                                                         shinydashboard::menuItem("Pathway2Clusters", tabName = "pathclu_tab")),
                                shinydashboard::menuItem("Drug Discovery", tabName = "dd_tab", icon = shiny::icon("capsules"))
    )
  ),
  
  shinydashboard::dashboardBody(
    shinydashboard::tabItems(
      shinydashboard::tabItem(tabName = "mat_sub_1",
                              shiny::fluidRow(
                                shiny::column(width = 6,
                                              shiny::div(style = "display: inline-block;",
                                                         shiny::uiOutput('jump2P1.1')
                                              ),
                                              shiny::div(style = "display: inline-block;",
                                                         shiny::uiOutput('jump2P1.2')
                                              )
                                ),
                                shiny::column(width = 6,
                                              shiny::uiOutput('jump2P2'), align = 'right'
                                )
                              ),
                              
                              shiny::fluidRow(
                                
                                shiny::column(width = 6,
                                              shiny::fluidRow(
                                                shiny::column(width = 12,
                                                              shiny::tags$hr(),
                                                              shinydashboard::box(
                                                                title = shiny::h2(shiny::span("Upload omics matrices", style = "font-weight: bold")),
                                                                width = 12, status = "primary", solidHeader = FALSE,
                                                                shiny::radioButtons(inputId = 'omics_example_opt',
                                                                                    label = shiny::h4(shiny::span('Load the example dataset', style = "font-weight: bold")),
                                                                                    choices = c('Yes', 'No'), selected = 'No'),
                                                                shiny::conditionalPanel(
                                                                  condition = 'input.omics_example_opt == "Yes"',
                                                                  shiny::checkboxGroupInput(inputId = 'omics_example_files',
                                                                                            label = shiny::h4(shiny::span('Select one or more omics matrices', style = "font-weight: bold")),
                                                                                            choices = c("BRCA_expr_HiSeq", "BRCA_proteome_CDAP",
                                                                                                        "BRCA_Methylation_Meth450", "BRCA_SCNA")
                                                                                            )
                                                                ),
                                                                shiny::conditionalPanel(
                                                                  condition = 'input.omics_example_opt == "No"',
                                                                  shiny::fileInput(inputId = "omics_files", label = shiny::h4(shiny::span('Select one or more files', style = "font-weight: bold")),
                                                                                   multiple = TRUE),
                                                                  shiny::uiOutput('omics_names'),
                                                                ),
                                                                shiny::tags$hr(),
                                                                shiny::actionButton('load_mat_button', 'Load')
                                                              )
                                                ),
                                                
                                                shiny::column(width = 12,
                                                              shiny::htmlOutput(outputId = "omics_sum_info")
                                                ),
                                                
                                                shiny::column(width = 12,
                                                              shinydashboard::box(
                                                                width = 12, status = "primary", solidHeader = FALSE,
                                                                title = shiny::h2(shiny::span("Upload drug files", style = "font-weight: bold")),
                                                                
                                                                shiny::radioButtons(inputId = 'drug_example_opt',
                                                                                    label = shiny::h4(shiny::span('Load the example dataset', style = "font-weight: bold")),
                                                                                    choices = c('Yes', 'No'), selected = 'No'),
                                                                shiny::conditionalPanel(
                                                                  condition = 'drug_example_opt == "No"',
                                                                  shiny::fileInput(inputId = "drug_files",
                                                                                   label = shiny::h4(shiny::span("Select one or more files", style = "font-weight: bold")),
                                                                                   multiple = TRUE),
                                                                  shiny::uiOutput('drug_names')
                                                                ),
                                                                
                                                                shiny::tags$hr(),
                                                                shiny::actionButton('load_drug_mat_button', 'Load')
                                                              )
                                                ),
                                                
                                                shiny::column(width = 12,
                                                              shiny::htmlOutput(outputId = "drug_sum_info")
                                                ),
                                                
                                                shiny::column(width = 12,
                                                              shinydashboard::box(
                                                                width = 12, status = "primary", solidHeader = FALSE,
                                                                title = shiny::h2(shiny::span("Upload annotation file", style = "font-weight: bold")),
                                                                shiny::radioButtons(inputId = 'anno_example_opt',
                                                                                    label = shiny::h4(shiny::span('Load the example annotation file', style = "font-weight: bold")),
                                                                                    choices = c('Yes', 'No'), selected = 'No'),
                                                                shiny::conditionalPanel(
                                                                  condition = 'input.anno_example_opt == "No"',
                                                                  shiny::fileInput(inputId = "anno_file", label = shiny::h4(shiny::span('Select a file', style = "font-weight: bold")),
                                                                                   multiple = FALSE)
                                                                ),
                                                                shiny::tags$hr(),
                                                                shiny::actionButton('load_anno_button', 'Load')
                                                              )
                                                ),
                                                
                                                shiny::column(width = 12,
                                                              shiny::htmlOutput(outputId = "anno_info")
                                                ),
                                                
                                                shiny::column(width = 12,
                                                              shiny::htmlOutput(outputId = "anno_info_extra")
                                                )
                                              )
                                ),
                                shiny::column(width = 6,
                                              shiny::tags$hr(),
                                              shinydashboard::infoBox(
                                                title = shiny::h3(shiny::span('Information', style = "font-weight: bold")),
                                                value = shiny::uiOutput('info_box_1'),
                                                subtitle = NULL,
                                                icon = shiny::icon("info"), color = "aqua", width = NULL,
                                                fill = FALSE)
                                )
                              )
      ),
      shinydashboard::tabItem(tabName = "mat_sub_2",
                              shiny::fluidRow(
                                shiny::column(width = 6,
                                              shiny::uiOutput('back2P1')),
                                shiny::column(width = 6,
                                              shiny::uiOutput('jump2P3'),  align = 'right'
                                )
                              ),
                              shiny::fluidRow(
                                shiny::column(width = 6,
                                              shiny::fluidRow(
                                                shiny::tags$hr(),
                                                shiny::column(width = 12,
                                                              #shiny::h2(shiny::span("Omics matrices Pre-processing", style = "font-weight: bold")),
                                                              shiny::uiOutput("process_omics_mat"),
                                                ),
                                                shiny::column(width = 12,
                                                              shiny::uiOutput('pro_matrices_sum_info')
                                                ),
                                                shiny::column(width = 12,
                                                              shiny::uiOutput('intersection')
                                                ),
                                                shiny::column(width = 12,
                                                              shiny::uiOutput('intersection_info')
                                                )
                                              )
                                ),
                                shiny::column(width = 6,
                                              shiny::tags$hr(),
                                              shiny::fluidRow(
                                                shiny::column(width = 12,
                                                              shinydashboard::infoBox(
                                                                title = shiny::h3(shiny::span('Information', style = "font-weight: bold")),
                                                                value = shiny::uiOutput('info_box_2'),
                                                                subtitle = NULL,
                                                                icon = shiny::icon("info"), color = "aqua", width = NULL,
                                                                fill = FALSE)
                                                ),
                                                shiny::column(width = 12,
                                                              shiny::uiOutput('download_proc_mat_box')
                                                )
                                              )
                                )
                              )
      ),
      
      shinydashboard::tabItem(tabName = 'net_sub_1',
                              shiny::fluidRow(
                                shiny::column(width = 6,
                                              shiny::uiOutput('back2P2')),
                                shiny::column(width = 6,
                                              shiny::uiOutput('jump2P4'),  align = 'right'
                                )
                              ),
                              shiny::fluidRow(
                                shiny::column(width = 6,
                                              shiny::tags$hr(),
                                              shiny::h2(shiny::span("Omics Network Inference", style = "font-weight: bold")),
                                              shiny::fluidRow(
                                                shiny::column(width = 12, 
                                                              shiny::uiOutput("omics_net_arguments")
                                                ),
                                                shiny::column(width = 12, 
                                                              shiny::uiOutput("omics_net_info")
                                                ),
                                                shiny::column(width = 12, 
                                                              shiny::h2(shiny::span("Drug Network Inference", style = "font-weight: bold")),
                                                              shiny::uiOutput("drug_net_arguments")
                                                ),
                                                shiny::column(width = 12, 
                                                              shiny::uiOutput("drug_net_info")
                                                )
                                              )
                                ),
                                shiny::column(width = 6,
                                              shiny::tags$hr(),
                                              shiny::fluidRow(
                                                shiny::column(width = 12,
                                                              shinydashboard::infoBox(
                                                                title = shiny::h3(shiny::span('Information', style = "font-weight: bold")),
                                                                value = shiny::uiOutput('info_box_3'),
                                                                subtitle = NULL,
                                                                icon = shiny::icon("info"), color = "aqua", width = NULL,
                                                                fill = FALSE)
                                                              
                                                ),
                                                shiny::column(width = 12,
                                                              shiny::uiOutput("download_net_box")
                                                )
                                              )
                                )
                              )
      ),
      shinydashboard::tabItem(tabName = 'net_sub_2',
                              shiny::fluidRow(
                                shiny::column(width = 6,
                                              shiny::uiOutput('back2P3')),
                                shiny::column(width = 6,
                                              shiny::div(style = "display: inline-block;",
                                                         shiny::uiOutput('jump2P5')
                                              ),
                                              shiny::div(style = "display: inline-block;",
                                                         shiny::uiOutput('jump2P1.2_from4')
                                              ), align = 'right'
                                )
                              ),
                              shiny::fluidRow(
                                shiny::column(width = 6,
                                              shiny::tags$hr(),
                                              shiny::h2(shiny::span("Omics Network Filtering", style = "font-weight: bold")),
                                              shiny::fluidRow(
                                                shiny::uiOutput("plot_net_box")
                                              )
                                ),
                                shiny::column(width = 6,
                                              shiny::tags$hr(),
                                              shiny::column(width = 12,
                                                            shinydashboard::infoBox(
                                                              title = shiny::h3(shiny::span('Information', style = "font-weight: bold")),
                                                              value = shiny::uiOutput('info_box_4'),
                                                              subtitle = NULL,
                                                              icon = shiny::icon("info"), color = "aqua", width = NULL,
                                                              fill = FALSE),
                                                            shiny::uiOutput("download_fnet_box")
                                              )
                                )
                              )
      ),
      
      shinydashboard::tabItem(tabName = "rwr_tab_1",
                              
                              shiny::fluidRow(
                                shiny::column(width = 6,
                                              shiny::uiOutput('back2P1_from1.1'),
                                              shiny::uiOutput('back2P4_from1.1')
                                ),
                                shiny::column(width = 6,
                                              shiny::uiOutput('jump2P5.1'),  align = 'right'
                                )
                              ),
                              
                              shiny::fluidRow(
                                shiny::column(width = 6,
                                              shiny::fluidRow(
                                                shiny::column(width = 12,
                                                              shiny::tags$hr(),
                                                              shinydashboard::box(
                                                                title = shiny::h2(shiny::span("Upload omics networks", style = "font-weight: bold")),
                                                                width = 12, status = "primary", solidHeader = FALSE,
                                                                
                                                                shiny::radioButtons(inputId = 'omicsNet_example_opt',
                                                                                    label = shiny::h4(shiny::span('Load the example dataset', style = "font-weight: bold")),
                                                                                    choices = c('Yes', 'No'), selected = 'No'),
                                                                shiny::conditionalPanel(
                                                                  condition = 'input.omicsNet_example_opt == "Yes"',
                                                                  shiny::checkboxGroupInput(inputId = 'omicsNet_example_files',
                                                                                            label = shiny::h4(shiny::span('Select one or more omic networks', style = "font-weight: bold")),
                                                                                            choices = c("BRCA_filt_coexpr_network", "BRCA_filt_prot_net",
                                                                                                        "BRCA_filtered_coMethy_network", "BRCA_filt_coamp_network.RDS",
                                                                                                        'BRCA_filt_codel_network.RDS'))
                                                                ),
                                                                shiny::conditionalPanel(
                                                                  condition = 'input.omicsNet_example_opt == "No"',
                                                                  shiny::fileInput("omics_net_files", label = shiny::h4(shiny::span( 'Select one or more files', style = "font-weight: bold")),
                                                                                   multiple = TRUE,
                                                                                   accept = c("text/csv", '.RDS', "text/comma-separated-values,text/plain", ".csv")
                                                                  ),
                                                                  shiny::uiOutput('omics_net_names')
                                                                ),
                                                                
                                                                shiny::hr(),
                                                                shiny::actionButton('load_omics_net_button', 'Load')
                                                              )
                                                ),
                                                shiny::column(width = 12,
                                                              shiny::htmlOutput(outputId = "omics_net_sum_info")
                                                ),
                                                
                                                shiny::column(width = 12,
                                                              shinydashboard::box(
                                                                title = shiny::h2(shiny::span("Upload drug networks", style = "font-weight: bold")),
                                                                width = 12, status = "primary", solidHeader = FALSE,
                                                                
                                                                shiny::radioButtons(inputId = 'drugNet_example_opt',
                                                                                    label = shiny::h4(shiny::span('Load the example dataset', style = "font-weight: bold")),
                                                                                    choices = c('Yes', 'No'), selected = 'No'),
                                                                shiny::conditionalPanel(
                                                                  condition = 'input.drugNet_example_opt == "No"',
                                                                  shiny::fileInput("drug_net_files", label = shiny::h4(shiny::span( 'Select one or more files', style = "font-weight: bold")),
                                                                                   multiple = TRUE,
                                                                                   accept = c("text/csv", '.RDS', "text/comma-separated-values,text/plain", ".csv")
                                                                  ),
                                                                  shiny::uiOutput('drug_net_names')
                                                                ),
                                                                shiny::actionButton('load_drug_net_button', 'Load')
                                                              )
                                                ),
                                                shiny::column(width = 12,
                                                              shiny::htmlOutput(outputId = "drug_net_sum_info")
                                                ),
                                                
                                                shiny::column(width = 12,
                                                              shinydashboard::box(
                                                                title = shiny::h2(shiny::span("Upload bipartite network", style = "font-weight: bold")),
                                                                width = 12, status = "primary", solidHeader = FALSE,
                                                                
                                                                shiny::radioButtons(inputId = 'bNet_example_opt',
                                                                                    label = shiny::h4(shiny::span('Load the example dataset', style = "font-weight: bold")),
                                                                                    choices = c('Yes', 'No'), selected = 'No'),
                                                                shiny::conditionalPanel(
                                                                  condition = 'input.bNet_example_opt == "No"',
                                                                  shiny::fileInput("bnet_file", label = shiny::h4(shiny::span( 'Select a file', style = "font-weight: bold")),
                                                                                   multiple = FALSE,
                                                                                   accept = c("text/csv", '.RDS', "text/comma-separated-values,text/plain", ".csv")
                                                                  )
                                                                ),
                                                                
                                                                shiny::actionButton('load_bnet_button', 'Load')
                                                              )
                                                ),
                                                shiny::column(width = 12,
                                                              shiny::htmlOutput(outputId = "bnet_sum_info")
                                                ),
                                                
                                                shiny::column(width = 12,
                                                              shinydashboard::box(
                                                                title = shiny::h2(shiny::span("Upload annotation file", style = "font-weight: bold")),
                                                                width = 12, status = "primary", solidHeader = FALSE,
                                                                
                                                                shiny::radioButtons(inputId = 'rwr_anno_example_opt',
                                                                                    label = shiny::h4(shiny::span('Load the example dataset', style = "font-weight: bold")),
                                                                                    choices = c('Yes', 'No'), selected = 'No'),
                                                                shiny::conditionalPanel(
                                                                  condition = 'input.rwr_anno_example_opt == "No"',
                                                                  shiny::fileInput(inputId = "anno_file1",
                                                                                   label = shiny::h4(shiny::span( 'Select a file', style = "font-weight: bold")),
                                                                                   multiple = FALSE)
                                                                ),
                                                                
                                                                shiny::actionButton('load_anno_button1', 'Load')
                                                              ),
                                                              shiny::htmlOutput(outputId = "anno_info1")
                                                )
                                                
                                              )
                                ),
                                shiny::column(width = 6,
                                              shiny::tags$hr(),
                                              shinydashboard::infoBox(
                                                title = shiny::h3(shiny::span('Information', style = "font-weight: bold")),
                                                value = shiny::uiOutput('info_box_5'),
                                                subtitle = NULL,
                                                icon = shiny::icon("info"), color = "aqua", width = NULL,
                                                fill = FALSE),
                                              shiny::htmlOutput(outputId = "anno_info_extra1"),
                                              shiny::htmlOutput(outputId = "anno_info_extra2")
                                )
                                
                              )
      ),
      
      shinydashboard::tabItem(tabName = "rwr_tab_2",
                              shiny::fluidRow(
                                shiny::column(width = 6,
                                              shiny::uiOutput('back2P4'),
                                              shiny::uiOutput('back2P1.1')),
                                shiny::column(width = 6,
                                              shiny::uiOutput('jump2P6'),  align = 'right'
                                )
                              ),
                              shiny::fluidRow(
                                shiny::column(width = 6,
                                              shiny::fluidRow(
                                                shiny::column(width=12,
                                                              shiny::tags$hr(),
                                                              shinydashboard::box(width = 12, status = "primary", solidHeader = FALSE,
                                                                                  title =  shiny::h2(shiny::span("Construction of Multiplex Omics Network", style = "font-weight: bold")),
                                                                                  shiny::radioButtons(inputId = 'weightMultiplex',
                                                                                                      label = shiny::h4(shiny::span('Do you want to create a weighted multiplex omics network?', style = "font-weight: bold")),
                                                                                                      choices = c("YES", "NO"), selected = "NO"),
                                                                                  shiny::tags$hr(),
                                                                                  shiny::conditionalPanel(
                                                                                    condition = 'input.weightMultiplex == "YES"',
                                                                                    shiny::radioButtons(inputId = 'pruneMultiplex',
                                                                                                        label = shiny::h4(shiny::span('Do you want to prune the multiplex omics network?', style = "font-weight: bold")),
                                                                                                        choices = c('YES', 'NO'), selected = 'NO'),
                                                                                    shiny::conditionalPanel(
                                                                                      condition = 'input.pruneMultiplex == "YES"',
                                                                                      shiny::numericInput(inputId = 'pruneMultiplex_th',
                                                                                                          label = shiny::h4(shiny::span('Select a threshold to prune the multiplex omics network', style = "font-weight: bold")),
                                                                                                          min = 0, max = 100, value = 50)
                                                                                    ),
                                                                                    shiny::tags$hr()
                                                                                  ),
                                                                                  shiny::actionButton(inputId = 'gen_multiplex_btn', label = 'Submit'))
                                                ),
                                                shiny::column(width=12,
                                                              shiny::htmlOutput(outputId = "multi_net_sum_info")
                                                ),
                                                shiny::column(width=12,
                                                              shinydashboard::box(width = 12, status = "primary", solidHeader = FALSE,
                                                                                  title =  shiny::h2(shiny::span("Construction of Multiplex Drug Network", style = "font-weight: bold")),
                                                                                  shiny::radioButtons(inputId = 'drug_weightMultiplex',
                                                                                                      label = shiny::h4(shiny::span('Do you want to create a weighted multiplex drug network?', style = "font-weight: bold")),
                                                                                                      choices = c("YES", "NO"), selected = "NO"),
                                                                                  shiny::tags$hr(),
                                                                                  shiny::conditionalPanel(
                                                                                    condition = 'input.drug_weightMultiplex == "YES"',
                                                                                    shiny::radioButtons(inputId = 'drug_pruneMultiplex',
                                                                                                        label = shiny::h4(shiny::span('Do you want to prune the multiplex drug network?', style = "font-weight: bold")),
                                                                                                        choices = c('YES', 'NO'), selected = 'NO'),
                                                                                    shiny::conditionalPanel(
                                                                                      condition = 'input.drug_pruneMultiplex == "YES"',
                                                                                      shiny::numericInput(inputId = 'drug_pruneMultiplex_th',
                                                                                                          label = shiny::h4(shiny::span('Select a threshold to prune the multiplex drug network', style = "font-weight: bold")),
                                                                                                          min = 0, max = 100, value = 50)
                                                                                    ),
                                                                                    shiny::tags$hr()
                                                                                  ),
                                                                                  shiny::actionButton(inputId = 'gen_multiplex_drug_btn', label = 'Submit'))
                                                ),
                                                shiny::column(width=12,
                                                              shiny::htmlOutput(outputId = "multi_drug_net_sum_info")
                                                ),
                                              )
                                ),
                                shiny::column(width = 6,
                                              shiny::tags$hr(),
                                              shinydashboard::infoBox(
                                                title = shiny::h3(shiny::span('Information', style = "font-weight: bold")),
                                                value = shiny::uiOutput('info_box_6'),
                                                subtitle = NULL,
                                                icon = shiny::icon("info"), color = "aqua", width = NULL,
                                                fill = FALSE)
                                )
                              )
      ),
      shinydashboard::tabItem(tabName = "rwr_tab_3",
                              shiny::fluidRow(
                                shiny::column(width = 6,
                                              shiny::uiOutput('back2P5')),
                                shiny::column(width = 6,
                                              shiny::uiOutput('jump2P7'),  align = 'right'
                                )
                              ),
                              
                              shiny::fluidRow(
                                
                                shiny::column(width = 6,
                                              shiny::tags$hr(),
                                              shinydashboard::box(width=12,
                                                                  status = "primary", solidHeader = FALSE,
                                                                  title = shiny::h2(shiny::span("Random Walk with Restart", style = "font-weight: bold")),
                                                                  shinydashboard::tabBox(width=12,
                                                                                         title = shiny::column(width = 6, shiny::uiOutput('download_rwr_mat')),
                                                                                         shiny::tabPanel(title = shiny::h3(shiny::span('Restart', style = "font-weight: bold")),
                                                                                                         shiny::fluidRow(
                                                                                                           shiny::column(width = 12,
                                                                                                                         shiny::sliderInput('restart', shiny::h4(shiny::span('Select the restarting parameter (r)', style = "font-weight: bold")),
                                                                                                                                            min = 0, max = 1, value = 0.7, step = 0.01)
                                                                                                           )
                                                                                                         ),
                                                                                                         shiny::tags$hr(),
                                                                                                         shiny::fluidRow(
                                                                                                           shiny::column(width = 12,
                                                                                                                         shiny::radioButtons('tao_opt',
                                                                                                                                             label = shiny::h4(shiny::HTML('<span style:"font-family = LM Roman 10"><b>Restarting Probabilities (Omics &tau;) </b></span>')),
                                                                                                                                             choices = c('Custumize restarting probabilities per layer', 'Use default restarting probabilities'),
                                                                                                                                             selected = 'Use default restarting probabilities'),
                                                                                                                         shiny::uiOutput(outputId = "tauBIO")
                                                                                                           )
                                                                                                         ),
                                                                                                         shiny::tags$hr(),
                                                                                                         shiny::fluidRow(
                                                                                                           shiny::column(width = 12,
                                                                                                                         shiny::radioButtons('drug_tao_opt',
                                                                                                                                             label = shiny::h4(shiny::HTML('<span style:"font-family = LM Roman 10"><b>Restarting Probabilities (Drug &tau;) </b></span>')),
                                                                                                                                             choices = c('Custumize restarting probabilities per layer', 'Use default restarting probabilities'),
                                                                                                                                             selected = 'Use default restarting probabilities'),
                                                                                                                         shiny::uiOutput(outputId = "tauDRUGS")
                                                                                                           )
                                                                                                         )
                                                                                         ),
                                                                                         shiny::tabPanel(title = shiny::h3(shiny::span('Jump', style = "font-weight: bold")),
                                                                                                         shiny::fluidRow(
                                                                                                           shiny::column(width = 12,
                                                                                                                         shiny::sliderInput('omics_delta',
                                                                                                                                            shiny::h4(shiny::HTML('<span style:"font-family = LM Roman 10"><b>Select the transition parameter (Omics &delta;) </b></span>')),
                                                                                                                                            min = 0, max = 1, value = 0.5, step = 0.01),
                                                                                                                         shiny::uiOutput(outputId = "drug_delta_slider")
                                                                                                           )
                                                                                                         ),
                                                                                                         
                                                                                                         shiny::tags$hr(),
                                                                                                         shiny::fluidRow(
                                                                                                           shiny::column(width = 12,
                                                                                                                         shiny::h4(shiny::span('Omics Transition layer matrix', style = "font-weight: bold")),
                                                                                                                         shiny::radioButtons(inputId = 'bioInf_transition',
                                                                                                                                             label = 'Do you want to create a biologically informed layer transition matrix?',
                                                                                                                                             choices = c('YES', 'NO'), selected = 'NO'),
                                                                                                                         shiny::conditionalPanel(
                                                                                                                           condition = 'input.bioInf_transition == "NO"',
                                                                                                                           shiny::uiOutput('omics_trans_mat')
                                                                                                                         )
                                                                                                           ),
                                                                                                           shiny::column(width = 12,
                                                                                                                         shiny::h4(shiny::span('Drug Transition layer matrix', style = "font-weight: bold")),
                                                                                                                         shiny::uiOutput('drugs_trans_mat')
                                                                                                           )
                                                                                                         ),
                                                                                                         
                                                                                                         shiny::tags$hr(),
                                                                                                         shiny::fluidRow(
                                                                                                           shiny::column(width = 12,
                                                                                                                         shiny::h4(shiny::span('RWR-NF option', style = "font-weight: bold")),
                                                                                                                         shiny::radioButtons('omics_jump_neigh', label = 'Connect nodes of different layers by neighborhood (omics multiplex)',
                                                                                                                                             
                                                                                                                                             c('YES', 'NO'), selected = 'NO'),
                                                                                                                         shiny::conditionalPanel(
                                                                                                                           condition = 'input.omics_jump_neigh == "YES"',
                                                                                                                           shiny::tags$hr(),
                                                                                                                           shiny::radioButtons('omics_weight_jump_neigh',
                                                                                                                                               shiny::h4(shiny::span('Weight inter-layers edges', style = "font-weight: bold")),
                                                                                                                                               c('YES', 'NO'), selected = 'NO')
                                                                                                                         )
                                                                                                           ),
                                                                                                           shiny::column(width = 12,
                                                                                                                         shiny::uiOutput('drug_jump_neigh_opt')
                                                                                                           )
                                                                                                         )
                                                                                         ),
                                                                                         shiny::tabPanel(title = shiny::h3(shiny::span('Transition', style = "font-weight: bold")),
                                                                                                         shiny::sliderInput('lambda', shiny::HTML("&lambda;"), min = 0, max = 1, value = 0.5, step = 0.01)
                                                                                         )
                                                                  ),
                                                                  shiny::fluidRow(
                                                                    shiny::column(width = 12,
                                                                                  shinydashboard::box(width = 12,
                                                                                                      shiny::numericInput('cores_rwr',
                                                                                                                          shiny::h4(shiny::span('Select the number of cores to run the RWR', style = "font-weight: bold")),
                                                                                                                          min = 1, max = 190, value = 1)
                                                                                  )
                                                                    )
                                                                  ),
                                                                  shiny::fluidRow(
                                                                    shiny::tags$hr(),
                                                                    shiny::column(width = 6, shiny::actionButton('rwr_button', 'Submit')),
                                                                    #shiny::column(width = 6, shiny::uiOutput('download_rwr_mat')),
                                                                  )
                                              )
                                ),
                                shiny::column(width = 6,
                                              shiny::tags$hr(),
                                              shinydashboard::infoBox(
                                                title = shiny::h3(shiny::span('Information', style = "font-weight: bold")),
                                                value = shiny::uiOutput('info_box_7'),
                                                subtitle = NULL,
                                                icon = shiny::icon("info"), color = "aqua", width = NULL,
                                                fill = FALSE)
                                )
                              )
      ),
      
      shinydashboard::tabItem(tabName = "dr_tab_1",
                              
                              shiny::fluidRow(
                                shiny::column(width = 6,
                                              shiny::uiOutput('back2P1_from1.2')
                                ),
                                shiny::column(width = 6,
                                              shiny::uiOutput('jump2P7_from1.2'),  align = 'right'
                                )
                              ),
                              
                              shiny::fluidRow(
                                shiny::column(width = 6,
                                              shiny::fluidRow(
                                                shiny::column(width = 12,
                                                              shiny::tags$hr(),
                                                              shinydashboard::box(
                                                                title = shiny::h2(shiny::span("Upload Similarity matrix", style = "font-weight: bold")),
                                                                width = 12, status = "primary", solidHeader = FALSE,
                                                                shiny::radioButtons(inputId = 'example_simMat', 
                                                                                    label = shiny::h4(shiny::span( 'Do you want load to select an example similarity matrix?', style = "font-weight: bold")),
                                                                                    choices = c('YES', 'NO'), selected = 'NO'),
                                                                
                                                                shiny::conditionalPanel(
                                                                  condition = "input.example_simMat == 'NO'",
                                                                  shiny::radioButtons(inputId = 'simMat_type', 
                                                                                      label = shiny::h4(shiny::span( 'Which type of similarity matrix do you want to load?', style = "font-weight: bold")),
                                                                                      choices = c('Similarity matrix', 'Embedded Similarity matrix', 'Other'), selected = 'Similarity matrix'),
                                                                  shiny::fileInput("simMat_file", label = shiny::h4(shiny::span( 'Select a file', style = "font-weight: bold")),
                                                                                   multiple = FALSE,
                                                                                   accept = c("text/csv", '.RDS', "text/comma-separated-values,text/plain", ".csv")
                                                                  ),
                                                                  shiny::uiOutput('simMat_name')
                                                                ),
                                                                shiny::conditionalPanel(
                                                                  condition = "input.example_simMat == 'YES'",
                                                                  shiny::selectInput(inputId = 'select_simMat', 
                                                                                     label = shiny::h4(shiny::span( 'Which type of similarity matrix do you want to load?', style = "font-weight: bold")),
                                                                                     choices = c('Embedded BRCA-5-omics FDA-drugs similarity matrix',
                                                                                                 'UMAP of Embedded BRCA-5-omics FDA-drugs similarity matrix')
                                                                  )
                                                                ),
                                                                
                                                                shiny::actionButton('load_simMat_button', 'Load')
                                                              )
                                                ),
                                                shiny::column(width = 12,
                                                              shiny::htmlOutput(outputId = "simMat_sum_info")
                                                ),
                                                shiny::column(width = 12,
                                                              shinydashboard::box(
                                                                title = shiny::h2(shiny::span("Upload annotation file", style = "font-weight: bold")),
                                                                width = 12, status = "primary", solidHeader = FALSE,
                                                                
                                                                shiny::radioButtons(inputId = 'example_simMat_anno', 
                                                                                    label = shiny::h4(shiny::span( 'Do you want to load an example annotation file?', style = "font-weight: bold")),
                                                                                    choices = c('YES', 'NO'), selected = 'NO'),
                                                                shiny::conditionalPanel(
                                                                  condition = "input.example_simMat_anno == 'NO'",
                                                                  shiny::fileInput(inputId = "anno_file2",
                                                                                   label = shiny::h4(shiny::span( 'Select a file', style = "font-weight: bold")),
                                                                                   multiple = FALSE)
                                                                ),
                                                                shiny::actionButton('load_anno_button2', 'Load')
                                                              ),
                                                              shiny::htmlOutput(outputId = "anno_info2")
                                                )
                                                
                                              )
                                ),
                                shiny::column(width = 6,
                                              shiny::tags$hr(),
                                              shinydashboard::infoBox(
                                                title = shiny::h3(shiny::span('Information', style = "font-weight: bold")),
                                                value = shiny::uiOutput('info_box_8'),
                                                subtitle = NULL,
                                                icon = shiny::icon("info"), color = "aqua", width = NULL,
                                                fill = FALSE),
                                              #shiny::htmlOutput(outputId = "anno_info_extra1"),
                                              #shiny::htmlOutput(outputId = "anno_info_extra2")
                                )
                                
                              )
      ),
      
      shinydashboard::tabItem(tabName = "dr_tab_2",
                              shiny::fluidRow(
                                shiny::column(width = 6,
                                              shiny::uiOutput('back2P6'),
                                              shiny::uiOutput('back2P1.2')),
                                shiny::column(width = 6,
                                              shiny::uiOutput('jump2P8'),  align = 'right'
                                )
                              ),
                              
                              shiny::fluidRow(
                                shiny::column(width = 5,
                                              shiny::fluidRow(
                                                shiny::column(width=12,
                                                              shiny::hr(),
                                                              shiny::uiOutput("dr_opt")
                                                ),
                                                shiny::column(width=12,
                                                              shiny::uiOutput('download_dr_box')
                                                )
                                              )
                                ),
                                
                                shiny::column(width = 7,
                                              shiny::fluidRow(
                                                shiny::column(width=12,
                                                              shiny::hr(),
                                                              shinydashboard::infoBox(
                                                                title = shiny::h3(shiny::span('Information', style = "font-weight: bold")),
                                                                value = shiny::uiOutput('info_box_9'),
                                                                subtitle = NULL,
                                                                icon = shiny::icon("info"), color = "aqua", width = NULL,
                                                                fill = FALSE)
                                                ),
                                              )
                                )
                              ),
                              
                              shiny::fluidRow(
                                shiny::column(width = 12,
                                              shinydashboard::tabBox(width = 12,
                                                                     shiny::tabPanel(title = shiny::h3(shiny::span('Plot', style = "font-weight: bold")), 
                                                                                     #shiny::uiOutput('dr_plot_box')
                                                                                     shiny::fluidRow(
                                                                                       shiny::column(width = 12, 
                                                                                                     shiny::uiOutput('dr_plot_box')
                                                                                       ),
                                                                                       shiny::column(width = 12,
                                                                                                     shiny::uiOutput(outputId = 'update_dr_plot_opt')
                                                                                       ),
                                                                                       shiny::column(width = 12,
                                                                                                     shiny::hr()
                                                                                       ),
                                                                                       shiny::column(width = 6,
                                                                                                     shiny::downloadButton(outputId = 'download_dr_plot', label = 'Download plot'),
                                                                                                     align = "left"
                                                                                       ),
                                                                                       shiny::column(width = 6,
                                                                                                     shiny::uiOutput(outputId = 'update_dr_plot_btn'),
                                                                                                     align = "right"
                                                                                       )
                                                                                       
                                                                                     )
                                                                      ), 
                                                                      shiny::tabPanel(title = shiny::h3(shiny::span('Density plot', style = "font-weight: bold")),
                                                                                      shiny::uiOutput('density_opt'),
                                                                                      shinycssloaders::withSpinner(plotly::plotlyOutput('density_plot', height = '600px'))
                                                                      )
                                                                     )
                                              )
                                ),
                                shiny::fluidRow(
                                  shiny::column(width = 6,
                                                shinydashboard::box(width = 12, collapsible = TRUE, 
                                                                    title = shiny::h3(shiny::span('Lower-dimensional table', style = "font-weight: bold")),
                                                                    DT::DTOutput('data_orig_table')),
                                                shiny::uiOutput('other_table_box')
                                  ),
                                  shiny::column(width = 6,
                                                shinydashboard::box(width = 12,
                                                                    title = shiny::h3(shiny::span('Selected table', style = "font-weight: bold")),
                                                                    DT::DTOutput('selected_data'))
                                  )
                                )
                              ),
                              
                              shinydashboard::tabItem(tabName = "cl_tab",
                                                      shiny::fluidRow(
                                                        shiny::column(width = 6,
                                                                      shiny::uiOutput('back2P7')),
                                                        shiny::column(width = 6,
                                                                      shiny::uiOutput('jump2P9'), align = 'right')
                                                      ),
                                                      
                                                      shiny::fluidRow(
                                                        shiny::column(width = 5,
                                                                      shiny::fluidRow(
                                                                        shiny::column(width=12,
                                                                                      shiny::hr(),
                                                                                      shiny::uiOutput("cl_opt"),
                                                                        ),
                                                                        shiny::column(width=12,
                                                                                      shiny::uiOutput('info_created_cl'),
                                                                        )
                                                                      )
                                                        ),
                                                        
                                                        shiny::column(width = 7,
                                                                      shiny::fluidRow(
                                                                        shiny::column(width=12,
                                                                                      shiny::hr(),
                                                                                      shinydashboard::infoBox(
                                                                                        title = shiny::h3(shiny::span('Information', style = "font-weight: bold")),
                                                                                        value = shiny::uiOutput('info_box_10'),
                                                                                        subtitle = NULL,
                                                                                        icon = shiny::icon("info"), color = "aqua", width = NULL,
                                                                                        fill = FALSE)
                                                                                      
                                                                        ),
                                                                        shiny::column(width=12,
                                                                                      shiny::conditionalPanel(
                                                                                        condition = "input.cluster_method == 'hclust'",
                                                                                        shiny::uiOutput(outputId = "plotHclust_box")                                                                        )
                                                                                      
                                                                        ),
                                                                        
                                                                      )
                                                                      
                                                        )
                                                      ),
                                                      shiny::fluidRow(
                                                        shiny::column(width=12,
                                                          shinydashboard::box(width = 12, 
                                                                              shiny::fluidRow(
                                                                                shiny::column(width = 12,
                                                                                              shiny::uiOutput(outputId = 'cl_plot_box'),
                                                                                              #shinycssloaders::withSpinner(plotly::plotlyOutput('cl_plot', height = '600px')),
                                                                                              shiny::uiOutput(outputId = 'manual_cl_plot_box')
                                                                                ),
                                                                                # shiny::column(width = 12,
                                                                                #               shiny::uiOutput(outputId = 'shape_opt'),
                                                                                #               shiny::uiOutput(outputId = 'manual_shape_opt')
                                                                                # ),
                                                                                shiny::column(width = 12,
                                                                                              shiny::hr()
                                                                                ),
                                                                                shiny::column(width = 6,
                                                                                              shiny::downloadButton(outputId = 'download_cl_plot', label = 'Download plot'),
                                                                                              align = "left"
                                                                                ),
                                                                                # shiny::column(width = 6,
                                                                                #               shiny::uiOutput(outputId = 'update_cl_plot_btn'),
                                                                                #               shiny::uiOutput(outputId = 'update_manual_cl_plot_btn'),
                                                                                #               align = "right"
                                                                                # )
                                                                              )
                                                          )
                                                        )
                                                      ),
                                                      shiny::fluidRow(
                                                        shiny::column(width=6,
                                                                      shinycssloaders::withSpinner(shiny::uiOutput('cluster_table_box')),
                                                                      shinycssloaders::withSpinner(shiny::uiOutput('manual_cluster_table_box'))
                                                        ),
                                                        shiny::column(width = 6,
                                                                      shinydashboard::box(width = 12,
                                                                                          title = shiny::h3(shiny::span('Selected cluster', style = "font-weight: bold")),
                                                                                          DT::DTOutput('cl_selected_data'),
                                                                                          shiny::hr(),
                                                                                          shiny::actionButton(inputId = 'add_cluster_btn', label = 'Submit')
                                                                                          )
                                                        )
                                                      )
                              ), 
                              shinydashboard::tabItem(tabName = "clupath_tab",
                                                      shiny::fluidRow(
                                                        shiny::column(width = 6,
                                                                      shiny::uiOutput('back2P8')),
                                                        shiny::column(width = 6,
                                                                      shiny::uiOutput('jump2P10'), align = 'right')
                                                      ),
                                                      
                                                      shiny::fluidRow(
                                                        shiny::column(width = 6,
                                                                      shiny::fluidRow(
                                                                        shiny::column(width=12,
                                                                                      shiny::hr(),
                                                                                      shinydashboard::box(title = shiny::h2(shiny::span("Pathway Enrichment Analysis", style = "font-weight: bold")), 
                                                                                                          width = 12, collapsible = TRUE, 
                                                                                                          solidHeader = FALSE, status = 'primary',
                                                                                                          shiny::uiOutput('cl_table_selector'),
                                                                                                          shiny::tags$hr(),
                                                                                                          shiny::actionButton(inputId = 'pea_btn', label = 'Submit')
                                                                                      )
                                                                        ),
                                                                        shiny::column(width=12,
                                                                                      shiny::hr(),
                                                                                      shiny::conditionalPanel(
                                                                                        condition = 'input.pea_btn >= 1',
                                                                                        shiny::uiOutput('clu_path_box')
                                                                                        
                                                                                      )
                                                                        )
                                                                      )
                                                        ),
                                                        shiny::column(width=6,
                                                                      shiny::fluidRow(
                                                                        shiny::column(width=12,
                                                                                      shiny::hr(),
                                                                                      shinydashboard::infoBox(
                                                                                        title = shiny::h3(shiny::span('Information', style = "font-weight: bold")),
                                                                                        value = shiny::uiOutput('info_box_11'),
                                                                                        subtitle = NULL,
                                                                                        icon = shiny::icon("info"), color = "aqua", width = NULL,
                                                                                        fill = FALSE)
                                                                        ),
                                                                        shiny::column(width=12,
                                                                                      shiny::uiOutput('info_created_pea'),
                                                                        )
                                                                        
                                                                      )
                                                        )
                                                      ),
                                                      
                                                      shiny::fluidRow(
                                                        shiny::column(width=12,
                                                                      shiny::conditionalPanel(
                                                                        condition = "input.showOne >= 1",
                                                                        shinydashboard::box(
                                                                          width = 12, 
                                                                          solidHeader = FALSE, status = 'primary',
                                                                          shinycssloaders::withSpinner(plotly::plotlyOutput(outputId = "TOPpathwayAnnotation", height = '600px' )), 
                                                                          
                                                                        )
                                                                      )            
                                                        )
                                                      ),
                                                      shiny::fluidRow(
                                                        shiny::column(width=12,
                                                                      shiny::conditionalPanel(
                                                                        condition = "input.showOne >= 1",
                                                                        shinydashboard::box(
                                                                          width = 12, 
                                                                          solidHeader = FALSE, status = 'primary',
                                                                          DT::DTOutput('TOPtable')
                                                                        )
                                                                        
                                                                      )            
                                                        )
                                                      )
                              ),
                              
                              shinydashboard::tabItem(tabName = "pathclu_tab",
                                                      shiny::fluidRow(
                                                        shiny::column(width = 6,
                                                                      shiny::uiOutput('back2P9')),
                                                        shiny::column(width = 6,
                                                                      shiny::uiOutput('jump2P11'), align ='right')
                                                      ),
                                                      
                                                      shiny::fluidRow(
                                                        shiny::column(width = 6,
                                                                      shiny::fluidRow(
                                                                        shiny::column(width=12,
                                                                                      shiny::hr(),
                                                                                      shinydashboard::box(title = shiny::h2(shiny::span("Pathway to Clusters Analysis", style = "font-weight: bold")), 
                                                                                                          width = 12, 
                                                                                                          solidHeader = FALSE, status = 'primary',
                                                                                                          shiny::uiOutput('pea_table_selector_2'),
                                                                                                          shiny::radioButtons('source',
                                                                                                                              'Select a database',
                                                                                                                              choices = c("CORUM", "GO:CC", "GO:MF", "HPA", "WP",
                                                                                                                                          "GO:BP", "KEGG", "REAC", "TF", "HP", "MIRNA"),
                                                                                                                              selected = character(0)
                                                                                                                              
                                                                                                          ),
                                                                                                          shiny::tags$hr(),
                                                                                                          shiny::actionButton(inputId = 'path2clust', label = 'Submit')
                                                                                      )
                                                                        )
                                                                      )
                                                        ),
                                                        shiny::column(width=6,
                                                                      shiny::fluidRow(
                                                                        shiny::column(width=12,
                                                                                      shiny::hr(),
                                                                                      shinydashboard::infoBox(
                                                                                        title = shiny::h3(shiny::span('Information', style = "font-weight: bold")),
                                                                                        value = shiny::uiOutput('info_box_12'),
                                                                                        subtitle = NULL,
                                                                                        icon = shiny::icon("info"), color = "aqua", width = NULL,
                                                                                        fill = FALSE)
                                                                        ),
                                                                        shiny::column(width=12,
                                                                                      shiny::conditionalPanel(
                                                                                        condition = 'input.path2clust >= 1',
                                                                                        shinydashboard::box(title = shiny::h2(shiny::span("Pathway to Clusters Plot", style = "font-weight: bold")), 
                                                                                                            width = 12, 
                                                                                                            solidHeader = FALSE, status = 'primary',
                                                                                                            shiny::uiOutput('pathway_selector'),
                                                                                                            shiny::tags$hr(),
                                                                                                            shiny::actionButton(inputId = "show" ,label = "Submit")
                                                                                        )
                                                                                      )
                                                                        )
                                                                      )
                                                        )
                                                      ),
                                                      
                                                      shiny::fluidRow(
                                                        shiny::column(width=12,
                                                                      shinydashboard::box(#title = shiny::h2(shiny::span("", style = "font-weight: bold")), 
                                                                        width = 12, 
                                                                        solidHeader = FALSE, status = 'primary',
                                                                        shinycssloaders::withSpinner(plotly::plotlyOutput(outputId = "pathwayPlot", height = '600px'))
                                                                      )
                                                        )
                                                      ),
                                                      
                                                      shiny::fluidRow(
                                                        shiny::column(width=12,
                                                                      shinydashboard::box(#title = shiny::h2(shiny::span("", style = "font-weight: bold")), 
                                                                        width = 12, 
                                                                        solidHeader = FALSE, status = 'primary',
                                                                        DT::DTOutput('Alltable') 
                                                                      )
                                                        )
                                                      )
                                                      
                              ),
                              
                              shinydashboard::tabItem(tabName = "dd_tab",
                                                      shiny::fluidRow(
                                                        shiny::column(width = 6,
                                                                      shiny::uiOutput('back2P10'))
                                                      ),
                                                      
                                                      shiny::fluidRow(
                                                        shiny::column(width = 6,
                                                                      shiny::fluidRow(
                                                                        shiny::column(width=12,
                                                                                      shiny::hr(),
                                                                                      shinydashboard::box(title = shiny::h2(shiny::span("Drug Discovery", style = "font-weight: bold")), 
                                                                                                          width = 12, 
                                                                                                          solidHeader = FALSE, status = 'primary',
                                                                                                          shiny::uiOutput('cl_table_selector2'),
                                                                                                          shiny::radioButtons(inputId = 'approach', label = 'Select an approach',
                                                                                                                              choices = c("drug",  "cluster"),
                                                                                                                              selected = "cluster"
                                                                                                          ),
                                                                                                          shiny::tags$hr(),
                                                                                                          shiny::actionButton('drug_button', 'Submit'),
                                                                                                          
                                                                                                          shiny::tags$hr(),
                                                                                                          shiny::conditionalPanel(
                                                                                                            condition = 'input.drug_button >= 1',
                                                                                                            shiny::uiOutput('drug_cluster_selector'),
                                                                                                            shiny::tags$hr(),
                                                                                                            shiny::actionButton("show_drug",
                                                                                                                                label = "Submit"
                                                                                                            )
                                                                                                          )
                                                                                      )
                                                                        )
                                                                      )
                                                        ),
                                                        shiny::column(width=6,
                                                                      shiny::fluidRow(
                                                                        shiny::column(width=12,
                                                                                      shiny::hr(),
                                                                                      shinydashboard::infoBox(
                                                                                        title = shiny::h3(shiny::span('Information', style = "font-weight: bold")),
                                                                                        value = shiny::uiOutput('info_box_13'),
                                                                                        subtitle = NULL,
                                                                                        icon = shiny::icon("info"), color = "aqua", width = NULL,
                                                                                        fill = FALSE)
                                                                        )
                                                                      )
                                                        )
                                                      ),
                                                      
                                                      shiny::fluidRow(
                                                        shiny::column(width=12,
                                                                      shiny::uiOutput(outputId = "drug_plot_box"),
                                                                      
                                                        )
                                                      ),
                                                      
                                                      
                                                      shiny::fluidRow(
                                                        shiny::column(width = 6, 
                                                                      shinydashboard::box(title = shiny::h2(shiny::span("DRUGS", style = "font-weight: bold")), 
                                                                                          width = 12, 
                                                                                          solidHeader = FALSE, status = 'primary',
                                                                                          DT::DTOutput('drug_table')
                                                                      )
                                                        ),
                                                        shiny::column(width = 6,
                                                                      shinydashboard::box(title = shiny::h2(shiny::span("GENES", style = "font-weight: bold")), 
                                                                                          width = 12, 
                                                                                          solidHeader = FALSE, status = 'primary',
                                                                                          DT::DTOutput('gene_table')
                                                                      )
                                                        )
                                                      )
                                                      
                              )
                              
      )
    ),
    
    shiny::tags$head(shiny::tags$style(shiny::HTML('* {font-family: "LM Roman 10"}; font-size: 1em;')))
  
)  
  
  server <- function(input, output, session) {
    
    ############################################################################################
    #                            DATA: LOADING DATA (INPUT)                       #
    ############################################################################################
    
    output$info_box_1 <- renderText({
      shiny::HTML("<br/> <span style='font-weight:normal;'>
            <p align='justify'>
            Welcome to the MiDNE Shiny app, an R package for network-based multi-omics and drug data integration.<br/>
            You can start the pipeline by uploading omics matrices and drug files on the following page,
            or click one of the two <b>buttons</b> in the top left-hand corner to proceed. </span> <br/>
            </p>
            <hr style='border-top: 1px solid white;'>

            <b>Upload Omic Matrices</b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            These matrices should have features on rows and the same or non-overlapping sets of samples on columns. </span> <br/>
            </p>
            <hr style='border-top: 1px solid white;'>

            <b>Use Example Dataset</b>: <br/>  <span style='font-weight:normal;'>
             <p align='justify'>
            Alternatively, you can upload the example dataset consisting of five omic matrices of BRCA patients,
            including transcriptomics, proteomics, epigenomics, and CNVs data. </span> </span>
            </p>
            <hr style='border-top: 1px solid white;'>

            <b>Upload Drugs File</b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            This step is not mandatory but essential for creating a drug network.
            The file should be a dataframe with one column containing the drug IDs, 
            and as many rows as the number of drugs. </span> <br/>
            </p>
            <hr style='border-top: 1px solid white;'>
            
            <b>Upload Annotation File</b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            This step is not mandatory but useful for enriching subsequent plots.
            The file should be a dataframe with a variable number of columns,
            where one column must contain the nodes IDs, the same ones present as column names in the uploaded files. </span> <br/>
            </p>
            <hr style='border-top: 1px solid white;'>

            <b>Use Example Annotation File</b>: <br/>  <span style='font-weight:normal;'>
            <p align='justify'>
            If you have selected the Example Dataset, you can use the related annotation file.
            Genes are classified based on different criteria retrieved from HGNC and COSMIC databases, while drugs are 
            annotated according to the related pathways and the FDA status (DrugBank). </span>
            </p>
            <hr style='border-top: 1px solid white;'>
            
            <span style='font-weight:normal;'>
            <p align='justify'>
            Accepted formats include .csv, .txt, and .RDS.  </span>
            </p> <hr style='border-top: 1px solid white;'>
            ")
      
    })
    
    
    ############## inputs loading #############
    
    ### Load omics matrices
    omics_files <- shiny::eventReactive(input$load_mat_button, {
      omicsFiles <- list()
      
      if (input$omics_example_opt == 'No'){
        inFiles <- input$omics_files
        if (is.null(inFiles))
          return(NULL)
        for (i in 1:nrow(inFiles)) {
          name <- tools::file_path_sans_ext(inFiles$name[i])
          input_name <- input[[paste0('name', i)]]
          new_name <- ifelse(nchar(input_name) == 0, name, input_name)
          
          ext <- tools::file_ext(inFiles$name[i])
          file <- switch(ext,
                         csv = vroom::vroom(inFiles$datapath[i], delim = ","),
                         tsv = vroom::vroom(inFiles$datapath[i], delim = "\t"),
                         RDS = readRDS(inFiles$datapath[i]),
                         validate("Invalid file; Please upload a .csv, .tsv or .RDS file")
          )
          file <- as.matrix(file)
          if (is.numeric(file)){
            sorted_file <- file[, sort(colnames(file))]
            omicsFiles[[new_name]] <- sorted_file
          }else{
            shinyalert::shinyalert("Type Error", "Uploaded Data is not a numeric matrix", closeOnClickOutside = TRUE, type = "error")
            returnValue()
          }
        }
        
      } else {
        input_list <- input$omics_example_files
        example_data <- list()
        for (file in input_list) {
          example_data[[ file ]] <- readRDS(paste0('MiDNE/DATA/data/biological/', file, '.RDS' ))
        }
        omicsFiles <- example_data
      }
      
      omicsFiles
    })
    
    
    observeEvent(input$load_mat_button, {
      shinyalert::shinyalert(
        title = "Wait",
        text = "Waiting for data loading",
        size = "xs",
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        type = "info",
        showConfirmButton = TRUE,
        confirmButtonText = "OK",
        confirmButtonCol = "#004192",
        showCancelButton = FALSE,
        imageUrl = "",
        animation = TRUE
      )
      
    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE
    )
    
    
    ## Load drug matrices 
    drug_files <- eventReactive(input$load_drug_mat_button, {
      drugFiles <- list()
      if (input$drug_example_opt == 'No'){
        
        inFiles <- input$drug_files
        if (is.null(inFiles))
          return(NULL)
        for (i in 1:nrow(inFiles)) {
          name <- tools::file_path_sans_ext(inFiles$name[i])
          input_name <- input[[paste0('drug_name', i)]]
          new_name <- ifelse(nchar(input_name) == 0, name, input_name)
          
          ext <- tools::file_ext(inFiles$name[i])
          file <- switch(ext,
                         csv = vroom::vroom(inFiles$datapath[i], delim = ","),
                         tsv = vroom::vroom(inFiles$datapath[i], delim = "\t"),
                         RDS = readRDS(inFiles$datapath[i]),
                         validate("Invalid file; Please upload a .csv, .tsv or .RDS file")
          )
          #file <- as.matrix(file)
          drugFiles[[new_name]] <- file
          
          # if (is.numeric(file)){
          #   sorted_file <- file[, sort(colnames(file))]
          #   drugFiles[[new_name]] <- sorted_file
          # }else{
          #   shinyalert::shinyalert("Type Error", "Uploaded Data is not a numeric matrix", closeOnClickOutside = TRUE, type = "error")
          #   returnValue()
          # }
        }
      } else {
        drugFiles$FDAdrugs <- readRDS('MiDNE/DATA/data/pharmacological/FDAdrugs.RDS' )
      }
      
      drugFiles
    })
    
    observeEvent(input$load_drug_mat_button, {
      shinyalert::shinyalert(
        title = "Wait",
        text = "Waiting for data loading",
        size = "xs",
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        type = "info",
        showConfirmButton = TRUE,
        confirmButtonText = "OK",
        confirmButtonCol = "#004192",
        showCancelButton = FALSE,
        imageUrl = "",
        animation = TRUE
      )
      
    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE
    )
    
    ### Load annotation file
    annotation <- reactiveVal(NULL)
    shiny::observeEvent(input$load_anno_button, {
      listFiles <- list()
      
      if (input$anno_example_opt == 'No'){
        inFiles <- input$anno_file
        if (is.null(inFiles)){
          return(NULL)
        }else{
          ext <- tools::file_ext(inFiles$name)
          file <- switch(ext,
                         csv = vroom::vroom(inFiles$datapath, delim = ","),
                         RDS = readRDS(inFiles$datapath),
                         validate("Invalid file; Please upload a .csv or .RDS file")
          )
          if (is.data.frame(file)){
            listFiles[['annotation']] <- file
            annotation(listFiles)
          }else {
            shinyalert::shinyalert("Type Error", "Uploaded Data is not a dataframe",closeOnClickOutside = TRUE, type = "error")
            returnValue()
          }
        }
      }else if (input$anno_example_opt == 'Yes'){
        listFiles[['annotation']] <- readRDS('MiDNE/DATA/annotation/all_genes_drugs_annotation.RDS' )
        annotation(listFiles)
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$load_anno_button, {
      shinyalert::shinyalert(
        title = "Wait",
        text = "Waiting for data loading",
        size = "xs",
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        type = "info",
        showConfirmButton = TRUE,
        confirmButtonText = "OK",
        confirmButtonCol = "#004192",
        showCancelButton = FALSE,
        imageUrl = "",
        animation = TRUE
      )
      
    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE
    )
    
    ############## set files names #############
    
    output$omics_names <- shiny::renderUI({
      omicsNames <- list()
      inFiles <- input$omics_files
      if (is.null(inFiles))
        return(NULL)
      for (i in 1:nrow(inFiles)) {
        name <- tools::file_path_sans_ext(inFiles$name[i])
        omicsNames[[i]] <- shiny::textInput(inputId = paste0('name', i),
                                            label = paste('Rename the ', name, ' file'),
                                            value = NULL)
      }
      omicsNames
    })
    
    output$drug_names <- renderUI({
      drugNames <- list()
      inFiles <- input$drug_files
      if (is.null(inFiles))
        return(NULL)
      for (i in 1:nrow(inFiles)) {
        name <- tools::file_path_sans_ext(inFiles$name[i])
        drugNames[[i]] <- shiny::textInput(inputId =paste0('drug_name', i), 
                                           label = paste('Enter new name for ', name, ' file'),
                                           value = NULL)
      }
      drugNames
    })
    
    
    
    ###
    
    output$jump2P2 <- renderUI({
      if (length(omics_files()) != 0){
        actionButton('jump2P2', label = 'Next', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }
    })
    
    observeEvent(input$jump2P2, {
      shinydashboard::updateTabItems(session, inputId = "tabs", selected = "mat_sub_2")
      shinyjs::runjs('$(".sidebar-menu .treeview").removeClass("active"); $("#mat_tab").closest(".treeview").addClass("active");')
    })
    
    ###
    output$jump2P1.1 <- renderUI({
      actionButton('jump2P1.1', label = 'Go to RWR', shiny::icon("paper-plane"),
                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
    })
    
    observeEvent(input$jump2P1.1, {
      shinydashboard::updateTabItems(session, inputId = "tabs", selected = "rwr_tab_1")
    })
    
    ###
    output$jump2P1.2 <- renderUI({
      actionButton('jump2P1.2', label = 'Go to DR', shiny::icon("paper-plane"),
                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
    })
    
    observeEvent(input$jump2P1.2, {
      shinydashboard::updateTabItems(session, inputId = "tabs", selected = "dr_tab_1")
    })
    
    
    ############################################################################################
    #                        DATA: PRE-PROCESSING (INPUT)                      #
    ############################################################################################
    
    shinyjs::addCssClass(selector = "a[data-value='mat_sub_2']", class = "inactiveLink")
    observe({
      if (length(omics_files()) != 0){
        shinyjs::removeCssClass(selector = "a[data-value='mat_sub_2']", class = "inactiveLink")
      }
    })
    
    ############
    output$info_box_2 <- renderText({
      shiny::HTML("<br/> <span style='font-weight:normal;'>
        <p align='justify'>
            In this section you can manage the omics matrices before inferring the networks. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <b> Omics Matrices Pre-processing </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            Here, you can decide for each omics matrix whether to remove samples with no detected features  and/or normalize them by column.
            When the <b> Submit </b> button is pressed, the omics matrices will be processed <u>all at once</u> </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <b> Omics Matrices Intersection </b>: <br/>  <span style='font-weight:normal;'>
            <p align='justify'>
            This box will be shown only if you have loaded more than one omics matrix and the <b> Submit </b> button have already been pressed.
            If you select <b> Yes</b>, only the overlapping sets of features will be kept.</span>
            </p> <hr style='border-top: 1px solid white;'>

            <span style='font-weight:normal;'>
            <p align='justify'>
            The processed matrices can be downloaded from the <b>download box</b> displayed below. </span>
            </p> <hr style='border-top: 1px solid white;'>
            ")
    })
    
    ############
    
    intersected_omics_mat <- shiny::reactiveValues()
    shiny::observeEvent(input$intersect_btn, {
      processed_mats <- shiny::reactiveValuesToList(proc_matrices)
      #shiny::isolate({
      if (input$intersect_opt == 'YES'){
        t_processed_mats <- lapply(processed_mats, t)
        intersected_by_features <- MoNETA::get_intersection_matrices(t_processed_mats)
        t_intersected_by_features <- lapply(intersected_by_features, t)
        
        print(dim(t_processed_mats[[1]]))
        print(dim(intersected_by_features[[1]]))
        print(dim(t_intersected_by_features[[1]]))
        
        intersected_omics_mat$matrices <- t_intersected_by_features
      }else {
        intersected_omics_mat$matrices <- processed_mats
      }
      
      shinyalert::shinyalert(
        title = "Success",
        text = 'Intersection step done! <br/>  Now press <b> "Next" </b> in the top right-hand corner to continue.',
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        type = "success",
        showConfirmButton = TRUE,
        confirmButtonText = "OK",
        confirmButtonCol = "#004192",
        showCancelButton = FALSE,
        imageUrl = "",
        animation = TRUE
      )
      
      #})
    })
    
    
    
    ############################################################################################
    #                        NETWORK : NETWORK INFERENCE (INPUT)                      #
    ############################################################################################
    
    shinyjs::addCssClass(selector = "a[data-value='net_sub_1']", class = "inactiveLink")
    observe({
      if (length(omics_files()) == 1 && length(reactiveValuesToList(proc_matrices)) != 0){
        shinyjs::removeCssClass(selector = "a[data-value='net_sub_1']", class = "inactiveLink")
      } else if (length(omics_files()) >= 1 && length(reactiveValuesToList(proc_matrices)) != 0 && shiny::isTruthy(input$intersect_btn)){
        shinyjs::removeCssClass(selector = "a[data-value='net_sub_1']", class = "inactiveLink")
      }else{
        return(NULL)
      }
      
    })
    
    ############
    output$info_box_3 <- renderText({
      shiny::HTML("<br/> <span style='font-weight:normal;'>
        <p align='justify'>
            In this section you can create the omics and drug networks. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>
            
            <span style='color: red;'>
            <p align='justify'>
            After clicking the button to infer a network, please wait for the process to finish before clicking the next one.
            If you click them one after the other, the second process will be queued, and the <b>Wait</b> message will be displayed accordingly. </span>
            </p> <hr style='border-top: 1px solid white;'> </span>

            <b> Omics Network Inference </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            You can choose between two methods for creating biological networks: correlation analysis (Spearman and Pearson) 
            from continuous data matrices, or co-occurrence networks using Fisher's exact test and post-hoc analysis. 
            For methylation data expressed as BETA values, you can set a threshold for data discretization. 
            In both cases, you can select a p-value correction method (FDR or Bonferroni).<br/>
            </p> <hr style='border-top: 1px solid white;'>
            
            <b> Drug Network Inference </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            The inferred drug network is an isolated graph.<br/>
            </p> <hr style='border-top: 1px solid white;'>
            
            <span style='font-weight:normal;'>
            <p align='justify'>
            The networks can be downloaded from the <b>download box</b> located on the left-hand side. </span>
            </p> <hr style='border-top: 1px solid white;'>
            ")
    })
    
    
    ############################################################################################
    #                        NETWORK : NETWORK FILTERING (INPUT)                      #
    ############################################################################################
    
    shinyjs::addCssClass(selector = "a[data-value='net_sub_2']", class = "inactiveLink")
    observe({
      if (length(shiny::reactiveValuesToList(gene_networks)) == length(omics_files())){
        shinyjs::removeCssClass(selector = "a[data-value='net_sub_2']", class = "inactiveLink")
      }else{
        return(NULL)
      }
    })
    
    output$info_box_4 <- renderText({
      shiny::HTML("<br/> <span style='font-weight:normal;'>
            <p align='justify'>
            Here, you can prune the edges of inferred networks by filtering the <b>weight</b> column. 
            If the network was generated using <u>correlation</u>, it is suggested to filter by selecting only the extremes of the range, 
            while if the <u>Fisher's exact test + post-hoc</u>  method was chosen, select only the weight interval above the desired threshold. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <b> Show button </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            When the <b>Show</b> button is pressed, the network and the related
            edge weight and node degree distributions (based on the chosen filtering threshold) will be displayed.
            This action does not trigger the update of the networks. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <b> Submit button </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            When the <b>Submit</b> button is pressed, the networks will be updated
            <u>all at once</u>, allowing you to proceed. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <span style='font-weight:normal;'>
            <p align='justify'>
            The filtered networks can be downloaded from the <b>download box</b> displayed below. </span>
            </p> <hr style='border-top: 1px solid white;'>
            ")
    })
    
    ############## omics network filtering #############
    
    filteredOmicsNetworks <- shiny::reactiveValues()
    shiny::observeEvent(
      eventExpr = {
        buttons <- paste0("update_net_", names(shiny::reactiveValuesToList(gene_networks)))
        list_of_buttons = NULL
        for(var in buttons) {
          list_of_buttons <- append(list_of_buttons, input[[var]])
        }
        req(list_of_buttons)
      },
      handlerExpr = {
        nets <- shiny::reactiveValuesToList(gene_networks)
        for (x in names(nets)){
          net <- nets[[x]]
          if (input[[paste0("range_option", x)]]){
            filt_net <- net[net$weight >= input[[paste0('omics_range', x)]][1] & net$weight <= input[[paste0('omics_range', x)]][2], ]
          }else{
            filt_net <- net[net$weight <= input[[paste0('omics_range', x)]][1] | net$weight >= input[[paste0('omics_range', x)]][2], ]
          }
          filteredOmicsNetworks[[paste0('f_', x)]] <- filt_net
          print('Filering network...done!')
          
          count_unet[[x]] <- 1
          #print(shiny::reactiveValuesToList(count_unet))
          net_update <- sum(unlist(shiny::reactiveValuesToList(count_unet)))
          
        }
        
        if (net_update == length(omics_files())){
          shinyalert::shinyalert(
            title = "Success",
            text = paste0(paste("<b>", x, "</b>"), " network updated. <br/>",
                          paste("<b> Network update status: </b>"), net_update, '/', length(omics_files()), "<br/>",
                          'Now press <b> "Next" </b> in the top right-hand corner to continue.'
            ),
            closeOnEsc = TRUE,
            closeOnClickOutside = TRUE,
            html = TRUE,
            type = "success",
            showConfirmButton = TRUE,
            confirmButtonText = "OK",
            confirmButtonCol = "#004192",
            showCancelButton = FALSE,
            imageUrl = "",
            animation = TRUE
          )
        } else{
          shinyalert::shinyalert(
            title = "Success",
            text = paste0(paste("<b>", x, "</b>"), " network updated. <br/>",
                          paste("<b> Network update status: </b>"), net_update, '/', length(omics_files())
            ),
            closeOnEsc = TRUE,
            closeOnClickOutside = TRUE,
            html = TRUE,
            type = "success",
            showConfirmButton = TRUE,
            confirmButtonText = "OK",
            confirmButtonCol = "#004192",
            showCancelButton = FALSE,
            imageUrl = "",
            animation = TRUE
          )
        }
        return(filt_net)
      },
      ignoreInit = T
    )
    
    
    ############################################################################################
    #                                 RWR: LOADING NETWORKS (INPUT)                            #
    ############################################################################################
    
    output$info_box_5 <- renderText({
      shiny::HTML("<br/> <span style='font-weight:normal;'>
            <p align='justify'>
            You can start the pipeline here by uploading omics networks,
            or click the <b> button </b> in the top left-hand corner to return back to the 'Data Uploading' page. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <b> Upload Omics/Drug Networks </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            The networks are undirected and are loaded as tables with three columns: source, target, and weight.
            The first two columns contain the nodes, while the last one contains the weights of the edges.
            The higher the weight, the stronger the interaction. A 'weight' column containing only '1' represents an unweighted network.</span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <b> Upload Bipartite Network </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            This is an unweighted, undirected and heterogeneous network as the source and target columns contain different nodes types, 
            genes and drugs respectivelly. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>
            
            <b> Upload Annotation File </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            This step is not mandatory but useful to enrich subsequent plots.
            The file should be a dataframe with a variable number of columns,
            where one column must contain the sample IDs representing the nodes in the uploaded files. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <span style='font-weight:normal;'>
            <p align='justify'>
            Accepted formats include .csv, .txt, and .RDS. </span>
             </p> <hr style='border-top: 1px solid white;'>
            ")
    })
    
    
    #################### load OMICS networks #####################
    
    loaded_omics_net_list <- shiny::eventReactive(input$load_omics_net_button, {
      req(input$omics_net_files)
      
      
      if (input$omicsNet_example_opt == 'No'){
        listFiles <- list()
        inFiles <- input$omics_net_files
        if (is.null(inFiles)){
          return(NULL)
        } else {
          for (i in 1:nrow(inFiles)) {
            name <- tools::file_path_sans_ext(inFiles$name[i])
            input_name <- input[[paste0('omics_net_name', i)]]
            new_name <- ifelse(nchar(input_name) == 0, name, input_name)
            
            ext <- tools::file_ext(inFiles$name[i])
            file <- switch(ext,
                           csv = vroom::vroom(inFiles$datapath[i], delim = ","),
                           RDS = readRDS(inFiles$datapath[i]),
                           validate("Invalid file; Please upload a .csv or .RDS file")
            )
            if (is.data.frame(file) & ncol(file) == 3){
              listFiles[[new_name]] <- file
            }else if (!is.data.frame(file) ){
              shinyalert::shinyalert("Type Error", "Uploaded Data is not a dataframe",closeOnClickOutside = TRUE, type = "error")
              returnValue()
            }else if (ncol(file) != 3){
              shinyalert::shinyalert("Column Error", "Uploaded Data has not 3 columns",closeOnClickOutside = TRUE, type = "error")
              returnValue()
            }
          }
        }
        
      } else {
        input_list <- input$omicsNet_example_files
        example_data <- list()
        for (file in input_list) {
          example_data[[ file ]] <- readRDS(paste0('MiDNE/DATA/networks/biological/', file, '.RDS' ))
        }
        listFiles <- example_data
      }
      return(listFiles)
    })
    
    
    observeEvent(input$load_omics_net_button, {
      shinyalert::shinyalert(
        title = "Wait",
        text = "Waiting for data loading",
        size = "xs",
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        type = "info",
        showConfirmButton = TRUE,
        confirmButtonText = "OK",
        confirmButtonCol = "#004192",
        showCancelButton = FALSE,
        imageUrl = "",
        animation = TRUE
      )
      
    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE
    )
    
    
    ############## set files names #############
    
    output$omics_net_names <- shiny::renderUI({
      inFiles <- input$omics_net_files
      if (is.null(inFiles)){
        return(NULL)
      }else{
        omicsNames <- list()
        for (i in 1:nrow(inFiles)) {
          name <- tools::file_path_sans_ext(inFiles$name[i])
          omicsNames[[i]] <- shiny::textInput(paste0('omics_net_name', i), paste('Rename the ', name, ' file'))
        }
        omicsNames
      }
    })
    
    
    #################### load DRUG networks #####################
    
    loaded_drug_net_list <- shiny::eventReactive(input$load_drug_net_button, {
      req(input$drug_net_files)
      
      if (input$drugNet_example_opt == 'No'){
        listFiles <- list()
        inFiles <- input$drug_net_files
        if (is.null(inFiles)){
          return(NULL)
        } else {
          for (i in 1:nrow(inFiles)) {
            name <- tools::file_path_sans_ext(inFiles$name[i])
            input_name <- input[[paste0('drug_net_name', i)]]
            new_name <- ifelse(nchar(input_name) == 0, name, input_name)
            
            ext <- tools::file_ext(inFiles$name[i])
            file <- switch(ext,
                           csv = vroom::vroom(inFiles$datapath[i], delim = ","),
                           RDS = readRDS(inFiles$datapath[i]),
                           validate("Invalid file; Please upload a .csv or .RDS file")
            )
            if (is.data.frame(file) & ncol(file) == 3){
              listFiles[[new_name]] <- file
            }else if (!is.data.frame(file) ){
              shinyalert::shinyalert("Type Error", "Uploaded Data is not a dataframe",closeOnClickOutside = TRUE, type = "error")
              returnValue()
            }else if (ncol(file) != 3){
              shinyalert::shinyalert("Column Error", "Uploaded Data has not 3 columns",closeOnClickOutside = TRUE, type = "error")
              returnValue()
            }
          }
        }
        
      } else {
        listFiles[['FDAdrugs_net']]  <- read_csv('MiDNE/DATA/networks/pharmacological/FDAdrugs_net.csv')
      }
      
      return(listFiles)
    })
    
    
    observeEvent(input$load_drug_net_button, {
      shinyalert::shinyalert(
        title = "Wait",
        text = "Waiting for data loading",
        size = "xs",
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        type = "info",
        showConfirmButton = TRUE,
        confirmButtonText = "OK",
        confirmButtonCol = "#004192",
        showCancelButton = FALSE,
        imageUrl = "",
        animation = TRUE
      )
      
    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE
    )
    
    
    ############## set files names #############
    
    output$drug_net_names <- shiny::renderUI({
      inFiles <- input$drug_net_files
      if (is.null(inFiles)){
        return(NULL)
      }else{
        drugNames <- list()
        for (i in 1:nrow(inFiles)) {
          name <- tools::file_path_sans_ext(inFiles$name[i])
          drugNames[[i]] <- shiny::textInput(paste0('drug_net_name', i), paste('Rename the ', name, ' file'))
        }
        drugNames
      }
    })
    
    
    
    #################### load BIPARTITE networks #####################
    
    bnetwork <- shiny::eventReactive(input$load_bnet_button, {
      req(input$bnet_file)
      
      if (input$bNet_example_opt == 'No'){
        inFiles <- input$bnet_file
        if (is.null(inFiles)){
          return(NULL)
        } else {
          name <- tools::file_path_sans_ext(inFiles$name)
          ext <- tools::file_ext(inFiles$name)
          file <- switch(ext,
                         csv = vroom::vroom(inFiles$datapath, delim = ","),
                         RDS = readRDS(inFiles$datapath),
                         validate("Invalid file; Please upload a .csv or .RDS file")
          )
          if (is.data.frame(file) & ncol(file) == 3){
            bnet <- file
          }else if (!is.data.frame(file) ){
            shinyalert::shinyalert("Type Error", "Uploaded Data is not a dataframe",closeOnClickOutside = TRUE, type = "error")
            returnValue()
          }else if (ncol(file) != 3){
            shinyalert::shinyalert("Column Error", "Uploaded Data has not 3 columns",closeOnClickOutside = TRUE, type = "error")
            returnValue()
          }
        }
        
      } else {
        bnet <- readRDS('MiDNE/DATA/networks/bipartite/FDA_active_DRUGBANK_bnet.RDS' )
      }
      
      return(bnet)
    })
    
    
    observeEvent(input$load_bnet_button, {
      shinyalert::shinyalert(
        title = "Wait",
        text = "Waiting for data loading",
        size = "xs",
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        type = "info",
        showConfirmButton = TRUE,
        confirmButtonText = "OK",
        confirmButtonCol = "#004192",
        showCancelButton = FALSE,
        imageUrl = "",
        animation = TRUE
      )
      
    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE
    )
    
    
    
    ### annotation1
    ### Load annotation file
    
    annotation1 <- reactiveVal(NULL)
    shiny::observeEvent(input$load_anno_button1, {
      req(input$anno_file1)
      
      if (input$rwr_anno_example_opt == 'No'){
        listFiles <- list()
        inFiles <- input$anno_file1
        if (is.null(inFiles)){
          return(NULL)
        } else {
          ext <- tools::file_ext(inFiles$name)
          new_name <- 'annotation1'
          file <- switch(ext,
                         csv = vroom::vroom(inFiles$datapath, delim = ","),
                         RDS = readRDS(inFiles$datapath),
                         validate("Invalid file; Please upload a .csv or .RDS file")
          )
          if (is.data.frame(file)){
            listFiles[[new_name]] <- file
          }else {
            shinyalert::shinyalert("Type Error", "Uploaded Data is not a dataframe",closeOnClickOutside = TRUE, type = "error")
            returnValue()
          }
        }
      } else {
        readRDS('MiDNE/DATA/annotation/all_genes_drugs_annotation.RDS')
      }
      
      annotation1(listFiles)
    })
    
    observeEvent(input$load_anno_button1, {
      shinyalert::shinyalert(
        title = "Wait",
        text = "Waiting for data loading",
        size = "xs",
        closeOnEsc = TRUE,
        closeOnClickOutside = TRUE,
        html = TRUE,
        type = "info",
        showConfirmButton = TRUE,
        confirmButtonText = "OK",
        confirmButtonCol = "#004192",
        showCancelButton = FALSE,
        imageUrl = "",
        animation = TRUE
      )
      
    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE
    )
    
    
    ############################################################################################
    #                                RWR: MULTIPLEX NETWORK (INPUT)                            #
    ############################################################################################
    
    shinyjs::addCssClass(selector = "a[data-value='rwr_tab_2']", class = "inactiveLink")
    observe({
      if (length(shiny::reactiveValuesToList(filteredOmicsNetworks)) == length(omics_files())){
        shinyjs::removeCssClass(selector = "a[data-value='rwr_tab_2']", class = "inactiveLink")
      }else{
        return(NULL)
      }
    })
    observe({
      if (length(loaded_omics_net_list()) != 0 ){
        shinyjs::removeCssClass(selector = "a[data-value='rwr_tab_2']", class = "inactiveLink")
      }else{
        return(NULL)
      }
    })
    
    ##############################################################################
    
    g_net_list <- shiny::reactive({
      if (length(shiny::reactiveValuesToList(filteredOmicsNetworks)) != 0) {
        return(shiny::reactiveValuesToList(filteredOmicsNetworks))
      }else{
        return(NULL)
      }
    })
    
    l_net_list <-  shiny::reactive({
      if (!is.null(loaded_omics_net_list())) {
        return(loaded_omics_net_list())
      } else {
        return(NULL)
      }
    })
    
    multiplex_network <- shiny::eventReactive(input$gen_multiplex_btn, {
      print('Generating the multiplex network...')
      net_list <- if (!is.null(g_net_list())) g_net_list() else l_net_list()
      
      if(input$weightMultiplex == 'YES'){
        multiplex <- MoNETA::create_multiplex(net_list, weighted = TRUE)
        multiplex_type <- 'Weighted'
        if(input$pruneMultiplex == 'YES'){
          if (!is.numeric(input$pruneMultiplex_th)) {
            shinyalert::shinyalert(title = 'Error', text = 'The threshold must be a real number.',closeOnClickOutside = TRUE,type = 'error')
            return(NULL)
          }
          multiplex <- MoNETA::prune_multiplex_network(multiplex, input$pruneMultiplex_th)
        }else{
          multiplex <- multiplex
        }
      }else{
        multiplex <- MoNETA::create_multiplex(net_list, weighted = FALSE)
        multiplex_type <- 'Unweighted'
      }
      
      shinyalert::shinyalert(title = 'Success', type = 'success',closeOnClickOutside = TRUE,
                             text = paste(multiplex_type,' Multiplex Omics Network created. <br/> Now press <b> "Next" </b> in the top left-hand corner to continue'), html = TRUE)
      return(multiplex)
    })
    
    
    ##############################################################################
    g_dnet_list <- shiny::reactive({
      if (length(shiny::reactiveValuesToList(drug_networks)) != 0) {
        return(shiny::reactiveValuesToList(drug_networks))
      }else{
        return(NULL)
      }
    })
    
    l_dnet_list <-  shiny::reactive({
      if (!is.null(loaded_drug_net_list())) {
        return(loaded_drug_net_list())
      } else {
        return(NULL)
      }
    })
    
    multiplex_drug_network <- shiny::eventReactive(input$gen_multiplex_drug_btn, {
      print('Generating the multiplex drug network...')
      net_list <- if (!is.null(g_dnet_list())) g_dnet_list() else l_dnet_list()
      
      if(input$drug_weightMultiplex == 'YES'){
        multiplex <- MoNETA::create_multiplex(net_list, weighted = TRUE)
        multiplex_type <- 'Weighted'
        if(input$drug_pruneMultiplex == 'YES'){
          if (!is.numeric(input$drug_pruneMultiplex_th)) {
            shinyalert::shinyalert(title = 'Error', text = 'The threshold must be a real number.',closeOnClickOutside = TRUE,type = 'error')
            return(NULL)
          }
          multiplex <- MoNETA::prune_multiplex_network(multiplex, input$drug_pruneMultiplex_th)
        }else{
          multiplex <- multiplex
        }
      }else{
        multiplex <- MoNETA::create_multiplex(net_list, weighted = FALSE)
        multiplex_type <- 'Unweighted'
      }
      
      shinyalert::shinyalert(title = 'Success', type = 'success',closeOnClickOutside = TRUE,
                             text = paste(multiplex_type,' Multiplex Drug Network created. <br/> Now press <b> "Next" </b> in the top left-hand corner to continue'), html = TRUE)
      return(multiplex)
    })
    
    
    ############################################################################################
    #                                   RWR: PARAMETERS (INPUT)                                #
    ############################################################################################
    
    shinyjs::addCssClass(selector = "a[data-value='rwr_tab_3']", class = "inactiveLink")
    observe({
      if (isTruthy(input$gen_multiplex_btn)){
        shinyjs::removeCssClass(selector = "a[data-value='rwr_tab_3']", class = "inactiveLink")
      }else{
        return(NULL)
      }
    })
    
    ######################## RWR-SIM ###########################
    
    RWR_output <- shiny::reactiveValues()
    shiny::observeEvent(input$rwr_button, {
      if (is.null(multiplex_network())){
        return(NULL)
      } else{
        
        omics_multiplex <- multiplex_network()
        drug_multiplex <- multiplex_drug_network() %>%  mutate_at(2:3, tolower)
        
        bnet <- bnetwork()
        colnames(bnet) <- c('source', 'target', 'weight')
        f_bnet <- bnet %>% filter(source %in% omics_multiplex$source, 
                                  tolower(target) %in% drug_multiplex$source)
        
        if (sum(is.na(omics_tau_list())) != 1) {
          if (sum(unlist(omics_tau_list())) != 1) {
            shinyalert::shinyalert("Error", "The sum of restarting probabilities per layer (&tau;) must be equal to 1.",closeOnClickOutside = TRUE, type = "error", html = TRUE)
            return(NULL)
          }
        }
        
        max_cores <- parallel::detectCores()
        if (!is.numeric(input$cores_rwr) | input$cores_rwr < 1 | input$cores_rwr > max_cores){
          shinyalert::shinyalert(title = 'Error', text = paste('Select a number of cores between 1 and', max_cores, '.'),closeOnClickOutside = TRUE, type = 'error')
          return(NULL)
        }
        
        
        if (input$bioInf_transition == 'YES') {
          net_list <- if (!is.null(g_net_list())) g_net_list() else l_net_list()
          omics_trans_mat <- MoNETA::create_layer_transition_matrix(net_list)
        }else {
          omics_trans_mat <- input$editable_trans_mat
        }
        
        drugs_trans_mat <- input$editable_drugs_trans_mat
        
        
        shinyalert::shinyalert(
          title = "Wait",
          text = "Waiting for RWR similarity matrix generation. <br/> Please check the <b> progressbar </b>  in the bottom right-hand corner.",
          closeOnEsc = TRUE,
          closeOnClickOutside = TRUE,
          html = TRUE,
          type = "info",
          showConfirmButton = TRUE,
          confirmButtonText = "OK",
          confirmButtonCol = "#004192",
          showCancelButton = FALSE,
          imageUrl = "",
          animation = TRUE
        )
        
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        
        progress$set(message = "MiDNE", detail = paste("Doing RWR"), value = 0)
        
        if (is.null(drug_multiplex)) {
          Summary <- list(
            'parameters' = c(
              '################################################################################',
              '#                             SUMMARY: RWR-M parameters                        #',
              '################################################################################',
              '', '',
              paste('- Number of layers:                     ', length(unique(omics_multiplex$EdgeType))), '',
              paste('- Restart probability (r):              ', input$restart), '',
              paste('- Cross-jumping probability (delta):    ', input$omics_delta), '',
              paste('- Restart probability per layer (tau):  ', paste(omics_tau_list(), collapse = ', ')), '',
              paste('- Weighted Multiplex Network:           ', input$weightMultiplex),'',
              paste('- Jumping Neighborhood permission:      ', input$omics_jump_neigh),'',
              paste('- Weighted jump_neigh edges:            ', input$omics_weight_jump_neigh),''
            ),
            'Transition layers matrix' = input$omics_trans_mat
          )
          
          
          rwr_simMat <- MoNETA::gen_sim_mat_M(
            
            network =  omics_multiplex, 
            restart = input$restart,
            delta = input$omics_delta, 
            tau = unlist(omics_tau_list()), 
            layer_transition = omics_trans_mat, 
            
            jump_neighborhood = ifelse(input$omics_jump_neigh == 'YES', TRUE, FALSE),  
            weighted_multiplex = ifelse(input$omics_weight_jump_neigh == 'YES', TRUE, FALSE),
            
            cores = input$cores_rwr)
          
          
        } else {
          
          Summary <- list(
            'parameters' = c(                   
              '################################################################################',
              '#                             SUMMARY: RWR-MH parameters                       #',
              '################################################################################',
              '', '', 
              paste('- Number of layers in multiplex omics network:                     ', length(unique(omics_multiplex[[1]]))),'',
              paste('- Number of layers in multiplex drugs network:                     ', length(unique(drug_multiplex[[1]]))),'',
              paste('- Restart probability (r):                                         ', input$restart),'',
              paste('- Cross-jumping probability among multiplexes (lambda):            ', input$lambda),'',
              paste('- Cross-jumping probability (delta1):                              ', input$omics_delta),'',
              paste('- Cross-jumping probability (delta2):                              ', input$drug_delta),'',
              paste('- Restart probability per layer (tau omics):                       ', paste(omics_tau_list(), collapse = ', ')),'',
              paste('- Restart probability per layer (tau drugs):                       ', paste(drugs_tau_list(), collapse = ', ')),'',
              paste('- Weighted Multiplex Omics Network:                                ', input$weightMultiplex),'',
              paste('- Weighted Multiplex Drugs Network:                                ', input$drug_weightMultiplex),'',
              paste('- Jumping Neighborhood permission in multiplex omics network:      ', input$omics_jump_neigh),'',
              paste('- Jumping Neighborhood permission in multiplex drugs network:      ', input$drug_jump_neigh),'',
              paste('- Weighted jump_neigh edges in multiplex omics network:            ', input$omics_weight_jump_neigh),'',
              paste('- Weighted jump_neigh edges in multiplex drugs network:            ', input$drug_weight_jump_neigh),''
            ), 
            
            'Transition layers matrix (omics network)' = omics_trans_mat,
            'Transition layers matrix (drugs network)' = drugs_trans_mat
          )
          
          fake_nodes <- drug_multiplex[stringr::str_detect(drug_multiplex[[3]], 'fk_'), 3][[1]]
          rwr_simMat_list <- gen_sim_mat_MH(
            
            network1 =  omics_multiplex, 
            network2 = drug_multiplex,
            network_B = f_bnet, 
            
            restart = input$restart, lambda = input$lambda,
            
            delta1 = input$omics_delta, 
            #delta2 = input$drug_delta, 
            delta2 = 0.5, 
            
            tau1 = unlist(omics_tau_list()), 
            tau2 = unlist(drugs_tau_list()),
            
            layer_transition_1 = omics_trans_mat, 
            layer_transition_2 = drugs_trans_mat, 
            
            jump_neighborhood_1 = ifelse(input$omics_jump_neigh == 'YES', TRUE, FALSE),  
            #jump_neighborhood_2 = ifelse(input$drug_jump_neigh == 'YES', TRUE, FALSE),
            jump_neighborhood_2 = FALSE,
            
            weighted_multiplex_1 = ifelse(input$omics_weight_jump_neigh == 'YES', TRUE, FALSE),
            #weighted_multiplex_2 = ifelse(input$drug_weight_jump_neigh == 'YES', TRUE, FALSE), 
            weighted_multiplex_2 = FALSE, 
            
            aggregation_method = 'Sum', get_completeRWRmat = FALSE,
            no_seed_nodes = fake_nodes,
            cores = input$cores_rwr)
          
          pre_rwr_simMat <- rwr_simMat_list$RWRMH_sim_mat
          
          nofk_rwr_simMat <- pre_rwr_simMat[!(rownames(pre_rwr_simMat) %in% fake_nodes),]
          rwr_simMat <-  t(t(nofk_rwr_simMat)/colSums(nofk_rwr_simMat))
         
        }
        
        progress$inc(0.8, detail = paste("Creating summary file!"))
        
        RWR_output$rwr_simMat <- rwr_simMat
        RWR_output$Summary <- Summary
        
        shinyalert::shinyalert(title = 'Success', type = 'success',closeOnClickOutside = TRUE,
                               text = 'Random Walk with Restart done! <br/> Now press <b> "Next" </b> in the top left-hand corner to continue', html = TRUE)
      }
    })
    
    ######################## download RWR-SIM mat ###########################
    
    output$download_rwr_mat <- shiny::renderUI({
      if (length(shiny::reactiveValuesToList(RWR_output)) == 0) {
        return(NULL)
      } else {
        shiny::downloadButton(
          outputId = 'download_rwr_btn',
          label = "Download",
          icon = shiny::icon("download")
        )
      }
    })
    
    output$download_rwr_btn <- shiny::downloadHandler(
      filename = function() {
        paste("MiDNEshiny_RWR_mat_", Sys.Date(), ".zip", sep = "")
      },
      content = function(file){
        temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
        dir.create(temp_directory)
        
        to_download <- shiny::reactiveValuesToList(RWR_output)
        for (obj in names(to_download)) {
          if (obj == 'Summary'){
            sink(file = glue::glue(temp_directory, "/{obj}.txt"), append =  TRUE)
            print(to_download[[obj]])
            sink(NULL)
          } else {
            file_name <- "RWR_sim_mat.csv"
            readr::write_csv(as.data.frame(to_download[[obj]]), file.path(temp_directory, file_name))
          }
        }
        zip::zip(
          zipfile = file,
          files = dir(temp_directory),
          root = temp_directory
        )
      }
    )
    
    
    ############################################################################################
    #                        DIMENSIONALITY REDUCTION: LOADING DATA (INPUT)                   #
    ############################################################################################
    
    loaded_simMat <- reactiveVal(NULL)
    shiny::observeEvent(input$load_simMat_button, {
      
      if (input$example_simMat == 'NO'){
        req(input$simMat_file)
        listFiles <- list()
        inFiles <- input$simMat_file
        if (is.null(inFiles)){
          return(NULL)
        } else {
          ext <- tools::file_ext(inFiles$name)
          new_name <- 'simMatFile'
          file <- switch(ext,
                         csv = vroom::vroom(inFiles$datapath, delim = ","),
                         RDS = readRDS(inFiles$datapath),
                         validate("Invalid file; Please upload a .csv or .RDS file")
          )
          if(ncol(file) > nrow(file)){
            rownames(file) <- paste0('component', 1:nrow(file))
          }
          
          if (is.matrix(file)){
            listFiles[[new_name]] <- file
          }else {
            shinyalert::shinyalert("Type Error", "Uploaded Data is not a matrix.", closeOnClickOutside = TRUE, type = "error")
            returnValue()
          }
        }
        loaded_simMat(listFiles)
        
        
      }else{
        listFiles <- list()
        if (input$select_simMat == 'Embedded BRCA-5-omics FDA-drugs similarity matrix'){
          sim_mat <- readRDS('MiDNE/DATA/similarity_matrices/emb_FDA_active_drug_5omics_uRWRMHmat.RDS')
          
        } else{
          sim_mat <- readRDS('MiDNE/DATA/similarity_matrices/umap_emb_FDA_active_drug_5omics_uRWRMHmat.RDS')
          rownames(sim_mat) <- paste0('component', 1:nrow(sim_mat))
        }
        listFiles[['simMatFile']] <- sim_mat
        loaded_simMat(listFiles)
      }
    })
    
    observeEvent(input$load_simMat_button, {
      if (isTruthy(input$simMat_file) & input$example_simMat == 'NO'){
        shinyalert::shinyalert(
          title = "Wait",
          text = "Waiting for data loading",
          size = "xs",
          closeOnEsc = TRUE,
          closeOnClickOutside = TRUE,
          html = TRUE,
          type = "info",
          showConfirmButton = TRUE,
          confirmButtonText = "OK",
          confirmButtonCol = "#004192",
          showCancelButton = FALSE,
          imageUrl = "",
          animation = TRUE
        )
      } else if (!isTruthy(input$simMat_file) & input$example_simMat == 'NO'){
        shinyalert::shinyalert(title = 'Error', text = 'Please select a file to be uploaded.', type = 'error',
                               closeOnEsc = TRUE, closeOnClickOutside = TRUE)
      }
    },  ignoreNULL = FALSE,  ignoreInit = TRUE )
    
    
    #### Load annotation file
    
    loaded_simMat_Anno <- reactiveVal(NULL)
    shiny::observeEvent(input$load_anno_button2, {
      
      if (input$example_simMat_anno == 'NO'){
      
        req(input$anno_file2)
        listFiles <- list()
        inFiles <- input$anno_file2
        if (is.null(inFiles)){
          return(NULL)
        } else {
          ext <- tools::file_ext(inFiles$name)
          new_name <- 'annotation2'
          file <- switch(ext,
                         csv = vroom::vroom(inFiles$datapath, delim = ","),
                         RDS = readRDS(inFiles$datapath),
                         validate("Invalid file; Please upload a .csv or .RDS file")
          )
          if (is.data.frame(file)){
            listFiles[[new_name]] <- file
          }else {
            shinyalert::shinyalert("Type Error", "Uploaded Data is not a dataframe.", closeOnClickOutside = TRUE, type = "error")
            returnValue()
          }
        }
        loaded_simMat_Anno(listFiles)
      }else{
        listFiles <- list()
        file <- readRDS('MiDNE/DATA/annotation/all_genes_drugs_annotation.RDS')
        listFiles[['annotation2']] <- file
        loaded_simMat_Anno(listFiles)
      }
    })
    
    observeEvent(input$load_anno_button2, {
      if(isTruthy(input$anno_file2) & input$example_simMat_anno == 'NO'){
        shinyalert::shinyalert(
          title = "Wait",
          text = "Waiting for data loading",
          size = "xs",
          closeOnEsc = TRUE,
          closeOnClickOutside = TRUE,
          html = TRUE,
          type = "info",
          showConfirmButton = TRUE,
          confirmButtonText = "OK",
          confirmButtonCol = "#004192",
          showCancelButton = FALSE,
          imageUrl = "",
          animation = TRUE
        )
      } else if(!isTruthy(input$anno_file2) & input$example_simMat_anno == 'NO'){
        shinyalert::shinyalert(title = 'Error', text = 'Please select a file to be uploaded.', type = 'error',
                               closeOnEsc = TRUE, closeOnClickOutside = TRUE)
      }
    },  ignoreNULL = FALSE,  ignoreInit = TRUE )
    
    
    
    
    
    
    ############################################################################################
    #                            DIMENSIONALITY REDUCTION: DR (INPUT)                          #
    ############################################################################################
    
    shinyjs::addCssClass(selector = "a[data-value='dr_tab_2']", class = "inactiveLink")
    observe({
      if (length(shiny::reactiveValuesToList(RWR_output)) != 0 | length(loaded_simMat()) != 0){
        shinyjs::removeCssClass(selector = "a[data-value='dr_tab_2']", class = "inactiveLink")
      }else{
        return(NULL)
      }
    })
    
    
    ######################## DR matrices ###########################
    
    dr_output <- shiny::reactiveValues()
    shiny::observeEvent(input$dr_btn, {
      if (length(shiny::reactiveValuesToList(RWR_output)) == 0 & length(loaded_simMat()) == 0){
        return(NULL)
      }else{
        #to_create(NULL)
        
        dr_input <- if (length(shiny::reactiveValuesToList(RWR_output)) != 0) shiny::reactiveValuesToList(RWR_output)$rwr_simMat else loaded_simMat()$simMatFile
        
        if (length(shiny::reactiveValuesToList(RWR_output)) != 0) {
          drType_input <- 'Similarity matrix'
        } else if (input$example_simMat == 'NO') {
          drType_input <- input$simMat_type
        } else {
          drType_input <-  if (input$select_simMat == 'Embedded BRCA-5-omics FDA-drugs similarity matrix') 'Embedded Similarity matrix' else 'Other'
        }
        
        
        if (drType_input == 'Similarity matrix'){ 
          rwr_mat <- dr_input
          req(input$dr_method)
          
          shinyalert::shinyalert(
            title = "Wait",
            text = "It might take a little time at the end. <br/> Please check the <b> progressbar </b> in the bottom right-hand corner.",
            closeOnEsc = TRUE,
            closeOnClickOutside = TRUE,
            html = TRUE,
            type = "info",
            showConfirmButton = TRUE,
            confirmButtonText = "OK",
            confirmButtonCol = "#004192",
            showCancelButton = FALSE,
            imageUrl = "",
            animation = TRUE
          )
          
          progress <- shiny::Progress$new()
          on.exit(progress$close())
          progress$set(message = "MiDNE", detail = paste("Doing Embedding"), value = 0)
          
          if (input$emb_opt == 'YES'){
            emb_mat <- MoNETA::get_embedding(matrix = rwr_mat, embedding_size = input$dim_emb, cores = input$cores_emb)
            dr_output$emb_mat <- emb_mat
          }else{
            emb_mat <- rwr_mat
          }
          
          if (input$dr_method == 'UMAP'){
            progress$inc(0.8, detail = paste("Doing UMAP"))
            dr_mat <- MoNETA::get_parallel_umap_embedding(matrix = emb_mat,
                                                          embedding_size = 2,
                                                          n_threads = input$threads)
          }else if (input$dr_method == 'PCA'){
            progress$inc(0.8, detail = paste("Doing PCA"))
            dr_mat <- MoNETA::get_pca_embedding(matrix = emb_mat, embedding_size = input$emb_size)
          }else if (input$dr_method == 'tSNE'){
            progress$inc(0.8, detail = paste("Doing tSNE"))
            dr_mat <- get_tsne_embedding(matrix = emb_mat,
                                         embedding_size = input$emb_size_tsne,
                                         perplexity = input$perplexity, max_iter = input$max_iter,
                                         num_threads = input$threads)
          }
          
          progress$inc(1, detail = paste("Done!"))
          dr_output$dr_mat <- dr_mat
          
          shinyalert::shinyalert(title = 'Success', type = 'success', closeOnClickOutside = TRUE,
                                 text = 'Dimensionality Reduction step done! <br/> Now press <b> "Next" </b> in the top left-hand corner to continue', html = TRUE)
          
        } else if (drType_input == 'Embedded Similarity matrix') {
          
          progress <- shiny::Progress$new()
          on.exit(progress$close())
          
          emb_mat <- dr_input
          if (input$dr_method == 'UMAP'){
            progress$set(message = "MiDNE", detail = paste("Doing UMAP"), value = 0)
            dr_mat <- MoNETA::get_parallel_umap_embedding(matrix = emb_mat,
                                                          embedding_size = 2,
                                                          n_threads = input$threads)
          }else if (input$dr_method == 'PCA'){
            progress$set(message = "MiDNE", detail = paste("Doing PCA"), value = 0)
            dr_mat <- MoNETA::get_pca_embedding(matrix = emb_mat, embedding_size = input$emb_size)
          }else if (input$dr_method == 'tSNE'){
            progress$set(message = "MiDNE", detail = paste("Doing tSNE"), value = 0)
            dr_mat <- get_tsne_embedding(matrix = emb_mat,
                                         embedding_size = input$emb_size_tsne,
                                         perplexity = input$perplexity, max_iter = input$max_iter,
                                         num_threads = input$threads)
          }
          
          progress$inc(1, detail = paste("Done!"))
          dr_output$dr_mat <- dr_mat
          
          shinyalert::shinyalert(title = 'Success', type = 'success', closeOnClickOutside = TRUE,
                                 text = 'Dimensionality Reduction step done! <br/> Now press <b> "Next" </b> in the top left-hand corner to continue', html = TRUE)
          
        }else if (drType_input == 'Other'){
          dr_mat <- dr_input  
          dr_output$dr_mat <- dr_mat
        }
        
      }
    })
    
    dr_plot <- shiny::reactive({
      x_lab <- '1'
      y_lab <- '2'
      if (length(shiny::reactiveValuesToList(RWR_output)) != 0 | (length(loaded_simMat()) != 0 & input$simMat_type == 'Similarity matrix')){
        method <- ifelse(input$dr_method == 'UMAP', 'UMAP of',
                         ifelse(input$dr_method == 'PCA', 'PCA of', 'tSNE of'))
        emb <- ifelse(input$emb_opt == 'YES', 'the Embedded RWR matrix', 'the RWR matrix')
        components <- 1:2
      } else if (length(loaded_simMat()) != 0 & input$simMat_type == 'Embedded Similarity matrix') {
        method <- ifelse(input$dr_method == 'UMAP', 'UMAP of',
                         ifelse(input$dr_method == 'PCA', 'PCA of', 'tSNE of'))
        emb <- 'the loaded Embedded RWR matrix'
        components <- 1:2
      } else if (length(loaded_simMat()) != 0 & input$simMat_type == 'Other'){
        method <- 'Projection'
        emb <- 'of the loaded low-dimensionality Similarity matrix'
        components <- c(input$select_c1, input$select_c2)
        x_lab <- input$select_c1
        y_lab <- input$select_c2
      }
      
      dr_mat <- shiny::reactiveValuesToList(dr_output)$dr_mat
      
      dr_plot <- MoNETA::plot_2D_matrix(coord = dr_mat[components,], nodes_anno = data.frame(id = colnames(dr_mat)),
                                        id_name = 'id', interactive = FALSE, wo_legend = FALSE) +
        ggplot2::ggtitle(paste(method, emb)) +
        ggplot2::xlab(paste0(input$dr_method, x_lab)) +  ggplot2::ylab( paste0(input$dr_method, y_lab)) +
        ggplot2::theme(text = ggplot2::element_text(family="LMRoman10"),
                       plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = 'bold'),
                       axis.title.x = ggplot2::element_text(size = 10),
                       axis.title.y = ggplot2::element_text(size = 10))
      return(dr_plot)
    })
    
    
    shiny::observeEvent(input$select_c1, {
      if (length(loaded_simMat()) != 0){
        updateSelectInput(inputId = 'select_c2', session, label = 'Select the second component',
                          choices = rownames(loaded_simMat()$simMatFile)[-which(rownames(loaded_simMat()$simMatFile) == input$select_c1)]
        )
      }
    })
    
    
    
    ############################################################################################
    #                                         CLUSTERING (INPUT)                               #
    ############################################################################################
    
    shinyjs::addCssClass(selector = "a[data-value='cl_tab']", class = "inactiveLink")
    observe({
     if (length(shiny::reactiveValuesToList(dr_output)) != 0){
       shinyjs::removeCssClass(selector = "a[data-value='cl_tab']", class = "inactiveLink")
     }else{
       return(NULL)
     }
    })
    
    
    showPlot <- reactiveVal(FALSE)
    
    anno <-  shiny::reactive({
      if (!is.null(annotation())) {
        return(annotation()$annotation)
      } else {
        return(NULL)
      }
    })
    
    anno1 <- shiny::reactive({
      if (!is.null(annotation1())) {
        return(annotation1()$annotation1)
      }else{
        return(NULL)
      }
    })
    
    anno2 <- shiny::reactive({
      if (!is.null(loaded_simMat_Anno())) {
        return(loaded_simMat_Anno()$annotation2)
      }else{
        return(NULL)
      }
    })
    
    to_update <-  shiny::reactiveVal(FALSE)
    to_create <-  shiny::reactiveVal(NULL)
    
    clusters <- shiny::reactiveValues()
    shiny::observeEvent(input$cl_btn, {
      if (length(shiny::reactiveValuesToList(dr_output)) == 0 ){
        return(NULL)
      }else{
        
        to_update(FALSE)
        
        dr_output <- shiny::reactiveValuesToList(dr_output)
        
        if (input$cl_on_emb == 'YES'){
          mat <- dr_output$emb_mat
        } else {
          mat <- dr_output$dr_mat
        }
        
        clusters$dr_mat <- dr_output$dr_mat
        t_mat <- t(mat)
        gPlot_mat <- dplyr::tibble("id" = rownames(t_mat), 'x' = t_mat[,1], 'y' = t_mat[,2])
        
        shinyalert::shinyalert(
          title = "Wait",
          text = paste0("It might take a little time at the end of the clustering analysis. <br/>",
                        "Please check the <b> progressbar </b>  in the bottom right-hand corner."
          ),
          closeOnEsc = TRUE,
          closeOnClickOutside = TRUE,
          html = TRUE,
          type = "info",
          showConfirmButton = TRUE,
          confirmButtonText = "OK",
          confirmButtonCol = "#004192",
          showCancelButton = FALSE,
          imageUrl = "",
          animation = TRUE
        )
        
        if (input$cluster_method == "kmeans"){
          #to_create(NULL)
          
          progress <- shiny::Progress$new()
          on.exit(progress$close())
          progress$set(message = "MiDNE", detail = paste("Doing kmeans"), value = 0)
          
          kmeans_clus <- kmeans(t_mat, input[[paste0(input$cluster_method, '_par')]])
          cluster <- gPlot_mat
          cluster$clust <- factor(kmeans_clus$cluster)
          
          progress$inc(1, detail = paste("Done!"))
          clusters[[paste0('kmeans_', input[[paste0(input$cluster_method, '_par')]])]]$cluster <- cluster
          
        } else if (input$cluster_method == "dbscan") {
          #to_create(NULL)
          progress <- shiny::Progress$new()
          on.exit(progress$close())
          progress$set(message = "MiDNE", detail = paste("Doing DBSCAN"), value = 0)
          
          db_clust <- fpc::dbscan(t_mat, input[[paste0(input$cluster_method, '_par')]])
          cluster <- gPlot_mat
          cluster$clust <- factor(db_clust$cluster)
          
          progress$inc(1, detail = paste("Done!"))
          clusters[[paste0('dbscan_', input[[paste0(input$cluster_method, '_par')]])]]$cluster <- cluster
          clusters[[paste0('dbscan_', input[[paste0(input$cluster_method, '_par')]])]]$num_clust <- length(unique(cluster$clust))
          
        } else if (input$cluster_method == "hclust")  {
          progress <- shiny::Progress$new()
          on.exit(progress$close())
          
          if (is.null(to_create())){
            progress$set(message = "MiDNE", detail = paste("Doing hclust"), value = 0)
            message('hclust')
            dendro <- t_mat %>%
              dist() %>%
              hclust()
            
            to_create(dendro)
            clusters$hclust_output <- dendro
            
          }else {
            dendro <- to_create()
            progress$set(message = "MiDNE", detail = paste("Doing Cut tree"), value = 0)
            message('cut tree')
            cutree_res <- cutree(dendro, h = input$hclust_par)
            num_clust <- length(unique(cutree_res))
            cluster <- gPlot_mat
            cluster$clust <- factor(cutree_res)
            clusters[[paste0('hclust_', input$hclust_par)]]$cluster <- cluster
            clusters[[paste0('hclust_', input$hclust_par)]]$num_clust <- num_clust
            
          }
          progress$inc(1, detail = paste("Done!"))
          
        #} else if (input$cluster_method == "annotation" | input$cluster_method == "annotation1" | input$cluster_method == "annotation2"){
        } else if (input$cluster_method == "manual"){
          # anno <- if ( length(annotation()) != 0) anno() else if (length(annotation1()) != 0) anno1() else anno2()
          # nodes <- colnames(mat)
          # id_col_check <- sapply(colnames(anno), FUN = function(x) sum(nodes %in% anno[[x]]))
          # id <- colnames(anno)[which.max(id_col_check)]
          # 
          # shape_col <- if (input$custom_anno_shape_selector != '-') input$custom_anno_shape_selector else NULL
          # shape_col_nm <- if (input$custom_anno_shape_selector != '-') 'shape' else NULL
          # 
          # f_anno <- anno[anno[[id]] %in% nodes, c(id, input$custom_anno_selector, shape_col)]
          # colnames(f_anno) <- c('id', 'clust', shape_col_nm)
          # 
          # cluster <- gPlot_mat %>% dplyr::inner_join(., f_anno, by = 'id')
          # print(cluster)
          
          clusters[[input$cluster_method]]$cluster <- gPlot_mat
          
          #}
        }
        
        shinyalert::shinyalert(
          title = "Success",
          text = paste0("Clustering Analysis done!"
          ),
          closeOnEsc = TRUE,
          closeOnClickOutside = TRUE,
          html = TRUE,
          type = "success",
          showConfirmButton = TRUE,
          confirmButtonText = "OK",
          confirmButtonCol = "#004192",
          showCancelButton = FALSE,
          imageUrl = "",
          animation = TRUE
        )
      }
    })
    
    
    shiny::observeEvent(input$add_cluster_btn, {
      
      if (!is.null(manual_cluster())) {
        selected_cluster <- manual_cluster() %>% 
          mutate(clust = rep(1, nrow(.))) %>% 
          select(-selected_)
        clusters$manual_cluster <- selected_cluster
        print(selected_cluster)
        print('added!!!')
        
        shinyalert::shinyalert(
          title = "Success",
          text = paste0("Cluster added!"
          ),
          closeOnEsc = TRUE,
          closeOnClickOutside = TRUE,
          html = TRUE,
          type = "success",
          showConfirmButton = TRUE,
          confirmButtonText = "OK",
          confirmButtonCol = "#004192",
          showCancelButton = FALSE,
          imageUrl = "",
          animation = TRUE
        )
      }
    })
    
    
    
    shiny::observeEvent(input$rm_cl_btn, {
      
      for (to_rm in input$dw_rm_cl_opt){
        clusters[[to_rm]] <- NULL
        
      }
      message(paste(input$dw_rm_cl_opt, ' removed, '))
    })
    
    
    output$h <-
      shiny::renderUI({
        req(to_create())
        if( !is.null(to_create())){
          tree <- as.dendrogram(to_create())
          shiny::sliderInput(inputId = "hclust_par", label = 'Select h to cut the dendrogram', min = 0,
                             max = round(base::attr(tree, "height"), digits = 2),
                             value = round(base::attr(tree, "height"), digits = 2), step = 0.01)
          }
      })
    
    title <-  shiny::eventReactive(input$cl_btn, {
      param <- if (isTruthy(input[[paste0(input$cluster_method, '_par')]])) input[[paste0(input$cluster_method, '_par')]] else NULL
      clusters_type <- shiny::reactiveValuesToList(clusters)[[paste0(input$cluster_method, '_', param)]]
      
      if (input$cluster_method == "kmeans"){
        return(paste0('Kmeans clustering (k: ', param, ')' ))
      } else if (input$cluster_method == "dbscan") {
        return(paste0('DBSCAN clustering (epsilon: ', param, ', #cluster: ', clusters_type$num_clust, ')' ))
      } else if (input$cluster_method == "hclust") {
        return(paste0('Hierarchical clustering (h: ', param, ', #cluster: ', clusters_type$num_clust, ')' ))
      } else {
        method <- ifelse(input$dr_method == 'UMAP', 'UMAP of',
                         ifelse(input$dr_method == 'PCA', 'PCA of', 'tSNE of'))
        emb <- ifelse(input$emb_opt == 'YES', 'the Embedded RWR matrix', 'the RWR matrix')
        return(paste(method, emb, '(color:', input$custom_anno_selector, ', shape:', input$custom_anno_shape_selector, ')'))
      }
    })
    
    
    cl_plot <- shiny::reactive({
      
      dr_mat <- shiny::reactiveValuesToList(clusters)$dr_mat
      param <- if (isTruthy(input[[paste0(input$cluster_method, '_par')]])) input[[paste0(input$cluster_method, '_par')]] else NULL
      if (!is.null(shiny::reactiveValuesToList(clusters)[[paste0(input$cluster_method, '_', param)]]$cluster)){
        cluster <- shiny::reactiveValuesToList(clusters)[[paste0(input$cluster_method, '_', param)]]$cluster
        clu_anno <- cluster %>% dplyr::arrange(clust)
        
      } else if (!is.null(shiny::reactiveValuesToList(clusters)$manual_cluster)) {
        cluster <- shiny::reactiveValuesToList(clusters)$manual_cluster
        dr_tab <- tibble( id = colnames(dr_mat), x = dr_mat[1,], y = dr_mat[2,], 
                          clust = 2)
        dr_tab[dr_tab$id %in% cluster$id, 'clust'] <- 1
        clu_anno <- dr_tab %>%  mutate_at(4, as.factor)
        
      } else{
        return(NULL)
      }
      
      cl_plot <- MoNETA::plot_2D_matrix(coord = dr_mat[1:2,], nodes_anno = clu_anno,
                                        id_name = 'id', id_anno_color = 'clust',
                                        interactive = FALSE, wo_legend = FALSE) +
        ggplot2::ggtitle(title()) +
        ggplot2::xlab(paste0(input$dr_method, '1')) +  ggplot2::ylab( paste0(input$dr_method, '2')) +
        ggplot2::theme(text = ggplot2::element_text(family="LMRoman10"),
                       plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = 'bold'),
                       axis.title.x = ggplot2::element_text(size = 10),
                       axis.title.y = ggplot2::element_text(size = 10))
    })
    
 
    
    ############################################################################################
    #                        ENRICHMENT ANALYSIS: CLUSTER2PATHWAY (INPUT)                     #
    ############################################################################################
    
    shinyjs::addCssClass(selector = "a[data-value='clupath_tab']", class = "inactiveLink")
    observe({
      if (length(shiny::reactiveValuesToList(clusters)) != 0){
        shinyjs::removeCssClass(selector = "a[data-value='clupath_tab']", class = "inactiveLink")
      }else{
        return(NULL)
      }
    })
    
    
    gProfiler_res <- shiny::reactiveValues()
    shiny::observeEvent(input$pea_btn, {
 
      if (!is.null(shiny::reactiveValuesToList(clusters)) #& is.null(shiny::reactiveValuesToList(gProfiler_res)[[input$select_cluList]])
      ){
        isolate({
          cl_table <- shiny::reactiveValuesToList(clusters)[[input$select_cluList]]
          
          if (input$select_cluList != 'manual_cluster'){
            cluster_anno <- cl_table$cluster %>% unstack(., id ~ clust)
          }else{
            cluster_anno <- list('1' = cl_table$id)
          }
          
          
          progress <- shiny::Progress$new()
          on.exit(progress$close())
          progress$set(message = "MiDNE", detail = paste("Doing Pathway Enrichment Analysis"), value = 0)
          message('Compute gProfiler...')
          enrich_analysis <- gprofiler2::gost(query = cluster_anno,
                                              organism = "hsapiens",
                                              ordered_query = FALSE,
                                              multi_query = FALSE
          )
          
          enrich_res <- enrich_analysis$result %>%
            select(query, term_name, term_size, query_size, intersection_size, source, p_value) %>%
            group_by(query) %>% dplyr::mutate_at('query', as.integer) %>% arrange(query)
          message('Done!')
          
          
          progress$inc(1, detail = paste("Done!"))
          
          gProfiler_res[[input$select_cluList]] <- enrich_res
          
        })
      }else {
        return(NULL)
      }
    })
    
    
    shiny::observeEvent(input$rm_pea_btn, {
      
      for (to_rm in input$dw_rm_pea_opt){
        gProfiler_res[[to_rm]] <- NULL
      }
      message(paste(input$dw_rm_pea_opt, ' removed, '))
    })
    
    
    topPath <- shiny::eventReactive(input$showOne, {
      isolate({
        if (is.null(input$anno_source) & is.null( shiny::reactiveValuesToList(gProfiler_res))){
          return(NULL)
        }else{
          
          pea_table <- shiny::reactiveValuesToList(gProfiler_res)[[input$select_peaTable]]
          top_enrich_path <- pea_table %>%
              dplyr::filter(source %in% input$anno_source) %>%
              dplyr::slice(which.min(p_value)) %>% ungroup()
          
          if (input$select_peaTable != 'manual_cluster'){
            top_anno <- shiny::reactiveValuesToList(clusters)[[input$select_peaTable]]$cluster %>%
              dplyr::mutate('query' = as.integer(clust)) %>%
              dplyr::left_join(., top_enrich_path, by = 'query', multiple = "all")
          }else{
            top_anno <- shiny::reactiveValuesToList(clusters)[[input$select_peaTable]] %>%
              dplyr::mutate('query' = as.integer(clust)) %>%
              dplyr::left_join(., top_enrich_path, by = 'query', multiple = "all")
          }
          
          
          
          
          
          return(list(top_anno = top_anno,
                      top_enrich_path = top_enrich_path))
        }
      })
    })
    
    ############################################################################################
    #                        ENRICHMENT ANALYSIS: PATHWAY2CLUSTERS (INPUT)                      #
    ############################################################################################
    shinyjs::addCssClass(selector = "a[data-value='pathclu_tab']", class = "inactiveLink")
    observe({
      if (length(shiny::reactiveValuesToList(gProfiler_res)) != 0){
        shinyjs::removeCssClass(selector = "a[data-value='pathclu_tab']", class = "inactiveLink")
      }else{
        return(NULL)
      }
    })
    
    
    path_annotation <- eventReactive(input$path2clust, {
      if (is.null(input$source) & is.null( shiny::reactiveValuesToList(gProfiler_res)))
        return(NULL)
      isolate({
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "MiDNE", detail = paste("Doing Pathway2Clusters"), value = 0)
        message('Running pathway2cluster function')
        
        pea_table <- shiny::reactiveValuesToList(gProfiler_res)[[input$select_peaTable2]]
        cl_table <- shiny::reactiveValuesToList(clusters)[[input$select_peaTable2]]
        
        if (input$select_peaTable != 'manual_cluster'){
          gene_cluster <- cl_table$cluster %>%
            mutate_at('clust', as.integer) %>%
            select(id, clust) %>% arrange(clust)
        }else{
          gene_cluster <- cl_table %>%
            dplyr::mutate_at('clust', as.integer) %>%
            select(id, clust) %>% arrange(clust)
        }
        
        path_cluster <- pea_table %>%
          filter(source == input$source) %>%
          select(query, term_name)
        
        message('DONE!')
        progress$inc(1, detail = paste("Done!"))
        return(list(path_cluster = path_cluster ,
                    gene_cluster = gene_cluster))
      })
      
    })
    
    
    ############################################################################################
    #                                   DRUG DISCOVERY (INPUT)                                 #
    ############################################################################################
    
    shinyjs::addCssClass(selector = "a[data-value='dd_tab']", class = "inactiveLink")
    observe({
      if (length(shiny::reactiveValuesToList(clusters)) != 0 & (!is.null(anno()) |  !is.null(anno1()) | !is.null(anno2()) ) ){
        anno <- if (!is.null(anno())) anno()  else if (!is.null(anno1())) anno1() else anno2()
        id_col_check <- sapply(colnames(anno), FUN = function(x) sum('drug' %in% anno[[x]]))
        type_col <- colnames(anno)[which.max(id_col_check)]
        
        drugs <- anno[anno[[type_col]] == 'drug',][['id']] 
        dr_mat <- shiny::reactiveValuesToList(dr_output)$dr_mat
        drugs_to_keep <- drugs[drugs %in%  colnames(dr_mat) ]
        
        if (length(drugs_to_keep) != 0) {
          shinyjs::removeCssClass(selector = "a[data-value='dd_tab']", class = "inactiveLink")
        }
      }else{
        return(NULL)
      }
    })
    
    
    drugs <- shiny::eventReactive(input$drug_button, {
      
      #if (is.null( multiplex_drug_network()) | is.null(loaded_simMat_Anno()$annotation2) ) {
      #  return(NULL)
      #}else{ 
      #if (!is.null( multiplex_drug_network() )) {
      #  drug_multiplex <- multiplex_drug_network()
      #  all_drugs <- sort(base::union(drug_multiplex[[2]], drug_multiplex[[3]]))
      #  drugs <- all_drugs[!stringr::str_detect(all_drugs, 'fk_')]
      #}else{
      message('Running drug approach')
      
      anno <- loaded_simMat_Anno()$annotation2 
      id_col_check <- sapply(colnames(anno), FUN = function(x) sum('drug' %in% anno[[x]]))
      type_col <- colnames(anno)[which.max(id_col_check)]
      
      drugs <- anno[anno[[type_col]] == 'drug',][['id']] 
      dr_mat <- shiny::reactiveValuesToList(dr_output)$dr_mat
      drugs_to_keep <- drugs[drugs %in%  colnames(dr_mat) ]
    
      #}
      return(drugs_to_keep)
      #}
    })
    
    
    ############################################################################################
    
    
    
    
    ############################################################################################
    #                            OMICS DATA : LOADING DATA (OUTPUT)                      #
    ############################################################################################
    
    ############## outputs loading info #############
    
    output$omics_sum_info <- shiny::renderUI({
      input$load_mat_button
      if (length(omics_files()) != 0){
        nsamples <- sapply(omics_files(), ncol)
        nfeatures <- sapply(omics_files(), nrow)
        shinydashboard::box(width = 12,
                            shiny::HTML(" "),
                            shiny::HTML(paste0(
                              paste("<b>", names(omics_files()), "</b> "), "loaded: ",
                              paste("<b>", nfeatures, "</b>"), " features and",
                              paste("<b>", nsamples, "</b>"), " samples <br/>")
                            )
        )
      } else {
        shiny::HTML("<br/>")
      }
    })
    
    output$drug_sum_info <- shiny::renderUI({
      input$load_drug_mat_button
      if (length(drug_files()) != 0){
        ncols <- sapply(drug_files(), ncol)
        nrows <- sapply(drug_files(), nrow)
        shinydashboard::box(width = 12,
                            shiny::HTML(" "),
                            shiny::HTML(paste0(
                              paste("<b>", names(drug_files()), "</b> "), "loaded: ",
                              paste("<b>", ncols, "</b>"), " columns and",
                              paste("<b>", nrows, "</b>"), " rows <br/>")
                            )
        )
      } else {
        shiny::HTML("<br/>")
      }
    })
    
    output$anno_info <- shiny::renderUI({
      input$load_anno_button
      if ( length(annotation()) != 0 )
      {
        ncolumns <- ncol(annotation()$annotation)
        nrows <- nrow(annotation()$annotation)
        shinydashboard::box(width = 12,
                            shiny::HTML(paste0(paste("<b>", 'Annotation file', "</b>"), " loaded: ",  paste("<b>", ncolumns, "</b>"), " columns and ", paste("<b>", nrows, "</b>"), " rows <br/>")),
                            shiny::HTML(paste0(paste("<b> Colnames </b>"), "<br/>")),
                            shiny::HTML(paste('&ensp;', colnames(annotation()$annotation), "<br/>"))
        )
      } else {
        shiny::HTML("<br/>")
      }
    })
    
    
    output$anno_info_extra <- shiny::renderUI({
      #input$load_anno_button
      if ( length(annotation()) != 0 & !is.null(omics_files())){
        genes <- lapply(omics_files(), rownames) %>% unlist() %>% unique()
        anno <- annotation()$annotation
        id_col_check <- sapply(colnames(anno), FUN = function(x) sum(genes %in% anno[[x]]))
        id <- colnames(anno)[which.max(id_col_check)]
        
        not_annotated_genes <- genes[which(!(genes %in% anno[[id]]))]
        extra_genes_in_anno <- anno[[id]][which(!(anno[[id]] %in% genes))]
        if (length(not_annotated_genes) == 0 & length(extra_genes_in_anno) == 0){
          return(NULL)
        }else {
          message1 <- NULL
          message2 <- NULL
          if (length(extra_genes_in_anno) != 0){
            message1 <- shiny::HTML(paste('&ensp;', paste("<b>", length(extra_genes_in_anno), "</b>"), " features are NOT PRESENT in the provided omics data collection", "<br/>"))
          }else if (length(not_annotated_genes) != 0){
            message2 <- shiny::HTML(paste('&ensp;', paste("<b>", length(not_annotated_genes), "</b>"), " features are NOT ANNOTATED in the provided annotation file", "<br/>"))
          }
          shinydashboard::box(width = 12,
                              shiny::HTML(paste("<b>", p('WARNING!', style = "color:red"), "</b>")), message1, message2
          )
        }
      } else {
        shiny::HTML("<br/>")
      }
    })
    
    ############################################################################################
    #                               OMICS DATA : PRE-PROCESSING (OUTPUT)                       #
    ############################################################################################
    
    ############## omics matrices pre-processing #############
 
    proc_matrices <- shiny::reactiveValues()
    
    # Reactive per definire i bottoni
    buttons <- shiny::reactive({
      if (length(omics_files()) != 0)
        paste0("process_mat_", names(omics_files()))
    })
    
    # ReactiveValues per tracciare lo stato di attivazione dei bottoni
    button_states <- reactiveValues()
    
    shiny::observeEvent(
      eventExpr = {
        req(buttons()) # Assicurati che i bottoni siano disponibili
        list_of_buttons = NULL
        for(var in buttons()) {
          list_of_buttons <- append(list_of_buttons, input[[var]])
        }
        list_of_buttons
      },
      handlerExpr = {
        # Verifica se almeno un bottone è stato attivato
        if (any(button_states$activated)) {
          print(any(button_states$activated))
          print(button_states$activated)
          for (x in names(omics_files())){
            processed_omics_mat <- omics_files()[[x]]
            if (input[[paste0("remove0col_", x)]] == 'YES')
              processed_omics_mat <- MoNETA::remove_zeros_cols(processed_omics_mat)
            if (input[[paste0("norm_", x)]] == 'YES')
              processed_omics_mat <- MoNETA::normalize_omics(processed_omics_mat)
            proc_matrices[[paste0('p_', x)]] <- processed_omics_mat
          }
          message('All matrices processed!')
          
          if (length(omics_files()) == 1){
            shinyalert::shinyalert(
              title = "Success",
              text = 'Processing step done! <br/> Now press <b> "Next" </b> in the top right-hand corner to continue.',
              closeOnEsc = TRUE,
              closeOnClickOutside = TRUE,
              html = TRUE,
              type = "success",
              showConfirmButton = TRUE,
              confirmButtonText = "OK",
              confirmButtonCol = "#004192",
              showCancelButton = FALSE,
              imageUrl = "",
              animation = TRUE
            )
          } else if (length(omics_files()) > 1){
            shinyalert::shinyalert(
              title = "Success",
              text = 'Processing step done! <br/>  Now move on to the <b> Intersection </b> step.',
              closeOnEsc = TRUE,
              closeOnClickOutside = TRUE,
              html = TRUE,
              type = "success",
              showConfirmButton = TRUE,
              confirmButtonText = "OK",
              confirmButtonCol = "#004192",
              showCancelButton = FALSE,
              imageUrl = "",
              animation = TRUE
            )
          }
        }
      },
      ignoreInit = TRUE
    )
    
    observe({
      req(buttons()) # Assicurati che i bottoni siano disponibili
      for (i in seq_along(buttons())) {
        button_states$activated[i] <- input[[buttons()[i]]]
      }
    })
    
    output$process_omics_mat <- shiny::renderUI({
      if (is.null(omics_files()))
        return(NULL)
      else{
        thetabs <- lapply(names(omics_files()), function(x) {
          
          tab <- shiny::tabPanel(title = x,
                                 shiny::fluidRow(
                                   shiny::column(width = 12,
                                                 shiny::radioButtons(inputId = paste0("remove0col_", x),
                                                                     label = shiny::h4(shiny::span("Do you want to remove zeros columns?", style = "font-weight: bold")),
                                                                     choices = c('YES', 'NO'), selected = 'NO'
                                                 ),
                                                 shiny::radioButtons(inputId = paste0("norm_", x),
                                                                     label =  shiny::h4(shiny::span("Do you want to normalize the omics matrix by column?", style = "font-weight: bold")),
                                                                     choices = c('YES', 'NO'), selected = 'NO'
                                                 )
                                   )
                                 ), shiny::hr(),
                                 shiny::fluidRow(
                                   shiny::column(width = 12, shiny::actionButton(paste0('process_mat_', x), "Submit"))
                                 )
          )
          return(tab)
        })
        thetabs$width <- 12
        do.call(shinydashboard::tabBox, thetabs)
      }
    })
    
    ############## Processed matrices info #############
    
    output$pro_matrices_sum_info <- shiny::renderUI({
      processed_mats <- shiny::reactiveValuesToList(proc_matrices)
      shiny::isolate({
        if (length(processed_mats) != 0){
          nsamples <- sapply(processed_mats, ncol)
          nfeatures <- sapply(processed_mats, nrow)
          shinydashboard::box(width = 12,
                              shiny::HTML(paste0(paste("<b>", names(processed_mats), "</b>"), " processed: ",
                                                 paste("<b>", nfeatures, "</b>"), " features and ", paste("<b>", nsamples, "</b>"), " samples <br/>")))
        } else {
          return(NULL)
        }
      })
    })
    
    ############## download processed omics matrices #############
    
    output$download_proc_mat_box <- shiny::renderUI({
      if (length(shiny::reactiveValuesToList(proc_matrices)) == 0) {
        return(NULL)
      } else {
        shinydashboard::box(width = 12, solidHeader = TRUE, collapsible = FALSE,
                            title = shiny::h3(shiny::span('Download', style = "font-weight: bold")),
                            shiny::checkboxGroupInput('download_proc_mat',
                                                      label = shiny::h4(shiny::span('Select one or more processed matrices', style = "font-weight: bold")),
                                                      choices = names(shiny::reactiveValuesToList(proc_matrices))
                            ),
                            
                            shiny::tags$hr(),
                            shiny::downloadButton(
                              outputId = 'download_proc_mat_btn',
                              label = "Download",
                              icon = shiny::icon("download")
                            )
        )
      }
    })
    
    output$download_proc_mat_btn <- shiny::downloadHandler(
      filename = function() {
        paste("MiDNEshiny_processed_mat_", Sys.Date(), ".zip", sep = "")
      },
      content = function(file){
        
        temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
        dir.create(temp_directory)
        
        to_download <- shiny::reactiveValuesToList(proc_matrices)[input$download_proc_mat]
        
        for (obj in names(to_download)) {
          file_name <- glue::glue("{obj}.csv")
          readr::write_csv(as.data.frame(to_download[[obj]]), file.path(temp_directory, file_name))
        }
        
        zip::zip(
          zipfile = file,
          files = dir(temp_directory),
          root = temp_directory
        )
      }
    )
    
    
    ############## not mandatory intersection #############
    
    output$intersection <- shiny::renderUI({
      loaded_mats <- omics_files()
      processed_mats <- shiny::reactiveValuesToList(proc_matrices)
      if (length(loaded_mats) > 1 & length(processed_mats) == length(loaded_mats))
        shinydashboard::box(width = 12,  status = "primary", solidHeader = FALSE,
                            shiny::h2(shiny::span("Omics matrices Intersection", style = "font-weight: bold")),
                            shiny::radioButtons(inputId = 'intersect_opt',
                                                label = shiny::h4(shiny::span('Do you want to intersect omics matrices by features?', style = "font-weight: bold")),
                                                choices = c('YES', 'NO'), selected = 'YES'
                            ),
                            shiny::hr(),
                            shiny::actionButton(inputId = 'intersect_btn', label = 'Submit')
        )
    })
    
    ############## intersection info #############
    
    output$intersection_info <- shiny::renderUI({
      int_mats <- shiny::reactiveValuesToList(intersected_omics_mat)$matrices
      shiny::isolate({
        if (length(int_mats) != 0){
          nsamples <- sapply(int_mats, ncol)
          nfeatures <- sapply(int_mats, nrow)
          shinydashboard::box(width = 12,
                              shiny::HTML(paste0(paste("<b>", names(int_mats), "</b>"), " processed: ",
                                                 paste("<b>", nfeatures, "</b>"), " features and ", paste("<b>", nsamples, "</b>"), " samples <br/>")))
        } else {
          return(NULL)
        }
      })
    })
    
    ###
    output$back2P1 <- renderUI({
      if (length(omics_files()) >= 1 ){
        actionButton('back2P1', label = 'Back', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$back2P1, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "mat_sub_1")
    })
    
    ###
    output$jump2P3 <- renderUI({
      if (length(omics_files()) == 1 && length(reactiveValuesToList(proc_matrices)) != 0){
        actionButton('jump2P3', label = 'Next', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      } else if (length(omics_files()) >= 1 && length(reactiveValuesToList(proc_matrices)) != 0 && isTruthy(input$intersect_btn)){
        actionButton('jump2P3', label = 'Next', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$jump2P3, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "net_sub_1")
    })
    
    
    ############################################################################################
    #                        NETWORK : NETWORK INFERENCE (OUTPUT)                      #
    ############################################################################################
    
    ############## output parameters for OMICS network inference  #############
    
    omics_network_list <- list(NULL)
    gene_networks <- shiny::reactiveValues()
    
    count_net <- shiny::reactiveValues()
    observe({
      if (length(shiny::reactiveValuesToList(intersected_omics_mat)$matrices) == 0){
        if (length(omics_files()) == 1){
          input_names <- names(shiny::reactiveValuesToList(proc_matrices))
        }else{
          return(NULL)
        }
      }else {
        input_names <- names(shiny::reactiveValuesToList(intersected_omics_mat)$matrices)
      }
      
      for (x in input_names){
        count_net[[x]] <- 0
      }
    })
    
    
    output$omics_net_arguments <- renderUI({
      if (length(shiny::reactiveValuesToList(intersected_omics_mat)$matrices) == 0){
        if (length(omics_files()) == 1){
          req(proc_matrices)
          mats <- shiny::reactiveValuesToList(proc_matrices)
        }else{
          return(NULL)
        }
      }else {
        req(intersected_omics_mat)
        mats <- shiny::reactiveValuesToList(intersected_omics_mat)$matrices
      }
      
      thetabs <- lapply(names(mats), function(x) {
        
        omics_network_list[[paste0('generate_omics_net_', x)]] <- observeEvent(input[[paste0('generate_omics_net_', x)]], {
          
          max_cores <- parallel::detectCores()
          
          if (!is.numeric(input[[paste0('cpu_', x)]]) | input[[paste0('cpu_', x)]] < 1 | input[[paste0('cpu_', x)]] > max_cores){
            shinyalert::shinyalert(title = 'Error', closeOnClickOutside = TRUE, text = paste('Select a number of cores between 1 and', max_cores, '.'), type = 'error')
            return(NULL)
          }
          if (input[[paste0('cpu_', x)]] == max_cores){
            shinyalert::shinyalert(title = 'Error', closeOnClickOutside = TRUE,
                                   text = tagList('You have selected the maximum number of cores available.'),
                                   type = 'error', html = TRUE
            )
            return(NULL)
          }
          
          shinyalert::shinyalert(
            title = "Wait",
            text = paste0("Waiting for the inference of", paste("<b>", x, "</b>"), "network. <br/>",
                          "Please check the <b> progressbar </b> in the bottom right-hand corner."
            ),
            closeOnEsc = TRUE,
            closeOnClickOutside = TRUE,
            html = TRUE,
            type = "info",
            showConfirmButton = TRUE,
            confirmButtonText = "OK",
            confirmButtonCol = "#004192",
            showCancelButton = FALSE,
            imageUrl = "",
            animation = TRUE
          )
          
          progress <- shiny::Progress$new()
          on.exit(progress$close())
          
          progress$set(message = "MiDNE", detail = paste("Inferring the network..."), value = 0)
          net <- omics_network_inference(omics_matrix = mats[[x]], 
                                         omics_type = input[[paste0('omics_type_', x)]], 
                                         inference_method = input[[paste0('inference_method_', x)]], 
                                         correction_method = input[[paste0('correction_method_', x)]], 
                                         th = input[[paste0('th_', x)]], 
                                         cpu = input[[paste0('cpu_', x)]])
          
          progress$inc(1, detail = paste("Network constructed!"))
          
          
          if (!is_tibble(net)){
            for (net_name in names(net)){
              print(net_name)
              print(net[[net_name]])
              gene_networks[[paste0(net_name, '_net')]] <- net[[net_name]]
            }
          }else{
            gene_networks[[paste0(x, '_net')]] <- net
          }
          
          print('Omics Network Inference...done!')
          
          count_net[[x]] <- 1
          #print(shiny::reactiveValuesToList(count_net))
          net_created <- sum(unlist(shiny::reactiveValuesToList(count_net)))
          
          
          if (net_created == length(omics_files())){
            shinyalert::shinyalert(
              title = "Success",
              text = paste0(paste("<b>", x, "</b>"), " network created.<br/>",
                            paste("<b> Network inferring status: </b>"), net_created, '/', length(omics_files()),"<br/>",
                            'Now press <b> "Next" </b> in the top right-hand corner to continue.'
              ),
              closeOnEsc = TRUE,
              closeOnClickOutside = TRUE,
              html = TRUE,
              type = "success",
              showConfirmButton = TRUE,
              confirmButtonText = "OK",
              confirmButtonCol = "#004192",
              showCancelButton = FALSE,
              imageUrl = "",
              animation = TRUE
            )
          } else {
            shinyalert::shinyalert(
              title = "Success",
              text = paste0(paste("<b>", x, "</b>"), " network created. <br/>",
                            paste("<b> Network inferring status: </b>"), net_created, '/', length(omics_files())
              ),
              closeOnEsc = TRUE,
              closeOnClickOutside = TRUE,
              html = TRUE,
              type = "success",
              showConfirmButton = TRUE,
              confirmButtonText = "OK",
              confirmButtonCol = "#004192",
              showCancelButton = FALSE,
              imageUrl = "",
              animation = TRUE
            )
          }
          
          return(net)
          
        })
        
        tab <- tabPanel(x, 
                        fluidRow(
                          shiny::column(width = 12, 
                                        selectInput(
                                          inputId = paste0('omics_type_', x),
                                          label = "Select a omics type ",
                                          choices =  c("Genomics", "Epigenomics", "Transcriptomics", 'Proteomics'),
                                        ),
                                        selectInput(
                                          inputId = paste0('inference_method_', x),
                                          label = "Select a network inference method",
                                          choices = list(
                                            "Association measures" = c("Pearson Correlation Coefficient", "Spearman Correlation Coefficient"),
                                            "Co-occurrence measures" = c("Fisher's exact test + post-hoc analysis", 'Other'))
                                        ),
                                        selectInput(
                                          inputId = paste0('correction_method_', x),
                                          label = "Select a network inference method",
                                          choices = c("bonferroni", 'fdr')
                                        ),
                                        numericInput(
                                          inputId = paste0('cpu_', x),
                                          label = "Select the number of cores",
                                          min = 1, max = 200, value = 1
                                        ),
                                        shiny::conditionalPanel(
                                          condition = paste('input.', paste0('omics_type_', x), '== "Epigenomics"'),
                                          sliderInput(
                                            inputId = paste0('th_', x),
                                            label = HTML("Select a threshold to discretized &beta;-values"), min = 0, max = 1, value = 0.3
                                          )
                                        )
                          ),
                          #shiny::column(width = 8,
                          #                shiny::uiOutput(paste0(x, "_spinner")),
                          #                shiny::uiOutput(paste0(x, '_custom_net')))
                        ), 
                        shiny::hr(),
                        shiny::fluidRow(
                          shiny::column(width = 6, shiny::actionButton(paste0('generate_omics_net_', x), paste0("Generate network ", x))),
                          #shiny::column(width = 6, shiny::uiOutput(paste0('update_omics_net_plot_', x)), align = 'right')
                        )
        )
        return(tab)
      })
      thetabs$width <- 12
      do.call(shinydashboard::tabBox, thetabs)
    })
    
    
    ############## output parameters for DRUG network inference  #############
    
    drugs_network_list <- list(NULL)
    drug_networks <- shiny::reactiveValues()
    
    d_count_net <- shiny::reactiveValues()
    observe({
      input_names <- names(drug_files())
      for (x in input_names){
        d_count_net[[x]] <- 0
      }
    })  
    
    
    output$drug_net_arguments <- renderUI({
      if (is.null(drug_files())){
        return(NULL)
      }else{ 
        
        drugtabs <- lapply(names(drug_files()), function(x) {
          
          drugs_network_list[[paste0('generate_drug_net_', x)]] <- observeEvent(input[[paste0('generate_drug_net_', x)]], {
            
            shinyalert::shinyalert(
              title = "Wait",
              text = paste0("Waiting for the inference of", paste("<b>", x, "</b>"), "network. <br/>",
                            "Please check the <b> progressbar </b> in the bottom right-hand corner."
              ),
              closeOnEsc = TRUE,
              closeOnClickOutside = TRUE,
              html = TRUE,
              type = "info",
              showConfirmButton = TRUE,
              confirmButtonText = "OK",
              confirmButtonCol = "#004192",
              showCancelButton = FALSE,
              imageUrl = "",
              animation = TRUE
            )
            
            if (input[[paste0('drug_type_', x)]] == "Fake Nodes"){
              dist <- NULL
              cpu <- NULL
            }else{
              dist <- input[[paste0('distance_measure_', x)]]
              cpu <- input[[paste0('cores_', x)]]
            }
            
            progress <- shiny::Progress$new()
            on.exit(progress$close())
            
            progress$set(message = "MiDNE", detail = paste("Inferring the network..."), value = 0)
            net <- drug_network_inference(drug_input = drug_files()[[x]], 
                                          drug_type = input[[paste0('drug_type_', x)]], 
                                          dist = dist, 
                                          cpu = cpu)
            
            progress$inc(1, detail = paste("Network constructed!"))
            drug_networks[[paste0(x, '_net')]] <- net
            print('Drug Network Inference...done!')
            
            d_count_net[[x]] <- 1
            #print(shiny::reactiveValuesToList(d_count_net))
            net_created <- sum(unlist(shiny::reactiveValuesToList(d_count_net)))
            
            
            if (net_created == length(drug_files())){
              shinyalert::shinyalert(
                title = "Success",
                text = paste0(paste("<b>", x, "</b>"), " network created.<br/>",
                              paste("<b> Network inferring status: </b>"), net_created, '/', length(drug_files()),"<br/>",
                              'Now press <b> "Next" </b> in the top right-hand corner to continue.'
                ),
                closeOnEsc = TRUE,
                closeOnClickOutside = TRUE,
                html = TRUE,
                type = "success",
                showConfirmButton = TRUE,
                confirmButtonText = "OK",
                confirmButtonCol = "#004192",
                showCancelButton = FALSE,
                imageUrl = "",
                animation = TRUE
              )
            } else {
              shinyalert::shinyalert(
                title = "Success",
                text = paste0(paste("<b>", x, "</b>"), " network created. <br/>",
                              paste("<b> Network inferring status: </b>"), net_created, '/', length(drug_files())
                ),
                closeOnEsc = TRUE,
                closeOnClickOutside = TRUE,
                html = TRUE,
                type = "success",
                showConfirmButton = TRUE,
                confirmButtonText = "OK",
                confirmButtonCol = "#004192",
                showCancelButton = FALSE,
                imageUrl = "",
                animation = TRUE
              )
            }
            return(net)
            
          })
          
          tab <- tabPanel(x, 
                          fluidRow(
                            shiny::column(width = 12, 
                                          selectInput(
                                            inputId = paste0('drug_type_', x),
                                            label = "Select a drug network type ",
                                            #choices =  c("Mechanism of action", "SMILES", "Fake Nodes"),
                                            choices =  "Fake Nodes",selected = "Fake Nodes"
                                          ),
                                          shiny::conditionalPanel(
                                            condition = paste0('input.drug_type_', x, '!="Fake Nodes"'),
                                            selectInput(
                                              inputId = paste0('distance_measure_', x),
                                              label = paste("Select a network inference method for ", x),
                                              choices = c("Levenshtein pairwise distance", "Hamming-Ipsen-Mikhailov distance", "Shortest Path")
                                            ),
                                            numericInput(
                                              inputId = paste0('cores_', x),
                                              label = "Select the number of cores",
                                              min = 1, max = 200, value = 1
                                            )
                                          )
                            )
                          ),
                          hr(),
                          fluidRow(
                            shiny::column(width = 6, actionButton(paste0('generate_drug_net_', x), "Generate network"))
                          )
          )
          return(tab)
        })
        drugtabs$width <- 12
        do.call(shinydashboard::tabBox, drugtabs)
      }
    })
    
    
    ############## info boxes  #############
    
    output$omics_net_info <- shiny::renderUI({
      if (length(reactiveValuesToList(gene_networks)) != 0){ 
        gene_net <- reactiveValuesToList(gene_networks)
        n_nodes <- sapply(gene_net, function(x) { length(base::union(x[[1]] ,x[[2]]))   })
        n_edges <- sapply(gene_net, nrow)
        
        shinydashboard::box(width = 12,
                            shiny::HTML(paste0(paste("<b>", names(gene_net), "</b>"), " loaded: ", 
                                               paste("<b>", n_nodes, "</b>"), " nodes and ", 
                                               paste("<b>", n_edges, "</b>"), " edges <br/>")))
      } else {
        shiny::HTML("<br/>")
      }
    })
    
    output$drug_net_info <- shiny::renderUI({
      if (length( reactiveValuesToList(drug_networks)) != 0){ 
        drug_net <- reactiveValuesToList(drug_networks)
        n_nodes <- sapply(drug_net, function(x) { length(base::union(x[[1]] ,x[[2]]))   })
        n_edges <- sapply(drug_net, nrow)
        
        shinydashboard::box(width = 12,
                            shiny::HTML(paste0(paste("<b>", names(drug_net), "</b>"), " loaded: ", 
                                               paste("<b>", n_nodes, "</b>"), " nodes and ", 
                                               paste("<b>", n_edges, "</b>"), " edges <br/>")))
        
      } else {
        shiny::HTML("<br/>")
      }
    })
    
    
    ############## download networks  #############
    
    output$download_net_box <- shiny::renderUI({
      if (length(shiny::reactiveValuesToList(gene_networks)) == 0 & length(shiny::reactiveValuesToList(drug_networks)) == 0) {
        return(NULL)
      }else{
        
        shinydashboard::box(shiny::h3(shiny::span('Download', style = "font-weight: bold")), solidHeader = TRUE, collapsible = FALSE, width = 12,
                            shiny::checkboxGroupInput('download_net',
                                                      label = shiny::h4(shiny::span('Select one or more networks', style = "font-weight: bold")),
                                                      choices = c( names(shiny::reactiveValuesToList(gene_networks)), names(shiny::reactiveValuesToList(drug_networks)) )
                            ),
                            shiny::tags$hr(),
                            shiny::downloadButton(
                              outputId = 'download_net_btn',
                              label = "Download",
                              icon = shiny::icon("download")
                            )
        )
      }
    })
    
    output$download_net_btn <- shiny::downloadHandler(
      filename = function() {
        paste("MiDNEshiny_inferred_net_", Sys.Date(), ".zip", sep = "")
      },
      content = function(file){
        
        temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
        dir.create(temp_directory)
        
        to_download <- c(shiny::reactiveValuesToList(gene_networks)[input$download_net], shiny::reactiveValuesToList(drug_networks)[input$download_net])
        
        for (obj in names(to_download)) {
          file_name <- glue::glue("{obj}.csv")
          readr::write_csv(as.data.frame(to_download[[obj]]), file.path(temp_directory, file_name))
        }
        
        zip::zip(
          zipfile = file,
          files = dir(temp_directory),
          root = temp_directory
        )
      }
    )
    
    
    ###
    output$back2P2 <- renderUI({
      if (length(omics_files()) == 1 && length(reactiveValuesToList(proc_matrices)) != 0){
        actionButton('back2P2', label = 'Back', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      } else if (length(omics_files()) >= 1 && length(reactiveValuesToList(proc_matrices)) != 0 && isTruthy(input$intersect_btn)){
        actionButton('back2P2', label = 'Back', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$back2P2, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "mat_sub_2")
    })
    
    ###
    output$jump2P4 <- renderUI({
      if (length(shiny::reactiveValuesToList(gene_networks)) == length(omics_files())){
        actionButton('jump2P4', label = 'Next', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$jump2P4, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "net_sub_2")
    })
    
    
    ############################################################################################
    #                               NETWORK : NETWORK FILTERING (OUTPUT)                       #
    ############################################################################################
    
    ############## Filtering networks box #############
    
    filteredOmicsNetworks_list <- list(NULL)
    filteredOmicsNetworks <- shiny::reactiveValues()
    
    
    count_unet <- shiny::reactiveValues()
    observe({
      input_names <- names(shiny::reactiveValuesToList(gene_networks))
      for (x in input_names){
        count_unet[[x]] <- 0
      }
    })
    
    
    output$plot_net_box <- shiny::renderUI({
      if (length(shiny::reactiveValuesToList(gene_networks)) == 0){
        return(NULL)
      }else
        nets <- shiny::reactiveValuesToList(gene_networks)
      
      thetabs <- lapply(names(nets), function(x) {
        
        output[[paste0(x, '_info')]] <-
          shiny::renderUI({
            input[[ paste0('button_', x)]]
            shiny::isolate({
              net <- nets[[x]]
              
              if (input[[paste0("range_option", x)]])
                fnet <- net[net$weight >= input[[paste0('omics_range', x)]][1] & net$weight <= input[[paste0('omics_range', x)]][2], ]
              else
                fnet <- net[net$weight <= input[[paste0('omics_range', x)]][1] | net$weight >= input[[paste0('omics_range', x)]][2], ]
              
              summary <- netSummary(dplyr::tibble(fnet))
              res <- c()
              for (i in 2:length(summary))
                res[i-1] <- (paste0(paste('<b>', names(summary)[i], '</b>'),
                                    ': <span style="font-family: LM Roman 10;"> ', toString(summary[[i]]), '</span> <br/>'))
              shiny::HTML(paste('<h4>', res, '</h4> '))
            })
          })
        
        filteredOmicsNetworks_list[[paste0('button_', x)]] <- shiny::observeEvent(input[[paste0('button_', x)]], {
          net <- nets[[x]]
          
          if (input[[paste0("range_option", x)]])
            fnet <- net[net$weight >= input[[paste0('omics_range', x)]][1] & net$weight <= input[[paste0('omics_range', x)]][2], ]
          else
            fnet <- net[net$weight <= input[[paste0('omics_range', x)]][1] | net$weight >= input[[paste0('omics_range', x)]][2], ]
          
          output[[paste0(x, "_spinner_weightDist")]] <- renderUI({
            shinycssloaders::withSpinner(plotly::plotlyOutput(paste0(x, "_weightDist")))
          })
          
          output[[paste0(x, "_weightDist")]] <- plotly::renderPlotly({
            #shiny::isolate({
            plotly::ggplotly(ggplot2::ggplot(fnet, ggplot2::aes(x= weight)) + ggplot2::geom_histogram() +
                               ggplot2::ggtitle('Weight Distribution') +
                               ggplot2::xlab("Weight") + ggplot2::ylab("Counts") +
                               ggplot2::theme(text = ggplot2::element_text(family="LM Roman 10"),
                                              plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = 'bold'),
                                              axis.title.x = ggplot2::element_text(size = 10),
                                              axis.title.y = ggplot2::element_text(size = 10))
            )
            #})
          })
          
          
          output[[paste0(x, "_spinner_nodeDegreeDist")]] <- renderUI({
            shinycssloaders::withSpinner(plotly::plotlyOutput(paste0(x, "_nodeDegreeDist")))
          })
          
          output[[paste0(x, "_nodeDegreeDist")]] <- plotly::renderPlotly({
            #shiny::isolate({
            fnet1 <- fnet
            fnet1$weight <- 1
            graph <- igraph::graph_from_data_frame(fnet1, directed = FALSE)
            degree_table <- data.frame(Degree = igraph::degree(graph))
            p <- plotly::ggplotly(ggplot2::ggplot(degree_table, ggplot2::aes(x = Degree)) + ggplot2::geom_histogram() +
                                    ggplot2::ggtitle('Nodes Degree Distribution') +
                                    ggplot2::ylab("Number of nodes") + ggplot2::xlab("Node degree") +
                                    ggplot2::theme(text = ggplot2::element_text(family="LM Roman 10"),
                                                   plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = 'bold'),
                                                   axis.title.x = ggplot2::element_text(size = 10),
                                                   axis.title.y = ggplot2::element_text(size = 10))
            )
            return(p)
            #})
          })
          
        })
        
        net <- nets[[x]]
        min <-  round(min(net$weight), digits = 2)
        max <-  round(max(net$weight), digits = 2)
        
        tab <- shiny::tabPanel(x,
                               shiny::fluidRow(
                                 shiny::column(width = 6,
                                               sliderInput(inputId = paste0("omics_range", x), "Select range:",
                                                           min = min, max = max, value = c(min, max)),
                                               shiny::radioButtons(inputId = paste0("range_option", x), "Do you want to select the inner range?",
                                                                   choices = c(TRUE, FALSE), selected = TRUE)
                                 ),
                                 shiny::column(width = 6,
                                               uiOutput(paste0(x, '_info')),
                                               shiny::tags$hr(style='border-top: 1px solid white;'),
                                               shiny::actionButton(inputId = paste0('button_', x), 'Show')
                                 )
                               ),
                               shiny::tags$hr(),
                               shiny::fluidRow(
                                 shiny::column(width = 6,
                                               shiny::uiOutput(paste0(x, "_spinner_weightDist"))
                                 ),
                                 shiny::column(width = 6,
                                               shiny::uiOutput(paste0(x, "_spinner_nodeDegreeDist"))
                                 )
                               ),
                               shiny::tags$hr(),
                               fluidRow(
                                 shiny::column(width = 12, actionButton(paste0('update_net_', x), 'Submit'), align = 'right')
                               )
                               
        )
        return(tab)
      })
      
      thetabs$width <- 12
      do.call(shinydashboard::tabBox, thetabs)
    })
    
    ############## Filtered networks info #############
    
    output$f_net_sum_info <- shiny::renderUI({
      f_nets <- shiny::reactiveValuesToList(filteredOmicsNetworks)
      shiny::isolate({
        if (length(f_nets) != 0){
          n_nodes <- sapply(f_nets, function(x){length(base::union(x[[1]], x[[2]]))})
          n_edges <- sapply(f_nets, nrow)
          shinydashboard::box(width = 4,
                              shiny::HTML(paste0(paste("<b>", names(f_nets), "</b>"), " filtered: ",
                                                 paste("<b>", n_nodes, "</b>"), " nodes and ",
                                                 paste("<b>", n_edges, "</b>"), " edges <br/>")))
        } else {
          return(NULL)
        }
      })
    })
    
    ############## download filtered networks  #############
    
    output$download_fnet_box <- shiny::renderUI({
      if (length(shiny::reactiveValuesToList(filteredOmicsNetworks)) == 0) {
        return(NULL)
      } else {
        shinydashboard::box(shiny::h3(shiny::span('Download', style = "font-weight: bold")), solidHeader = TRUE, collapsible = FALSE, width = 12,
                            shiny::checkboxGroupInput('download_fnet',
                                                      label = shiny::h4(shiny::span('Select one or more filtered networks', style = "font-weight: bold")),
                                                      choices = names(shiny::reactiveValuesToList(filteredOmicsNetworks))
                            ),
                            shiny::tags$hr(),
                            shiny::downloadButton(
                              outputId = 'download_fnet_btn',
                              label = "Download",
                              icon = shiny::icon("download")
                            )
        )
      }
    })
    
    output$download_fnet_btn <- shiny::downloadHandler(
      filename = function() {
        paste("MiDNEshiny_filtered_net_", Sys.Date(), ".zip", sep = "")
      },
      content = function(file){
        
        temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
        dir.create(temp_directory)
        
        to_download <- shiny::reactiveValuesToList(filteredOmicsNetworks)[input$download_fnet]
        
        for (obj in names(to_download)) {
          file_name <- glue::glue("{obj}.csv")
          readr::write_csv(as.data.frame(to_download[[obj]]), file.path(temp_directory, file_name))
        }
        
        zip::zip(
          zipfile = file,
          files = dir(temp_directory),
          root = temp_directory
        )
      }
    )
    
    ###
    output$back2P3 <- renderUI({
      if (length(shiny::reactiveValuesToList(gene_networks)) == length(omics_files())){
        actionButton('back2P3', label = 'Back', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$back2P3, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "net_sub_1")
    })
    
    ###
    output$jump2P5 <- renderUI({
      if (length(shiny::reactiveValuesToList(filteredOmicsNetworks)) != 0){
        actionButton('jump2P5', label = 'Go to RWR multiplex', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$jump2P5, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "rwr_tab_2")
    })
    
    
    ###
    output$jump2P1.2_from4 <- renderUI({
      if (length(shiny::reactiveValuesToList(filteredOmicsNetworks)) != 0){
        actionButton('jump2P1.2_from4', label = 'Go to Network loading', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$jump2P1.2_from4, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "rwr_tab_1")
    })
    
    
    ############################################################################################
    #                             RWR: LOADING NETWORKS (OUTPUT)                          #
    ############################################################################################
    
    ############## outputs OMICS loading info #############
    
    output$omics_net_sum_info <- shiny::renderUI({
      if (length(loaded_omics_net_list()) != 0){
        n_nodes <- sapply(loaded_omics_net_list(), function(x){length(unique(c(x[[1]], x[[2]])))})
        n_edges <- sapply(loaded_omics_net_list(), nrow)
        shinydashboard::box(width = 12,
                            shiny::HTML(paste0(paste("<b>", names(loaded_omics_net_list()), "</b>"), " loaded: ",  paste("<b>", n_nodes, "</b>"), " nodes and ", paste("<b>", n_edges, "</b>"), " edges <br/>")))
      } else {
        shiny::HTML("<br/>")
      }
    })
    
    ############## outputs DRUGS loading info #############
    
    output$drug_net_sum_info <- shiny::renderUI({
      if (length(loaded_drug_net_list()) != 0){
        n_nodes <- sapply(loaded_drug_net_list(), function(x){length(base::union(x[[1]], x[[2]]))})
        n_edges <- sapply(loaded_drug_net_list(), nrow)
        shinydashboard::box(width = 12,
                            shiny::HTML(paste0(paste("<b>", names(loaded_drug_net_list()), "</b>"), " loaded: ",  paste("<b>", n_nodes, "</b>"), " nodes and ", paste("<b>", n_edges, "</b>"), " edges <br/>")))
      } else {
        shiny::HTML("<br/>")
      }
    })
    
    ############## bnet  info #############
    
    output$bnet_sum_info <- shiny::renderUI({
      if ( shiny::isTruthy(bnetwork()))
      {
        n_genes <- length(unique(bnetwork()[[1]]))
        n_drugs <- length(unique(bnetwork()[[2]]))
        n_edges <- nrow(bnetwork())
        shinydashboard::box(width = 12,
                            shiny::HTML(paste0(paste("<b> Bipartite Network </b>"), " loaded: ",  
                                               paste("<b>", n_genes, "</b>"), " genes, ", 
                                               paste("<b>", n_drugs, "</b>"), " drugs and ", 
                                               paste("<b>", n_edges, "</b>"), " edges <br/>")))
      }else {
        shiny::HTML("<br/>")
      }
    })
    
    
    ############## annotation1  info #############
    
    output$anno_info1 <- shiny::renderUI({
      input$load_anno_button1
      if ( shiny::isTruthy(annotation1()))
      {
        ncolumns <- ncol(annotation1()$annotation1)
        nrows <- nrow(annotation1()$annotation1)
        shinydashboard::box(width = 12,
                            shiny::HTML(paste0(paste("<b>", 'Annotation file', "</b>"), " loaded: ",  paste("<b>", ncolumns, "</b>"), " columns and ", paste("<b>", nrows, "</b>"), " rows <br/>")),
                            shiny::HTML(paste0(paste("<b> Colnames </b>"), "<br/>")),
                            shiny::HTML(paste('&ensp;', colnames(annotation1()$annotation1), "<br/>"))
        )
      } else {
        shiny::HTML("<br/>")
      }
    })
    
    output$anno_info_extra1 <- shiny::renderUI({
      if ( shiny::isTruthy(annotation1()) && length(omics_files()) != 0){
        genes <- lapply(omics_files(), rownames) %>% unlist() %>% unique()
        anno <- annotation1()$annotation1
        id_col_check <- sapply(colnames(anno), FUN = function(x) sum(genes %in% anno[[x]]))
        id <- colnames(anno)[which.max(id_col_check)]
        not_annotated_genes <- genes[which(!(genes %in% anno[[id]]))]
        extra_genes_in_anno <- anno[[id]][which(!(anno[[id]] %in% genes))]
        if (length(not_annotated_genes) == 0 & length(extra_genes_in_anno) == 0){
          return(NULL)
        }else {
          message1 <- NULL
          message2 <- NULL
          if (length(extra_genes_in_anno) != 0){
            message1 <- shiny::HTML(paste('&ensp;', paste("<b>", length(extra_genes_in_anno), "</b>"), " features are not present in the omics data collection", "<br/>"))
          }else if (length(not_annotated_genes) != 0){
            message2 <- shiny::HTML(paste('&ensp;', paste("<b>", length(not_annotated_samples), "</b>"), " features are not annotated in the provided annotation file", "<br/>"))
          }
          shinydashboard::box(width = 12,
                              shiny::HTML(paste("<b>", p('WARNING!', style = "color:red"), "</b>")), message1, message2
          )
        }
      } else {
        shiny::HTML("<br/>")
      }
    })
    
    # output$anno_info_extra2 <- shiny::renderUI({
    #   if ( shiny::isTruthy(annotation2()) && length(loaded_omics_net_list()) != 0){
    #     samples <- lapply(loaded_omics_net_list(), function(x) {base::union(x[,1], x[,2])}) %>% unlist() %>% unique()
    #     
    #     anno <- annotation2()$annotation2
    #     id_col_check <- sapply(colnames(anno), FUN = function(x) sum(samples %in% anno[[x]]))
    #     id <- colnames(anno)[which.max(id_col_check)]
    #     
    #     not_annotated_samples <- samples[which(!(samples %in% anno[[id]]))]
    #     extra_samples_in_anno <- anno[[id]][which(!(anno[[id]] %in% samples))]
    #     if (length(not_annotated_samples) == 0 & length(extra_samples_in_anno) == 0){
    #       return(NULL)
    #     }else {
    #       message1 <- NULL
    #       message2 <- NULL
    #       if (length(extra_samples_in_anno) != 0){
    #         message1 <- shiny::HTML(paste('&ensp;', paste("<b>", length(extra_samples_in_anno), "</b>"), " samples are not present in the omics data collection", "<br/>"))
    #       }else if (length(not_annotated_samples) != 0){
    #         message2 <- shiny::HTML(paste('&ensp;', paste("<b>", length(not_annotated_samples), "</b>"), " samples are not annotated in the provided annotation file", "<br/>"))
    #       }
    #       shinydashboard::box(width = 12,
    #                           shiny::HTML(paste("<b>", p('WARNING!', style = "color:red"), "</b>")), message1, message2
    #       )
    #     }
    #   } else {
    #     shiny::HTML("<br/>")
    #   }
    # })
    
    ###
    output$back2P1_from1.1 <- renderUI({
      actionButton('back2P1_from1.1', label = 'Go to Data', shiny::icon("paper-plane"),
                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
    })
    
    observeEvent(input$back2P1_from1.1, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "mat_sub_1")
    })
    
    ###
    output$back2P4_from1.1 <- renderUI({
      if (length(shiny::reactiveValuesToList(filteredOmicsNetworks)) == length(omics_files())){
        actionButton('back2P4_from1.1', label = 'Back', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$back2P4_from1.1, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "net_sub_2")
    })
    
    
    ###
    output$jump2P5.1 <- renderUI({
      if (shiny::isTruthy(loaded_omics_net_list()) | shiny::isTruthy(loaded_drug_net_list()) | shiny::isTruthy(bnetwork()) ){
        actionButton('jump2P5.1', label = 'Next', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$jump2P5.1, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "rwr_tab_2")
    })
    
    
    ############################################################################################
    #                               RWR: MULTIPLEX NETWORK (OUTPUT)                            #
    ############################################################################################
    
    output$info_box_6 <- renderText({
      HTML("<br/> <span style='font-weight:normal;'>
            <p align='justify'>
            In this section of the app, you can create a multiplex network,
            which is a multilayered network with as many layers as the number of networks inferred in the previous steps. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <b> Weighted Multiplex </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
             If you choose to create a <u>weighted</u> multiplex network, the edge weights of each
             network will be retained and contribute to the node visit probability
             in the next phase of the pipeline. Subsequently,
             you can filter the multiplex network by selecting a threshold.  </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <b> Unweighted Multiplex </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
             If you decide not to weight the multiplex network, the weight column will be filled with <b>1</b>,
             and consequently, each <span style='color:red;'><b>intra-layer</b> </span> edge will have an equal probability of being traversed by the random walker. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>
            ")
    })
    
    ############## output multiplex omics net info #############
    
    output$multi_net_sum_info <- shiny::renderUI({
      input$gen_multiplex_btn
      shiny::isolate({
        multiplex <- multiplex_network()
        if (!is.null(multiplex)){
          multiplex_type <- ifelse(input$weightMultiplex == 'YES', 'Weighted', 'Unweighted')
          n_layers <- multiplex %>%  dplyr::select(EdgeType) %>% unique(.) %>% nrow(.)
          n_nodes <- base::union(multiplex$source, multiplex$target) %>% length(.)
          n_edges <- multiplex %>% nrow(.)
          shinydashboard::box(width = 12,
                              shiny::HTML(paste0(paste("<b>",multiplex_type, "Multiplex network:</b>"), "<br/>")),
                              shiny::HTML(paste('&ensp;', "<b>", n_layers, "</b>", " layers, <br/>")),
                              shiny::HTML(paste('&ensp;',"<b>", n_nodes, "</b>", " nodes <br/>")),
                              shiny::HTML(paste('&ensp;',"<b>", n_edges, "</b>", "edges <br/>"))
          )
        } else {
          shiny::HTML("<br/>")
        }
      })
    })
    
    ############## output multiplex drug net info #############
    
    output$multi_drug_net_sum_info <- shiny::renderUI({
      input$gen_multiplex_drug_btn
      shiny::isolate({
        multiplex <- multiplex_drug_network()
        if (!is.null(multiplex)){
          multiplex_type <- ifelse(input$drug_weightMultiplex == 'YES', 'Weighted', 'Unweighted')
          n_layers <- multiplex %>%  dplyr::select(EdgeType) %>% unique(.) %>% nrow(.)
          n_nodes <- base::union(multiplex[[2]], multiplex[[3]]) %>% length(.)
          n_edges <- multiplex %>% nrow(.)
          shinydashboard::box(width = 12,
                              shiny::HTML(paste0(paste("<b>", multiplex_type, "Multiplex network:</b>"), "<br/>")),
                              shiny::HTML(paste('&ensp;', "<b>", n_layers, "</b>", " layers, <br/>")),
                              shiny::HTML(paste('&ensp;',"<b>", n_nodes, "</b>", " nodes <br/>")),
                              shiny::HTML(paste('&ensp;',"<b>", n_edges, "</b>", "edges <br/>"))
          )
        } else {
          shiny::HTML("<br/>")
        }
      })
    })
    
    ###
    output$back2P4 <- renderUI({
      if (length(shiny::reactiveValuesToList(filteredOmicsNetworks)) == length(omics_files())){
        actionButton('back2P4', label = 'Back', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    output$back2P1.1 <- renderUI({
      if (!is.null(loaded_omics_net_list()) | shiny::isTruthy(loaded_drug_net_list()) | shiny::isTruthy(bnetwork()) ){
        actionButton('back2P1.1', label = 'Back', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$back2P4, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "net_sub_2")
    })
    
    observeEvent(input$back2P1.1, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "rwr_tab_1")
    })
    
    ###
    output$jump2P6 <- renderUI({
      if (isTruthy(input$gen_multiplex_btn)){
        actionButton('jump2P6', label = 'Next', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")}else{
                       return(NULL)
                     }
    })
    
    observeEvent(input$jump2P6, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "rwr_tab_3")
    })
    
    ############################################################################################
    #                                    RWR: PARAMETERS (OUTPUT)                              #
    ############################################################################################
    
    output$info_box_7 <- renderText({
      HTML("<br/> <span style='font-weight:normal;'>
            <p align='justify'>
            Here, you can integrate multi-omics network through the Random Walk with Restart (RWR) algorithm,
            that simulates the traversal of an imaginary particle, called random walker. <br/>
            <hr style='border-top: 1px solid white;'>

            The random walker can:<br/>
            i) walk within a layer of the multiplex through <b><span style='color:red;'>intra-layer</span></b> edges; <br/>
            ii) return back to the <span style='color:red;'>seed</span>  (<b>Restart</b> panel); <br/>
            iii) jump to another node in a different layer of the same network (<b>Transition</b> panel). <br/>
            iv) jump to another node of a different network (<b>Transition</b> panel). <br/>
            </p> </span> <hr style='border-top: 1px solid white;'>

            <b> Restart panel</b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
                 The parameter <b>r</b> represents the probability for the random walker
                 to return back to the seed, while the vector parameter <b>&tau;</b> measures
                 the probability of restarting in the seed of each layer. If you <u>use
                 default restarting probabilities</u>, all layers will have the same probability of being visited during the restart. <br/>
            </p> </span>
             <hr style='border-top: 1px solid white;'>

            <b>Jump panel</b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
                The parameter <b>&delta;</b> represents the general probability for the random walker
                to jump between layers of that network. <br/>
                &ensp; Actually, it is possible to tune the jumping probability
                between all possible pairs of layers by modifing the <b>Transition Layer Matrix</b>.
                If the <u>biologically informed</u> strategy is selected, the random walker can transition between
                omics layers with probability proportional to the number of shared relationships. <br/>
                &ensp; By default, the <b><span style='color:red;'>inter-layers</span></b> edges connect nodes representing the same entity across layers, but
                if you select the <b>RWR-NF</b> option, the edge set will also include their neighborhood across layers.
                Only if RWR-NF is selected, it is possible to weight these edges.
                For further details, please refer to the following paper
                <span style= 'color: blue;'> <u>  <a href='https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04029-3' target='_blank'> (link)</a></u></span>.<br/>
             </p> </span>
             <hr style='border-top: 1px solid white;'>

              <b> Transition panel</b>: <br/> <span style='font-weight:normal;'>
              <p align='justify'>
                 The parameter <b>&lambda;</b> represents the probability for the random walker
                 to transition from the omics network to the drug network, and viceversa. <br/>
             </p> </span>
             <hr style='border-top: 1px solid white;'>


             <span style='font-weight:normal;'>
             <p align='justify'>
             The RWR similarity matrix and the parameter summary file can be downloaded by cliking the <b>download button</b> displayed in the top right-hand corner of the left box. </span>
             </p> <hr style='border-top: 1px solid white;'>
            ")
    })
    
    ################# Omics parameters ################
    
    #### Taus
    output$tauBIO <- shiny::renderUI({
      if (input$tao_opt == 'Custumize restarting probabilities per layer'){
        net_list <- if (!is.null(g_net_list())) g_net_list() else l_net_list()
        numInput <- length(net_list)
        lapply(1:numInput, function(i) {
          shiny::numericInput(
            inputId = paste0('omics_tau', i),
            label = shiny::HTML("&tau;", i, "(", names(net_list)[i] ,")"),
            min = 0,
            max = 1,
            value = 0,
            step = 0.1)
        })
      }else{
        return(NULL)
      }
    })
    
    omics_tau_list <- shiny::reactive({
      if (input$tao_opt != 'Custumize restarting probabilities per layer'){
        return(NA)
      }else{
        numInput <- if (!is.null(g_net_list())) length(g_net_list()) else length(l_net_list())
        lapply(1:numInput, function(i) {
          input[[paste0("omics_tau", i)]]
        })
      }
    })
    
    #### Layer Transition matrix
    observe({
      net_list <- if (!is.null(g_net_list())) g_net_list() else l_net_list()
      m <- shiny::reactive({
        L1 <- length(net_list)
        mat <- matrix(1, ncol = L1 , nrow = L1)
        mat <- mat/(L1-1)
        diag(mat) <- 0
        colnames(mat) <- rownames(mat) <- names(net_list)
        return(mat)
      })
      output$omics_trans_mat <- shiny::renderUI({
        div(
          shinyMatrix::matrixInput(inputId = "editable_trans_mat", value = m(), class = "numeric")
        )
      })
    })
    
    ################# Drug parameters ################
    #### delta 2 --> cross-jumping probability per layer in Drugs multiplex
    output$drug_delta_slider <- renderUI({
      if (length(loaded_drug_net_list()) == 0 & length(shiny::reactiveValuesToList(drug_networks)) == 0)
        return(NULL)
      shiny::sliderInput('drug_delta',
                         shiny::h4(shiny::HTML('<span style:"font-family = LM Roman 10"><b>Select the transition parameter (Drug &delta;) </b></span>')),
                         min = 0, max = 1, value = 0.5, step = 0.01)
    })
    
    #### Taus
    output$tauDRUGS <- renderUI({
      if (length(loaded_drug_net_list()) == 0 & length(shiny::reactiveValuesToList(drug_networks)) == 0){
        return(NULL)
      }else {
        dnet <- if (!is.null(l_dnet_list())) l_dnet_list() else g_dnet_list()
        if (input$drug_tao_opt == 'Custumize restarting probabilities per layer'){
          numInput <- length(dnet)
          lapply(1:numInput, function(i) {
            numericInput(
              inputId = paste0('drugs_tau', i),
              label = shiny::HTML("&tau;", i, "(", names(dnet)[i] ,")"),
              min = 0,
              max = 1,
              value = 0,
              step = 0.1)
          })
        }else{
          return(NULL)
        }
      }
    })
    
    drugs_tau_list <- reactive({
      if (input$drug_tao_opt != 'Custumize restarting probabilities per layer'){
        return(NA)
      }else{
        dnet <- if (!is.null(l_dnet_list())) l_dnet_list() else g_dnet_list()
        numInput <- length(dnet)
        lapply(1:numInput, function(i) {
          input[[paste0("drugs_tau", i)]]
        })
      }
    })
    
    #### Layer Transition matrix
    
    observe({
      if (!is.null(g_dnet_list()) | !is.null(l_dnet_list())) {
        net_list <- if (!is.null(g_dnet_list())) g_dnet_list() else l_dnet_list()
        dm <- reactive({
          L2 <- length(net_list)
          dmat <- matrix(1, ncol = L2, nrow = L2)
          dmat <- dmat/(L2-1)
          diag(dmat) <- 0
          colnames(dmat) <- rownames(dmat) <- names(net_list)
          return(dmat)
        })
        output$drugs_trans_mat <- renderUI({
          div(
            shinyMatrix::matrixInput(inputId = "editable_drugs_trans_mat", value = dm(), class = "numeric")
          )
        })
      }
    })
    
    #### Jump cross prob in drugs multiplex
    output$drug_jump_neigh_opt <- renderUI({
      if (length(loaded_drug_net_list()) == 0 & length(shiny::reactiveValuesToList(drug_networks)) == 0)
        return(NULL)
      output <- tagList()
      output[[1]] <-
        shiny::radioButtons('drug_jump_neigh', 'Connect nodes of different layers by neighborhood (drugs multiplex)', c('YES', 'NO'), selected = 'NO')
      output[[2]] <-
        shiny::conditionalPanel(
          condition = 'input.drug_jump_neigh == "YES"',
          shiny::tags$hr(),
          shiny::radioButtons('drug_weight_jump_neigh',
                              shiny::h4(shiny::span('Weight inter-layers edges', style = "font-weight: bold")),
                              c('YES', 'NO'), selected = 'NO')
        )
      return(output)
    })
    
    
    #####################
    output$back2P5 <- renderUI({
      
      if (isTruthy(input$gen_multiplex_btn)){
        actionButton('back2P5', label = 'Back', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$back2P5, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "rwr_tab_2")
    })
    
    ###
    output$jump2P7 <- renderUI({
      if (length(shiny::reactiveValuesToList(RWR_output)) != 0){
        actionButton('jump2P7', label = 'Next', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$jump2P7, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "dr_tab_2")
    })
    
    
    ############################################################################################
    #                       DIMENSIONALITY REDUCTION: LOADING DATA (OUTPUT)                    #
    ############################################################################################
    
    output$info_box_8 <- renderText({
      shiny::HTML("<br/> <span style='font-weight:normal;'>
            <p align='justify'>
            You can start the pipeline here by uploading a similarity matrix, an embedded similarity matrix or 
            a lower-dimensional matrix (created via UMAP, PCA, tSNE, etc...),
            or click the <b> button </b> in the top left-hand corner to return back to the 'Data Uploading' page. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <b>Upload Similarity Matrix</b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            In this section, you can upload either i) a similarity matrix (the output of RWR),
            a square matrix with features on both rows and columns, or ii) an embedded similarity matrix 
            (the output of denoising phase or a dimensionality reduction algorithm) with as many columns
            as features and as many rows as selected latent dimensions. </span> <br/>
            </p>
            <hr style='border-top: 1px solid white;'>

            <b>Use Example Dataset</b>: <br/>  <span style='font-weight:normal;'>
             <p align='justify'>
            You can also upload example matrices obtained after conducting RWR-MH on a 5-layer biological network 
            (transcriptomics, proteomics, epigenomics, and CNVs data from BRCA paients)
            and a drug network consisting of pharmacologically active and FDA-approved drugs deposited on DrugBank. </span> </span>
            </p>
            <hr style='border-top: 1px solid white;'>

            <b>Upload Annotation File</b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            This step is not mandatory but useful for enriching subsequent plots.
            The file should be a dataframe with a variable number of columns,
            where one column must contain the nodes IDs, the same ones present as column names in the uploaded files. </span> <br/>
            </p>
            <hr style='border-top: 1px solid white;'>

            <b>Use Example Annotation File</b>: <br/>  <span style='font-weight:normal;'>
            <p align='justify'>
            In the annotation file, genes are classified based on difseveral criteria retrieved from HGNC and COSMIC databases, 
            while drugs are annotated according to the related pathways and the FDA status (DrugBank). </span>
            </p>
            <hr style='border-top: 1px solid white;'>
            
            <span style='font-weight:normal;'>
            <p align='justify'>
            Accepted formats include .csv, .txt, and .RDS.  </span>
            </p> <hr style='border-top: 1px solid white;'>
            ")
      
    })
    
    
    output$simMat_sum_info <- shiny::renderUI({
      input$load_simMat_button
      if ( length(loaded_simMat()) != 0 )
      {
        ncolumns <- ncol(loaded_simMat()$simMatFile)
        nrows <- nrow(loaded_simMat()$simMatFile)
        shinydashboard::box(width = 12,
                            shiny::HTML(paste0(paste("<b>", 'Similarity matrix', "</b>"), " loaded: ",  
                                               paste("<b>", ncolumns, "</b>"), " columns and ", paste("<b>", nrows, "</b>"), " rows <br/>")
                            )
        )
      } else {
        shiny::HTML("<br/>")
      }
    })
    
    
    output$anno_info2 <- shiny::renderUI({
      input$load_anno_button2
      if ( length(loaded_simMat_Anno()) != 0 )
      {
        ncolumns <- ncol(loaded_simMat_Anno()$annotation2)
        nrows <- nrow(loaded_simMat_Anno()$annotation2)
        shinydashboard::box(width = 12,
                            shiny::HTML(paste0(paste("<b>", 'Annotation file', "</b>"), " loaded: ",  paste("<b>", ncolumns, "</b>"), " columns and ", paste("<b>", nrows, "</b>"), " rows <br/>")),
                            shiny::HTML(paste0(paste("<b> Colnames </b>"), "<br/>")),
                            shiny::HTML(paste('&ensp;', colnames(loaded_simMat_Anno()$annotation2), "<br/>"))
        )
      } else {
        shiny::HTML("<br/>")
      }
    })
    
    
    
    output$back2P1_from1.2 <- renderUI({
        actionButton('back2P1_from1.2', label = 'Go to Data', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
    })
    observeEvent(input$back2P1_from1.2, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "mat_sub_1")
    })
    
    ###
    output$jump2P7_from1.2 <- renderUI({
      if (length(loaded_simMat()) != 0 | length(loaded_simMat_Anno()) != 0){
        actionButton('jump2P7_from1.2', label = 'Next', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    observeEvent(input$jump2P7_from1.2, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "dr_tab_2")
    })
    
    
    
    ############################################################################################
    #                            DIMENSIONALITY REDUCTION: DR (OUTPUT)                         #
    ############################################################################################
    
    output$info_box_9 <- renderText({
      HTML("<br/> <span style='font-weight:normal;'>
              <p align='justify'>
              In this section of the app, you can denoise and decrease the dimensionality of the RWR similarity matrix
              in order to facilitate further analyses. </span> <br/>
              </p> <hr style='border-top: 1px solid white;'>

              <b> Denoising via Embedding </b>: <br/> <span style='font-weight:normal;'>
              <p align='justify'>
               This step is not mandatory and could take a lot of time.
               Starting from a (nxn)-matrix, you will get a (dxn)-matrix, where d&lt;&lt;n.
               For further details, please refer to the following paper
               <span style= 'color: blue;'> <u>  <a href='https://www.nature.com/articles/s41598-021-87987-1' target='_blank'> (link)</a></u></span>.<br/>
              </p> </span> <hr style='border-top: 1px solid white;'>

              <b> Dimensionality reduction </b>: <br/> <span style='font-weight:normal;'>
              <p align='justify'>
               Whether or not you decide to denoise the data, it is necessary
               to reduce the dimensionality of the RWR similarity matrix to visualize
               the multi-omics associations in a 2D space. You can choose <b>UMAP</b>, <b>PCA</b> or <b>tSNE</b> method.</span> <br/>
              </p> <hr style='border-top: 1px solid white;'>

              <b> Output </b>: <br/> <span style='font-weight:normal;'>
              <p align='justify'>
               <u>Plot</u>: an interactive plot linked to the two tables below.
               Row selection in the left table will increment the size of the corresponding points on the plot. <br/>
               Selecting a window of points on the plot will populate the rigth table.  <br/>
               <u>Density Plot</u>: an interactive plot useful for studying the density distribution of all/specific nodes. </span> <br/>
              </p> <hr style='border-top: 1px solid white;'>

              <span style='font-weight:normal;'>
              <p align='justify'>
              The denoised and the dimensionality-reducted matrices can be downloaded from the <b>download box</b> located on the left-hand side,
              while the plot can be downloaded by clicking on the <b>download plot</b> button or on the <i class='fas fa-camera'></i> icon
              in the interactive plot. </span>
              </p> <hr style='border-top: 1px solid white;'>
              ")
    })
    
    
    ###
    
    output$dr_opt <- shiny::renderUI({
      
      output$emb_opt_box <- renderUI({
        if (length(shiny::reactiveValuesToList(RWR_output)) != 0 | (input$example_simMat == 'NO' & length(loaded_simMat()) != 0 & input$simMat_type == 'Similarity matrix')){
          output <- tagList()
          output[[1]] <-
            shiny::radioButtons(inputId = 'emb_opt',
                                label = shiny::h4(shiny::span('Do you want to denoise the RWR-mat via embedding?', style = "font-weight: bold")),
                                choices = c('YES', 'NO'), selected = 'NO')
          output[[2]] <-
            shiny::conditionalPanel(
              condition = 'input.emb_opt== "YES"',
              shiny::numericInput(inputId = 'dim_emb',
                                  label = shiny::h4(shiny::span('Select the latent dimensions of the embedding',style = "font-weight: bold")),
                                  min = 1, max = 100, value = 64),
              shiny::numericInput(inputId = 'cores_emb',
                                  label = shiny::h4(shiny::span('Select the number of CPU for the embedding function', style = "font-weight: bold")),
                                  min = 1, max = 100, value = 1)
            )
          return(output)
        } else{
          return(NULL)
        }
      })
      
      output$dr_opt_box <- renderUI({
        if (length(shiny::reactiveValuesToList(RWR_output)) != 0 | (length(loaded_simMat()) != 0 & (
          (input$example_simMat == 'NO' & input$simMat_type != 'Other') | 
          (input$example_simMat == 'YES' & startsWith(input$select_simMat, 'E'))
        )
        )){
          output <- tagList()
          output[[1]] <-
            shiny::selectInput(inputId = 'dr_method',
                               label = shiny::h4(shiny::span('Select a dimensionality reduction method', style = "font-weight: bold")),
                               choices = c('PCA', 'UMAP', 'tSNE'), selected = 'UMAP')
          output[[2]] <-
            shiny::conditionalPanel(
              condition = 'input.dr_method == "PCA"',
              shiny::numericInput(inputId = 'emb_size',
                                  label = shiny::h4(shiny::span('Select the SIZE of the output embedding',style = "font-weight: bold")),
                                  value = 64),
            )
          output[[3]] <-
            shiny::conditionalPanel(
              condition = 'input.dr_method == "tSNE"',
              shiny::numericInput(inputId = 'emb_size_tsne',
                                  label = shiny::h4(shiny::span('Select the SIZE of the output embedding',style = "font-weight: bold")),
                                  min = 1, max = 3, value = 2),
              shiny::numericInput(inputId = 'max_iter',
                                  label = shiny::h4(shiny::span('Select the number of iterations',style = "font-weight: bold")),
                                  value = 20000),
              shiny::numericInput(inputId = 'perplexity',
                                  label = shiny::h4(shiny::span('Select the PERPLEXITY parameter',style = "font-weight: bold")),
                                  min = 1, max = 100, value = 5)
            )
          output[[4]] <-
            shiny::conditionalPanel(
              condition = 'input.dr_method != "PCA"',
              shiny::numericInput(inputId = 'threads',
                                  label = shiny::h4(shiny::span('Select the number of CORES', style = "font-weight: bold")),
                                  value = 1)
            )
          return(output)
        }else{
          return(NULL)
        }
      })
      
      output$componentSelector <- renderUI({
        if (length(loaded_simMat()) != 0 & (input$simMat_type == 'Other' | startsWith(input$select_simMat, 'U')) ){
          output <- tagList()
          output[[1]] <-
            shiny::selectInput(inputId = 'select_c1', 
                               label = 'Select the first component', 
                               choices = rownames(loaded_simMat()$simMatFile),
                               selected = NULL)
          output[[2]] <-
            shiny::selectInput(inputId = 'select_c2', 
                               label = 'Select the second component', 
                               choices = rownames(loaded_simMat()$simMatFile)
            )
          
          return(output)
        } else{
          return(NULL)
        }
      })
      
      
      shinydashboard::box(width = 12, status = 'primary', solidHeader = FALSE, collapsible = FALSE,
                          title = shiny::h2(shiny::span("Dimensionality Reduction", style = "font-weight: bold")),
                          shiny::uiOutput('emb_opt_box'),
                          shiny::hr(),
                          shiny::uiOutput('dr_opt_box'),
                          shiny::uiOutput('componentSelector'),
                          shiny::hr(),
                          shiny::actionButton("dr_btn", "Submit")
      )
    })
    
    ##############################################################################
    #----------------------------------------------------------------------------
    output$update_dr_plot_opt <- shiny::renderUI({
      #shiny::isolate({
        if ( !is.null(anno()) |  !is.null(anno1()) | !is.null(anno2())) {
          anno <- if ( !is.null(anno()) ) anno() else if ( !is.null(anno1()) ) anno1() else anno2()
          output <- tagList()
          output[[1]] <-
            shiny::selectInput(inputId ='dr_plot_color',
                               label = 'Select the color of the points',
                               choices = colnames(anno), selected = NULL)
          output[[2]] <-
            shiny::selectInput(inputId = 'dr_plot_shape',
                               label = 'Select the shape of the points',
                               choices = c('-', colnames(anno)), selected = '-')
          output
        }else{
          return(NULL)
        }
      #})
    })
    
    
    output$update_dr_plot_btn <- shiny::renderUI({
      #shiny::isolate({
      if  ( !is.null(anno()) |  !is.null(anno1()) | !is.null(anno2())) {
        shiny::actionButton(inputId = 'update_dr_plot_btn', label = 'Update plot')
      }else{
        return(NULL)
      }
      #})
    })
    
    
    data_orig_list <- shiny::reactive({
      if (length(shiny::reactiveValuesToList(dr_output)) != 0){
        t_mat <- t( shiny::reactiveValuesToList(dr_output)$dr_mat)
        
        
        
        shared_input <- dplyr::tibble("id" = rownames(t_mat), 'x' = t_mat[,1], 'y' = t_mat[,2])
        
        data_shared <- crosstalk::SharedData$new(shared_input)
        color_base <- "dodgerblue3"
        color_select <- "red"
        var_x <- "x"
        var_y <- "y"
        
        return(list('data_shared' = data_shared,
                    'var_x' = var_x,
                    'var_y' = var_y,
                    'color_base' = color_base,
                    'color_select' = color_select))
      }else{
        return(NULL)
      }
    })
    
    #----------------------------------------------------------------------------
    output$dr_plot_box <- shiny::renderUI({
      req(input$dr_btn)
      shiny::isolate({
        if (length(shiny::reactiveValuesToList(dr_output)) == 0){
          return(NULL)
        } else { 
          if (length(shiny::reactiveValuesToList(RWR_output)) != 0 | (length(loaded_simMat()) != 0 & input$simMat_type == 'Similarity matrix')){
            method <- ifelse(input$dr_method == 'UMAP', 'UMAP of',
                             ifelse(input$dr_method == 'PCA', 'PCA of', 'tSNE of'))
            emb <- ifelse(input$emb_opt == 'YES', 'the Embedded RWR matrix', 'the RWR matrix')
            components <- 1:2
          } else if (length(loaded_simMat()) != 0 & input$simMat_type == 'Embedded Similarity matrix') {
            method <- ifelse(input$dr_method == 'UMAP', 'UMAP of',
                             ifelse(input$dr_method == 'PCA', 'PCA of', 'tSNE of'))
            emb <- 'the loaded Embedded RWR matrix'
            components <- 1:2
          } else if (length(loaded_simMat()) != 0 & input$simMat_type == 'Other'){
            method <- 'Projection'
            emb <- 'of the loaded low-dimensionality Similarity matrix'
            components <- c(input$select_c1, input$select_c2)
            x_lab <- input$select_c1
            y_lab <- input$select_c2
          }
          
          dr_mat <- shiny::reactiveValuesToList(dr_output)$dr_mat
          
          
          output$dr_plot1 <- plotly::renderPlotly({
            
            var_x <- data_orig_list()$var_x
            var_y <- data_orig_list()$var_y
            color_select <- data_orig_list()$color_select
            color_base<- data_orig_list()$color_base
            data_orig <- data_orig_list()$data_shared$origData()
            
            dr_plot1 <-
              plotly::plot_ly(x = data_orig[[var_x]], y = data_orig[[var_y]], source = 'dr_plot1') %>%
              plotly::add_markers(
                color = I(color_base),
                selected = list(
                  marker = list(color = color_select)
                ),
                text = ~paste("x: ", data_orig[[var_x]], "\n",
                              "y: ", data_orig[[var_y]], "\n",
                              "nodeID: ", data_orig$id),
                hoverinfo = "text"
              ) %>%
              plotly::layout(
                dragmode = "select",
                font = list(family = "LMRoman10", color = "black"),
                title = list(text = paste(method, emb), font = list(family = "LMRoman10", color = "black", size = 20)),
                xaxis = list(title = list(text = paste0(input$dr_method, '1'))),
                yaxis = list(title = list(text =  paste0(input$dr_method, '2')))
              ) #%>%
            #plotly::config(
            #    toImageButtonOptions = list(filename = paste('MiDNEshiny', gsub(x = method, replacement = '_', pattern = ' '),
            #                                                 gsub(x = emb, replacement = '_', pattern = ' '), 'plot', sep='_'),
            #                                format = "png", width = 1800, height = 800)
            #)
            
            #}
            
            s <- input$data_orig_table_rows_selected
            print(s)
            if (length(s)) {
              dr_plot1 <- dr_plot1 %>%
                plotly::add_trace(data = data_orig[s, , drop = FALSE],
                                  x = ~x, y = ~y, type = 'scatter',
                                  mode = "markers", color = 'red',
                                  marker = list(size = 20)) %>%
                plotly::layout(showlegend = FALSE)
            }
            
            dr_plot1 <- dr_plot1 %>%
              plotly::ggplotly(source = "dr_plot1") %>%
              plotly::event_register("plotly_brushing") %>%
              plotly::event_register("plotly_relayout") %>%
              plotly::event_register("plotly_restyle")
            #})
            return(dr_plot1)
          })
          
          
          
          #---------------------------- NEW --------------------------------------
          shiny::observeEvent(input$update_dr_plot_btn, {
            
            output$dr_plot_spinner <- renderUI({
              shinycssloaders::withSpinner(visNetwork::visNetworkOutput('dr_plot1'))
            })
            
            output$dr_plot1 <- visNetwork::renderVisNetwork({
              shiny::isolate({
                
                anno <- if ( !is.null(anno()) ) anno() else if ( !is.null(anno1()) ) anno1() else anno2()
                color <- if (input$dr_plot_color != '-') input$dr_plot_color else NA
                shape_col <- if (input$dr_plot_shape != '-') input$dr_plot_shape else NA
                
                var_x <- data_orig_list()$var_x
                var_y <- data_orig_list()$var_y
                #color_select <- data_orig_list()$color_select
                #color_base<- data_orig_list()$color_base
                data_orig <- data_orig_list()$data_shared$origData()
                
                
                nodes <- data_orig$id
                id_col_check <- sapply(colnames(anno), FUN = function(x) sum(nodes %in% anno[[x]]))
                id <- colnames(anno)[which.max(id_col_check)]
                f_anno <- anno[anno[[id]] %in% nodes, ]
                nodes_anno <- f_anno %>% dplyr::as_tibble(.) %>% dplyr::select(id, everything())
                
                data_orig <- data_orig %>% dplyr::left_join(nodes_anno, by = "id")
              
                if (!is.na(color) & !is.na(shape_col)) {
                  dr_plot1 <-
                    plotly::plot_ly(x = data_orig[[var_x]], y = data_orig[[var_y]], source = 'dr_plot1',
                                    type = "scatter", color = data_orig[[color]], 
                                    mode = "markers", symbol = data_orig[[shape_col]], 
                                    marker = list(size = 5)
                    ) %>%
                    plotly::add_markers(
                      colors = "Set2",
                      text = ~paste("x: ", data_orig[[var_x]], "\n",
                                    "y: ", data_orig[[var_y]], "\n",
                                    "nodeID: ", data_orig[['id']], "\n",
                                    "color: ", data_orig[[color]], "\n",
                                    "shape: ", data_orig[[shape_col]]
                      ),
                      hoverinfo = "text"
                    ) %>% 
                    plotly::layout(
                      showlegend = TRUE,
                      dragmode = "select",
                      font = list(family = "LMRoman10", color = "black"),
                      title = list(text = paste(method, emb), font = list(family = "LMRoman10", color = "black", size = 20)),
                      xaxis = list(title = list(text = paste0(input$dr_method, '1'))),
                      yaxis = list(title = list(text =  paste0(input$dr_method, '2')))
                    ) #%>%
                  #plotly::config(
                  #    toImageButtonOptions = list(filename = paste('MiDNEshiny', gsub(x = method, replacement = '_', pattern = ' '),
                  #                                                 gsub(x = emb, replacement = '_', pattern = ' '), 'plot', sep='_'),
                  #                                format = "png", width = 1800, height = 800)
                  #)
                  
                  
                  
                } else if (!is.na(color) & is.na(shape_col)) {
                  dr_plot1 <-
                    plotly::plot_ly(x = data_orig[[var_x]], y = data_orig[[var_y]], source = 'dr_plot1',
                                    type = "scatter", color = data_orig[[color]], 
                                    mode = "markers", 
                                    marker = list(size = 5)
                    ) %>%
                    plotly::add_markers(
                      colors = "Set2",
                      text = ~paste("x: ", data_orig[[var_x]], "\n",
                                    "y: ", data_orig[[var_y]], "\n",
                                    "nodeID: ", data_orig[['id']], "\n",
                                    "color: ", data_orig[[color]]
                      ),
                      hoverinfo = "text"
                    )
                }
                
                dr_plot1 <- dr_plot1 %>%
                  plotly::layout(
                    showlegend = TRUE,
                    legend = list(
                      orientation = 'h',
                      y = -0.3,
                      font = list(size = 10)
                    ),
                    margin = list(b = 100),
                    dragmode = "select",
                    font = list(family = "LMRoman10", color = "black"),
                    title = list(text = paste(method, emb), font = list(family = "LMRoman10", color = "black", size = 20)),
                    xaxis = list(title = list(text = paste0(input$dr_method, '1'))),
                    yaxis = list(title = list(text =  paste0(input$dr_method, '2')))
                  ) #%>%
                #plotly::config(
                #    toImageButtonOptions = list(filename = paste('MiDNEshiny', gsub(x = method, replacement = '_', pattern = ' '),
                #                                                 gsub(x = emb, replacement = '_', pattern = ' '), 'plot', sep='_'),
                #                                format = "png", width = 1800, height = 800)
                #)
                
                
                s <- input$data_orig_table_rows_selected
                print(paste('selected updated: ', s))
                if (length(s)) {
                  dr_plot1 <- dr_plot1 %>%
                    plotly::add_trace(data = data_orig[s, , drop = FALSE],
                                      x = ~x, y = ~y, type = 'scatter',
                                      mode = "markers", color = 'red',
                                      marker = list(size = 20)) #%>%
                    #plotly::layout(showlegend = FALSE)
                }
                
                dr_plot1 <- dr_plot1 %>%
                  plotly::ggplotly(source = "dr_plot1") %>%
                  plotly::event_register("plotly_brushing") %>%
                  plotly::event_register("plotly_relayout") %>%
                  plotly::event_register("plotly_restyle")
                
                
                return(dr_plot1)
              })
            })
          })
          
          
          #-----------------------------------------------------------------------
          
          
          # listen to the brushing event and draw a
          # rect shape that mimics the brush
          observe({
            brush <- plotly::event_data("plotly_brushing", source = 'dr_plot1')
            #print(paste('brush:', brush))
            
            # if the brush is undefined, remove all shapes and exit
            if (is.null(brush)) {
              plotly::plotlyProxy("dr_plot1", session) %>%
                plotly::plotlyProxyInvoke("relayout", list(shapes = NULL))
              return()
            }
            
            # mimc the brush as a rect shape
            brush_rect <- list(
              type = "rect",
              x0 = brush$x[1],
              x1 = brush$x[2],
              y0 = brush$y[1],
              y1 = brush$y[2],
              fillcolor = NA,
              line = list(
                color = "black",
                dash = "dot",
                width = 1
              )
            )
            #print(paste('brush_rect:', brush_rect))
            
            # draw the rect shape and turn off brush coloring
            # imposed by plotly.js
            plotly::plotlyProxy("dr_plot1", session) %>%
              plotly::plotlyProxyInvoke("relayout", list(shapes = list(brush_rect))) %>%
              plotly::plotlyProxyInvoke("restyle", "selectedpoints", list(list()))
          })
          
          # A reactive value that tracks the dimensions of the brush
          brush <- reactiveVal()
          observe({
            evt <- plotly::event_data(event = "plotly_relayout", source = 'dr_plot1')
            #print(paste('evt:', evt))
            
            val <- if (!is.null(evt$shapes)) {
              evt$shapes
            } else if (!is.null(evt[["shapes[0].x0"]])) {
              list(
                x0 = evt[["shapes[0].x0"]],
                x1 = evt[["shapes[0].x1"]],
                y0 = evt[["shapes[0].y0"]],
                y1 = evt[["shapes[0].y1"]]
              )
            }
            brush(val)
          })
          
          # map the brush limits to a data selection
          observe({
            # if brush isn't active, no selection is active
            if (is.null(brush())) {
              data_orig_list()$data_shared$selection(FALSE)
              return()
            }
            
            data_orig <- data_orig_list()$data_shared$origData()
            var_x <- data_orig_list()$var_x
            var_y <- data_orig_list()$var_y
            
            selection <- data.table::between(data_orig[[var_x]], brush()$x0, brush()$x1) &
              data.table::between(data_orig[[var_y]], brush()$y0, brush()$y1)
            #print(paste('selection:', selection))
            
            data_orig_list()$data_shared$selection(selection)
          })
          
          # update the marker colors
          observe({
            dat <- data_orig_list()$data_shared$data(withSelection = TRUE)
            color_select <- data_orig_list()$color_select
            color_base<- data_orig_list()$color_base
            
            color <- dplyr::if_else(dat$selected_, color_select, color_base)
            plotly::plotlyProxy(outputId = "dr_plot1", session) %>%
              plotly::plotlyProxyInvoke("restyle", "marker.color", list(color), 0)
          })
          
          
          # display the selected data
          output$selected_data <- DT::renderDT({
            if (length(shiny::reactiveValuesToList(dr_output)) != 0){
              dat <- data_orig_list()$data_shared$data(withSelection = TRUE)
              filtered_tab <- rstatix::filter(dat, selected_)
              
              DT::datatable(filtered_tab,
                            extensions = 'Buttons',
                            options = list(dom = "Blfrtip",
                                           buttons = list("copy",
                                                          list(extend = 'csv',   filename =  paste("MiDNE", input$cluster_method, 'table',sep = "_")),
                                                          list(extend = 'excel', filename =  paste("MiDNE", input$cluster_method, 'table',sep = "_"))
                                           ),
                                           paging = TRUE,
                                           pageLength = 10,
                                           scrollX = TRUE,
                                           scrollY = TRUE,
                                           autoWidth = TRUE
                            ),
                            selection = 'multiple',
                            filter = 'top',
                            rownames = FALSE
              )
            }
          })
          
          return( shinycssloaders::withSpinner(plotly::plotlyOutput('dr_plot1', height = '600px')) )
          
        }
      })
    })
    
    #----------------------------------------------------------------------------
    
    # display the selected points bigger
    output$data_orig_table <- DT::renderDT({
      if (length(shiny::reactiveValuesToList(dr_output)) != 0) {
        
        dr_table <- data_orig_list()$data_shared$origData()
        
        DT::datatable(dr_table,
                      extensions = 'Buttons',
                      options = list(dom = "Blfrtip",
                                     buttons = list("copy",
                                                    list(extend = 'csv',   filename =  paste("MiDNE", 'dimensionality reducted', 'table',sep = "_")),
                                                    list(extend = 'excel', filename =  paste("MiDNE", 'dimensionality reducted', 'table',sep = "_"))
                                     ),
                                     paging = TRUE,
                                     pageLength = 10,
                                     scrollX = TRUE,
                                     scrollY = TRUE,
                                     autoWidth = TRUE
                      ),
                      selection = 'multiple',
                      filter = 'top',
                      rownames = FALSE
        )
      }else{
        return(NULL)
      }
    })
    
    
    
    #### Create a density plot in the same tabBox oh the dr_plot
    output$density_opt <- shiny::renderUI({
      req(input$dr_btn)
      shiny::isolate({
        observeEvent(list(anno(), anno1(), anno2()), {
          anno_table <- if (!is.null(anno())) anno() else if (!is.null(anno1())) anno1() else if (!is.null(anno2())) anno2() else NULL
          if (!is.null(anno_table)) {
            nodes <- c('gene', 'drug')
            type_col_check <- sapply(colnames(anno_table), FUN = function(x) sum(nodes %in% anno_table[[x]]))
            type_col <- colnames(anno_table)[which.max(type_col_check)]
            
            dr_tab <- shiny::reactiveValuesToList(dr_output)$dr_mat
            anno <- anno_table %>% filter(id %in% colnames(dr_tab))
            
            if ('gene' %in% anno[[type_col]] || 'drug' %in% anno[[type_col]]) {
              options <- c('All', unique(anno[[type_col]]))
              shiny::updateSelectInput(session, "density_on_points", choices = options, selected = 'All')
            }
          }
        })
        tagList(
          shiny::selectInput(inputId = 'density_on_points',
                             label = shiny::h5(shiny::span('Which points do you want to perform the estimation on?', style = "font-weight: bold")),
                             choices = 'All', selected = 'All'),
          shiny::actionButton(inputId = 'density_plot_btn', label = 'Submit')
        )
      })
    })
    
    output$density_plot <- plotly::renderPlotly({
      req(input$density_plot_btn)
      isolate({
        if (length(shiny::reactiveValuesToList(dr_output)) == 0){
          return(NULL)
        } else { 
          
          t_mat <- t( shiny::reactiveValuesToList(dr_output)$dr_mat)
          dr_tab <- dplyr::tibble("id" = rownames(t_mat), 'x' = t_mat[,1], 'y' = t_mat[,2])
          
          if (!is.null(anno()) | !is.null(anno1()) | !is.null(anno2())) {
            pre_anno <- if (!is.null(anno())) anno()  else if (!is.null(anno1())) anno1() else anno2()
            print(pre_anno)
            anno <- pre_anno %>%  filter(id %in% dr_tab$id)
            print(anno)
            nodes <- c('gene', 'drug')
            type_col_check <- sapply(colnames(anno), FUN = function(x) sum(nodes %in% anno[[x]]))
            print(type_col_check)
            idx <- ifelse(max(type_col_check) > 0, which.max(type_col_check), NA)
            print(idx)
            
            if ( !is.na(idx)) {
              colnames(anno)[idx] <- 'type'
              gPlot_data <- dr_tab %>%
                dplyr::left_join(anno, by = "id")
            }else {
              gPlot_data <- dr_tab %>%
                dplyr::mutate(type = "node")
            }
            
          } else {
            gPlot_data <- dr_tab %>%
              dplyr::mutate(type = "node")
          }
          
          if(input$density_on_points == 'All') {
            dr_mat_density <- dr_tab
          } else if (input$density_on_points == 'drug'){
            to_keep <- gPlot_data %>% filter(type == 'drug') %>% select(id)
            dr_mat_density <- dr_tab %>%  filter(id %in% to_keep[['id']])
          }else{
            to_keep <- gPlot_data %>% filter(type == 'gene') %>% select(id)
            dr_mat_density <- dr_tab %>%  filter(id %in% to_keep[['id']])
          }
          
          d_plot <-
            plotly::ggplotly(ggplot2::ggplot(dr_mat_density, aes(x=x, y=y) ) +
                               ggplot2::stat_density_2d(mapping = aes(fill= ..density..),  geom = "raster", contour = FALSE) +
                               ggplot2::geom_point(data = gPlot_data, 
                                                   aes(x=x, y=y, color = type, label = id), alpha=1, size = 0.5) +
                               theme_bw()
            )
          
          print('plot created')
          return(d_plot)
        }
      })
      
    })
    
    
    
    #### ADD other information 
    output$other_table_box <- shiny::renderUI({
    
      if (shiny::isTruthy(bnetwork())){
        
        output$networks <- DT::renderDT({
          if (length(shiny::reactiveValuesToList(dr_output)) != 0) {
           
            DT::datatable(bnetwork(),
                          extensions = 'Buttons',
                          options = list(
                                         search = list(regex = TRUE, caseInsensitive = TRUE),
                                         dom = "Blfrtip",
                                         buttons = list("copy",
                                                        list(extend = 'csv',   filename =  paste("MiDNE", input$select_data_orig, 'table',sep = "_")),
                                                        list(extend = 'excel', filename =  paste("MiDNE", input$select_data_orig, 'table',sep = "_"))
                                         ),
                                         paging = TRUE,
                                         pageLength = 10,
                                         scrollX = TRUE,
                                         scrollY = TRUE,
                                         autoWidth = TRUE
                          ),
                          selection = 'multiple',
                          filter = 'top',
                          rownames = FALSE
            )
          }
        })
        
        return(shinydashboard::box(width = 12, collapsible = TRUE,
                                   title = shiny::h3(shiny::span('Network Table', style = "font-weight: bold")),
                                   #shiny::uiOutput('select_data_orig'),
                                   DT::DTOutput('networks')))
      }else{
        return(NULL)
      }
    })
    
    
    #----------------------------------------------------------------------------

    
    ############## download dr matrices  #############
    
    output$download_dr_box <- shiny::renderUI({
      if (length(shiny::reactiveValuesToList(dr_output)) == 0) {
        return(NULL)
      } else {
        shinydashboard::box(shiny::h3(shiny::span('Download', style = "font-weight: bold")), solidHeader = TRUE, collapsible = FALSE, width = 12,
                            shiny::checkboxGroupInput('download_dr',
                                                      label = shiny::h4(shiny::span('Select one or more low-dimensional RWR matrices', style = 'font-weight: bold')),
                                                      choices = names(shiny::reactiveValuesToList(dr_output))
                            ),
                            shiny::tags$hr(),
                            shiny::downloadButton(
                              outputId = 'download_dr_btn',
                              label = "Download",
                              icon = shiny::icon("download")
                            )
        )
      }
    })
    
    output$download_dr_btn <- shiny::downloadHandler(
      filename = function() {
        paste("MiDNEshiny_dimensionalityReduction_", Sys.Date(), ".zip", sep = "")
      },
      content = function(file){
        
        temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
        dir.create(temp_directory)
        
        to_download <- shiny::reactiveValuesToList(dr_output)[input$download_dr]
        
        for (obj in names(to_download)) {
          file_name <- if (obj == 'dr_mat') glue::glue("{input$dr_method}_{obj}.csv") else glue::glue("{obj}.csv")
          readr::write_csv(as.data.frame(to_download[[obj]]), file.path(temp_directory, file_name))
        }
        
        zip::zip(
          zipfile = file,
          files = dir(temp_directory),
          root = temp_directory
        )
      }
    )
    
    ############## download plot cluster in 2D #############
    
    output$download_dr_plot <- shiny::downloadHandler(
      filename = function() { paste0('MiDNEshiny_', input$dr_method, '_plot.png', sep='') },
      content = function(file) {
        ggplot2::ggsave(filename = file, plot = dr_plot(), device = "png", width = 10, height = 5, limitsize = FALSE)
      }
    )
    
    ###
    output$back2P6 <- renderUI({
      if (length(shiny::reactiveValuesToList(RWR_output)) != 0){
        actionButton('back2P6', label = 'Back', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$back2P6, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "rwr_tab_3")
    })
    
    ###
    output$back2P1.2 <- renderUI({
      actionButton('back2P1.2', label = 'Back', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
    })
    
    observeEvent(input$back2P1.2, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "dr_tab_1")
    })
    
    ###
    output$jump2P8 <- renderUI({
      if (length(shiny::reactiveValuesToList(dr_output)) != 0){
        actionButton('jump2P8', label = 'Next', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$jump2P8, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "cl_tab")
    })
    
    
    
    ############################################################################################
    #                                     CLUSTERING (OUTPUT)                                  #
    ############################################################################################
    
    output$info_box_10 <- renderText({
      HTML("<br/> <span style='font-weight:normal;'>
            <p align='justify'>
            Here, you can conduct a clustering analysis on the dimensionality-reducted RWR
            similarity matrix. You can choose different clustering algorithms. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <b> Clustering algorithms </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
             You can generate clusters starting from an embedded matrix, if you 
             decided to denoise the RWR similarity matrix, or from the dimensionality-reducted matrix
             via UMAP/PCA/tSNE. <br/>
             &ensp; If you select <span style='color:steelblue4;'><b>kmeans</b></span>, you just have  to 
             set the number of desired clusters (<b>k</b>).<br/>
             &ensp; If you select <span style='color:steelblue4;'><b>dbscan</b></span>, you just have to 
             set the <b>&epsilon;</b> parameter.<br/>
             &ensp; If you select <span style='color:steelblue4;'><b>hclust</b></span>, you just have to 
             follow these steps: <br/>
             &ensp;&ensp;&ensp;i) <u>press</u> the <b>Submit</b> button; <br/>
             &ensp;&ensp;&ensp;ii) <u>look at</u> the <b>dendrogram</b>; <br/>
             &ensp;&ensp;&ensp;iii) <u>move</u> the slider to select a threshold <b> h </b> to cut the dendrogram; <br/>
             &ensp;&ensp;&ensp;iv) <u>click</u> the <b> 'Submit' </b> button again. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>
            
            <b> Manual method </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
             You  can interactively select a region of points on the plot, and subsequently populate
             the table in the bottom right. If you press the <b>Update</b> button in the bottom left of the table,
             a single cluster of points will be generated for analysis in the subsequent stages. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>
            
            ")
    })
    
    
    ############## output cluster box #############
    
    output$cl_opt <- shiny::renderUI({
      
      shinydashboard::box(width = 12, status = 'primary', solidHeader = FALSE,
                          title = shiny::h2(shiny::span("Clustering Analysis", style = "font-weight: bold")),
                          shiny::conditionalPanel(
                            condition = 'input.emb_opt == "YES"',
                            shiny::radioButtons(inputId = 'cl_on_emb',
                                                label = shiny::h4(shiny::span('Do you want to cluster the embedded matrix?',style = "font-weight: bold")),
                                                choices = c('YES', 'NO'), selected = 'NO'),
                          ),
                          shiny::selectInput(inputId = 'cluster_method',
                                             label = shiny::h4(shiny::span('Select a clustering method',style = "font-weight: bold")),
                                             choices = list(
                                               'algorithms' = c("kmeans", "dbscan" , "hclust"),
                                               'customize' = list("manual")  #anno_list[!sapply(anno_list, is.null)]
                                             ),
                                             selected = "kmeans"
                          ),
                         
                          shiny::conditionalPanel(
                            condition = "input.cluster_method == 'kmeans'",
                            shiny::numericInput(inputId = "kmeans_par",
                                                label = shiny::h4(shiny::span("Number of clusters for kmeans",style = "font-weight: bold")),
                                                min = 2, max= 500, value = 50
                            )
                          ),
                          shiny::conditionalPanel(
                            condition = "input.cluster_method == 'dbscan'",
                            shiny::sliderInput(inputId = "dbscan_par",
                                               label = shiny::h4(shiny::span("epsilon",style = "font-weight: bold")),
                                               min = 0, max = 1,value = 0.1, step = 0.01
                            )
                          ),
                          shiny::conditionalPanel(
                            condition = "input.cluster_method == 'hclust'",
                            shiny::uiOutput('h'),
                            
                          ),
                          
                          shiny::hr(),
                          shiny::actionButton("cl_btn", "Submit")
      )
    })
    
    ############## plot cluster in 2D #############
  
    ##### SHARED OBJ
    cl_data_orig_list <- shiny::reactive({
      isolate({
        
      if (length(shiny::reactiveValuesToList(clusters)$manual) != 0 ){
    
          clusters <- shiny::reactiveValuesToList(clusters)
          clu_anno <- clusters$manual$cluster
          
          data_shared <- crosstalk::SharedData$new(clu_anno)
          color_base <- "dodgerblue3"
          color_select <- "red"
          var_x <- "x"
          var_y <- "y"
          
          return(list('data_shared' = data_shared,
                      'var_x' = var_x,
                      'var_y' = var_y,
                      'color_base' = color_base,
                      'color_select' = color_select))
        
      }else{
        return(NULL)
      }
      
    })
    })
    
    output$cl_plot_box <- shiny::renderUI({
      req(input$cl_btn)
      shiny::isolate({
        if (length(shiny::reactiveValuesToList(clusters)) != 0){
          param <- if (isTruthy(input[[paste0(input$cluster_method, '_par')]])) input[[paste0(input$cluster_method, '_par')]] else NULL
          if (is.null(shiny::reactiveValuesToList(clusters)[[paste0(input$cluster_method, '_', param)]]$cluster)){
            return(NULL)
          }else{
            clusters <- shiny::reactiveValuesToList(clusters)
            var_x <- 'x'
            var_y <- 'y'
            clu_anno <- clusters[[paste0(input$cluster_method, '_', param)]]$cluster %>% dplyr::arrange(clust)
            
            output$cl_plot <- plotly::renderPlotly({
              #isolate({
              cl_plot <- clu_anno %>% 
                plotly::plot_ly(x = ~x, y = ~y, source = 'cl_plot',
                                type = "scatter", color = ~clust,
                                mode = "markers", 
                                marker = list(size = 5)
                ) %>%
                plotly::add_markers(
                  colors = "Set2",
                  text = ~paste("x: ", clu_anno$x, "\n",
                                "y: ", clu_anno$y, "\n",
                                "nodeID: ",  clu_anno$id, "\n",
                                "color: ",  clu_anno$clust
                  ),
                  hoverinfo = "text"
                ) %>%
                plotly::layout(
                  dragmode = "select",
                  font = list(family = "LMRoman10", color = "black"),
                  title = list(text = title(), font = list(family = "LMRoman10", color = "black", size = 20)),
                  xaxis = list(title = list(text = paste0(input$dr_method, '1'))),
                  yaxis = list(title = list(text =  paste0(input$dr_method, '2')))
                ) 
              
              cl_plot <- cl_plot %>%
                plotly::ggplotly(source = "cl_plot") %>%
                plotly::event_register("plotly_brushing") %>%
                plotly::event_register("plotly_relayout") %>%
                plotly::event_register("plotly_restyle")
              
              to_update(TRUE)
              
              return(cl_plot)
              #})
            })
            
            return(shinycssloaders::withSpinner(plotly::plotlyOutput('cl_plot', height = '600px')))
          }
        }
      })
    })
    
  
    manual_cluster <- shiny::reactiveVal()
    output$manual_cl_plot_box <- shiny::renderUI({
      req(input$cl_btn)
      if (showPlot()){
      isolate({
        
      if (length(shiny::reactiveValuesToList(clusters)$manual) != 0 ) {
        
        var_x <- 'x'
        var_y <- 'y'
        
        clu_anno <- cl_data_orig_list()$data_shared$origData()
        
        output$manual_cl_plot <- plotly::renderPlotly({
          #if (showPlot()){
          isolate({
            manual_cl_plot <- clu_anno %>% 
              plotly::plot_ly(x = ~x, y = ~y, source = 'manual_cl_plot',
                              type = "scatter", color = 'dodgerblue3',
                              mode = "markers", 
                              marker = list(size = 5)
              ) %>%
              plotly::add_markers(
                #colors = "Set2",
                text = ~paste("x: ", clu_anno$x, "\n",
                              "y: ", clu_anno$y, "\n",
                              "nodeID: ",  clu_anno$id
                ),
                hoverinfo = "text"
              ) %>%
              plotly::layout(
                dragmode = "select",
                font = list(family = "LMRoman10", color = "black"),
                title = list(text = title(), font = list(family = "LMRoman10", color = "black", size = 20)),
                xaxis = list(title = list(text = paste0(input$dr_method, '1'))),
                yaxis = list(title = list(text =  paste0(input$dr_method, '2')))
              ) #%>%
            #plotly::config(
            #    toImageButtonOptions = list(filename = paste('MiDNEshiny', gsub(x = method, replacement = '_', pattern = ' '),
            #                                                 gsub(x = emb, replacement = '_', pattern = ' '), 'plot', sep='_'),
            #                                format = "png", width = 1800, height = 800)
            #)
            
            #}
            
            s <- input$cl_data_orig_table_rows_selected
            if (length(s)) {
              manual_cl_plot <- manual_cl_plot %>%
                plotly::add_trace(data =  clu_anno[s, , drop = FALSE],
                                  x = ~x, y = ~y, type = 'scatter',
                                  mode = "markers", color = 'red',
                                  marker = list(size = 20)) %>%
                plotly::layout(showlegend = FALSE)
            }
            
            manual_cl_plot <- manual_cl_plot %>%
              plotly::ggplotly(source = "manual_cl_plot") %>%
              plotly::event_register("plotly_brushing") %>%
              plotly::event_register("plotly_relayout") %>%
              plotly::event_register("plotly_restyle")
            
            to_update(TRUE)
            
            return(manual_cl_plot)
          })
          #}
        })
        
      
        ######################## BRUSHING ################################
        
        observe({
          brush <- plotly::event_data("plotly_brushing", source = 'manual_cl_plot')
          #print(paste('brush:', brush))
          
          # if the brush is undefined, remove all shapes and exit
          if (is.null(brush)) {
            plotly::plotlyProxy("manual_cl_plot", session) %>%
              plotly::plotlyProxyInvoke("relayout", list(shapes = NULL))
            return()
          }
          
          # mimc the brush as a rect shape
          brush_rect <- list(
            type = "rect",
            x0 = brush$x[1],
            x1 = brush$x[2],
            y0 = brush$y[1],
            y1 = brush$y[2],
            fillcolor = NA,
            line = list(
              color = "black",
              dash = "dot",
              width = 1
            )
          )
          #print(paste('brush_rect:', brush_rect))
          
          # draw the rect shape and turn off brush coloring
          # imposed by plotly.js
          plotly::plotlyProxy("manual_cl_plot", session) %>%
            plotly::plotlyProxyInvoke("relayout", list(shapes = list(brush_rect))) %>%
            plotly::plotlyProxyInvoke("restyle", "selectedpoints", list(list()))
        })
        
        # A reactive value that tracks the dimensions of the brush
        brush <- reactiveVal()
        observe({
          evt <- plotly::event_data(event = "plotly_relayout", source = 'manual_cl_plot')
          #print(paste('evt:', evt))
          
          val <- if (!is.null(evt$shapes)) {
            evt$shapes
          } else if (!is.null(evt[["shapes[0].x0"]])) {
            list(
              x0 = evt[["shapes[0].x0"]],
              x1 = evt[["shapes[0].x1"]],
              y0 = evt[["shapes[0].y0"]],
              y1 = evt[["shapes[0].y1"]]
            )
          }
          brush(val)
        })
        
        
        # map the brush limits to a data selection
        observe({
          # if brush isn't active, no selection is active
          if (is.null(brush())) {
            cl_data_orig_list()$data_shared$selection(FALSE)
            return()
          }
          
          data_orig <- cl_data_orig_list()$data_shared$origData()
          var_x <- cl_data_orig_list()$var_x
          var_y <- cl_data_orig_list()$var_y
          
          selection <- data.table::between(data_orig[[var_x]], brush()$x0, brush()$x1) &
            data.table::between(data_orig[[var_y]], brush()$y0, brush()$y1)
          #print(paste('selection:', selection))
          
          cl_data_orig_list()$data_shared$selection(selection)
        })
        
        # update the marker colors
        observe({
          dat <- cl_data_orig_list()$data_shared$data(withSelection = TRUE)
          color_select <- cl_data_orig_list()$color_select
          color_base<- cl_data_orig_list()$color_base
          
          color <- dplyr::if_else(dat$selected_, color_select, color_base)
          plotly::plotlyProxy(outputId = "manual_cl_plot", session) %>%
            plotly::plotlyProxyInvoke("restyle", "marker.color", list(color), 0)
        })
        
        # display the selected data
        output$cl_selected_data <- DT::renderDT({
          
            dat <- cl_data_orig_list()$data_shared$data(withSelection = TRUE)
            
            filtered_tab <- rstatix::filter(dat, selected_)
            manual_cluster(filtered_tab)
            
            print(filtered_tab)
            
            DT::datatable(filtered_tab,
                          extensions = 'Buttons',
                          options = list(dom = "Blfrtip",
                                         buttons = list("copy",
                                                        list(extend = 'csv',   filename =  paste("MiDNE", input$cluster_method, 'table',sep = "_")),
                                                        list(extend = 'excel', filename =  paste("MiDNE", input$cluster_method, 'table',sep = "_"))
                                         ),
                                         paging = TRUE,
                                         pageLength = 10,
                                         scrollX = TRUE,
                                         scrollY = TRUE,
                                         autoWidth = TRUE
                          ),
                          selection = 'multiple',
                          filter = 'top',
                          rownames = FALSE
            )
            
        })
        
        return(
          
           shinycssloaders::withSpinner(
             plotly::plotlyOutput('manual_cl_plot', height = '600px')
             )
           
          )
      } else{  
        shiny::tagList()
      }
        
      })
      }
    })
    
    
    observeEvent(input$cl_btn, {
      if (input$cluster_method != 'manual') {
        showPlot(FALSE) 
      } else {
        showPlot(TRUE) 
      }
    })
    
   
          
    ### dendrogram
    output$plotHclust_box <- shiny::renderUI({
      
      if (!is.null(shiny::reactiveValuesToList(clusters)$hclust_output)){
        tree <- shiny::reactiveValuesToList(clusters)$hclust_output 
        
        output$plotHclust <- shiny::renderPlot({
          req(input$cl_btn)
          
          shiny::isolate({
            #if (!isTruthy(input$hclust_par)) {
            if (is.null(shiny::reactiveValuesToList(clusters)[[paste0('hclust_', input$hclust_par)]]$num_clust)) {
              tree %>% plot(labels=FALSE, sub='', xlab='')
              
              shinyalert::shinyalert(
                title = "Next step",
                text = paste0("i) Look at the <b> dendrogram </b> that was just displayed; <br/>
                                          ii) <b> move the slider </b> to select a threshold <b> 'h' </b> to cut the dendrogram; <br/>
                                          iii) click the <b> 'Submit' </b> button again."
                ),
                closeOnEsc = TRUE,
                closeOnClickOutside = TRUE,
                html = TRUE,
                type = "info",
                showConfirmButton = TRUE,
                confirmButtonText = "OK",
                confirmButtonCol = "#004192",
                showCancelButton = FALSE,
                imageUrl = "",
                animation = TRUE
              )
              
            } else{
              h_cluster <- shiny::reactiveValuesToList(clusters)[[paste0('hclust_', input$hclust_par)]]
              num_clust <- h_cluster$num_clust
              
              
              
              the_bars <- h_cluster$cluster %>% 
                dplyr::select(-c(x,y)) %>%
                tibble::column_to_rownames(var = "id") %>%
                dplyr::mutate_if(is.character, as.factor) %>%
                dplyr::mutate_if(is.factor, as.numeric)
              
              
              par(mar = c(10, 3, 3, 4) + 0.1,
                  xpd = NA) # allow content to go into outer margin
              
              tree_plot <- tree %>% stats::as.dendrogram() %>% 
                #dendextend::set("labels_col", k = num_clust) %>%
                dendextend::set("labels", rep('', length(tree$labels))) %>%
                dendextend::set("branches_k_color", k = num_clust) %>%
                plot(horiz=FALSE, axes=TRUE)
              
              final_tree_plot <- tree_plot +
                abline(h = input$hclust_par, lty = 'dotdash') +
                dendextend::colored_bars(colors = the_bars, dend = stats::as.dendrogram(tree),
                                         rowLabels = colnames(the_bars), add = TRUE,
                                         y_shift = -1
                )
              return(final_tree_plot)
            }
          })
        })
        shinydashboard::box(width = 12,
                            shinycssloaders::withSpinner(shiny::plotOutput(outputId = "plotHclust"))
                            )
      }else {
        return(NULL)
      }
    })
    
    ### table
    
    output$cluster_table_box <- shiny::renderUI({
      cl_method <- input$cluster_method
      param <- if (isTruthy(input[[paste0(cl_method, '_par')]])) input[[paste0(cl_method, '_par')]] else NULL
      if (!is.null(shiny::reactiveValuesToList(clusters)[[paste0(cl_method, '_', param)]]$cluster) ){
        
        output$cl_table <- DT::renderDT(server = FALSE, {
          input$cl_btn
          shiny::isolate({
            #print(shiny::reactiveValuesToList(clusters)[[paste0(input$cluster_method, '_', param)]]$cluster)
           
            cl_table <- shiny::reactiveValuesToList(clusters)[[paste0(input$cluster_method, '_', param)]]$cluster #%>% dplyr::arrange(clust)
            
            DT::datatable(cl_table,
                          extensions = 'Buttons',
                          options = list(dom = "Blfrtip",
                                         buttons = list("copy",
                                                        list(extend = 'csv',   filename =  paste("MiDNE", cl_method, param, 'table',sep = "_")),
                                                        list(extend = 'excel', filename =  paste("MiDNE", cl_method, param, 'table',sep = "_"))
                                         ),
                                         paging = TRUE,
                                         pageLength = 10,
                                         scrollX = TRUE,
                                         scrollY = TRUE,
                                         autoWidth = TRUE
                          ),
                          selection = 'multiple',
                          filter = 'top',
                          rownames = FALSE
            )
          })
        
        })
        shinydashboard::box(width = 12,
                            title = shiny::h3(shiny::span('Clustering table', style = "font-weight: bold")),
                            DT::DTOutput('cl_table'))
      } else{
        return(NULL)
      }
    })
    
    
    output$manual_cluster_table_box <- shiny::renderUI({
      if (!is.null(shiny::reactiveValuesToList(clusters)$manual)){
        
        output$cl_data_orig_table <- DT::renderDT(server = FALSE, {
          input$cl_btn
          shiny::isolate({
            cl_table <- cl_data_orig_list()$data_shared$origData()
            DT::datatable(cl_table,
                          extensions = 'Buttons',
                          options = list(dom = "Blfrtip",
                                         buttons = list("copy",
                                                        list(extend = 'csv',   filename =  paste("MiDNE lower-Dimensional table",sep = "_")),
                                                        list(extend = 'excel', filename =  paste("MiDNE lower-Dimensional table",sep = "_"))
                                         ),
                                         paging = TRUE,
                                         pageLength = 10,
                                         scrollX = TRUE,
                                         scrollY = TRUE,
                                         autoWidth = TRUE
                          ),
                          selection = 'multiple',
                          filter = 'top',
                          rownames = FALSE
            )
          })
        
        })
        shinydashboard::box(width = 12,
                            title = shiny::h3(shiny::span('Manual clustering table', style = "font-weight: bold")),
                            DT::DTOutput('cl_data_orig_table'))
      } else{
        return(NULL)
      }
    })
    
    
    ############## info created cluster tables #############
    
    output$info_created_cl <- renderUI({
      if (!is.null(shiny::reactiveValuesToList(clusters)) ) {
        cl_list <- shiny::reactiveValuesToList(clusters)
        
        cl_tables_names <- c()
        for (name in names(cl_list)){
          if (name != 'dr_mat' & !is.null(cl_list[[name]])) {
            cl_tables_names <- c(cl_tables_names, name)
          }
        }
        
        shinydashboard::box(width = 12, title = shiny::h3(shiny::span("Download or Remove box", style = "font-weight: bold")),
                            #shiny::HTML(" "),
                            shiny::checkboxGroupInput(inputId = 'dw_rm_cl_opt', label = 'Download or remove clustering tables',
                                                      choices = cl_tables_names, selected = NULL),
                            shiny::tags$hr(),
                            shiny::actionButton(inputId = 'rm_cl_btn', label = 'Remove'),
                            shiny::downloadButton(outputId = 'dw_cl_btn', 
                                                  label = 'Download', 
                                                  icon = shiny::icon('download')),
                            #shiny::HTML(paste0(
                            #  paste("<b>", cl_tables, "</b> "), "created <br/>")
                            #)
        )
      }
    })
    
    output$dw_cl_btn <- shiny::downloadHandler(
      filename = function() {
        paste("MiDNEshiny_clusterTable_", Sys.Date(), ".zip", sep = "")
      },
      content = function(file){
        
        temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
        dir.create(temp_directory)
        
        to_download <- shiny::reactiveValuesToList(clusters)[input$dw_rm_cl_opt]
        
        for (obj in names(to_download)) {
          file_name <- glue::glue("{obj}.csv")
          readr::write_csv(as.data.frame(to_download[[obj]]), file.path(temp_directory, file_name))
        }
        
        zip::zip(
          zipfile = file,
          files = dir(temp_directory),
          root = temp_directory
        )
      }
    )
    
    
    
    ############## download plot cluster in 2D #############
    
    output$download_cl_plot <- shiny::downloadHandler(
      filename = function() { paste('MiDNEshiny', #input$dr_method, 
                                    input$cluster_method, 'plot.png', sep='_') },
      content = function(file) {
        ggplot2::ggsave(filename = file, plot = cl_plot(), device = "png", width = 10, height = 5, limitsize = FALSE)
      }
    )
    
    ###
    output$back2P7 <- renderUI({
      if (length(shiny::reactiveValuesToList(dr_output)) != 0){
        actionButton('back2P7', label = 'Back', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$back2P7, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "dr_tab_2")
    })
    
    ###
    output$jump2P9 <- renderUI({
      if (length(shiny::reactiveValuesToList(clusters)) != 0){
        actionButton('jump2P9', label = 'Next', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$jump2P9, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "clupath_tab")
    })
    
    
    ############################################################################################
    #                        ENRICHMENT ANALYSIS: CLUSTER2PATHWAY (OUTPUT)                     #
    ############################################################################################
    
    output$info_box_11 <- renderText({
      HTML("<br/> <span style='font-weight:normal;'>
              <p align='justify'>
              Here, you can conduct the Pathway Enrichment Analysis (PEA) for a particular list of clusters obtained in the previous step. <br/>
              Then, you can perform the <b> Cluster2Pathway Analysis </b>, by filtering the PEA table based on one or more databases. 
              The output will be a table and a plot in which each cluster is associated with the most statistically significant pathway/function.
              </span> <br/>
              </p> <hr style='border-top: 1px solid white;'>
              ")
    })
  
    
    output$clu_path_box <- shiny::renderUI({
      
      if  (!is.null(shiny::reactiveValuesToList(gProfiler_res))) {
        shinydashboard::box(title = shiny::h2(shiny::span("Cluster to Pathway Analysis", style = "font-weight: bold")), 
                            solidHeader = FALSE, status = 'primary', width = 12, 
                            shiny::uiOutput('pea_table_selector'),
                            shiny::checkboxGroupInput('anno_source',
                                                      'Select one or more source terms as filters',
                                                      choices = list("CORUM", "GO:CC", "GO:MF", "HPA", "WP",
                                                                     "GO:BP", "KEGG", "REAC", "TF", "HP", "MIRNA")
                            ),
                            shiny::tags$hr(),
                            shiny::actionButton("showOne", label = "ShowOne")
                            
        )
      }else{
        return(NULL)
      }
    })
    
    
    output$cl_table_selector <- shiny::renderUI({
      if (!is.null(shiny::reactiveValuesToList(clusters))){
        
        clusters <- shiny::reactiveValuesToList(clusters)
        cl_choices <- c()
        for (name in names(clusters)){
          if (!is.null(clusters[[name]])) {
            cl_choices <- c(cl_choices, name)
          }
        }
        
        shiny::selectInput(inputId = 'select_cluList', label = 'Select a cluster table',
                           choices = cl_choices[!(cl_choices %in% c('dr_mat', 'hclust_output', 'manual'))], selected = NULL)
      }
    })
    
    output$pea_table_selector <- shiny::renderUI({
      if (!is.null(shiny::reactiveValuesToList(gProfiler_res))){
        pea_list <- shiny::reactiveValuesToList(gProfiler_res)
        pea_choices <- c()
        for (name in names(pea_list)){
          if (!is.null(pea_list[[name]])) {
            pea_choices <- c(pea_choices, name)
          }
        }
        shiny::selectInput(inputId = 'select_peaTable', label = 'Select a pea table',
                           choices = pea_choices, selected = NULL)
      }
    })
    
    output$TOPpathwayAnnotation <- plotly::renderPlotly({
      input$showOne
      isolate({
        #if (length(shiny::reactiveValuesToList(clusters)) == 0){
        #  return(NULL)
        #}else{
        
        dr_mat <- shiny::reactiveValuesToList(clusters)$dr_mat[1:2,]
        
        if (input$select_peaTable != 'manual_cluster'){
          top_anno <- topPath()$top_anno
        }else{
          pre_top_anno <- topPath()$top_anno
          to_add <- colnames(dr_mat)[! colnames(dr_mat) %in% pre_top_anno$id]
          
          t_mat <- t(dr_mat)
          dr_table <-dplyr::tibble("id" = rownames(t_mat), 'x' = t_mat[,1], 'y' = t_mat[,2])
          
          top_anno <- dplyr::full_join(pre_top_anno, dr_table, by = 'id')
          top_anno[ is.na(top_anno$clust), 'clust'] <- 2
        }
        
        MoNETA::plot_2D_matrix(coord = dr_mat,
                               nodes_anno = top_anno,
                               id_name = 'id', id_anno_color = 'term_name', id_anno_shape = 'clust',
                               interactive = TRUE, wo_legend = TRUE,
                               title = paste('The most significant pathway per cluster (', input$select_peaTable, ')'))
        #}
      })
    })
    
    output$TOPtable <-  DT::renderDT(server = FALSE,{
      input$showOne
      req(topPath())
      DT::datatable(topPath()$top_enrich_path,
                    extensions = 'Buttons',
                    options = list(paging = TRUE,    ## paginate the output
                                   pageLength = 15,  ## number of rows to output for each page
                                   scrollX = TRUE,   ## enable scrolling on X axis
                                   scrollY = TRUE,   ## enable scrolling on Y axis
                                   autoWidth = TRUE,
                                   dom = "Blfrtip",
                                   buttons = list("copy", list(
                                     extend = "collection",
                                     buttons = c("csv", "excel"),
                                     text = "Download"
                                   )
                                   )
                    ),
                    selection = 'multiple', ## enable selection of a single row
                    filter = 'top',        ## include column filters at the bottom
                    rownames = FALSE
      )
    })
    
    
    output$info_created_pea <- renderUI({
      if (!is.null(shiny::reactiveValuesToList(gProfiler_res)) ) {
        pea_list <- shiny::reactiveValuesToList(gProfiler_res)
        
        pea_tables_names <- c()
        for (name in names(pea_list)){
          if (!is.null(pea_list[[name]])) {
            pea_tables_names <- c(pea_tables_names, name)
          }
        }
        
        shinydashboard::box(width = 12, title = shiny::h3(shiny::span("Download or Remove box", style = "font-weight: bold")),
                            #shiny::HTML(" "),
                            shiny::checkboxGroupInput(inputId = 'dw_rm_pea_opt', label = 'Download or remove pea tables',
                                                      choices = pea_tables_names, selected = NULL),
                            shiny::tags$hr(),
                            shiny::actionButton(inputId = 'rm_pea_btn', label = 'Remove'),
                            shiny::downloadButton(outputId = 'dw_pea_btn', 
                                                  label = 'Download', 
                                                  icon = shiny::icon('download'))
        )
      }
    })
    
    output$dw_pea_btn <- shiny::downloadHandler(
      filename = function() {
        paste("MiDNEshiny_peaTable_", Sys.Date(), ".zip", sep = "")
      },
      content = function(file){
        
        temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
        dir.create(temp_directory)
        
        to_download <- shiny::reactiveValuesToList(gProfiler_res)[input$dw_rm_pea_opt]
        
        for (obj in names(to_download)) {
          file_name <- glue::glue("{obj}.csv")
          readr::write_csv(as.data.frame(to_download[[obj]]), file.path(temp_directory, file_name))
        }
        
        zip::zip(
          zipfile = file,
          files = dir(temp_directory),
          root = temp_directory
        )
      }
    )
    
    
    ###
    output$back2P8 <- renderUI({
      if (length(shiny::reactiveValuesToList(clusters)) != 0){
        actionButton('back2P8', label = 'Back', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$back2P8, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "cl_tab")
    })
    
    ###
    output$jump2P10 <- renderUI({
      if (length(shiny::reactiveValuesToList(gProfiler_res)) != 0){
        actionButton('jump2P10', label = 'Next', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$jump2P10, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "pathclu_tab")
    })
    
    
    
    
    ############################################################################################
    #                        ENRICHMENT ANALYSIS: PATHWAY2CLUSTERS (OUTPUT)                     #
    ############################################################################################
    
    output$info_box_12 <- renderText({
      HTML("<br/> <span style='font-weight:normal;'>
              <p align='justify'>
              After conducting the pathway enrichment analysis,
              here you can verify the enrichment specificity of a particular pathway across all defined clusters. 
              Select a single database to filter the enrichment table, and then select a pathway/functionality to test its specificity. </span> <br/>
              </p> <hr style='border-top: 1px solid white;'>
              ")
    })
    
    output$pea_table_selector_2 <- shiny::renderUI({
      if (!is.null(shiny::reactiveValuesToList(gProfiler_res))){
        pea_list <- shiny::reactiveValuesToList(gProfiler_res)
        pea_choices <- c()
        for (name in names(pea_list)){
          if (!is.null(pea_list[[name]])) {
            pea_choices <- c(pea_choices, name)
          }
        }
        shiny::selectInput(inputId = 'select_peaTable2', label = 'Select a pea table',
                           choices = pea_choices, selected = NULL)
      }
    })
    
    output$pathway_selector <- shiny::renderUI({
      input$path2clust
      req(path_annotation())
      if (is.null(path_annotation()))
        return(NULL)
      isolate({
        selectInput('pathway_selector',
                    'Select one pathway to highlight in the plot',
                    choices = sort(unique(path_annotation()$path_cluster$term_name))
        )
      })
    })
    
    
    output$pathwayPlot <- plotly::renderPlotly({
      input$show
      if (length(shiny::reactiveValuesToList(clusters)) == 0) {
        return(NULL)}
      isolate({
        enriched_clusters <- path_annotation()$path_cluster %>% filter(term_name == input$pathway_selector) %>% dplyr::pull(., query)
        annotation_table <- path_annotation()$gene_cluster %>% mutate(selected_pathway = ifelse(!(clust %in% enriched_clusters), 'FALSE', 'TRUE'))
        
        clusters <- shiny::reactiveValuesToList(clusters)
        MoNETA::plot_2D_matrix(coord = clusters$dr_mat[1:2,], nodes_anno = annotation_table,
                               id_name = 'id', id_anno_color = 'selected_pathway',
                               interactive = TRUE, wo_legend = TRUE,
                               title = paste('Enriched clusters for "', input$pathway_selector, '" (', input$select_peaTable2, ')'))
      })
    })
    
    
    
    output$Alltable <-  DT::renderDT(server = FALSE, {
      req(input$show)
      if (is.null(shiny::reactiveValuesToList(gProfiler_res) ))
        return(NULL)
      isolate({
        gprof_res <- shiny::reactiveValuesToList(gProfiler_res)[[input$select_peaTable2]]
        onedb <- filter(gprof_res, source == input$source)
        DT::datatable(onedb,
                      extensions = 'Buttons',
                      options = list(paging = TRUE,
                                     pageLength = 10,
                                     scrollX = TRUE,
                                     scrollY = TRUE,
                                     autoWidth = TRUE,
                                     dom = "Blfrtip",
                                     buttons = list("copy", list(
                                       extend = "collection",
                                       buttons = c("csv", "excel"),
                                       text = "Download"
                                     )
                                     )
                      ),
                      selection = 'multiple',
                      filter = 'top',
                      rownames = FALSE
        )
      })
    })
    
    
    ###
    output$back2P9 <- renderUI({
      if (length(shiny::reactiveValuesToList(gProfiler_res)) != 0){
        actionButton('back2P9', label = 'Back', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$back2P9, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "clupath_tab")
    })
    
    ###
    output$jump2P11 <- renderUI({
      if (length(shiny::reactiveValuesToList(clusters)) != 0 & (!is.null(anno()) |  !is.null(anno1()) | !is.null(anno2()) ) ){
        anno <- if (!is.null(anno())) anno()  else if (!is.null(anno1())) anno1() else anno2()
        id_col_check <- sapply(colnames(anno), FUN = function(x) sum('drug' %in% anno[[x]]))
        type_col <- colnames(anno)[which.max(id_col_check)]
        
        drugs <- anno[anno[[type_col]] == 'drug',][['id']] 
        dr_mat <- shiny::reactiveValuesToList(dr_output)$dr_mat
        drugs_to_keep <- drugs[drugs %in%  colnames(dr_mat) ]
        
        if (length(drugs_to_keep) != 0) {
          actionButton('jump2P11', label = 'Next', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
        }
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$jump2P11, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "dd_tab")
    })
    
    
    
    ############################################################################################
    #                                  DRUG DISCOVERY (OUTPUT)                                #
    ############################################################################################
    
    output$info_box_13 <- renderText({
      HTML("<br/> <span style='font-weight:normal;'>
              <p align='justify'>
              Here, you can explore the neighborhood of drugs among the clusters generated with a specific clustering algorithm. </span> <br/>
              </p> <hr style='border-top: 1px solid white;'>

              <b> Drug-based approach </b>: <br/> <span style='font-weight:normal;'>
              <p align='justify'>
               If you are interesting in a specific drug, you can select it from the dropdown  menu and highlight 
               the cluster to which the drug belongs on the plot. Additionally, 
               two separate tables will display the genes and drugs in its neighborhood. </span>.<br/>
              </p> </span> <hr style='border-top: 1px solid white;'>

              <b> Cluster-based approach </b>: <br/> <span style='font-weight:normal;'>
              <p align='justify'>
               Here, you can explore only the clusters that contain at least one drug, and for the selected cluster, 
               retrieve the genes and the drugs in two different tables.</span> <br/>
              </p> <hr style='border-top: 1px solid white;'>
              ")
    })
    
    
    output$cl_table_selector2 <- shiny::renderUI({
      if (!is.null(shiny::reactiveValuesToList(clusters))){
        clusters <- shiny::reactiveValuesToList(clusters)
        cl_choices <- c()
        for (name in names(clusters)){
          if (!is.null(clusters[[name]])) {
            cl_choices <- c(cl_choices, name)
          }
        }
        shiny::selectInput(inputId = 'select_cluList2', label = 'Select a cluster table',
                           choices = cl_choices[!(cl_choices %in% c('dr_mat', 'hclust_output', 'manual', 'manual_cluster'))], selected = NULL)
      }
    })
    
    
    output$drug_cluster_selector <- renderUI({
      input$drug_button
      isolate({
        if (length(drugs()) == 0 & length(shiny::reactiveValuesToList(clusters)) == 0){
          return(NULL)} 
        
        if (input$approach == 'drug'){
          shiny::selectInput('drug_selector',
                             'Select one drug to highlight in the plot',
                             choices = drugs(), selected = NULL)
        } else {
          drug_cluster <- shiny::reactiveValuesToList(clusters)[[input$select_cluList2]]$cluster %>% 
            dplyr::mutate_if(is.factor, as.numeric) %>%
            filter(id %in% drugs()) %>% select(clust) %>% unique(.) %>% sort(.)
          shiny::selectInput('cluster_selector',
                             'Select one cluster to highlight in the plot',
                             choices = drug_cluster[[1]], selected = NULL)
        }
      })
    })
    
   
    output$drug_plot_box <- shiny::renderUI({
      input$show_drug
     
      isolate({
        if (length(shiny::reactiveValuesToList(clusters)) == 0 | length(drugs()) == 0){
          return(NULL)
        }else {
          output$drug_plot <- plotly::renderPlotly({
            isolate({
              
              clusters <- shiny::reactiveValuesToList(clusters)
              clu_anno <- clusters[[input$select_cluList2]]$cluster %>% arrange(clust)
              
              if (input$approach == 'cluster') {
                clu_anno <- clu_anno %>% mutate(clu_col = ifelse(clust == as.integer(input$cluster_selector), 'cluster', 'others'))
                clu_anno <- clu_anno %>% mutate(shape = ifelse(id %in% drugs(), 'drug', 'gene'))
                print(clu_anno)
                
                p <- MoNETA::plot_2D_matrix(coord = clusters$dr_mat[1:2,], nodes_anno = clu_anno,
                                            id_name = 'id', id_anno_color = 'clu_col', id_anno_shape = 'shape',
                                            title =  paste('Cluster2Drugs:', 
                                                           paste0('"',  input$cluster_selector, '"'), 
                                                           'cluster highlighted ',
                                                           paste0('(',  input$select_cluList2, ')')
                                            ),
                                            interactive = TRUE, wo_legend = FALSE)
                
                
              } else if  (input$approach == 'drug') {
                id_cluster <- clu_anno %>% filter(id == input$drug_selector) %>% select(clust)
                clu_anno <- clu_anno %>% mutate(clu_col = ifelse(clust == id_cluster[[1]], 'cluster', 'others'))
                clu_anno <- clu_anno %>% mutate(shape = ifelse(id %in% drugs(), 'drug', 'gene'))
                
                p <- MoNETA::plot_2D_matrix(coord =  clusters$dr_mat[1:2,], nodes_anno = clu_anno,
                                            id_name = 'id', id_anno_color = 'clu_col', id_anno_shape = 'shape',
                                            title = paste('Drug2Cluster:', 
                                                          paste0('"',  input$drug_selector, '"'), 
                                                          'drug found in the', id_cluster[[1]], 'cluster ',
                                                          paste0('(',  input$select_cluList2, ')')
                                            ),
                                            interactive = TRUE, wo_legend = FALSE)
              }
              return(p)
              
            })
          })
          
          return(shinydashboard::box(width = 12, #height = 500,
                                     shinycssloaders::withSpinner(plotly::plotlyOutput('drug_plot', height = '600px')),
                                     #shiny::hr(),
                                     #shiny::downloadButton(outputId = 'download_drug_plot', label = 'Download plot')
          )
          )
        }
        
      })
    })
    
    
    output$drug_table <- DT::renderDT(server = FALSE, {
      input$show_drug
      isolate({
        if (is.null(shiny::reactiveValuesToList(clusters)) | length(drugs()) == 0) {
          return(NULL)
        } else{
          req(input$cluster_selector)
          cluster <- shiny::reactiveValuesToList(clusters)[[input$select_cluList2]]$cluster %>% dplyr::mutate_if(is.factor, as.integer)
          
          if (input$approach == 'cluster'){
            
            onedb <- dplyr::filter(cluster , clust == as.integer(input$cluster_selector)) %>%
              filter(id %in% drugs())
            DT::datatable(onedb,
                          extensions = 'Buttons',
                          options = list(paging = TRUE,
                                         pageLength = 15,
                                         scrollX = TRUE,
                                         scrollY = TRUE,
                                         autoWidth = TRUE,
                                         dom = "Blfrtip",
                                         buttons = list("copy", list(
                                           extend = "collection",
                                           buttons = c("csv", "excel"),
                                           text = "Download"
                                         )
                                         )
                          ),
                          
                          selection = 'multiple',
                          filter = 'top',
                          rownames = FALSE
            )
          } else if (input$approach == 'drug') {
            id_cluster <- cluster %>% filter(id == input$drug_selector) %>% select(clust)
            onedb <- dplyr::filter(cluster , clust == id_cluster[[1]]) %>%
              filter(id %in% drugs())
            DT::datatable(onedb,
                          options = list(paging = TRUE,
                                         pageLength = 15,
                                         scrollX = TRUE,
                                         scrollY = TRUE,
                                         autoWidth = TRUE,dom = "Blfrtip",
                                         buttons = list("copy", list(
                                           extend = "collection",
                                           buttons = c("csv", "excel"),
                                           text = "Download"
                                         )
                                         )
                          ),
                          extensions = 'Buttons',
                          selection = 'multiple',
                          filter = 'top',
                          rownames = FALSE
            )
          }
        }
      })
    })
    
    
    output$gene_table <- DT::renderDT(server = FALSE, {
      req(input$show_drug)
      isolate({
        if (is.null(shiny::reactiveValuesToList(clusters)[[input$select_cluList2]]$cluster) & length(drugs()) == 0) {
          return(NULL)
        }else{
          req(input$cluster_selector)
          cluster <- shiny::reactiveValuesToList(clusters)[[input$select_cluList2]]$cluster %>% dplyr::mutate_if(is.factor, as.integer)
          
          if (input$approach == 'cluster') {
            onedb <- dplyr::filter(cluster , clust == as.integer(input$cluster_selector)) %>%
              filter(!(id %in% drugs()))
            DT::datatable(onedb,
                          extensions = 'Buttons',
                          options = list(paging = TRUE,
                                         pageLength = 15,
                                         scrollX = TRUE,
                                         scrollY = TRUE,
                                         autoWidth = TRUE,
                                         dom = "Blfrtip",
                                         buttons = list("copy", list(
                                           extend = "collection",
                                           buttons = c("csv", "excel"),
                                           text = "Download"
                                         )
                                         )
                          ),
                          selection = 'multiple',
                          filter = 'top',
                          rownames = FALSE
            )
          } else if (input$approach == 'drug') {
            id_cluster <- cluster %>% filter(id == input$drug_selector) %>% select(clust)
            onedb <- dplyr::filter(cluster , clust == id_cluster[[1]]) %>%
              filter(!(id %in% drugs()))
            DT::datatable(onedb,
                          options = list(paging = TRUE,
                                         pageLength = 15,
                                         scrollX = TRUE,
                                         scrollY = TRUE,
                                         autoWidth = TRUE,
                                         dom = "Blfrtip",
                                         buttons = list("copy", list(
                                           extend = "collection",
                                           buttons = c("csv", "excel"),
                                           text = "Download"
                                         )
                                         )
                          ),
                          extensions = 'Buttons',
                          selection = 'multiple',
                          filter = 'top',
                          rownames = FALSE
            )
          }
        }
      })
    })
    
    ###
    output$back2P10 <- renderUI({
      if (length(shiny::reactiveValuesToList(gProfiler_res)) != 0){
        actionButton('back2P10', label = 'Back', shiny::icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
      }else{
        return(NULL)
      }
    })
    
    observeEvent(input$back2P10, {
      shinydashboard::updateTabItems(session, inputId = "tabs",
                                     selected = "pathclu_tab")
    })
    
  }
