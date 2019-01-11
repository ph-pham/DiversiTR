sideMenu <- sidebarMenu(
  id = "sideTabs",
  menuItem(
    selected = T,
    text = "About",
    icon = icon("info"),
    tabName = "aboutTab"
  ),
  convertMenuItem(menuItem(
    text = "Upload RDS file",
    icon = icon("upload"),
    fileInput(
      "RDSfile",
      NULL,
      multiple = T,
      accept = c("rds",
                 ".rds"),
      buttonLabel = "Choose File(s)",
      placeholder = "No file(s) selected"
    ),
    tabName = "uploadRDStab"),
    tabName = "uploadRDStab"
  ),
  convertMenuItem(menuItem(
    text = "Upload text files",
    icon = icon("upload"),
    selectizeInput(
      "source",
      'Source ?',
      choices = list("ClonotypeR", "MiXCR"),
      options = list(onInitialize = I('function() { this.setValue(""); }'))
    ),
    checkboxGroupInput(
      "chain",
      "chain ?",
      choiceNames = c("alpha", "beta"),
      choiceValues = c("A", "B")
    ),
    fileInput(
      "sInfofile",
      NULL,
      multiple = T,
      accept = c(
        "txt/tsv",        
        "text/tabulation-separated-values,text/plain",
        ".txt",
        ".tsv"
      ),
      buttonLabel = "Choose samples info File(s)",
      placeholder = "No file(s) selected"
    ),
    conditionalPanel(
      condition = "output.canUpload",
      fileInput(
        "samplefiles",
        NULL,
        multiple = T,
        accept = c(
          "txt/tsv",        
          "text/tabulation-separated-values,text/plain",
          ".txt",
          ".tsv"
        ),
        buttonLabel = "Choose sample File(s)",
        placeholder = "No file(s) selected"
      )
    ), tabName = "uploadTXTtab"
  ), tabName = "uploadTXTtab"),
  menuItemOutput("showDataTab"),
  #On peut ajouter les autres attribus d'un repseqexperiment
  #menuItemOutput : diversity in DT
  menuItemOutput("singleSampleTab"),
  menuItemOutput("multipleSampleTab"),
  shinyjs::useShinyjs(),
  menuItemOutput("downloadRDS")
)