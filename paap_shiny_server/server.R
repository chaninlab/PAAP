library(protr)
library(seqinr)
library(randomForest)


shinyServer(function(input, output, session) {
 
  mod <- readRDS("model.rds")

  observe({
    
    shinyjs::hide("downloadData") # Hide download button before input submission
    if(input$submitbutton>0)
      shinyjs::show("downloadData") # Show download button after input submission
    
    FASTADATA <- ''
    fastaexample <- '>AHTP1
AHEPVK
>AHTP2
AQTQSL
>AHTP3
ESIINF
>AHTP4
GGVIPN
>nonAHTP1
QVQADR
>nonAHTP2
RDKYFL
>nonAHTP3
RGDVIL
>nonAHTP4
RGGGLE
'
    
    if(input$addlink>0) {
      isolate({
        FASTADATA <- fastaexample
        updateTextInput(session, inputId = "Sequence", value = FASTADATA)
      })
    }
  })
  
  datasetInput <- reactive({
    
    inFile <- input$file1 
    inTextbox <- input$Sequence

    if (is.null(inTextbox)) {
      return("Please insert/upload sequence in FASTA format")
    } else {
      if (is.null(inFile)) {
        # Read data from text box
        x <- inTextbox
        write.fasta(sequence = x, names = names(x),
                    nbchar = 80, file.out = "text.fasta")
        x <- readFASTA("text.fasta")
        
        # Feature extraction for Testing set
        xtest <- read.fasta('example.fasta', seqtype="AA", as.string = TRUE)###read data
        xtest2 <- xtest[(sapply(xtest, protcheck))]###check special symbol
        m2 = length(xtest2)
        paactest <- matrix(nrow = m2, ncol = 23)
        for(i in 1:m2){ 
          paactest[i, ]= extractPAAC(xtest2[[i]][1],lambda = 3, w = 0.1 , props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"))
        }
        data <- data.frame(paactest)
        
        # Predicting unknown sequences
        results <- data.frame(Prediction= predict(mod,data))
        
        print(results)
      } 
      else {  
        # Read data from uploaded file
        x <- readFASTA(inFile$datapath)
        
        # Feature extraction for Testing set
        xtest <- read.fasta('example.fasta', seqtype="AA", as.string = TRUE)###read data
        xtest2 <- xtest[(sapply(xtest, protcheck))]###check special symbol
        m2 = length(xtest2)
        paactest <- matrix(nrow = m2, ncol = 23)
        for(i in 1:m2){ 
          paactest[i, ]= extractPAAC(xtest2[[i]][1],lambda = 3, w = 0.1 , props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"))
        }
        data <- data.frame(paactest)
        
        # Predicting unknown sequences
        results <- data.frame(Prediction= predict(mod,data))
        
        print(results)
      }
    }
  })
  
  output$contents <- renderPrint({
    if (input$submitbutton>0) { 
      isolate(datasetInput()) 
    } else {
      return("Server is ready for prediction.")
    }
  })
  
  output$downloadData <- downloadHandler(
    filename = function() { paste('predicted_results', '.csv', sep='') },
    content = function(file) {
      write.csv(datasetInput(), file, row.names=FALSE)
    })
  
})
