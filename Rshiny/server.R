server <- function(input, output, session) {
  # test.generator <- reactive({
  #   if(input$image_ratio == "0.5"){
  #     
  #   }
  #   gcde.mnist.x34image.cnn.test()
  # })
  
  output$mnist_output <- renderPlot({
  if(input$image_ratio == "0.5"){
    test.generator = gcde.mnist.xhalfimage.cnn.test()
    if(input$start_check){
      # scaling process when input == 0.5
      x_train <- array(mnist$train$x[, 1:28, 1:14], dim = c(60000, 28,14,1))
      x_train <- (x_train/255-0.5)/0.5
      y_train <- array((mnist$train$x/255 - 0.5)/0.5, dim = c(60000, 28, 28, 1))
      x_test <- array(mnist$test$x[, 1:28, 1:14], dim = c(10000, 28,14,1))
      x_test <- (x_test/255 - 0.5)/0.5
      y_test <- array((mnist$test$x/255 - 0.5)/0.5, dim = c(10000, 28, 28, 1))
      

      # load the pretrained weight
      load_model_weights_hdf5(test.generator$generator, "trained/mnist_xhalfimage_cnn_primal_generator.h5")
      withProgress(message = 'Please wait a second. :)',
                   plot.whiletraining.mnist.xhalfimage.cnn(test.generator$generator.t, "pic19", dual = "fenchel",train.x=x_test, pic.num = input$fig_num+1,test_num=as.numeric(input$num_check)),
                   min = 0,
                   max = 1,
                   value = 10
      )
    }else{
      # Noting to do
    }  
  } 
  else{ # image_ratio == "0.75"
    test.generator = gcde.mnist.x34image.cnn.test()
    if(input$start_check){
      # scaling process when input ==0.75
      x_train <- array(mnist$train$x[, 1:28, 1:28], dim = c(60000, 28,28,1))
      x_train <- (x_train/255-0.5)/0.5 # scaling
      x_train[1:60000, 15:28, 15:28,1] <- 0
      y_train <- array((mnist$train$x[,15:28,15:28]/255 - 0.5)/0.5, dim = c(60000, 14, 14, 1))
      x_test <- array(mnist$test$x[, 1:28, 1:28], dim = c(10000, 28,28,1))
      x_test[1:10000, 15:28,15:28,1] <- 0
      x_test <- (x_test/255 - 0.5)/0.5
      y_test <- array((mnist$test$x[, 15:28, 15:28]/255 - 0.5)/0.5, dim = c(10000, 14, 14, 1))
      

      # load the pretrained weight
      load_model_weights_hdf5(test.generator$generator, "trained/mnist_x34image_cnn_primal_generator.h5")
      withProgress(message = 'Please wait a second. :)',
                   plot.whiletraining.mnist.x34image.cnn(test.generator$generator.t, "pic19", dual = "fenchel",train.x=x_test, pic.num = input$fig_num+1,test_num=as.numeric(input$num_check)),
                   min = 0,
                   max = 1,
                   value = 10
      )
    }else{
      # Noting to do
    }    
  }
    
  

  # load_model_weights_hdf5(test.generator()$generator, "trained/mnist_x34image_cnn_primal_generator.h5")
  # #print(test.generator)
  # withProgress(message = 'Please wait a seconds. :)',plot.whiletraining.mnist.x34image.cnn(test.generator()$generator.t, "pic19", dual = "fenchel", pic.num = 5))
    
  })
}
