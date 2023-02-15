# For virtual env
# install.packages("tensorflow")
library(tensorflow)
# #install_tensorflow()
library(reticulate)
# conda_create("tensorflow_200")
# install_tensorflow(envname="tensorflow_200")
#use_condaenv(condaenv = "tensorflow_200",required = TRUE)
#install_keras(envname="tensorflow_200")
# py_config()


#------------------------------------------------------------------------------#
##################################   SET UP  ##########################################
library(keras)
library(tensorflow)
library(xtable)
library(progress)#progress bar
library(parallel)#Used in generate table
library(pracma)#gramschimidt for orthonormal matrix
library(abind)
library(tidyverse)
library(shiny)
library(shinydashboard)

#Set to use double-precision floating point
is_keras_available()
tf$keras$backend$set_floatx('float32')
#set to only allocate GPU memory when needed.
if(length(tf$config$experimental$list_physical_devices("GPU")))
  tf$config$experimental$set_memory_growth(tf$config$experimental$list_physical_devices("GPU")[[1]], TRUE)


#------------------------------------------------------------------------------#
# data setting
load("mnist.Rdata")
# str(mnist)
# dim(mnist$train$x)
# image(mnist$train$x[1,,])
# image(mnist$train$x[2,,])
# image(mnist$train$x[3,,])
# image(mnist$train$x[4,,])

# mnist$train$y[1:4]

#x is the leftbottom 1/4 image. Y is the rest of the image.
# x_train <- array(mnist$train$x[, 1:28, 1:28], dim = c(60000, 28,28,1))
# x_train <- (x_train/255-0.5)/0.5 # scaling
# x_train[1:60000, 15:28, 15:28,1] <- 0
# # image(x_train[1,,,1]) # 3/4 ??
# 
# 
# y_train <- array((mnist$train$x[,15:28,15:28]/255 - 0.5)/0.5, dim = c(60000, 14, 14, 1))
# # image(y_train[1,,,1]) # 1/4 ??
# 
# 
# x_test <- array(mnist$test$x[, 1:28, 1:28], dim = c(10000, 28,28,1))
# x_test[1:10000, 15:28,15:28,1] <- 0
# x_test <- (x_test/255 - 0.5)/0.5
# y_test <- array((mnist$test$x[, 15:28, 15:28]/255 - 0.5)/0.5, dim = c(10000, 14, 14, 1))

#------------------------------------------------------------------------------#
# 3/4 diagnostic function 
plot.whiletraining.mnist.x34image.cnn <- function(generator.t, ite, train.x = x_test, 
                                                  orig.image = mnist$test$x[1:62,,],
                                                  labels = mnist$train$y, padding = "black", dual = "fenchel",
                                                  pic.num = 19,test_num){
  #x11(width = 12, height = 8, title = paste0("iteration = ", ite))
  test.id <- c(4, 3, 2, 19, 5, 16, 12, 1, 62, 8)
  test.id  <- test.id[test_num + 1]
  real.plot <- function(){
    par(mfrow = c(10,pic.num +3), mar = c(0,0,0,0))
    for (i in test.id) {
      plot.new()
      plot.image <- matrix(ifelse(padding == "white", 255, -255), 32, 32)
      plot.image[3:30,3:30] <- train.x[i, , ,1] * 255
      
      n <- 1
      dim.eta  <- 100
      
      ori.im <- matrix(ifelse(padding == "white", 255, 0), 32, 32)
      ori.im[3:30,3:30] <- orig.image[i,,]
      image(t(ori.im)[, ncol(ori.im):1], axes = FALSE, 
            col = gray.colors(255, start = 0, end = 1))
      rect(0,0,1,1, border = "white")
      #text(0,0.95, labels[i], cex=1, pos=4, col="red", offset=0)
      for(j in 1:pic.num){
        t <- as.array(predict(generator.t,list(array(train.x[i, , , ], dim = c(1,28, 28,1)), 
                                               matrix(rnorm(n*dim.eta), nrow = n, ncol = dim.eta)), 
                              training = FALSE))
        plot.image[17:30, 17:30] <- t[1, 1:14, 1:14,1]*255
        if(j == 1){
          plot.image[17:30, 17:30] <- -255
        }
        image(t(plot.image)[, ncol(plot.image):1], axes = FALSE, col = gray.colors(255, start = 0, end = 1))
        lines(x = c(0, 0.5, 0.5, 1, 1, 0), y = c(0, 0, 0.5, 0.5, 1,1), col = "gray")
        #text(0,0.95, labels[i], cex=1, pos=4, col="red", offset=0)
      }
      plot.new()
    }
  }
  real.plot()

}

# 1/2 diagnostic function
plot.whiletraining.mnist.xhalfimage.cnn <- function(generator.t, ite, train.x = x_test, 
                                                    orig.image = mnist$test$x[1:62,,],
                                                    labels = mnist$train$y, padding = "black", dual = "fenchel",
                                                    pic.num = 19,test_num){
  #x11(width = 12, height = 8, title = paste0("iteration = ", ite))
  test.id <- c(4, 3, 2, 19, 5, 16, 12, 1, 62, 8)
  test.id  <- test.id[test_num + 1]
  
  real.plot <- function(){
    par(mfrow = c(10,pic.num + 3), mar = c(0,0,0,0))
    for (i in test.id) {
      plot.new()
      plot.image <- matrix(ifelse(padding == "white", 255, -255), 32, 32)
      plot.image[3:30,3:16] <- train.x[i, , ,1] * 255
      
      n <- 1
      dim.eta  <- 100
      
      ori.im <- matrix(ifelse(padding == "white", 255, 0), 32, 32)
      ori.im[3:30,3:30] <- orig.image[i,,]
      image(t(ori.im)[, ncol(ori.im):1], axes = FALSE, 
            col = gray.colors(255, start = 0, end = 1))
      rect(0,0,1,1, border = "white")
      #text(0,0.95, labels[i], cex=1, pos=4, col="red", offset=0)
      for(j in 1:pic.num){
        t <- as.array(predict(generator.t,
                              list(array(train.x[i, , , ], dim = c(1,28, 14,1)), 
                                   matrix(rnorm(n*dim.eta), nrow = n, ncol = dim.eta)), 
                              training = FALSE))
        plot.image[3:30, 17:30] <- t[1, 1:28, 15:28,1]*255
        if(j == 1){
          plot.image[3:30, 17:30] <- -255
        }
        image(t(plot.image)[, ncol(plot.image):1], axes = FALSE, col = gray.colors(255, start = 0, end = 1))
        rect(0,0,0.5,1, border = "gray")
        #text(0,0.95, labels[i], cex=1, pos=4, col="red", offset=0)
      }
      plot.new()
    }
  }
  #dir.create("training_plot/xhalfimage")
  real.plot()
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# only predict ; not for train!
gcde.mnist.x34image.cnn.test = function(disc.learning.rate = 0.0001,
                                        gene.learning.rate = 0.0001,
                                        kl.dual = FALSE,
                                        v.label = "new",
                                        dual.form = NULL,
                                        dim.eta= 100,
                                        lambda.mse= 0,
                                        lambda.kl = 20 ){
  
  discriminator_input_x <- layer_input(shape = c(28, 28,1))
  discriminator_input_y <- layer_input(shape = c(14, 14, 1))
  
  disc.x.cnn <- layer_conv_2d(discriminator_input_x, 64, kernel_size = 5, strides = 2,
                              padding = "same")%>%
    layer_batch_normalization(trainable = TRUE)%>%
    layer_activation_leaky_relu(alpha = 0.2)%>%
    layer_conv_2d(128, kernel_size = 4, strides = 1, padding = "same")%>%
    layer_batch_normalization(trainable = TRUE)%>%
    layer_activation_leaky_relu(alpha = 0.2)%>%
    layer_flatten()
  disc.y.cnn <- layer_conv_2d(discriminator_input_y, 32, kernel_size = 5, strides = 2,
                              padding = "same")%>%
    layer_batch_normalization(trainable = TRUE)%>%
    layer_activation_leaky_relu(alpha = 0.2)%>%
    layer_conv_2d(64, kernel_size = 5, strides = 2, padding = "same")%>%
    layer_batch_normalization(trainable = TRUE)%>%
    layer_activation_leaky_relu(alpha = 0.2)%>%
    layer_flatten()
  disc.xy <- layer_concatenate(list(disc.x.cnn, disc.y.cnn), axis = -1L)%>%
    #layer_conv_2d(128, kernel_size = 7, strides = 1, padding = "valid")%>%
    #layer_batch_normalization(trainable = TRUE)%>%
    #layer_activation_leaky_relu(alpha = 0.2)%>%
    #layer_flatten()%>%
    # 64 -64
    layer_dense(128)%>%
    layer_batch_normalization(trainable = TRUE)%>%
    layer_activation_leaky_relu(alpha = 0.2)%>%
    layer_dense(64)%>%
    layer_batch_normalization(trainable = TRUE)%>%
    layer_activation_leaky_relu(alpha = 0.2)
  
  disc.output <- layer_dense(disc.xy, 1, activation = NULL)
  
  
  
  discriminator <- keras_model(inputs = c(discriminator_input_x, discriminator_input_y), 
                               outputs = disc.output)
  
  discriminator.structure <- capture.output(summary(discriminator))
  discriminator.optimizer <- tf$keras$optimizers$Adam(learning_rate = disc.learning.rate, beta_1 = 0.5)
  
  
  #Note that y_pred is actually V. V is the label of {-1,1} indicating True or generated Y.
  if(kl.dual == FALSE)
    discri.loss <- function(y_true, y_pred){
      #Add orthonormal constraints here
      if(useddr & (lambda.ortho.loss != 0)){
        weightmatrix <- reducer$weights[[1]]
        btb <- tf$matmul(tf$transpose(weightmatrix), weightmatrix)
        loss.orthonormal <- 
          tf$norm(btb - tf$linalg$diag(k_constant(1, dtype = "float32", shape = dr.nodes.num)))
      } else loss.orthonormal <- tf$constant(0, dtype = "float32")
      
      k_mean(tf$math$log(k_constant(1, dtype = "float32") + tf$math$exp(-y_true*y_pred) ))+ 
        k_constant(lambda.ortho.loss , dtype = "float32") * loss.orthonormal
    }
  
  #Note that y_true is v labels. V is 1 if it is from true data. V is -1 if it is from generated(fake).
  
  if(kl.dual == TRUE){
    if(v.label == "old")
      stop("Use new labels in dual form, please.")
    
    if(v.label == "new")
      if(dual.form == "fenchel"){
        message("Fenchel dual used.")
        discri.loss <- function(y_true, y_pred){
          #Fenchel dual
          - k_mean(y_pred * (y_true + 1) * (1/2)  - tf$math$exp(y_pred) * (y_true - 1) * (-1/2) )
        }
      } else if(dual.form == "donsker"){
        message("Donsker dual used.")
        discri.loss <- function(y_true, y_pred){
          #Donsker dual
          - k_mean(y_pred * (y_true + 1) * (1/2) * 2) + 
            tf$math$log(k_mean(tf$math$exp(y_pred) * (y_true - 1) * (-1/2)) * 2  )
        }
      } else stop("Either fenchel or donsker dual may be used in dual form.")
  }
  
  
  discriminator %>% compile(
    optimizer = discriminator.optimizer,
    loss = discri.loss
  )
  #*********************************   END of Define discriminator  ********************************
  
  
  
  ###################################   Define GENERATOR   ##########################################
  # Set discriminator weights to non-trainable
  # (will only apply to the `generator` model)
  freeze_weights(discriminator)
  #freeze_weights(discriminator)
  
  gene.input.x              <- layer_input(shape = c(28, 28, 1))
  gene.input.eta            <- layer_input(shape = c(dim.eta))
  
  #gene.conv <- layer_conv_2d(gene.input.x, 32, kernel_size = 5, strides = 2,
  #                           padding = "same")%>%
  #    layer_batch_normalization(trainable = TRUE)%>%
  #    layer_activation_leaky_relu(alpha = 0.2)%>%
  #    layer_conv_2d(64, kernel_size = 4, strides = 2, padding = "valid")%>%
  #    layer_batch_normalization(trainable = TRUE)%>%
  #    layer_activation_leaky_relu(alpha = 0.2)%>%
  #    layer_flatten()
  gene.conv <- layer_flatten(gene.input.x)%>%
    #64
    layer_dense(128)%>%
    layer_batch_normalization(trainable = TRUE)%>%
    layer_activation_leaky_relu(alpha = 0.2)
  gene.combined <- layer_concatenate(list(gene.conv, gene.input.eta))%>%
    layer_reshape(target_shape = c(1,1,228))
  gene.output.t <- layer_conv_2d_transpose(gene.combined, 256, kernel_size = c(7,7), strides = c(1,1), 
                                           padding = "valid")%>%
    layer_batch_normalization(trainable = TRUE)%>%
    layer_activation_leaky_relu(alpha = 0.2)%>%
    layer_conv_2d_transpose(128, kernel_size = c(5,5), strides = c(2,2), padding = "same")%>%
    layer_batch_normalization(trainable = TRUE)%>%
    layer_activation_leaky_relu(alpha = 0.2)%>%
    layer_conv_2d_transpose(1, kernel_size = c(5,5), strides = c(1,1), padding = "same")%>%
    layer_activation(activation = "tanh")
  
  
  
  #calculate the log density ratio
  gene.output.d <- discriminator(list(gene.input.x, gene.output.t))
  
  #estimate the conditional expectation
  generator.t <- keras_model(inputs = c(gene.input.x, gene.input.eta), outputs = c(gene.output.t))
  
  
  #t is T(X,eta). d is D(T(X,eta)). g is g(X).
  
  generator <- keras_model(inputs = c(gene.input.x, gene.input.eta), 
                           outputs = c(gene.output.t, gene.output.d))
  
  generator.structure <- capture.output(summary(generator))
  #summary(generator)
  
  
  #generator.optimizer <- optimizer_adam( 
  #    lr = gene.learning.rate, 
  #    decay = gene.lr.decay)
  
  generator.optimizer <- tf$keras$optimizers$Adam(learning_rate = gene.learning.rate, beta_1 = 0.5)
  
  loss.zero <- function(y_true, y_pred){
    #Add orthonormal loss if use ddr
    if(useddr & (lambda.ortho.loss != 0)){
      weightmatrix <- reducer$weights[[1]]
      btb <- tf$matmul(tf$transpose(weightmatrix), weightmatrix)
      loss.orthonormal <- 
        tf$norm(btb - tf$linalg$diag(k_constant(1, dtype = "float32", shape = dr.nodes.num)))
    } else loss.orthonormal <- tf$constant(0, dtype = "float32")
    k_constant(lambda.ortho.loss , dtype = "float32") * loss.orthonormal
    #k_constant(0, dtype = "double")
  }
  
  #The kl loss part for both primal and dual form are the same
  if(v.label == "new"){
    loss.kl <- function(y_true, y_pred){
      k_mean(y_pred)}
  }
  
  if(v.label == "old"){
    if(kl.dual == TRUE) stop("If you want to use dual form, please use new label.")
    loss.kl <- function(y_true, y_pred){
      k_mean(-y_pred)
    }
  }
  
  if(lambda.mse == 0){
    loss.mse <- function(y_true, y_pred)
      k_constant(0, dtype = "float32")
  } else {
    loss.mse <- function(y_true, y_pred){
      k_mean(tf$math$squared_difference(y_true, y_pred)) * dim.y
    }
  }
  
  generator %>% compile(
    optimizer = generator.optimizer, loss = list(loss.zero, loss.kl),
    loss_weights = list(1, lambda.kl)
  )
  return(list(generator = generator, generator.t = generator.t,
              discriminator = discriminator)
  )
}

# 1/2
test.generator = gcde.mnist.x34image.cnn.test()
test.generator$generator
test.generator = gcde.mnist.xhalfimage.cnn.test()
test.generator$generator
gcde.mnist.xhalfimage.cnn.test = function(disc.learning.rate = 0.0001,
                                        gene.learning.rate = 0.0001,
                                        kl.dual = FALSE,
                                        v.label = "new",
                                        dual.form = NULL,
                                        dim.eta= 100,
                                        lambda.mse= 0,
                                        lambda.kl = 20 ){
  
  discriminator_input_x <- layer_input(shape = c(28, 14,1))
  discriminator_input_y <- layer_input(shape = c(28, 28, 1))
  
  disc.x.cnn <- layer_conv_2d(discriminator_input_x, 32, kernel_size = 5, strides = 2,
                              padding = "same")%>%
    layer_batch_normalization(trainable = TRUE)%>%
    layer_activation_leaky_relu(alpha = 0.2)%>%
    layer_conv_2d(64, kernel_size = 4, strides = 1, padding = "same")%>%
    layer_batch_normalization(trainable = TRUE)%>%
    layer_activation_leaky_relu(alpha = 0.2)%>%
    layer_flatten()
  disc.y.cnn <- layer_conv_2d(discriminator_input_y, 64, kernel_size = 5, strides = 2,
                              padding = "same")%>%
    layer_batch_normalization(trainable = TRUE)%>%
    layer_activation_leaky_relu(alpha = 0.2)%>%
    layer_conv_2d(128, kernel_size = 5, strides = 2, padding = "same")%>%
    layer_batch_normalization(trainable = TRUE)%>%
    layer_activation_leaky_relu(alpha = 0.2)%>%
    layer_flatten()
  disc.xy <- layer_concatenate(list(disc.x.cnn, disc.y.cnn), axis = -1L)%>%
    #layer_conv_2d(128, kernel_size = 7, strides = 1, padding = "valid")%>%
    #layer_batch_normalization(trainable = TRUE)%>%
    #layer_activation_leaky_relu(alpha = 0.2)%>%
    #layer_flatten()%>%
    # 64 -64
    layer_dense(256)%>%
    layer_batch_normalization(trainable = TRUE)%>%
    layer_activation_leaky_relu(alpha = 0.2)%>%
    layer_dense(128)%>%
    layer_batch_normalization(trainable = TRUE)%>%
    layer_activation_leaky_relu(alpha = 0.2)
  
  disc.output <- layer_dense(disc.xy, 1, activation = NULL)
  
  
  
  discriminator <- keras_model(inputs = c(discriminator_input_x, discriminator_input_y), 
                               outputs = disc.output)
  
  discriminator.structure <- capture.output(summary(discriminator))
  #summary(discriminator)
  
  # To stabilize training, we use learning rate decay
  # and gradient clipping (by value) in the optimizer.
  #discriminator.optimizer <- optimizer_adam( 
  #    lr = disc.learning.rate, 
  #    decay = disc.lr.decay)
  discriminator.optimizer <- tf$keras$optimizers$Adam(learning_rate = disc.learning.rate, beta_1 = 0.5)
  
  
  #Note that y_pred is actually V. V is the label of {-1,1} indicating True or generated Y.
  if(kl.dual == FALSE)
    discri.loss <- function(y_true, y_pred){
      #Add orthonormal constraints here
      if(useddr & (lambda.ortho.loss != 0)){
        weightmatrix <- reducer$weights[[1]]
        btb <- tf$matmul(tf$transpose(weightmatrix), weightmatrix)
        loss.orthonormal <- 
          tf$norm(btb - tf$linalg$diag(k_constant(1, dtype = "float32", shape = dr.nodes.num)))
      } else loss.orthonormal <- tf$constant(0, dtype = "float32")
      
      k_mean(tf$math$log(k_constant(1, dtype = "float32") + tf$math$exp(-y_true*y_pred) ))+ 
        k_constant(lambda.ortho.loss , dtype = "float32") * loss.orthonormal
    }
  
  #Note that y_true is v labels. V is 1 if it is from true data. V is -1 if it is from generated(fake).
  if(kl.dual == TRUE){
    if(v.label == "old")
      stop("Use new labels in dual form, please.")
    
    if(v.label == "new")
      if(dual.form == "fenchel"){
        message("Fenchel dual used.")
        discri.loss <- function(y_true, y_pred){
          #Fenchel dual
          - k_mean(y_pred * (y_true + 1) * (1/2)  - tf$math$exp(y_pred) * (y_true - 1) * (-1/2) )
        }
      } else if(dual.form == "donsker"){
        message("Donsker dual used.")
        discri.loss <- function(y_true, y_pred){
          #Donsker dual
          - k_mean(y_pred * (y_true + 1) * (1/2) * 2) + 
            tf$math$log(k_mean(tf$math$exp(y_pred) * (y_true - 1) * (-1/2)) * 2  )
        }
      } else stop("Either fenchel or donsker dual may be used in dual form.")
  }
  
  
  discriminator %>% compile(
    optimizer = discriminator.optimizer,
    loss = discri.loss
  )
  #*********************************   END of Define discriminator  ********************************
  
  
  
  ###################################   Define GENERATOR   ##########################################
  # Set discriminator weights to non-trainable
  # (will only apply to the `generator` model)
  freeze_weights(discriminator)
  #freeze_weights(discriminator)
  
  gene.input.x              <- layer_input(shape = c(28, 14, 1))
  gene.input.eta            <- layer_input(shape = c(dim.eta))
  
  #gene.conv <- layer_conv_2d(gene.input.x, 32, kernel_size = 5, strides = 2,
  #                           padding = "same")%>%
  #    layer_batch_normalization(trainable = TRUE)%>%
  #    layer_activation_leaky_relu(alpha = 0.2)%>%
  #    layer_conv_2d(64, kernel_size = 4, strides = 2, padding = "valid")%>%
  #    layer_batch_normalization(trainable = TRUE)%>%
  #    layer_activation_leaky_relu(alpha = 0.2)%>%
  #    layer_flatten()
  gene.conv <- layer_flatten(gene.input.x)%>%
    #64
    layer_dense(128)%>%
    layer_batch_normalization(trainable = TRUE)%>%
    layer_activation_leaky_relu(alpha = 0.2)
  gene.combined <- layer_concatenate(list(gene.conv, gene.input.eta))%>%
    layer_reshape(target_shape = c(1,1,228))
  gene.output.t <- layer_conv_2d_transpose(gene.combined, 256, kernel_size = c(7,7), strides = c(1,1), 
                                           padding = "valid")%>%
    layer_batch_normalization(trainable = TRUE)%>%
    layer_activation_leaky_relu(alpha = 0.2)%>%
    layer_conv_2d_transpose(128, kernel_size = c(5,5), strides = c(2,2), padding = "same")%>%
    layer_batch_normalization(trainable = TRUE)%>%
    layer_activation_leaky_relu(alpha = 0.2)%>%
    layer_conv_2d_transpose(1, kernel_size = c(5,5), strides = c(2,2), padding = "same")%>%
    layer_activation(activation = "tanh")
  
  
  
  #calculate the log density ratio
  gene.output.d <- discriminator(list(gene.input.x, gene.output.t))
  
  #estimate the conditional expectation
  generator.t <- keras_model(inputs = c(gene.input.x, gene.input.eta), outputs = c(gene.output.t))
  
  
  #t is T(X,eta). d is D(T(X,eta)). g is g(X).
  
  generator <- keras_model(inputs = c(gene.input.x, gene.input.eta), 
                           outputs = c(gene.output.t, gene.output.d))
  
  generator.structure <- capture.output(summary(generator))
  #summary(generator)
  
  
  #generator.optimizer <- optimizer_adam( 
  #    lr = gene.learning.rate, 
  #    decay = gene.lr.decay)
  
  generator.optimizer <- tf$keras$optimizers$Adam(learning_rate = gene.learning.rate, beta_1 = 0.5)
  
  loss.zero <- function(y_true, y_pred){
    #Add orthonormal loss if use ddr
    if(useddr & (lambda.ortho.loss != 0)){
      weightmatrix <- reducer$weights[[1]]
      btb <- tf$matmul(tf$transpose(weightmatrix), weightmatrix)
      loss.orthonormal <- 
        tf$norm(btb - tf$linalg$diag(k_constant(1, dtype = "float32", shape = dr.nodes.num)))
    } else loss.orthonormal <- tf$constant(0, dtype = "float32")
    k_constant(lambda.ortho.loss , dtype = "float32") * loss.orthonormal
    #k_constant(0, dtype = "double")
  }
  
  #The kl loss part for both primal and dual form are the same
  if(v.label == "new"){
    loss.kl <- function(y_true, y_pred){
      k_mean(y_pred)}
  }
  
  if(v.label == "old"){
    if(kl.dual == TRUE) stop("If you want to use dual form, please use new label.")
    loss.kl <- function(y_true, y_pred){
      k_mean(-y_pred)
    }
  }
  
  if(lambda.mse == 0){
    loss.mse <- function(y_true, y_pred)
      k_constant(0, dtype = "float32")
  } else {
    loss.mse <- function(y_true, y_pred){
      k_mean(tf$math$squared_difference(y_true, y_pred)) * dim.y
    }
  }
  
  generator %>% compile(
    optimizer = generator.optimizer, loss = list(loss.zero, loss.kl),
    loss_weights = list(1, lambda.kl)
  )
  
  return(list(generator = generator, generator.t = generator.t,
              discriminator = discriminator)
  )
}

#------------------------------------------------------------------------------#
#test.generator = gcde.mnist.x34image.cnn.test()
# load the trained weight.
# load_model_weights_hdf5(test.generator$generator, "trained/mnist_x34image_cnn_primal_generator.h5")
# plotting
# plot.whiletraining.mnist.x34image.cnn(test.generator$generator.t, "pic19", dual = "fenchel", pic.num = 5)
#------------------------------------------------------------------------------#












