
convertMenuItem <- function(mi, tabName) {
  mi$children[[1]]$attribs['data-toggle']="tab"
  mi$children[[1]]$attribs['data-value'] = tabName
  if(length(mi$attribs$class)>0 && mi$attribs$class=="treeview"){
    mi$attribs$class=NULL
  }
  mi
}


header <-dashboardHeader(titleWidth = 250, title = em("GCDS simulation"))
sidebar<-dashboardSidebar(
  sidebarUserPanel("Seungjun & Minseok",subtitle= p("Simulation with MNIST example!")),
  sidebarMenu(id = "menu",
              wellPanel(style = "background-color: #170c45",  
                        radioButtons("image_ratio", "Input image size(MNIST)",choices = c("1/2 patches"=1/2,"3/4 patches"=3/4)
                                  ),
                        ),
              wellPanel(style = "background-color: #170c45",  
                        sliderInput("fig_num", "# of predictions for 1 input image patch",min=1,max=10,value=5
                                  ),
                        ),
              wellPanel(style = "background-color: #170c45",  
                        checkboxGroupInput("num_check", "what numbers would you like to display?",
                                           c("0"=0,"1"=1,"2"=2,"3"=3,"4"=4,"5"=5,"6"=6
                                             ,"7"=7,"8"=8,"9"=9),
                                           selected=(c("0","1","2","3","4")),
                                           inline=TRUE
                                           )
                        ),
              wellPanel(style = "background-color: #170c45",  
                        checkboxInput("start_check", "Prediction start!"
                        ),
              )
  )
)

body<- dashboardBody( 
  tabBox(width = 12, 
         
         tabPanel("Paper",
                  h3("A Deep Generative Approach to Conditional Sampling(2022) "),
                  p("Xingyu Zhou et al."),
                  br(),
                  br(),
                  strong("Abstract"),
                  p("We propose a deep generative approach to sampling from a conditional distribution based on a unified formulation of conditional distribution and generalized nonparametric regression function using the noise-outsourcing lemma. The proposed approach aims at learning a conditional generator, so that a random sample from the target conditional distribution can be obtained by transforming a sample drawn from a reference distribution. The conditional generator is estimated nonparametrically with neural networks by matching appropriate joint distributions using the Kullback-Liebler divergence. An appealing aspect of our method is that it allows either of or both the predictor and the response to be high-dimensional and can handle both continuous and discrete type predictors and responses. We show that the proposed method is consistent in the sense that the conditional generator converges in distribution to the underlying conditional distribution under mild conditions. Our numerical experiments with simulated and benchmark image data validate the proposed method and demonstrate that it outperforms several existing conditional density estimation methods. Supplementary materials for this article are available online."),
                  br(),
                  a(href = "https://www.tandfonline.com/doi/full/10.1080/01621459.2021.2016424",
                    "▶ More specific infomation is here. "),
                  br(),
                  a(href = "https://github.com/Seungjun1588/expri/tree/main/Rshiny",
                    "▶ R coodes for this shiny apps are here.") 

                  ),
         tabPanel("Structure of Generator",
                  tags$img(src = "gen_image.png",height = 800,witdh = 800)),
         tabPanel("Structure of Discriminator",
                  tags$img(src = "disc_image.png",height = 800,witdh = 800)),
         tabPanel("Prediction",
                  plotOutput("mnist_output"))
  ),

)

dashboardPage(header,sidebar,body,skin = "purple") # skin: