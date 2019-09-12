datdir = "tcga_refined_data" ## Data for analysise
sparsematrixdir = "tcga_sparse_matrix"
ECTDdir = "tcga_ECTD_matrix"
resultdir = "tcga_results"
biodir = "tcga_Biological_Evaluation"

## Graphical Measures
dat.normal.v <- 828804.5 # Volumn of Adjacency Matrix
dat.cancer.v <- 732223.2 # Volumn of Adjacency Matrix

# global clustering coefficient
dat.normal.glasso.ECTD.globalcc <- 100 
dat.cancer.glasso.ECTD.globalcc <- 100
dat.normal.MoorePenrose.ECTD.globalcc <- 100 
dat.cancer.MoorePenrose.ECTD.globalcc <- 100
dat.normal.cov.globalcc <- 92.15904
dat.cancer.cov.globalcc <- 92.92874

# multiplot
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}