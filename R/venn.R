library(VennDiagram)
library(RColorBrewer)
library(tools)
library(here)


loadSets <- function(file_dir) {
  set_list <- list()
  files <- list.files(file_dir)
  print(files)
  for (f in files) {
    set_list[[tools::file_path_sans_ext(f)]] <- read.csv(file.path(file_dir, f), header = F)[, 1]
  }
  set_list
}

drawVenn <- function(file_dir, out_file_name) {
  set_list <- loadSets(file_dir)
  myCol <- brewer.pal(length(set_list), "Pastel2")[seq(1, length(set_list))]
  
  print(paste(length(set_list), " detected"))
  nameAdj <- switch(length(set_list),
                    NULL,
                    list(c(0, 0), c(0, 0)),
                    list(c(0, 0), c(0, 0), c(0, 0)),
                    list(c(0,0), c(0,0) , c(0,0) , c(0,0)),
                    list(c(0.6,1) , c(0,0) , c(0,0) , c(1,1) , c(1,0)),
                    )
  
  venn.diagram(
    x = set_list,
    category.names = names(set_list),
    filename = out_file_name,
    # output=TRUE,
    # 
    # # Output features
    imagetype="png" ,
    height = 1920 ,
    width = 1920 ,
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    # lty = "dashed",
    col = myCol,
    fill = myCol,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "text",
    cat.just=nameAdj,
    
    ext.pos = 30,
    ext.dist = -0.05,
    ext.length = 0.85,
    ext.line.lwd = 2,
    ext.line.lty = "dashed"
    # cat.pos = c(-27, 27, 135),
    # cat.dist = c(0.055, 0.055, 0.085),
    # cat.fontfamily = "sans",
    # rotation = 1
  )
}

drawVenn("../../DOX/results/1_degs/deg_venn/data/CC/down", 
         "../../DOX/results/1_degs/deg_venn/venn/CC/down.png")
drawVenn("../../DOX/results/1_degs/deg_venn/data/CC/up", 
         "../../DOX/results/1_degs/deg_venn/venn/CC/up.png")

drawVenn("../../DOX/results/1_degs/deg_venn/data/CO/down", 
         "../../DOX/results/1_degs/deg_venn/venn/CO/down.png")
drawVenn("../../DOX/results/1_degs/deg_venn/data/CO/up", 
         "../../DOX/results/1_degs/deg_venn/venn/CO/up.png")

drawVenn("../../DOX/results/1_degs/deg_venn/data/MC/down", 
         "../../DOX/results/1_degs/deg_venn/venn/MC/down.png")
drawVenn("../../DOX/results/1_degs/deg_venn/data/MC/up", 
         "../../DOX/results/1_degs/deg_venn/venn/MC/up.png")

drawVenn("../../DOX/results/1_degs/deg_venn/data/HT/down", 
         "../../DOX/results/1_degs/deg_venn/venn/HT/down.png")
drawVenn("../../DOX/results/1_degs/deg_venn/data/HT/up", 
         "../../DOX/results/1_degs/deg_venn/venn/HT/up.png")