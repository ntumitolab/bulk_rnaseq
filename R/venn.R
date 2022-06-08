library(VennDiagram)
library(RColorBrewer)
library(tools)



loadSets <- function(file_dir) {
  set_list <- list()
  files <- list.files(file_dir)
  print(files)
  for (f in files) {
    set_list[[tools::file_path_sans_ext(f)]] <- read.csv(file.path(file_dir, f), header = F)[, 1]
  }
  set_listㄛ
}

drawVenn <- function(file_dir, out_file_name) {
  set_list <- loadSets(file_dir)
  myCol <- brewer.pal(length(set_list), "Pastel2")
  
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

drawVenn("../../OneDrive/桌面/Research/DOX/degs/ipsc_ups", 
         "../../OneDrive/桌面/Research/DOX/results/venn/ipsc_up_venn.png")
drawVenn("../../OneDrive/桌面/Research/DOX/degs/ipsc_downs", 
         "../../OneDrive/桌面/Research/DOX/results/venn/ipsc_down_venn.png")

drawVenn("../../OneDrive/桌面/Research/DOX/degs/cell_type_ups", 
         "../../OneDrive/桌面/Research/DOX/results/venn/cell_type_ups_venn.png")

drawVenn("../../OneDrive/桌面/Research/DOX/degs/cell_type_downs", 
         "../../OneDrive/桌面/Research/DOX/results/venn/cell_type_downs_venn.png")



drawVenn("../../OneDrive/桌面/Research/DOX/degs/dosage_ups", 
         "../../OneDrive/桌面/Research/DOX/results/venn/dosage_ups_venn.png")

drawVenn("../../OneDrive/桌面/Research/DOX/degs/dosage_downs", 
         "../../OneDrive/桌面/Research/DOX/results/venn/dosage_downs_venn.png")