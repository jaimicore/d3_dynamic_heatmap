#####################
## Load R packages ##
#####################
required.libraries <- c("amap",         ## Distance calculation
                        "RColorBrewer"  ## Color palette
                        )
for (lib in required.libraries) {
  if (!require(lib, character.only=TRUE)) {
    install.packages(lib)
    suppressPackageStartupMessages(library(lib, character.only=TRUE))
  }
}


##################################################################################################
## Read arguments from command line
##
## Some variables are mandatory: If they are not declared from the command line the program ends
###################################################################################################
message("; Reading arguments from command line")
args <- commandArgs(trailingOnly=TRUE)
if (length(args >= 1)) {
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}


if (!exists("input.tab.file")) {
  stop("Missing mandatory argument (Input table): input.tab.file ")
  
} 

if (!exists("page.title")) {
  page.title <- "D3 dynamic heatmap"
  
} 

if (!exists("d3.base.file")) {
  d3.base.file <- file.path("../HTML/D3/d3.min.js")
  
} 

if (!exists("results.dir")) {
  results.dir <- file.path("results", "D3_heatmap")
  
} 

if (!exists("heatmap.template.file")) {
  heatmap.template.file <- file.path("HTML/Template/dynamic_heatmap_d3.html")
  
} 

if (!exists("color.classes")) {
  color.classes <- 9
  
} 

if (!exists("color.palette")) {
  color.palette <- "YlOrRd"
  
}
color.classes <- as.numeric(color.classes)


###########
## Debug ##
###########
# input.tab.file <- "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/D3_heatmap/Demo_files/AUC_NWDsignificantScore_heatmap_compare.txt"
# input.tab.file <- "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/D3_heatmap/Demo_files/RNAseq_example.tab"
# results.dir <- "/home/jamondra/Documents/PostDoc/Mathelier_lab/Projects/D3_heatmap/demo_results"


###########################
## Create results folder ##
###########################
message("; Creating result folders")
folders.to.create <- list()
folders.to.create[["tables"]] <- file.path(results.dir, "tables")
folders.to.create[["html"]] <- file.path(results.dir, "html")
thrash <- sapply(folders.to.create, dir.create, recursive = TRUE, showWarnings = FALSE)


######################
## Read input table ##
######################
input.tab <- read.table(input.tab.file, sep = "\t", header = TRUE)

## Check if the matrix have column and row names
## Otherwise, assign generic names
## Columns
if( is.null(colnames(input.tab)) ){
  
  generic.col.names <- paste("column",
                             1:ncol(input.tab),
                             sep = "_")
  colnames(input.tab) <- generic.col.names
}

## Rows
if( is.null(rownames(input.tab)) ){
  
  generic.row.names <- paste("row",
                             1:nrow(input.tab),
                             sep = "_")
  rownames(input.tab) <- generic.row.names
}


## For testing, this should be adapted
input.tab <- abs(input.tab)


########################
## Cluster the matrix ##
########################
## Run the hierarchical clustering with four methods (average + complete + single + ward). 
## Save the order of the nodes
comp.order.list.cols <- list()
comp.order.list.rows <- list()

thrash <- sapply(c("average", "complete", "single", "ward.D"), function(m){
  
  message("; clustering input table - Method ", m)
  
  ######################
  ## Cluster the rows ##
  ######################
  row.dist.tab <- Dist(input.tab, method = 'pearson')
  row.order <- hclust(row.dist.tab, method = m)[[3]]
  
  #########################
  ## Cluster the columns ##
  #########################
  input.tab.t <- t(input.tab)
  col.dist.tab <- Dist(input.tab.t, method = 'pearson')
  col.order <- hclust(col.dist.tab, method = m)[[3]]
  
  
  ## Export the order in a  list
  comp.order.list.rows[[m]] <<- paste(row.order, collapse = ",")
  comp.order.list.cols[[m]] <<- paste(col.order, collapse = ",")

})
## Add default values
comp.order.list.rows[["default"]] <- paste(1:nrow(input.tab), collapse = ",")
comp.order.list.cols[["default"]] <- paste(1:ncol(input.tab), collapse = ",")


###############################################
## Convert and export the input table in TSV ##
###############################################
## NOTE: TSV is the format read by D3

tsv.tab <- NULL
for(j in 1:dim(input.tab)[1]){
  for(i in 1:dim(input.tab)[2]){
    tsv.tab <<- rbind(tsv.tab, matrix(c(j,i, as.numeric(input.tab[j,i])), nrow = 1))
  }
}
colnames(tsv.tab) <- c("Row", "Col", "Value")
tsv.tab <- as.data.frame(tsv.tab)
tsv.tab$Value <- round(tsv.tab$Value, digits = 5)
tsv.tab.file <- file.path(folders.to.create[["tables"]], "D3_dynamic_heatmap.tsv")
message("; Exporting input table in TSV format: ", tsv.tab.file)
write.table(tsv.tab, file = tsv.tab.file, sep = "\t", quote = FALSE, row.names = FALSE)


###################################
## Legend, colors and dimensions ##
###################################
message("; Generating color palette")
## Color palette (from RcolorBrewer)
rgb.palette <- colorRampPalette(brewer.pal(color.classes, color.palette), space="Lab")
white <- "#FFFFFF"
## Ten color classes: white + 9 colors from the palette
color.codes <- append(white,rgb.palette(9))
color.codes.rev <- rev(append(white,rgb.palette(9)))
color.codes <- paste("'", color.codes, "'",collapse=",")
color.codes.rev <- paste("'", color.codes.rev, "'",collapse=",")


## Domain numbers
domain <- seq(from = min(input.tab),
              to = max(input.tab),
              length.out = 9
              )
domain <- round(domain, digits = 5)
domain <- paste(domain, collapse=",")


## Legend numbers
## NOTE: Data legend is longer (+1) then the Domain
legend <- rev(seq(from = min(input.tab),
              to = max(input.tab),
              length.out = 10))

legend <- round(legend, digits = 5)
legend <- paste(legend, collapse=",")


## Heatmap width
cell.size <- 30
heatmap.width <- 350 + (cell.size * ncol(input.tab))

## Margin
labels.tab <- c(colnames(input.tab), rownames(input.tab))
margin.heatmap <- (max(as.vector(sapply(labels.tab, nchar))) + 2) * 10

## Col and Row labels
col.lab <- paste(paste("'", colnames(input.tab), "'", sep = ""), collapse = ",")
row.lab <- paste(paste("'", rownames(input.tab), "'", sep = ""), collapse = ",")


#############################################################
## Fill the HTML template                                  ##
## Substitute the words marked in the template by the data ##
#############################################################
message("; Filling the template")
## Load the template
html.template <- readLines(heatmap.template.file)

## Insert title
html.template <- gsub(html.template, pattern = "--Page_title--", replacement = page.title)

## Insert D3 base
## NOTE: it must be the relative path to the html document
html.template <- gsub(html.template, pattern = "--d3--", replacement = d3.base.file)

## Column and row number
html.template <- gsub(html.template, pattern = "--c_numb--", replacement = ncol(input.tab))
html.template <- gsub(html.template, pattern = "--r_numb--", replacement = nrow(input.tab))

## Column order (default + clustered)
html.template <- gsub(html.template, pattern = "--col_default_nb--", replacement = comp.order.list.cols[["default"]])
html.template <- gsub(html.template, pattern = "--comp_average_c_number--", replacement = comp.order.list.cols[["average"]])
html.template <- gsub(html.template, pattern = "--comp_complete_c_number--", replacement = comp.order.list.cols[["complete"]])
html.template <- gsub(html.template, pattern = "--comp_single_c_number--", replacement = comp.order.list.cols[["single"]])
html.template <- gsub(html.template, pattern = "--comp_ward_c_number--", replacement = comp.order.list.cols[["ward.D"]])

## Row order (default + clustered)
html.template <- gsub(html.template, pattern = "--row_default_nb--", replacement = comp.order.list.rows[["default"]])
html.template <- gsub(html.template, pattern = "--comp_average_r_number--", replacement = comp.order.list.rows[["average"]])
html.template <- gsub(html.template, pattern = "--comp_complete_r_number--", replacement = comp.order.list.rows[["complete"]])
html.template <- gsub(html.template, pattern = "--comp_single_r_number--", replacement = comp.order.list.rows[["single"]])
html.template <- gsub(html.template, pattern = "--comp_ward_r_number--", replacement = comp.order.list.rows[["ward.D"]])

## Column and row labels
html.template <- gsub(html.template, pattern = "--col_labels--", replacement = col.lab)
html.template <- gsub(html.template, pattern = "--row_labels--", replacement = row.lab)

## TSV path file
tsv.relative.path <- file.path("tables", basename(tsv.tab.file))
html.template <- gsub(html.template, pattern = "--file--", replacement = tsv.relative.path)

## Color hexa codes
html.template <- gsub(html.template, pattern = "--color_scale--", replacement = color.codes)
html.template <- gsub(html.template, pattern = "--color_scale_rev--", replacement = color.codes.rev)

## Insert domain numbers
html.template <- gsub(html.template, pattern = "--domain--", replacement = domain)

## Insert legend numbers
html.template <- gsub(html.template, pattern = "--data_legend--", replacement = legend)

## Set heatmap width
html.template <- gsub(html.template, pattern = "--heatmap_width--", replacement = heatmap.width)

## Set heatmap margins
html.template <- gsub(html.template, pattern = "--left--", replacement = margin.heatmap)

## Set cell size
html.template <- gsub(html.template, pattern = "--cell_size--", replacement = cell.size)


##############################
## Export HTML with heatmap ##
##############################
message("; Exporting the heatmap")
heatmap.html.file <- file.path(results.dir, "D3_dynamic_heatmap.html")
write(html.template, file = heatmap.html.file)
