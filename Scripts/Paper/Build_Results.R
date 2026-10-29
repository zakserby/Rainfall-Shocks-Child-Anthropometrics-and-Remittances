#===============================================================================
# MS Thesis
# Results document
#===============================================================================

# This script puts every table and figure in the output folder into one LaTeX file and also a PDF file. We use TinyTex. Note that appendix and miscellaneous figures are not included
# Must be run after the summary stats, the regressions and regression stats do files and r scripts, since it only reads their output. 00_Master.do runs it last

#===============================================================================
# 1). Preliminary
#===============================================================================

rm(list=ls())
# Close a log left open by an earlier run that stopped partway, same as Contdid_regs.R
while (sink.number() > 0) sink()

# Load packages. tinytex compiles the LaTeX file
options(repos = c(CRAN = "https://cloud.r-project.org"))
want <- c("tinytex")
need <- want[!(want %in% installed.packages()[,"Package"])]
if (length(need)) install.packages(need)
lapply(want, function(i) require(i, character.only=TRUE))
rm(want, need)

# Directories, same as the other R scripts: the project root is one folder up from Scripts
dir <- list()
dir$root <- dirname(getwd())
# If the project root is not one folder up from the working folder, go up from the working folder until the folder that has Scripts/00_Master.do
if (!dir.exists(paste(dir$root,"/Data",sep=""))) {
  dir$root <- getwd()
  while (!file.exists(paste(dir$root,"/Scripts/00_Master.do",sep="")) && dirname(dir$root) != dir$root) dir$root <- dirname(dir$root)
  if (!file.exists(paste(dir$root,"/Scripts/00_Master.do",sep=""))) stop("Project folder not found. Set the working directory to the Scripts folder (open MS Thesis.Rproj) or run 00_Master.do")
}
dir$paper <- paste(dir$root,"/Paper",sep="")
dir$tables <- paste(dir$root,"/Output/Tables",sep="")
dir$figures <- paste(dir$root,"/Output/Figures",sep="")
dir$output_logs <- paste(dir$root,"/Output/Logs",sep="")
dir.create(dir$paper, recursive = TRUE, showWarnings = FALSE)
for (s in c("Main","Miscellaneous")) {
  dir.create(paste(dir$output_logs,"/",s,sep=""), recursive = TRUE, showWarnings = FALSE)
}

# Everything printed below is also saved to this log
sink(paste(dir$output_logs,"/Main/build_results.log",sep=""), split = TRUE)

# TinyTeX is installed once if there is no LaTeX on the computer. The usual download link was down in October 2026, so if it fails we take the same
# bundle from TinyTeX's GitHub releases and unpack it where tinytex looks for it
if (!tinytex::is_tinytex() && Sys.which("pdflatex") == "") {
  ok <- tryCatch({tinytex::install_tinytex(); TRUE}, error = function(e) FALSE)
  if (!ok) {
    if (!("jsonlite" %in% installed.packages()[,"Package"])) install.packages("jsonlite")
    os <- ifelse(Sys.info()[["sysname"]] == "Darwin", "darwin", "linux")
    rel <- jsonlite::fromJSON("https://api.github.com/repos/rstudio/tinytex-releases/releases/latest")
    url <- grep(paste0("TinyTeX-1-", os, "-.*\\.tar\\.(xz|gz)$"), rel$assets$browser_download_url, value = TRUE)[1]
    file <- file.path(tempdir(), basename(url))
    download.file(url, file, mode = "wb")
    home <- ifelse(os == "darwin", "~/Library", "~")
    untar(file, exdir = path.expand(home))
  }
}
# Put the TinyTeX programs on the path for this session, so pdflatex is found
if (tinytex::is_tinytex()) {
  bin <- list.dirs(file.path(tinytex::tinytex_root(), "bin"), recursive = FALSE)[1]
  Sys.setenv(PATH = paste(bin, Sys.getenv("PATH"), sep = .Platform$path.sep))
}

#===============================================================================
# 2). Tables
#===============================================================================

# Each table file is its own LaTeX document, so we keep only what is between \begin{document} and \end{document}
# Meaning this function keeps only the middle of each file: the tables themselves, without the setup and ending
table_body <- function(file) {
  lines <- readLines(file, warn = FALSE)
  # This finds which line number holds \begin{document} and similar the end one with \end{document}
  start <- grep("^\\\\begin\\{document\\}", lines)
  end <- grep("^\\\\end\\{document\\}", lines)
  lines[(start + 1):(end - 1)]
}

#===============================================================================
# 3). Figures
#===============================================================================

# NOTE: The main purpose of this is to write the title above each figure in Results.pdf, built from the image’s file name
# First defining so lookups that will be used in the function
outcome <- c(haz = "Height-for-Age", waz = "Weight-for-Age", bmiz = "BMI-for-Age")
group <- c(boys = "Boys", girls = "Girls", rural = "Rural Children", urban = "Urban Children", all3 = "Children in All Three Waves")
# Defining the function
figure_title <- function(name) {
  # Breaking apart the file name
  parts <- strsplit(sub("\\.png$", "", name), "_")[[1]]
  # Finding the outcome and group and giving them the appropropriate names
  v <- outcome[parts[parts %in% names(outcome)]]
  g <- group[parts[parts %in% names(group)]]
  # This finds the figure type, whatever’s left after removing the outcome and subgroup pieces (i.e. event study)
  base <- paste(parts[!(parts %in% c(names(outcome), names(group)))], collapse = "_")
  # This is basicaly constructing the title
  family <- gsub("_", " ", base)
  if (base == "trends") family <- paste0("Mean ", v, " by Wave, Treated and Control, Binary DID")
  if (base == "event_study") family <- paste0("Event Study, ", v, ", Binary DID")
  if (base == "event_study_ddd") family <- paste0("Triple Difference Event Study, ", v, ", Binary DID")
  if (base == "contdid_binscatter") family <- paste0("Change in ", v, " by 2009 Rainfall Deficit, Continuous DID")
  if (base == "contdid_att_clustered") family <- paste0("ATT by 2009 Rainfall Deficit, ", v, ", Continuous DID")
  if (base == "contdid_acrt_clustered") family <- paste0("ACRT by 2009 Rainfall Deficit, ", v, ", Continuous DID")
  if (base == "contdid_spi_att_clustered") family <- paste0("ATT by SPI Drought Dose, ", v, ", Continuous DID")
  if (base == "drought_map_spi") family <- "2009 Rainfall by Municipality, SPI Classes"
  if (base == "drought_map_binary") family <- "2009 Drought by Municipality, Binary DID Treatment"
  # This joins the type title and the subgroup with a comma
  paste(c(family, g), collapse = ", ")
}

# Creating one figure, the title on top and the figure below it, in that order
# The function takes two inputs: the figures subfolder, e.g. "Main" and the image file, e.g. "event_study_waz.png"
figure_block <- function(sub, name) {
  c("\\begin{figure}[H]", "\\centering",
    paste0("\\caption{", figure_title(name), "}"),
    paste0("\\includegraphics[width=0.8\\textwidth]{", name, "}"), "\\end{figure}", "")
}

# This function decides the order of the figures in Results.pdf. This function sorts them by:
# 1) Family aka figure type: drought maps first, then event studies, triple difference event studies, trends, and the continuous DID figures
# 2) Outcome within each family: HAZ, WAZ, BMIZ
# 3) Subgroup: the plain version before subgroup versions
figure_order <- function(files) {
  fam <- c("event_study_", "trends_", "contdid_binscatter", "contdid_att_clustered", "contdid_acrt_clustered", "contdid_spi_att_clustered")
  # Each file gets the position of its family in fam: startsWith checks which family its name begins with, which()[1] takes that family's position, and sapply repeats this for every file, e.g. event_study_waz.png is 1
  f_rank <- sapply(files, function(x) which(startsWith(x, fam))[1])
  f_rank[startsWith(files, "event_study_ddd")] <- 1.5
  f_rank[startsWith(files, "drought_map_")] <- 0
  v_rank <- sapply(files, function(x) which(c(grepl("_haz", x), grepl("_waz", x), grepl("_bmiz", x)))[1])
  # order() sorts by family, then outcome, then name length, then alphabetically, and files[...] puts the names in that order
  files[order(f_rank, v_rank, nchar(files), files)]
}
# list.files lists the .png files in Output/Figures/Main, and figure_order sorts them into the order they appear in Results.pdf
main_figs <- figure_order(list.files(paste(dir$figures,"/Main",sep=""), pattern = "\\.png$"))

#===============================================================================
# 4). Write Results.tex
#===============================================================================

# This builds the whole Results.tex file as one list of text lines and saves it
# It’s the step that puts everything together by implementing the document setup, adds a title page, inserts the middles of the three table files and the Main figures, and closes the document
# The summary stats, the regression results, and the appendix coded with table_body are simply section headings
tex <- c("\\documentclass[12pt]{article}", "\\usepackage[margin=1in]{geometry}", "\\usepackage{booktabs}", "\\usepackage{longtable}",
         "\\usepackage{caption}", "\\usepackage{graphicx}", "\\usepackage{float}",
         "\\DeclareCaptionLabelSeparator{emdash}{---}", "\\captionsetup{font=normalsize,labelfont=sc,textfont=sc,labelsep=emdash}", "\\captionsetup[figure]{position=top,justification=centering}",
         "\\def\\sym#1{\\ifmmode^{#1}\\else\\(^{#1}\\)\\fi}",
         "\\graphicspath{{../Output/Figures/Main/}{Output/Figures/Main/}{./}}",
         "\\title{Results First Draft}", "\\author{Zachary Serby}", "\\date{\\today}",
         "\\begin{document}", "\\maketitle", "\\thispagestyle{empty}", "\\clearpage", "",
         "\\subsection*{Summary Statistics}", "",
         table_body(paste(dir$tables,"/Main/summary_stats_did_table.tex",sep="")), "",
         "\\clearpage", "\\subsection*{Regression Results}", "",
         table_body(paste(dir$tables,"/Main/regression_tables_main.tex",sep="")), "",
         "\\clearpage", unlist(lapply(main_figs, function(f) figure_block("Main", f))),
         "\\clearpage", "\\section*{Appendix}", "",
         table_body(paste(dir$tables,"/Appendix/regression_tables_appendix.tex",sep="")), "",
         "\\end{document}")
# Writes every line in tex and saves it 
writeLines(tex, paste(dir$paper,"/Results.tex",sep=""))
cat("Results.tex written with", length(main_figs), "main figures\n")

#===============================================================================
# 5). Compile Results.pdf
#===============================================================================

# This chunk turns Results.tex into Results.pdf, and cleans up a bit
setwd(dir$paper)
# Compiles PDF
tinytex::latexmk("Results.tex", engine = "pdflatex", clean = FALSE)
if (file.exists("Results.log")) file.rename("Results.log", paste(dir$output_logs,"/Miscellaneous/Results_latex.log",sep=""))
for (f in c("Results.aux", "Results.out", "Results.fls", "Results.fdb_latexmk", "Results.toc")) {
  if (file.exists(f)) file.remove(f)
}
cat("Results.pdf written to", dir$paper, "\n")
sink()

# The end
