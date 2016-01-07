#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")

parser$add_argument("-q", "--quietly", action="store_false",
                    dest="verbose", help="Print little output")

parser$add_argument("-p", "--points", type="integer", default=5,
                    help="minimum number of iRT peptides to make regression [default %(default)s]",
                    metavar="number")

parser$add_argument("-s", "--suffix", type="character", default="_irt_corrected",
                    help="Suffix to add before ending [default %(default)s]")

parser$add_argument("irt_table", nargs=1, 
                    help="iRT table with peptide\\tiRT")

parser$add_argument("input_file", nargs=1, 
                    help="File/pattern library to correct iRT times in, must be .csv ending")

parser$add_argument("output_dir", nargs=1, 
                    help="Output directory")



args <- parser$parse_args()

input_file <- args$input_file
output_dir <- args$output_dir
suffix <- args$suffix
points <- args$points
irt_table <- args$irt_table

#input_file = "../../results/2016-01-04/spectrast_results/umpire_34_single_irt_cons_openswath.csv"
#irt_table = "~/projects/microSWATH_library/results/2015-12-05/iRT_microSWATH.txt"

#output_dir = "./results/2016-01-04/final_libraries"
if (output_dir != "." || output_dir != "./") {
  dir.create(output_dir, recursive = T)
}
  
irt.dataset <- read.table(irt_table, quote="\"", comment.char="")
names(irt.dataset) = c("Peptide.Sequence", "iRT")

#fixing irt in the libraries
# irt.dataset  = normalized %>% filter(toRemove != 1 ) %>% group_by(Peptide.Sequence) %>% 
#   summarise(meanIRT = round(mean(iRT, na.rm=T)*100,4), sd = sd(iRT), meanRT = mean(apex, na.rm=T)*60) %>% 
#   arrange(meanIRT) %>% dplyr::select(Peptide.Sequence, meanIRT, meanRT)

#input_path = "./results/2016-01-04/spectrast_results"
#output_dir = "./results/2016-01-04/final_libraries"

filesToProcess = dir(path=dirname(input_file), pattern = basename(input_file), recursive=F, full.names = T)

for (i in 1:length(filesToProcess)) {
  input_file = filesToProcess[i]
  lib.raw <- read.delim(input_file)
  tmp.dataset = data.frame(rt = lib.raw$Tr_recalibrated)
  tmp = cbind(irt.dataset, data.frame (rt = lib.raw[match(irt.dataset$Peptide.Sequence, lib.raw$PeptideSequence), "Tr_recalibrated"]))
  stopifnot(nrow(na.omit(tmp)) > points)
  
  lm.fit = lm(iRT~rt, data = na.omit(tmp))
  message(paste("correcting", input_file))
  summary = summary(lm.fit)
  if (summary$r.squared < 0.9) {
    message(paste("WARNING: cannot correct", input_file))
    next
  }
  
  lib.raw$Tr_recalibrated = predict(lm.fit, tmp.dataset)
  output_file = sub(pattern = ".csv", x = basename(input_file), replacement = paste(suffix, ".csv", sep = ""))
  tmp.file_path = paste(output_dir, output_file, sep = "/")
  write.table(x = lib.raw, file = tmp.file_path, row.names = F, col.names = T, quote = F, sep = "\t")
}









