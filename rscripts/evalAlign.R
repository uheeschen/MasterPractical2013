#######################################################################
#File:        evalAlign.R
#Description: perform pairwise comparison of sequence search methods
#######################################################################

library(ggplot2)

loadBlastTSV <- function(file) {
  data <- read.csv(file=file, sep='\t', header=FALSE, comment.char='#')
  names(data) <- c('query', 
                       'hit', 
                       'identity', 
                       'algn_len', 
                       'mismatch', 
                       'gap_open',
                       'q_start',
                       'q_end',
                       's_start', 
                       's_end', 
                       'eval',
                       'score')
  return(data)
}

###
# configure test case
###
input.filenames <- c('blastp_BCKDHA.tsv', 
                 'psiblast-2_iterations_eval_0.002-refseq_BCKDHA_protein.fasta.tsv',
                 'psiblast-2_iterations_eval_10e-10-refseq_BCKDHA_protein.fasta.tsv',
                 'psiblast-10_iterations_eval_0.002-refseq_BCKDHA_protein.fasta.tsv',
                 'psiblast-10_iterations_eval_10e-10-refseq_BCKDHA_protein.fasta.tsv')
data.titles <- c('Blast of BCKDHA',
                 'PSI-Blast of BCKDHA [iter. 2, eval. 0.002]',
                 'PSI-Blast of BCKDHA [iter. 2, eval. 10e-10]',
                 'PSI-Blast of BCKDHA [iter. 10, eval. 0.002]',
                 'PSI-Blast of BCKDHA [iter. 10, eval. 10e-10]')
input.files <- file.path(
  '/home/wei/git/MasterPractical2013/results/task02/01-seq-search/search-results',
  #'/mnt/home/student/weish/master-practical-2013/task02/01-seq-search/results',
  input.filenames)

###
# evaluation: proteion classification
###
for (f in 1:length(input.files)) {
  filename <- input.filenames[f]
  data <- loadBlastTSV(input.files[f])
  data <- ggplot(data, aes(x=eval), aes_string(x='E-value'))
  #plot distribution of E-values
  data + geom_histogram(binwidth=0.1) + scale_y_log10()
  ggsave(paste(
    'evalues-distribution', 
    filename, 
    '.pdf', 
    sep='_'))
}