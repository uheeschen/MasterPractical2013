#######################################################################
#File:        evalAlign.R
#Description: perform pairwise comparison of sequence search methods
#######################################################################

library(ggplot2)

library(RSQLite)

loadHHBlitsTSV <- function(file) {
  data <- read.csv(file=file, sep='\t', header=TRUE)
  data$eval <- data$e_value
  return(data)
}

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

#get pdb hits from sequence search results
getPDBEntries <- function(data)
{
  d <- data[grep('pdb\\|pdb', data$hit),]
  d$chain <- gsub(pattern='^.*(\\w)$', replacement='\\1', d$hit)
  d$pdb <- gsub(pattern='^.*(\\w{4})_\\w$', replacement='\\1', d$hit)
  d$pdb_chain <- paste(d$pdb, d$chain, sep='_')
  return(d)
}

#get uniprot hits from sequence search results
getUniProtEntries <- function(data)
{
  d <- data[grep("^tr|sp\\|", data$hit), ]
  d$uniprot <- gsub(pattern="^(tr|sp)\\|(\\w*)\\|.*$", replacement="\\2", d$hit)
  return(d)
}

#load parsable file: dir.cla.scop.txt
loadSCOPClassification <- function(file)
{
  data <- read.csv(file=file, sep='\t', comment.char='#')
  names(data) <- c('domain',
                   'pdb',
                   'chain',
                   'sccs',
                   'sunid',
                   'fullstr')
  data$pdb_chain <- paste(data$pdb, substr(data$chain,1,1), sep='_')
  return(data)
}

getSCOPFold <- function(scopData)
{
  gsub(pattern="^(\\w\\.\\d*)\\..*$", replacement="\\1", scopData$sccs)
}

getUniprotGO <- function(data, database.path)
{
  conn <- dbConnect(dbDriver('SQLite'), dbname=database.path)
  uniprot_acs <- unique(data$uniprot)
  #create temperory table
  dbGetQuery(conn, "CREATE TEMPORARY TABLE temp_uniprot_ac(uniprot_ac CHAR(6));")
  dbSendPreparedQuery(
    conn, 
    "INSERT INTO temp_uniprot_ac(uniprot_ac) VALUES(:uniprot_ac)",
    bind.data=data.frame(uniprot_ac=data$uniprot))
  result <- dbGetQuery(
    conn, 
    "SELECT count(t.uniprot_ac) as count, go_term as go
    FROM temp_uniprot_ac t,uniprot_go_mapping m 
    WHERE t.uniprot_ac = m.uniprot_ac GROUP BY go_term ORDER BY count")
  dbSendQuery(conn, "DROP TABLE temp_uniprot_ac;")
  dbDisconnect(conn)
  return(result)
}

###
# configure test cases for reference sequences
###
querys <- c('BCKDHA', 'BCKDHB', 'DBT', 'DLD')
data.path <- 
  '/home/wei/git/MasterPractical2013/results/task02/01-seq-search/search-results'
color.palette <- 'Set2'
database.path <- '/home/wei/idmapping.sqlite3'
image.type <- '.png'

#load SCOP classification
scop <- loadSCOPClassification(
  file='/home/wei/git/MasterPractical2013/data/SCOP/dir.cla.scop.txt_1.75')

for (query in querys) {
  #configure input files
  input.filenames <- c(paste('blastp_', query, '.tsv', sep=''),
      paste('psiblast-2_iterations_eval_0.002-refseq_', query, '_protein.fasta.tsv', sep=''),
      paste('psiblast-2_iterations_eval_10e-10-refseq_', query, '_protein.fasta.tsv', sep=''),
      paste('psiblast-10_iterations_eval_0.002-refseq_', query, '_protein.fasta.tsv', sep=''),
      paste('psiblast-10_iterations_eval_10e-10-refseq_', query, '_protein.fasta.tsv', sep=''),
      paste('hhblits_refseq_', query, '_protein.fasta.hhr.tsv', sep=''))
  data.names <- c('blast', 
                  'psiblast(iter. 2, e-val. 0.002)', 
                  'psiblast(iter. 2, e-val. 10e-10)',
                  'psiblast(iter. 10, e-val. 0.002)',
                  'psiblast(iter. 10, e-val. 10e-10)',
                  'hhblits')
  data.titles <- c('Blast',
                   'PSI-Blast[iter. 2, eval. 0.002]',
                   'PSI-Blast[iter. 2, eval. 10e-10]',
                   'PSI-Blast[iter. 10, eval. 0.002]',
                   'PSI-Blast[iter. 10, eval. 10e-10]',
                   'HHBlits')
  
  #load data
  data <- list()
  evals <- c()
  methods <- c()
  identities <- c()
  for (index in 1:length(input.filenames))
  {
    data.name <- data.names[index]
    print(data.name)
    input.path <- file.path(data.path, input.filenames[index])
    if (data.name != 'hhblits') 
    {
      frame <- loadBlastTSV(input.path)
    } else {
      frame <- loadHHBlitsTSV(input.path)
    }
    data[[ data.name ]] <- frame
    n <- length(frame$eval)
    evals <- c(evals, frame$eval)
    methods <- c(methods, rep(x=data.name, n))
    identities <- c(identities, frame$identity)
  }
  DATA <- data.frame(evalue=evals, identity=identities, method=methods)

###
# evaluation: e-value distribution
###
  PLOT <- ggplot(DATA, aes(x=evalue))
  PLOT + geom_density(aes(colour=factor(method), fill=factor(method)), alpha=.7) +
    scale_x_log10() + scale_alpha(range=c(0, 1)) +
    ggtitle(paste('E-value distribtution (', query, ')', sep='')) +
    xlab('E-value') + ylab('Density')+
    scale_colour_brewer(palette=color.palette) +
    scale_fill_brewer(palette=color.palette)
  ggsave(paste("e-value-distribution_", query, image.type, sep=''), width=8.3, height=6.8, dpi=100)
###
# evaluation: identity distribution
###  
  PLOT <- ggplot(DATA, aes(x=identity))
  PLOT + geom_density(aes(colour=factor(method), fill=factor(method)), alpha=.7) +
    scale_alpha(range=c(0, 1)) +
    ggtitle(paste('Identity distribtution (', query, ')', sep='')) +
    xlab('Identity') + ylab('Density')+
    scale_colour_brewer(palette=color.palette) +
    scale_fill_brewer(palette=color.palette)
  ggsave(paste("identity_distribution_", query,image.type, sep=''), width=8.3, height=6.8, dpi=100)

###
# evaluation: intersection curve
###
  thresholds <- c(0, 1e-100, 1e-90, 1e-80, 1e-70, 
             1e-60, 1e-50, 1e-40, 1e-30, 1e-20, 
             1e-10, 1, 10)
  for (index1 in 1:length(data.names))
  {
    dataset1 <- data[[index1]]
    hits1 <- dataset1$hit
    meth.name <- data.names[index1]
    threshold <- c()
    quality <- c()
    method <- c()
    
    for (index2 in 1:length(data.names))
    {
      if (index2 == index1) {
        next();
      }
      dataset2 <- data[[index2]]
      hits2 <- dataset2$hit
      curmeth <- data.names[index2]
      for (thr in thresholds)
      {
        intersection <- length( 
          intersect(
            hits1[ dataset1$eval <= thr ],
            hits2[ dataset2$eval <= thr ]) )
        avg.intersect <- 0.5 * (
          intersection / length(hits1) + 
          intersection / length(hits2) )
        
        threshold <- c(threshold, thr)
        method <- c(method, curmeth)
        quality <- c(quality, avg.intersect)
      }
    }
    
    performance <- data.frame(threshold, quality, method)
    PLOT <- ggplot(performance, aes(x=threshold, y=quality, colour=method))
    PLOT + geom_line(size=1) + geom_point() + scale_x_log10() +
      ggtitle( 
        paste('Relative intersections (comparison to ',meth.name,')', sep='') ) +
      xlab('E-value') + ylab('intersecting results') +
      scale_colour_brewer(palette=color.palette) + ylim(0,1)
    ggsave( paste("intersection_to_", meth.name,"_", query,image.type, sep=''), width=8.3, height=6.8, dpi=100)
  }
  
###
# evaluation: protein structure classification
###
  for (methId in 1:length(data.names))
  {
    method <- data.names[methId]
    dataset <- getPDBEntries(data[[ method ]])
    subsetSCOP <- scop[ scop$pdb_chain %in% dataset$pdb_chain, ]
    d <- data.frame(folds=getSCOPFold(subsetSCOP))
    if (length(d$folds) == 0)
    {
      next();
    }
    d$freq <- rep(0, length(d$folds))
    ftable <- table(d$folds)
    for (fold in names(ftable))
    {
      d$freq[ d$folds == fold ] <- ftable[ fold ]
    }
    PLOT <- ggplot(d, aes(x=reorder(folds, freq)))
    PLOT + geom_histogram() + 
      ggtitle("Histogram of fold classes from annotated pdb hits") +
      xlab('fold class')
    ggsave( paste("SCOP_histogram_", method, "_", query, image.type, sep=''), width=8.3, height=6.8, dpi=100 );
  }
  
###
# evaluation: Gene Ontology
###
#   common_gos <- list()
  for (methId in 1:length(data.names))
  {
    method <- data.names[methId]
    dataset <- getUniProtEntries(data[[ method ]])
    if (length(dataset$uniprot) == 0)
    {
      next();
    }
    result <- getUniprotGO(dataset, database.path=database.path)
    total.count <- length(dataset$uniprot)
    #result$percent <- result$count / length(dataset$uniprot)
    result$percent <- result$count / total.count
    #result <- result[ result$count > 0.05 * total.count, ]
    result <- result[ order(result$count, decreasing=TRUE)[1:5], ]
    PLOT <- ggplot(result, aes(x=reorder(go, count), y=percent))
    PLOT + geom_bar(stat='identity') + coord_flip() + 
      ggtitle("Histogram of top-5 go terms from annotated uniprot hits") +
      xlab('GO term') + ylab('Frequency') + ylim(0,1)
    ggsave( paste("GO_histogram_", method, "_", query, ".png", sep=''), width=8.3, height=6.8, dpi=100 );
    
#     go_terms <- data.frame(threshold = thresholds, count=rep(0, length(thresholds)))
#     for (thr in thresholds)
#     {
#       subset <- dataset[ dataset$eval <= thr, ]
#       result <- getUniprotGO(dataset, database.path=database.path)
#       go_terms$count[ go_terms$threshold == thr ] <- length(unique(result$go))
#     }
#     common_gos[method] <- go_terms
  }
  
  save.image(file=paste('Workspace_', query, '.RData', sep=''), compress="xz")
}
