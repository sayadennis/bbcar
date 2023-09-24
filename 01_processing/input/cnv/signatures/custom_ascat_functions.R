library(stringr)

#' Custom version of readAlleleCountFiles for alleles files
#'
#' The readAlleleCountFiles function needed to be customized since the original version did not work with the alleleCounter output that I was working with. 
#' @param allelecounts_fn An object of class characterstring. Path to the allele counts file to be read.
#' @param minCounts An object of class "int". Minimum threshold for allele counts. The rows that do not meet this minimum are removed.
#' @return Returns an object of class "data.frame". Allele counts file read as a data.frame.
#' @examples
#' tissue_input_data <- custom_readAlleleCountFiles(
#'     "/projects/b1131/saya/bbcar/data/02b_cnv/signatures/01_alleleCounts/1009_tissue_alleleFrequencies.txt",
#'     )
custom_readAlleleCountFiles = function(allelecounts_fn, minCounts) {
    if (file.exists(allelecounts_fn)) {
        data = data.frame(data.table::fread(allelecounts_fn, sep='\t', showProgress=F, header=T), stringsAsFactors=F)
        data = data[data[,7]>=minCounts,]
        rownames(data) = paste0(data[,1], '_', data[,2])
    } else {
        data=data.frame(data.table())
    }
    return(data)
}

#' Custom function to read the alleles file
#'
#' This function is customized to add options to add/remove the "chr" string from the chromosome column values. ("chr1" vs. "1")
#' @param prefix An object of class characterstring. Prefix to the allele count file, up to the chromosome number. 
#' @param suffix An object of class characterstring. Suffix to the allele count file i.e. what comes after the chromosome number. 
#' @param chrom_names An object of class list (?). List of chromosome names, like c(1:22, "XX")
#' @param add_chr_string An object of class boolean. Indicate preference to add "chr" to values of the chrom column.
#' @param remove_chr_string An object of class boolean. Indicate preference to remove "chr" from values of the chrom column.
#' @return Returns an object of class "data.frame". Allele counts file read as a data.frame.
#' @examples
#' # Add some code illustrating how to use the function
readAllelesFiles = function(prefix, suffix, chrom_names, add_chr_string=FALSE, remove_chr_string=FALSE) {
    files = paste0(prefix, chrom_names,suffix)
    files = files[sapply(files,function(x) file.exists(x) && file.info(x)$size>0)]
    stopifnot(length(files)>0)
    data=do.call(rbind,lapply(files,function(x) {
        tmp=data.frame(data.table::fread(x,sep='\t',showProgress=F,header=T))
        tmp=tmp[!is.na(tmp[,2] & !is.na(tmp[,3])),]
        tmp=tmp[!duplicated(tmp[,1]),]
        tmp$chromosome=gsub(paste0(prefix,'(',paste(chrom_names,collapse='|'),')',suffix),'\\1',x)
        if (add_chr_string) tmp$chromosome=paste0('chr',tmp$chromosome)
        if (remove_chr_string) tmp$chromosome=gsub('^chr','',tmp$chromosome)
        tmp=tmp[,c(4,1:3)]
        rownames(tmp)=paste0(tmp[,1],'_',tmp[,2])
        return(tmp)
  }))
  stopifnot(nrow(data)>0)
  return(data)
}

#' Custom function to generate and write the BAF and LogR files.
#'
#' This function is customized to work with my alleleCounts output, which seems a bit different than what ASCAT assumes. This could be due to the fact that I am using an alleleCounter from a container, since the ASCAT built-in could not be run on an HPC.
#' @param tissue_name An object of class characterstring. Name of the tissue e.g. "1009_tissue"
#' @param germline_name An object of class characterstring. Name of the germline e.g. "1009_germline"
#' @param dout An object of class characterstring. Output directory WITHOUT trailing "/"
#' @param tissueAlleleCountsFile.prefix An object of class characterstring.
#' @param normalAlleleCountsFile.prefix An object of class characterstring.
#' @param alleles.prefix An object of class characterstring. Prefix used in readAllelesFiles() call. 
#' @param gender An object of class characterstring. Expected values are one of c('XX','XY')
#' @param genomeVersion An object of class characterstring. Expected values are one of c('hg19','hg38')
#' @param chrom_names An object of class list (?). List of chromosome names, like c(1:22, "XX")
#' @param minCounts An object of class int. The minimum count threshold to set for reading allele counts files.
#' @param BED_file An object of class characterstring. The interval file for exome sequencing.
#' @param probloci_file An object of class characterstring. Problematic loci to be removed.
#' @return Returns None. Results are written to files. 
#' @examples
#' # Add some code illustrating how to use the function
custom_getBAFsAndLogRs = function(
        tissue_name, germline_name, dout, 
        tissueAlleleCountsFile.prefix, normalAlleleCountsFile.prefix, 
        alleles.prefix, gender, genomeVersion, 
        chrom_names=c(1:22), minCounts=20, BED_file=NA, probloci_file=NA, seed=as.integer(Sys.time())
    ) {
    set.seed(seed)
    stopifnot(gender %in% c('XX','XY'))
    stopifnot(genomeVersion %in% c('hg19','hg38'))

    samplename = str_split(tissue_name, '_')[[1]][1]

    if (
        file.exists(paste0(din, "/", tissue_name,"_alleleFrequencies.txt")) 
        & file.exists(paste0(din, "/", germline_name,"_alleleFrequencies.txt"))
    ) {
        tissueLogR_file=paste0(dout, "/tissue_normal/", tissue_name,"_tissueLogR.txt")
        tissueBAF_file=paste0(dout, "/tissue_normal/", tissue_name,"_tissueBAF.txt")
        normalLogR_file=paste0(dout, "/tissue_normal/", germline_name,"_germlineLogR.txt")
        normalBAF_file=paste0(dout, "/tissue_normal/", germline_name,"_germlineBAF.txt")
        # Load data, only keep SNPs with enough coverage
        tissue_input_data = custom_readAlleleCountFiles(paste0(tissueAlleleCountsFile.prefix, ".txt"), minCounts=1)
        normal_input_data = custom_readAlleleCountFiles(paste0(normalAlleleCountsFile.prefix, ".txt"), minCounts=1)
        allele_data = readAllelesFiles(alleles.prefix, ".txt", chrom_names)
        # Synchronise DFs
        rownames(allele_data) = paste("chr", rownames(allele_data), sep = "") # rownames are "1_10583" rather than "chr1_10583" by default
        allele_data$chromosome = paste0("chr", allele_data$chromosome) # chromosome column is by default "1" not "chr1" 
        matched_data = Reduce(intersect, list(rownames(tissue_input_data), rownames(normal_input_data), rownames(allele_data)))
        # tissue_input_data = tissue_input_data[rownames(tissue_input_data) %in% matched_data,]
        # normal_input_data = normal_input_data[rownames(normal_input_data) %in% matched_data,]
        # allele_data = allele_data[rownames(allele_data) %in% matched_data,]
        tissue_input_data=tissue_input_data[matched_data,]
        normal_input_data=normal_input_data[matched_data,]
        allele_data=allele_data[matched_data,]
        rm(matched_data)
        # If a probloci file is provided, remove those
        if (!is.na(probloci_file)) {
            stopifnot(file.exists(probloci_file) && file.info(probloci_file)$size>0)
            probloci=data.frame(data.table::fread(probloci_file,sep='\t',showProgress=F,header=T),stringsAsFactors=F)
            probloci=paste0('chr', probloci[,1],'_',probloci[,2]) # probloci=paste0(gsub('^chr','',probloci[,1]),'_',probloci[,2])
            probloci=which(rownames(tissue_input_data) %in% probloci)
            if (length(probloci)>0) {
                tissue_input_data = tissue_input_data[-probloci,]
                normal_input_data = normal_input_data[-probloci,]
                allele_data = allele_data[-probloci,]
            } else {
                warning('The probloci did not remove any SNPs, it might be worth checking the data.')
            }
            rm(probloci)
        }
        stopifnot(isTRUE(all.equal(allele_data[,1],tissue_input_data[,1]))
                  && isTRUE(all.equal(allele_data[,1],normal_input_data[,1]))
        )
        stopifnot(isTRUE(all.equal(allele_data[,2],tissue_input_data[,2]))
                  && isTRUE(all.equal(allele_data[,2],normal_input_data[,2]))
        )
        # tissue_input_data = tissue_input_data[,3:6]
        # normal_input_data = normal_input_data[,3:6]
        # If a BED is provided, only look at SNPs within those intervals
        if (!is.na(BED_file)) {
            stopifnot(file.exists(BED_file) && file.info(BED_file)$size>0)
            BED=read_BED_as_table(BED_file)[,1:3]
            colnames(BED)=c('chr','start','end')
            # BED$chr=gsub('^chr','',BED$chr)
            BED$start=BED$start+1; BED$end=BED$end+1 # Start is 0-based in BED files
            BED=BED[BED$chr %in% paste0("chr", chrom_names),]
            if (nrow(BED)==0) stop('Major issue with BED file, please double-check its content')
            requireNamespace("GenomicRanges")
            requireNamespace("IRanges")
            overlaps=GenomicRanges::findOverlaps(
                GenomicRanges::GRanges(seqnames=BED$chr,ranges=IRanges::IRanges(start=BED$start,end=BED$end)),
                GenomicRanges::GRanges(
                    seqnames=allele_data$chromosome,
                    ranges=IRanges::IRanges(start=allele_data$position,end=allele_data$position)
                )
            )
            if (length(overlaps)>0) {
                tissue_input_data=tissue_input_data[unique(overlaps@to),]
                normal_input_data=normal_input_data[unique(overlaps@to),]
                allele_data=allele_data[unique(overlaps@to),]
            } else {
                print(head(allele_data))
                print(head(BED))
                stop('The overlap between the BED file and loci is empty. Data must be checked!')
            }
            rm(BED,overlaps)
        }
        # Obtain depth for both alleles for tissue and normal
        # len = nrow(allele_data)
        tissue_input_data$REF=allele_data[,3] # tissue_input_data[cbind(1:len,allele_data[,3])]
        tissue_input_data$ALT=tissue_input_data$Good_depth # allele_data[,4] # tissue_input_data[cbind(1:len,allele_data[,4])]
        normal_input_data$REF=allele_data[,3] # normal_input_data[cbind(1:len,allele_data[,3])]
        normal_input_data$ALT=normal_input_data$Good_depth # allele_data[,4] # normal_input_data[cbind(1:len,allele_data[,4])]
        # Make sure that ALT+REF fit with minimal counts
        TO_KEEP=which(
            (tissue_input_data$REF+tissue_input_data$ALT>=1) & 
            (normal_input_data$REF+normal_input_data$ALT>=minCounts)
        )
        stopifnot(length(TO_KEEP)>0)
        allele_data=allele_data[TO_KEEP,]
        tissue_input_data=tissue_input_data[TO_KEEP,]
        normal_input_data=normal_input_data[TO_KEEP,]
        rm(TO_KEEP)
        # Prepare allele counts to derive BAF and logR
        len = nrow(allele_data)
        mutCount1 = tissue_input_data$REF
        mutCount2 = tissue_input_data$ALT
        totalTissue = mutCount1 + mutCount2
        normCount1 = normal_input_data$REF
        normCount2 = normal_input_data$ALT
        totalNormal = normCount1 + normCount2
        rm(tissue_input_data,normal_input_data)
        normalBAF = vector(length=len, mode="numeric")
        tissueBAF = vector(length=len, mode="numeric")
        normalLogR = vector(length=len, mode="numeric")
        tissueLogR = vector(length=len, mode="numeric")
        # Output raw (=unmirrored) BAF from some downstream analyses (e.g. refphase)
        normalBAF_unmirrored=normCount2/totalNormal
        tissueBAF_unmirrored=mutCount2/totalTissue
        normal.BAF_unmirrored = data.frame(Chromosome=allele_data$chromosome, Position=allele_data$position, baf=normalBAF_unmirrored, ID=rownames(allele_data), row.names=4, stringsAsFactors=F)
        tissue.BAF_unmirrored = data.frame(Chromosome=allele_data$chromosome, Position=allele_data$position, baf=tissueBAF_unmirrored, ID=rownames(allele_data), row.names=4, stringsAsFactors=F)
        colnames(tissue.BAF_unmirrored)[3]=tissue_name
        colnames(normal.BAF_unmirrored)[3]=germline_name
        write.table(tissue.BAF_unmirrored,file=gsub('\\.txt$','_rawBAF.txt',tissueBAF_file), row.names=T, quote=F, sep="\t", col.names=NA)
        write.table(normal.BAF_unmirrored,file=gsub('\\.txt$','_rawBAF.txt',normalBAF_file), row.names=T, quote=F, sep="\t", col.names=NA)
        rm(normalBAF_unmirrored,tissueBAF_unmirrored,normal.BAF_unmirrored,tissue.BAF_unmirrored)
        # Randomise A and B alleles
        selector = round(runif(len))
        normalBAF[which(selector==0)] = normCount1[which(selector==0)] / totalNormal[which(selector==0)]
        normalBAF[which(selector==1)] = normCount2[which(selector==1)] / totalNormal[which(selector==1)]
        tissueBAF[which(selector==0)] = mutCount1[which(selector==0)] / totalTissue[which(selector==0)]
        tissueBAF[which(selector==1)] = mutCount2[which(selector==1)] / totalTissue[which(selector==1)]
        # Normalise tissueLogR to normalLogR
        tissueLogR = totalTissue/totalNormal
        tissueLogR = log2(tissueLogR/mean(tissueLogR, na.rm=T))
        rm(selector)
        # For males, chrX needs to be adjusted as logR baseline will be 0 because of T/N ratio
        if (gender=='XY') {
            # PAR1 and PAR2 information should be a mix of chrX and chrY so we should expect 1+1 (1 from X and 1 from Y).
            # nonPAR should be X-specific and baseline is 1+0 so logR needs to be decreased according to gamma parameter (ascat.runAscat)
            if (genomeVersion=='hg19') {
                nonPAR=c(2699521,154931043)
            } else if (genomeVersion=='hg38') {
                nonPAR=c(2781480,155701382)
            }
            nonPAR=which(allele_data$chromosome %in% c('X','chrX') & allele_data$position>=nonPAR[1] & allele_data$position<=nonPAR[2])
            tissueLogR[nonPAR]=tissueLogR[nonPAR]-1
        }
        # Create the output data.frames
        tissue.LogR = data.frame(
            Chromosome=allele_data$chromosome, Position=allele_data$position, logr=tissueLogR, ID=rownames(allele_data), 
            row.names=4, stringsAsFactors=F)
        tissue.BAF = data.frame(
            Chromosome=allele_data$chromosome, Position=allele_data$position, baf=tissueBAF, ID=rownames(allele_data), 
            row.names=4, stringsAsFactors=F)
        normal.LogR = data.frame(
            Chromosome=allele_data$chromosome, Position=allele_data$position, logr=normalLogR, ID=rownames(allele_data), 
            row.names=4, stringsAsFactors=F)
        normal.BAF = data.frame(
            Chromosome=allele_data$chromosome, Position=allele_data$position, baf=normalBAF, ID=rownames(allele_data), 
            row.names=4, stringsAsFactors=F)
        colnames(tissue.LogR)[3]=tissue_name
        colnames(tissue.BAF)[3]=tissue_name
        colnames(normal.LogR)[3]=germline_name
        colnames(normal.BAF)[3]=germline_name
        # Save data.frames to disk
        write.table(tissue.LogR,file=tissueLogR_file, row.names=T, quote=F, sep="\t", col.names=NA)
        write.table(tissue.BAF,file=tissueBAF_file, row.names=T, quote=F, sep="\t", col.names=NA)
        write.table(normal.LogR,file=normalLogR_file, row.names=T, quote=F, sep="\t", col.names=NA)
    }
}

#' Custom read.table function that skips all lines starting with "@"
#'
#' @param int_filename An object of class characterstring. Path to the interval file to be read.
#' @return Returns an object of class "data.frame". Interval BED file read as a data.frame.
#' @examples
#' BED <- read_BED_as_table(
#'     "/projects/b1131/saya/bbcar/data/02b_cnv/signatures/01_alleleCounts/1009_tissue_alleleFrequencies.txt",
#'     ) 
read_BED_as_table <- function(int_filename) {
    # Step 1: Count the number of rows starting with "@"
    num_skip <- 0
    con <- file(int_filename, "r")
    while (length(line <- readLines(con, n = 1)) > 0) {
      if (substr(line, 1, 1) == "@") {
        num_skip <- num_skip + 1
      } else {
        break
      }
    }
    close(con)
    
    # Step 2: Read the file using read.table() and skip the top num_skip rows
    data <- read.table(int_filename, sep = "\t", header = FALSE, skip = num_skip)

    return(data)
}

