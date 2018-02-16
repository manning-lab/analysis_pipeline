library(argparser)
library(TopmedPipeline)
library(SeqVarTools)
library(dplyr)
sessionInfo()

argp <- arg_parser("LocusZoom plots")
argp <- add_argument(argp, "config", help="path to config file")
argp <- add_argument(argp, "--segment", help="row in locus_file to plot", default=1, type="integer")
argp <- add_argument(argp, "--assoc_file", help="Association File", type="character")
argp <- add_argument(argp, "--locus_file", help="Locus File", type="character")
argp <- add_argument(argp, "--gds_file", help="GDS File", type="character")
argp <- add_argument(argp, "--variantIDcol", help="column in association file with variantID matching GDS file; NA for missing", default=NA, type="character")


argv <- parse_args(argp)
config <- readConfig(argv$config)
segment <- argv$segment

gdsfile <- argv$gds_file
print(gdsfile)

variantIDcol <- argv$variantIDcol

print(argv)
required <- c()
optional <- c("flanking_region"=500,
              "gds_file"=NA,
              "genome_build"="hg19",
              "ld_sample_include"=NA,
              "locus_type"="variant",
              "out_prefix"="locuszoom",
              "track_file"=NA,
              "track_file_type"="window",
              "track_label"="",
              "track_threshold"=5e-8)
config <- setConfigDefaults(config, required, optional)
print(config)

stopifnot(config["locus_type"] %in% c("variant", "region"))

# read selected locus
locus <- read.table(argv$locus_file, header=TRUE, as.is=TRUE)[segment,]
stopifnot(all(c("chr", "pop") %in% names(locus)))
print(locus)

# population for LD
pop <- toupper(locus$pop)
stopifnot(pop %in% c("TOPMED", "AFR", "AMR", "ASN", "EUR", "EAS", "SAS"))

## get association test results
var.chr <- locus$chr
#assocfile <- insertChromString(config["assoc_file"], var.chr)
#assoc <- getobj(assocfile)
assoc <- read.table(argv$assoc_file,header=T,sep=",",as.is=T)
head(assoc)

if(variantIDcol == "NA") {
  # variantID is not in the association file, need to add
  gds <- seqOpen(gdsfile)
  variant.id.all <- seqGetData(gds,"variant.id")
  chr.all <- seqGetData(gds,"chromosome")
  pos.all <- seqGetData(gds,"position")
  allele.all <- seqGetData(gds,"allele")
  seqClose(gds)

  small.ref <- data.frame(variantID=variant.id.all, chr=chr.all,pos=pos.all,allele=allele.all)
  print(head(small.ref))

  print(dim(assoc))
  assoc <- merge(assoc,small.ref,by=c("chr","pos"))
  print(dim(assoc))
} else {
assoc$variantID <- assoc[,variantIDcol]
}

head(assoc)


if (config["locus_type"] == "variant") {
    #stopifnot("variantID" %in% names(locus))
    markername <- locus$MarkerID

    variant <- assoc$variantID[which(assoc$pos==strsplit(markername,":")[[1]][2])]
    flank <- as.numeric(config["flanking_region"]) * 1000
    var.pos <- assoc$pos[assoc$variantID == variant]
    start <- var.pos - flank
    end <- var.pos + flank

    lz.name <- paste0("chr", var.chr, ":", var.pos)
    ld.region <- paste0("--refsnp \"", lz.name, "\"", " --flank ", config["flanking_region"], "kb")
    prefix <- paste0(config["out_prefix"], "_var", variant, "_ld_", pop)
    maf <- assoc$MAF[assoc$variantID == variant]
    title <- paste(lz.name, "- MAF:", formatC(maf, digits=3))

} else if (config["locus_type"] == "region") {
    stopifnot(all(c("start", "end") %in% names(locus)))
    start <- locus$start
    end <- locus$end

    ld.region <- paste("--chr", var.chr, "--start", start, "--end", end)
    prefix <- paste0(config["out_prefix"], "_ld_", pop)
    title <- ""
}
head(assoc)

## construct METAL-format file
print(var.chr)
print(start)
print(end)

assoc <- assoc %>%
    filter(chr == var.chr, pos > start, pos < end) %>%
    select(variantID, chr, pos, ends_with("pval"))
names(assoc)[4] <- "pval"
assoc <- assoc %>%
    filter(pval < 0.1)
dim(assoc)
summary(assoc)
head(assoc)

assoc.filename <- tempfile(pattern="assoc",tmpdir="/home")
writeMETAL(assoc, file=assoc.filename)

# LD
if (pop != "TOPMED") {
    ld.cmd <- paste("--pop", pop, "--source 1000G_Nov2014")
    ld.title <- paste("LD: 1000G", pop)
} else {
    if (!is.na(config["ld_sample_include"])) {
        sample.id <- getobj(config["ld_sample_include"])
    } else {
        sample.id <- NULL
    }
    if (config["locus_type"] == "variant") {
        ref.var <- variant
    } else {
        ref.var <- filter(assoc, pval == min(pval))$variantID
    }
    #gdsfile <- insertChromString(config["gds_file"], var.chr)
    ld <- calculateLD(gdsfile, variant.id=assoc$variantID, ref.var=ref.var, sample.id=sample.id)
    ld.filename <- tempfile(pattern="ld",tmpdir="/home")

    writeLD(assoc, ld, ref.var, file=ld.filename)

    ld.cmd <- paste("--ld", ld.filename)
    ld.title <- "LD: TOPMed"
}
title <- if (title == "") ld.title else paste(ld.title, title, sep=" - ")

## construct BED track file
if (!is.na(config["track_file"])) {
    trackfile <- insertChromString(config["track_file"], var.chr)
    track <- getAssoc(trackfile, config["track_file_type"]) %>%
        filter(pval < as.numeric(config["track_threshold"]))
    track.filename <- tempfile()
    writeBED(track, file=track.filename, track.label=config["track_label"])
    track.cmd <- paste("--bed-tracks", track.filename)
} else {
    track.cmd <- ""
}

command <- paste("/usr/local/locuszoom-standalone/bin/locuszoom",
                 "theme=publication",
                 "--cache None",
                 "--no-date",
                 "--plotonly",
                 "--gene-table gencode",
                 "--build", config["genome_build"],
                 "--chr", var.chr,
                 "--metal", assoc.filename,
                 track.cmd,
                 ld.cmd,
                 ld.region,
                 "--prefix ", prefix,
                 paste0("title=\"", title, "\""),
                 paste0("signifLine=\"", -log10(5e-8), "\" signifLineColor=\"gray\" signifLineWidth=\"2\""),
                 "ylab=\"-log10(p-value) from single variant test\"")

cat(paste(command, "\n"))
system(command)

#unlink(assoc.filename)
#if (exists("track.filename")) unlink(track.filename)
#if (exists("ld.filename")) unlink(ld.filename)
