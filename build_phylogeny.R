
if(!require("optparse", character.only=T,quietly = T, warn.conflicts = F)){
  install.packages("optparse",repos = "http://cran.us.r-project.org")
  library("optparse", character.only=T,quietly = T, warn.conflicts = F)
}
#----------------------------------
# Input options
#----------------------------------
option_list = list(
  make_option(c("-i", "--donor_id"), action="store", default='Patient', type='character', help="Patient/donor ID to add to names of output files"),
  make_option(c("-v", "--input_nv"), action="store", default=NULL, type='character', help="Input NV matrix (rows are variants, columns are samples)"),
  make_option(c("-r", "--input_nr"), action="store", default=NULL, type='character', help="Input NR matrix (rows are variants, columns are samples)"),
  make_option(c("-c", "--cgpvaf_output"), action="store", default=NULL, type='character', help="CGPVaf output file, instead of NR/NV matrices - can be multiple files, i.e. indel and snv data for the same donor (comma-separated)"),
  make_option(c("-o", "--output_dir"), action="store", default="", type='character', help="Output directory for files"),
  make_option(c("-b", "--beta_binom_shared"), action="store", default=T, type='logical', help="Only run beta-binomial filter on shared mutations. If FALSE, run on all mutations, before germline/depth filtering"),
  make_option(c("-n", "--ncores"), action="store", default=1, type='numeric', help="Number of cores to use for the beta-binomial step"),
  make_option(c("--normal_flt"), action="store", default='PDv37is', type='character', help="Name of the dummy normal to exclude from cgpVAF output"),
  make_option(c("--snv_rho"), action="store", default=0.1, type='numeric', help="Rho value threshold for SNVs"),
  make_option(c("--indel_rho"), action="store", default=0.15, type='numeric', help="Rho value threshold for indels"),
  make_option(c("--min_cov"), action="store", default=10, type='numeric', help="Lower threshold for mean coverage across variant site"),
  make_option(c("--max_cov"), action="store", default=500, type='numeric', help="Upper threshold for mean coverage across variant site"),
  make_option(c("--only_snvs"), action="store", default=T, type='logical', help="If indel file is provided, only use SNVs to construct the tree (indels will still be mapped to branches)"),
  make_option(c("--split_trees"), action="store", default=T, type='logical', help="If both indels and SNVs are provided, plot trees separately for each."),
  make_option(c("--keep_ancestral"), action="store", default=F, type='logical', help="Keep an ancestral branch in the phylogeny for mutation mapping"),
  make_option(c("-x","--exclude_samples"), action="store", default=NULL, type='character', help="Option to manually exclude certain samples from the analysis, separate with a comma"),
  make_option(c("--cnv_samples"), action="store", default=NULL, type='character', help="Samples with CNVs, exclude from germline/depth-based filtering, separate with a comma"),
  make_option(c("--vaf_absent"), action="store", default=0.1, type='numeric', help="VAF threshold (autosomal) below which a variant is absent"),
  make_option(c("--vaf_present"), action="store", default=0.3, type='numeric', help="VAF threshold (autosomal) above which a variant is present"),
  make_option(c("-m", "--mixmodel"), action="store", default=F, type='logical', help="Use a binomial mixture model to filter out non-clonal samples?"),
  make_option(c("-t", "--tree_mut_pval"), action="store", default=0.01, type='numeric', help="Pval threshold for treemut's mutation assignment"),
  make_option(c("-g", "--genotype_conv_prob"), action="store", default=F, type='logical', help="Use a binomial mixture model to filter out non-clonal samples?"),
  make_option(c("-p", "--min_pval_for_true_somatic"), action="store", default=0.05, type='numeric', help="Pval threshold for somatic presence if generating a probabilistic genotype matrix"),
  make_option(c("--min_variant_reads_shared"), action="store", default=2, type='numeric', help="Minimum variant reads used in generating a probabilistic genotype matrix"),
  make_option(c("--min_vaf_shared"), action="store", default=2, type='numeric', help="Minimum VAF used in generating a probabilistic genotype matrix"),
  make_option(c("--create_multi_tree"), action="store", default=T, type='logical', help="Convert dichotomous tree from MPBoot to polytomous tree"),
  make_option(c("--mpboot_path"), action="store", default="", type='character', help="Path to MPBoot executable"),
  make_option(c("--germline_cutoff"), action="store", default=-5, type='numeric', help="Log10 of germline qval cutoff"),
  make_option(c("--genomeFile"), action="store", default="/nfs/cancer_ref01/Homo_sapiens/37/genome.fa", type='character', help="Reference genome fasta for plotting mutational spectra"),
  make_option(c("--plot_spectra"), action="store", default=F, type='logical', help="Plot mutational spectra?"),
  make_option(c("--max_muts_plot"), action="store", default=5000, type='numeric', help="Maximum number of SNVs to plot in mutational spectra"),
  make_option(c("--lowVAF_filter"), action="store", default=0, type='numeric', help="Minimum VAF threshold to filter out subclonal variants. Disabled by default."),
  make_option(c("--lowVAF_filter_positive_samples"), action="store", default=0, type='numeric', help="Read number to apply exact binomial filter for samples with more than given number of reads. Disabled by default."),
  make_option(c("--VAF_treshold_mixmodel"), action="store", default=0.3, type='numeric', help="VAF threshold for the mixture modelling step to consider a sample clonal")
)
opt = parse_args(OptionParser(option_list=option_list, add_help_option=T))

print(opt)

dp_pos=opt$lowVAF_filter_positive_samples
ncores=opt$ncores
lowVAF_threshold=opt$lowVAF_filter
normal_flt=opt$normal_flt
snv_rho=opt$snv_rho
genomeFile=opt$genomeFile
plot_spectra=opt$plot_spectra
VAF_treshold=opt$VAF_treshold_mixmodel
indel_rho=opt$indel_rho
min_cov=opt$min_cov
max_cov=opt$max_cov
output_dir=opt$output_dir
only_snvs=opt$only_snvs
germline_cutoff=opt$germline_cutoff
if(is.null(opt$exclude_samples)) {samples_exclude=NULL} else {samples_exclude=unlist(strsplit(x=opt$exclude_samples,split = ","))}
if(is.null(opt$cnv_samples)) {samples_with_CNVs=NULL} else {samples_with_CNVs=unlist(strsplit(x=opt$cnv_samples,split = ","))}
if(is.null(opt$cgpvaf_output)) {cgpvaf_paths=NULL} else {cgpvaf_paths=unlist(strsplit(x=opt$cgpvaf_output,split = ","))}
keep_ancestral=opt$keep_ancestral
patient_ID=opt$donor_id
output_dir=opt$output_dir
nv_path=opt$input_nv
nr_path=opt$input_nr
max_muts_plot=opt$max_muts_plot
VAF_present=opt$vaf_present
VAF_absent=opt$vaf_absent
mixmodel=opt$mixmodel
split_trees=opt$split_trees
genotype_conv_prob=opt$genotype_conv_prob
min_pval_for_true_somatic_SHARED = opt$min_pval_for_true_somatic
min_variant_reads_SHARED=opt$min_variant_reads_shared
min_vaf_SHARED=opt$min_vaf_shared
tree_mut_pval=opt$tree_mut_pval
beta_binom_shared=opt$b
create_multi_tree=opt$create_multi_tree
path_to_mpboot=opt$mpboot_path

#----------------------------------
# Load packages (install if they are not installed yet)
#----------------------------------
options(stringsAsFactors = F)
cran_packages=c("ggplot2","ape","seqinr","stringr","data.table","tidyr","dplyr","VGAM","MASS","devtools")
bioconductor_packages=c("Rsamtools","GenomicRanges")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}
if (!require("BiocManager", quietly = T, warn.conflicts = F))
  install.packages("BiocManager")
for(package in bioconductor_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    BiocManager::install(as.character(package))
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}
if(!require("treemut", character.only=T,quietly = T, warn.conflicts = F)){
  install_git("https://github.com/NickWilliamsSanger/treemut")
  library("treemut",character.only=T,quietly = T, warn.conflicts = F)
}

#----------------------------------
# Functions
#----------------------------------

plot_spectrum = function(bed,save,add_to_title="",genomeFile = "/nfs/cancer_ref01/Homo_sapiens/37/genome.fa"){
  mutations=as.data.frame(bed)
  colnames(mutations) = c("chr","pos","ref","mut")
  mutations$pos=as.numeric(mutations$pos)
  mutations = mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")) & mutations$chr %in% c(1:22,"X","Y"),]
  mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$chr, IRanges(as.numeric(mutations$pos)-1, 
                                                                                     as.numeric(mutations$pos)+1))))
  # 2. Annotating the mutation from the pyrimidine base
  ntcomp = c(T="A",G="C",C="G",A="T")
  mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
  mutations$trinuc_ref_py = mutations$trinuc_ref
  for (j in 1:nrow(mutations)) {
    if (mutations$ref[j] %in% c("A","G")) { # Purine base
      mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
      mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
    }
  }
  
  # 3. Counting subs
  freqs = table(paste(mutations$sub,paste(substr(mutations$trinuc_ref_py,1,1),substr(mutations$trinuc_ref_py,3,3),sep="-"),sep=","))
  sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
  freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
  
  xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
  
  if(!is.null(save)) pdf(save,width=12,height=4)
  if(is.null(save)) dev.new(width=12,height=4)
  colvec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
  y = freqs_full; maxy = max(y)
  h = barplot(y, las=2, col=colvec, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab="Number mutations", main=paste0("Number of mutations: ",sum(freqs_full), add_to_title))
  for (j in 1:length(sub_vec)) {
    xpos = h[c((j-1)*16+1,j*16)]
    rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colvec[j*16])
    text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
  }    
  if(!is.null(save)) dev.off()
  
}

exact.binomial=function(gender,NV,NR,cutoff=-5,qval_return=F){
  # Function to filter out germline variants based on unmatched
  # variant calls of multiple samples from same individual (aggregate coverage
  # ideally >150 or so, but will work with less). NV is matrix of reads supporting 
  # variants and NR the matrix with total depth (samples as columns, mutations rows, 
  # with rownames as chr_pos_ref_alt or equivalent). Returns a logical vector, 
  # TRUE if mutation is likely to be germline.
  
  XY_chromosomal = grepl("X|Y",rownames(NR))
  autosomal = !XY_chromosomal
  
  if(gender=="female"){
    NV_vec = rowSums(NV)
    NR_vec = rowSums(NR)
    pval = rep(1,length(NV_vec))
    for (n in 1:length(NV_vec)){
      if(NR_vec[n]>0){
        pval[n] = binom.test(x=NV_vec[n],
                             n=NR_vec[n],
                             p=0.5,alt='less')$p.value
      }
    }
  }
  # For male, split test in autosomal and XY chromosomal part
  if(gender=="male"){
    pval=rep(1,nrow(NV))
    NV_vec = rowSums(NV)[autosomal]
    NR_vec = rowSums(NR)[autosomal]
    pval_auto = rep(1,sum(autosomal))
    pval_XY = rep(1,sum(XY_chromosomal))
    
    for (n in 1:sum(autosomal)){
      if(NR_vec[n]>0){
        pval_auto[n] = binom.test(x=NV_vec[n],
                                  n=NR_vec[n],
                                  p=0.5,alt='less')$p.value
      }
    }
    NV_vec = rowSums(NV)[XY_chromosomal]
    NR_vec = rowSums(NR)[XY_chromosomal]
    for (n in 1:sum(XY_chromosomal)){
      if(NR_vec[n]>0){
        pval_XY[n] = binom.test(x=NV_vec[n],
                                n=NR_vec[n],
                                p=0.95,alt='less')$p.value
      }
    }
    pval[autosomal]=pval_auto
    pval[XY_chromosomal]=pval_XY
  }
  qval = p.adjust(pval,method="BH")
  if(qval_return){
    return(qval)
  }else{
    germline = log10(qval)>cutoff
    return(germline)
  }
}

estimateRho_gridml = function(NV_vec,NR_vec) {
  # Function to estimate maximum likelihood value of rho for beta-binomial
  rhovec = 10^seq(-6,-0.05,by=0.05) # rho will be bounded within 1e-6 and 0.89
  mu=sum(NV_vec)/sum(NR_vec)
  ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=NV_vec, size=NR_vec, rho=rhoj, prob=mu, log=T)))
  return(rhovec[ll==max(ll)][1])
}

beta.binom.filter = function(NR,NV){
  # Function to apply beta-binomial filter for artefacts. Works best on sets of
  # clonal samples (ideally >10 or so). As before, takes NV and NR as input. 
  # Optionally calculates pvalue of likelihood beta-binomial with estimated rho
  # fits better than binomial. This was supposed to protect against low-depth variants,
  # but use with caution. Returns logical vector with good variants = TRUE
  
  rho_est = pval = rep(NA,nrow(NR))
  for (k in 1:nrow(NR)){
    rho_est[k]=estimateRho_gridml(NV_vec = as.numeric(NV[k,]),
                                  NR_vec=as.numeric(NR[k,]))
  }
  return(rho_est)
}

dbinomtrunc = function(x, size, prob, minx=4) {
  dbinom(x, size, prob) / pbinom(minx-0.1, size, prob, lower.tail=F)
}

estep = function(x,size,p.vector,prop.vector,ncomp, mode){
  ## p.vector = vector of probabilities for the individual components
  ## prop.vector = vector of proportions for the individual components
  ## ncomp = number of components
  p.mat_estep = matrix(0,ncol=ncomp,nrow=length(x))
  for (i in 1:ncomp){
    if(mode=="Truncated") p.mat_estep[,i]=prop.vector[i]*dbinomtrunc(x,size,prob=p.vector[i])
    if(mode=="Full") p.mat_estep[,i]=prop.vector[i]*dbinom(x,size,prob=p.vector[i])
  }
  norm = rowSums(p.mat_estep) ## normalise the probabilities
  p.mat_estep = p.mat_estep/norm
  LL = sum(log(norm)) ## log-likelihood
  
  ## classification of observations to specific components (too crude?)
  which_clust = rep(1,length(x))
  if(ncomp>1){
    which_clust = apply(p.mat_estep, 1, which.max)
  }
  
  list("posterior"=p.mat_estep,
       "LL"=LL,
       "Which_cluster"=which_clust)
}

mstep = function(x,size,e.step){
  # estimate proportions
  prop.vector_temp = colMeans(e.step$posterior)
  # estimate probabilities
  p.vector_temp = colSums(x/size*e.step$posterior) / colSums(e.step$posterior)
  
  list("prop"=prop.vector_temp,
       "p"=p.vector_temp)   
}

em.algo = function(x,size,prop.vector_inits,p.vector_inits,maxit=5000,tol=1e-6,nclust,binom_mode){
  ## prop.vector_inits =  initial values for the mixture proportions
  ## p.vector_inits =  initial values for the probabilities 
  
  # Initiate EM
  flag = 0
  e.step = estep(x,size,p.vector = p.vector_inits,prop.vector = prop.vector_inits,ncomp=nclust,mode=binom_mode)
  m.step = mstep(x,size,e.step)
  prop_cur = m.step[["prop"]]
  p_cur = m.step[["p"]]
  cur.LL = e.step[["LL"]]
  LL.vector = e.step[["LL"]]
  
  # Iterate between expectation and maximisation steps
  for (i in 2:maxit){
    e.step = estep(x,size,p.vector = p_cur,prop.vector = prop_cur,ncomp=nclust,mode=binom_mode)
    m.step = mstep(x,size,e.step)
    prop_new = m.step[["prop"]]
    p_new = m.step[["p"]]
    
    LL.vector = c(LL.vector,e.step[["LL"]])
    LL.diff = abs((cur.LL - e.step[["LL"]]))
    which_clust = e.step[["Which_cluster"]]
    # Stop iteration if the difference between the current and new log-likelihood is less than a tolerance level
    if(LL.diff < tol){ flag = 1; break}
    
    # Otherwise continue iteration
    prop_cur = prop_new; p_cur = p_new; cur.LL = e.step[["LL"]]
    
  }
  if(!flag) warning("Didn’t converge\n")
  
  BIC = log(length(x))*nclust*2-2*cur.LL
  AIC = 4*nclust-2*cur.LL
  list("LL"=LL.vector,
       "prop"=prop_cur,
       "p"=p_cur,
       "BIC"=BIC,
       "AIC"=AIC,
       "n"=nclust,
       "Which_cluster"=which_clust)
}

binom_mix = function(x,size,nrange=1:3,criterion="BIC",maxit=5000,tol=1e-6, mode="Full"){
  ## Perform the EM algorithm for different numbers of components
  ## Select best fit using the Bayesian Information Criterion (BIC) 
  ## or the Akaike information criterion (AIC)
  i=1
  results = list()
  BIC_vec = c()
  AIC_vec = c()
  
  for (n in nrange){
    ## Initialise EM algorithm with values from kmeans clustering
    init = kmeans(x/size,n)
    prop_init = init$size/length(x)
    p_init = init$centers
    
    results[[i]] = em.algo(x,size,prop.vector_inits = prop_init,p.vector_inits=p_init,nclust=n,maxit,tol,binom_mode=mode)
    BIC_vec = c(BIC_vec,results[[i]]$BIC)
    AIC_vec = c(AIC_vec,results[[i]]$AIC)
    i=i+1
  }
  if (criterion=="BIC"){
    results[[which.min(BIC_vec)]]$BIC_vec=BIC_vec
    return(results[[which.min(BIC_vec)]])
  }
  if (criterion=="AIC"){
    return(results[[which.min(AIC_vec)]])
  }
}

binom_pval_matrix = function(NV,NR,gender,qval_return=F) {
  NR_nonzero=NR
  NR_nonzero[NR_nonzero==0]=1
  pval_mat <- matrix(0, nrow = nrow(NV), ncol = ncol(NV))
  rownames(pval_mat)=rownames(NV)
  colnames(pval_mat)=colnames(NV)
  if(gender == "male") {
    for(i in 1:nrow(NV)) {
      for (j in 1:ncol(NV)) {
        if (!grepl("X|Y",rownames(NV)[1])) {pval_mat[i,j] <- binom.test(NV[i,j], NR_nonzero[i,j], p = 0.5, alternative = "less")$p.value}
        else {pval_mat[i,j] <- binom.test(NV[i,j], NR_nonzero[i,j], p = 0.95, alternative = "less")$p.value}
      }
    }
  } else if(gender == "female") {
    for(i in 1:nrow(NV)) {
      for (j in 1:ncol(NV)) {
        pval_mat[i,j] <- binom.test(NV[i,j], NR_nonzero[i,j], p = 0.5, alternative = "less")$p.value
      }
    }
  }
  if(qval_return){
    qval_mat=matrix(p.adjust(as.vector(pval_mat), method='BH'),ncol=ncol(pval_mat))
    rownames(qval_mat)=rownames(NV)
    colnames(qval_mat)=colnames(NV)
    return(qval_mat)
  }else{
    return(pval_mat)
  }
}

apply_mix_model=function(NV,NR,plot=T,prop_cutoff=0.15){
  peak_VAF=rep(0,ncol(NV))
  names(peak_VAF)=colnames(NV)
  autosomal=!grepl("X|Y",rownames(NV))
  for(s in colnames(NV)){
    muts_include=NV[,s]>3&autosomal
    NV_vec=NV[muts_include,s]
    NR_vec=NR[muts_include,s]
    res=binom_mix(NV_vec,NR_vec,mode="Truncated",nrange=1:3)
    saveRDS(res,paste0(output_dir,s,"_binom_mix.Rdata"))
    
    if(plot){
      pdf(paste0(output_dir,s,"_binom_mix.pdf"))
      p=hist(NV_vec/NR_vec,breaks=20,xlim=c(0,1),col='gray',freq=F,xlab="Variant Allele Frequency",
             main=paste0(s,", (n=",length(NV_vec),")"))
      cols=c("red","blue","green","magenta","cyan")
      
      y_coord=max(p$density)-0.5
      y_intv=y_coord/5
      
      for (i in 1:res$n){
        depth=rpois(n=5000,lambda=median(NR_vec))
        sim_NV=unlist(lapply(depth,rbinom,n=1,prob=res$p[i]))
        sim_VAF=sim_NV/depth
        sim_VAF=sim_VAF[sim_NV>3]
        dens=density(sim_VAF)
        lines(x=dens$x,y=res$prop[i]*dens$y,lwd=2,lty='dashed',col=cols[i])
        y_coord=y_coord-y_intv/2
        text(y=y_coord,x=0.9,label=paste0("p1: ",round(res$p[i],digits=2)))
        segments(lwd=2,lty='dashed',col=cols[i],y0=y_coord+y_intv/4,x0=0.85,x1=0.95)
      }
      dev.off()
    }
    peak_VAF[s]=max(res$p[res$prop>prop_cutoff])
  }
  return(peak_VAF)
}

add_ancestral_outgroup=function(tree,outgroup_name="Ancestral"){
  #This function adds the ancestral tip at the end
  tmp=tree$edge
  N=length(tree$tip.label)
  newroot=N+2
  renamedroot=N+3
  ancestral_tip=N+1
  tmp=ifelse(tmp>N,tmp+2,tmp)
  
  tree$edge=rbind(c(newroot,renamedroot),tmp,c(newroot,ancestral_tip))
  tree$edge.length=c(0,tree$edge.length,0)
  
  tree$tip.label=c(tree$tip.label,outgroup_name)
  tree$Nnode=tree$Nnode+1
  mode(tree$Nnode)="integer"
  mode(tree$edge)="integer"
  return(tree)
}

low_vaf_in_pos_samples = function(NR, NV, gender, define_pos = 3, qval_return = F) {
  pval=rep(0,nrow(NR))
  if(gender == "male") {
    for(n in 1:nrow(NR)) {
      NV_vec=NV[n,]
      NR_vec=NR[n,]
      if(any(NV_vec >= define_pos)){
        NV_vec_pos=NV_vec[which(NV_vec >= define_pos)]
        NR_vec_pos=NR_vec[which(NV_vec >= define_pos)]
        if (grepl("X|Y",rownames(NR)[n])) {
          pval[n]=binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.95, alt = "less")$p.value
        } else {
          pval[n]=binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.5, alt = "less")$p.value
        }
      }
    }
  } else if(gender == "female") {
    for(n in 1:nrow(NR)) {
      NV_vec=NV[n,]
      NR_vec=NR[n,]
      if(any(NV_vec >= define_pos)){
        NV_vec_pos=NV_vec[which(NV_vec >= define_pos)]
        NR_vec_pos=NR_vec[which(NV_vec >= define_pos)]
        pval[n]=binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.5, alt = "less")$p.value
      }
    }
  }
  if(qval_return){
    return(p.adjust(pval,method="BH"))
  }else{
    return(pval)
  }
}

#----------------------------------
# Read in data
#----------------------------------
print("Reading in data...")

if(!is.null(cgpvaf_paths)){
  if(length(cgpvaf_paths)==1){
    data = fread(cgpvaf_paths,header=T,data.table=F)
    Muts = paste(data$Chrom,data$Pos,data$Ref,data$Alt,sep="_")
    NR = data[,grepl("DEP",colnames(data))&!grepl(paste(c(normal_flt,samples_exclude),collapse="|"),colnames(data))]
    NV = data[,grepl("MTR",colnames(data))&!grepl(paste(c(normal_flt,samples_exclude),collapse="|"),colnames(data))]
    rownames(NV)=rownames(NR)=Muts
    samples=colnames(NR)=colnames(NV)=gsub("_DEP","",colnames(NR))
  }else{
    NR=NV=Muts=c()
    for(n in 1:length(cgpvaf_paths)){
      data = fread(cgpvaf_paths[n],header=T,data.table=F)
      Muts = c(Muts,paste(data$Chrom,data$Pos,data$Ref,data$Alt,sep="_"))
      NR = rbind(NR,data[,grepl("DEP",colnames(data))&!grepl(paste(c(normal_flt,samples_exclude),collapse="|"),colnames(data))])
      NV = rbind(NV,data[,grepl("MTR",colnames(data))&!grepl(paste(c(normal_flt,samples_exclude),collapse="|"),colnames(data))])
    }
    rownames(NV)=rownames(NR)=Muts
    samples=colnames(NR)=colnames(NV)=gsub("_DEP","",colnames(NR))
  }
}else{    
  if(!is.null(nr_path)&!is.null(nv_path)){
    NR = fread(nr_path,data.table=F)
    rownames(NR)=NR[,1]
    NR=NR[,-1]
    
    NV = fread(nv_path,data.table=F)
    rownames(NV)=NV[,1]
    NV=NV[,-1]
    
    samples=colnames(NV)
    Muts=rownames(NV)
  }else{
    print("Please provide either NV and NR files or a path to CGPVaf output")
    break
  }
}

Muts_coord=matrix(ncol=4,unlist(strsplit(Muts,split="_")),byrow = T)
if(all(nchar(Muts_coord[,3])==1&nchar(Muts_coord[,4]))==1){
  mut_id="snv"
} else{
  if(all(nchar(Muts_coord[,3])>1|nchar(Muts_coord[,4])>1)){
    mut_id="indel"
  } else{
    mut_id="both"
  }
}
print(paste0("Mutations in data:", mut_id))

XY_chromosomal = grepl("X|Y",Muts)
autosomal = !XY_chromosomal
xy_depth=mean(rowMeans(NR[XY_chromosomal,]))
autosomal_depth=mean(rowMeans(NR[autosomal,]))

gender='male'
if(xy_depth>0.8*autosomal_depth) gender='female'

noCNVs=!samples%in%samples_with_CNVs

#----------------------------------
# Filtering
#----------------------------------
if(output_dir!="") system(paste0("mkdir -p ",output_dir))
print("Starting filtering...")

filter_df=as.data.frame(matrix(ncol=4,unlist(strsplit(rownames(NV),split="_")),byrow = T))
rownames(filter_df)=rownames(NV)
colnames(filter_df)=c("Chr","Pos","Ref","Alt")

filter_df$Mean_Depth=rowMeans(NR[,noCNVs])
# Filter out variant sites with high and low depth across samples
if(gender=='male'){
  filter_df$Depth_filter = (rowMeans(NR[,noCNVs])>min_cov&rowMeans(NR[,noCNVs])<max_cov&autosomal)|
    (rowMeans(NR[,noCNVs])>(min_cov/2)&rowMeans(NR[,noCNVs])<(max_cov/2)&XY_chromosomal)
}else{
  filter_df$Depth_filter = rowMeans(NR)>min_cov&rowMeans(NR)<max_cov
}

# Filter out variants likely to be germline
germline_qval=exact.binomial(gender=gender,NV=NV[,noCNVs],NR=NR[,noCNVs],qval_return=T) 
filter_df$Germline_qval=germline_qval
filter_df$Germline=as.numeric(log10(germline_qval)<germline_cutoff)


if(lowVAF_threshold>0){
  NR_nonzero=NR
  NR_nonzero[NR_nonzero==0]=1
  VAF=NV/NR_nonzero
  filter_df$lowVAF=rowSums(VAF>lowVAF_threshold)>0
}

if(beta_binom_shared){
  print("Running beta-binomial on shared mutations...")
  
  if(lowVAF_threshold>0){
    NR_flt=NR[filter_df$Germline&
                filter_df$Depth_filter&
                filter_df$lowVAF,]
    NV_flt=NV[filter_df$Germline&
                filter_df$Depth_filter&
                filter_df$lowVAF,]
  }else{
    NR_flt=NR[filter_df$Germline&
                filter_df$Depth_filter,]
    NV_flt=NV[filter_df$Germline&
                filter_df$Depth_filter,]
  }
  
  NR_flt_nonzero=NR_flt
  NR_flt_nonzero[NR_flt_nonzero==0]=1
  
  # Find shared variants and run beta-binomial filter  
  shared_muts=rownames(NV_flt)[rowSums(NV_flt>0)>1]
  
  if(ncores>1){
    rho_est=unlist(mclapply(shared_muts,function(x){
      estimateRho_gridml(NR_vec=as.numeric(NR_flt_nonzero[x,]),NV_vec=as.numeric(NV_flt[x,]))
    },mc.cores=ncores))
  }else{
    rho_est = beta.binom.filter(NR=NR_flt_nonzero[shared_muts,],NV=NV_flt[shared_muts,])
  }
  
  filter_df$Beta_binomial=filter_df$Rho=NA
  filter_df[shared_muts,"Rho"]=rho_est
  filter_df[shared_muts,"Beta_binomial"]=1
  
  if(mut_id=="snv")flt_rho=rho_est<snv_rho
  if(mut_id=="indel")flt_rho=rho_est<indel_rho
  if(mut_id=="both"){
    Muts_coord=matrix(ncol=4,unlist(strsplit(shared_muts,split="_")),byrow = T)
    is.indel=nchar(Muts_coord[,3])>1|nchar(Muts_coord[,4])>1
    flt_rho=(rho_est<indel_rho&is.indel)|(rho_est<snv_rho&!is.indel)
  }
  rho_filtered_out = shared_muts[flt_rho]
  filter_df[rho_filtered_out,"Beta_binomial"]=0
  
  NR_filtered = NR_flt[!rownames(NR_flt)%in%rho_filtered_out,]
  NV_filtered = NV_flt[!rownames(NV_flt)%in%rho_filtered_out,]
  
  if(dp_pos>0){
    pval_vec=low_vaf_in_pos_samples(NR=NR_filtered,NV=NV_filtered,gender=gender,define_pos=dp_pos)
    filter_df$lowVAF_in_pos_samples=NA
    filter_df[rownames(NR_filtered),"lowVAF_in_pos_samples"]=pval_vec>1e-3
  }
  
}else{
  print("Running beta-binomial on ALL mutations...")
  
  if(ncores>1){
    rho_est=unlist(mclapply(1:nrow(NR),function(x){
      estimateRho_gridml(NR_vec=as.numeric(NR[x,]),NV_vec=as.numeric(NV[x,]))
    },mc.cores=ncores))
  }else{
    rho_est=beta.binom.filter(NR=NR, NV=NV)
  }
  
  filter_df$Rho=rho_est
  if(mut_id=="snv")filter_df$Beta_binomial=as.numeric(rho_est>snv_rho&!is.na(rho_est))
  if(mut_id=="indel")filter_df$Beta_binomial=as.numeric(rho_est>indel_rho&!is.na(rho_est))
  if(mut_id=="both"){
    is.indel=nchar(filter_df$Ref)>1|nchar(filter_df$Alt)>1
    filter_df$Beta_binomial=as.numeric(((rho_est>indel_rho&is.indel)|(rho_est>snv_rho&!is.indel))&!is.na(rho_est))
  }
  
  
  filter_names=c("Depth_filter","Germline","Beta_binomial")
  if(lowVAF_threshold>0){
    filter_names=c(filter_names,"lowVAF")
  }
  if(dp_pos>0){
    pval_vec=low_vaf_in_pos_samples(NR=NR,NV=NV,gender=gender,define_pos=dp_pos)
    filter_df$lowVAF_in_pos_samples_pval=pval
    filter_df$lowVAF_in_pos_samples=pval_vec>1e-3
    filter_names=c(filter_names,"lowVAF_in_pos_samples")
  }
  NV_filtered=NV[rowSums(filter_df[,filter_names])==length(filter_names),]
  NR_filtered=NR[rowSums(filter_df[,filter_names])==length(filter_names),]
}

write.table(NR_filtered,paste0(output_dir,patient_ID,"_",mut_id,"_NR_filtered_all.txt"))
write.table(NV_filtered,paste0(output_dir,patient_ID,"_",mut_id,"_NV_filtered_all.txt"))
write.table(filter_df,paste0(output_dir,patient_ID,"_",mut_id,"_filtering_all.txt"))
#----------------------------------
# Mix_model decompose
#----------------------------------
if(mixmodel){
  print("Running mixture modelling...")
  peak_VAF=apply_mix_model(NV=NV_filtered,NR=NR_filtered)
  clonal=peak_VAF>VAF_treshold
  NV_filtered=NV_filtered[,clonal]
  NR_filtered=NR_filtered[,clonal]
}

#----------------------------------
# Plot mutational spectra
#----------------------------------

if(plot_spectra){
  print("Plotting mutational spectra...")
  filter_snvs_df=filter_df[nchar(filter_df$Ref)==1&nchar(filter_df$Alt)==1,]
  if(sum(!filter_snvs_df$Germline)<max_muts_plot){
    plot_spectrum(filter_snvs_df[filter_snvs_df$Germline==0,1:4], save=paste0(output_dir,patient_ID,"_",mut_id,"_germline_spectrum.pdf"),genomeFile = genomeFile)
  }else{
    subset=sample(which(filter_snvs_df$Germline==0),max_muts_plot,replace = F)
    plot_spectrum(filter_snvs_df[subset,1:4], save=paste0(output_dir,patient_ID,"_",mut_id,"_germline_spectrum.pdf"),genomeFile = genomeFile)
  }
  
  if(sum(!filter_snvs_df$Depth_filter)<max_muts_plot){
    plot_spectrum(filter_snvs_df[filter_snvs_df$Depth_filter==0,1:4], save=paste0(output_dir,patient_ID,"_",mut_id,"_highlow_depth_spectrum.pdf"),genomeFile = genomeFile)
  }else{
    subset=sample(which(filter_snvs_df$Depth_filter==0),max_muts_plot,replace = F)
    plot_spectrum(filter_snvs_df[subset,1:4], save=paste0(output_dir,patient_ID,"_",mut_id,"_highlow_depth_spectrum.pdf"),genomeFile = genomeFile)
  }
  
  if(sum(filter_snvs_df$Beta_binomial==0&!is.na(filter_snvs_df$Beta_binomial))<max_muts_plot){
    plot_spectrum(filter_snvs_df[filter_snvs_df$Beta_binomial==0&!is.na(filter_snvs_df$Beta_binomial),1:4], save=paste0(output_dir,patient_ID,"_",mut_id,"_bbinomial_spectrum.pdf"),genomeFile = genomeFile)
  }else{
    subset=sample(which(filter_snvs_df$Beta_binomial==0&!is.na(filter_snvs_df$Beta_binomial)),max_muts_plot,replace = F)
    plot_spectrum(filter_snvs_df[subset,1:4], save=paste0(output_dir,patient_ID,"_",mut_id,"_bbinomial_spectrum.pdf"),genomeFile = genomeFile)
  }
  
  if(nrow(NV_filtered)<max_muts_plot){
    plot_spectrum(filter_snvs_df[rownames(NV_filtered),1:4], save=paste0(output_dir,patient_ID,"_",mut_id,"_filtered_spectrum.pdf"),genomeFile = genomeFile)
  }else{
    subset=sample(rownames(NV_filtered),max_muts_plot,replace = F)
    plot_spectrum(filter_snvs_df[subset,1:4], save=paste0(output_dir,patient_ID,"_",mut_id,"_filtered_spectrum.pdf"),genomeFile = genomeFile)
  }
}

#----------------------------------
# Discretize genotype matrix and create fasta file for MPBoot
#----------------------------------
print("Constructing a fasta file...")

NR_flt_nonzero=NR_filtered
NR_flt_nonzero[NR_flt_nonzero==0]=1
XY_chromosomal=grepl("X|Y",rownames(NR_filtered))
autosomal=!XY_chromosomal

if(genotype_conv_prob){
  pval_matrix=binom_pval_matrix(NR=NR_filtered, NV=NV_filtered,gender=gender)
  if(!is.na(min_variant_reads_SHARED)) {min_variant_reads_mat <- NV_filtered >= min_variant_reads_SHARED} else {min_variant_reads_mat=1}
  if(!is.na(min_pval_for_true_somatic_SHARED)) {min_pval_for_true_somatic_mat <- pval_matrix > min_pval_for_true_somatic_SHARED} else {min_pval_for_true_somatic_mat=1}
  if(!is.na(min_vaf_SHARED[1]) & gender=="female") {
    min_vaf_mat <- NV_filtered/NR_flt_nonzero>min_vaf_SHARED[1]
  } else if(!is.na(min_vaf_SHARED) & gender=="male") {
    min_vaf_mat=matrix(0,ncol=ncol(NV_filtered),nrow=nrow(NV_filtered))
    min_vaf_mat[XY_chromosomal,]=NV_filtered[XY_chromosomal,]/NR_flt_nonzero[XY_chromosomal,] > min_vaf_SHARED[2]
    min_vaf_mat[autosomal,]=NV_filtered[autosomal,]/NR_flt_nonzero[!autosomal,] > min_vaf_SHARED[1]
  } else {min_vaf_mat=1}
  genotype_bin = min_variant_reads_mat * min_pval_for_true_somatic_mat * min_vaf_mat
  #Select the "not sure" samples by setting genotype to 0.5.  THIS IS THE ONLY SLIGHTLY OPAQUE BIT OF THIS FUNCTION - SET EMPIRICALLY FROM EXPERIMENTATION.
  genotype_bin[NV_filtered > 0 & pval_matrix > 0.01 & genotype_bin != 1] <- 0.5 #If have any mutant reads, set as "?" as long as p-value > 0.01
  genotype_bin[NV_filtered >= 3 & pval_matrix > 0.001 & genotype_bin != 1] <- 0.5 #If have high numbers of mutant reads, should set as "?" even if incompatible p-value (may be biased sequencing)
  genotype_bin[(NV_filtered == 0) & (pval_matrix > 0.05)] <- 0.5 #Essentially if inadequate depth to exclude mutation, even if no variant reads
  write.table(pval_matrix,paste0(output_dir,patient_ID,"_",mut_id,"_filtered_binom_pval_mat.txt"))
}else{
  genotype_bin=as.matrix(NV_filtered/NR_flt_nonzero)
  if(gender=="male"){
    genotype_bin[autosomal,][genotype_bin[autosomal,]<VAF_absent]=0
    genotype_bin[autosomal,][genotype_bin[autosomal,]>=VAF_present]=1
    genotype_bin[XY_chromosomal,][genotype_bin[XY_chromosomal,]<(2*VAF_absent)]=0
    genotype_bin[XY_chromosomal,][genotype_bin[XY_chromosomal,]>=(2*VAF_present)]=1
    genotype_bin[genotype_bin>0&genotype_bin<1]=0.5
  }
  if(gender=="female"){
    genotype_bin[genotype_bin<VAF_absent]=0
    genotype_bin[genotype_bin>=VAF_present]=1
    genotype_bin[genotype_bin>0&genotype_bin<1]=0.5
  }
}
present_vars_full=rowSums(genotype_bin>0)>0

if(only_snvs){
  Muts_coord=matrix(ncol=4,unlist(strsplit(rownames(genotype_bin),split="_")),byrow = T)
  is.indel=nchar(Muts_coord[,3])>1|nchar(Muts_coord[,4])>1
  genotype_bin=genotype_bin[!is.indel,]
}

#Create dummy fasta consisting of As (WT) and Ts (Mutant)
Ref = rep("A",nrow(genotype_bin))
Alt = rep("T",nrow(genotype_bin))
dna_strings = list()
dna_strings[1]=paste(Ref,sep="",collapse="") #Ancestral sample
for (n in 1:ncol(genotype_bin)){
  Mutations = Ref
  Mutations[genotype_bin[,n]==0.5] = '?'
  Mutations[genotype_bin[,n]==1] = Alt[genotype_bin[,n]==1]
  dna_string = paste(Mutations,sep="",collapse="")
  dna_strings[n+1]=dna_string
}

names(dna_strings)=c("Ancestral",colnames(genotype_bin))
require(seqinr)
write.fasta(dna_strings, names=names(dna_strings),paste0(output_dir,patient_ID,"_",mut_id,"_for_MPBoot.fa"))

#----------------------------------
# Build tree with MPBoot
#----------------------------------
print("Building a tree...")

system(paste0(path_to_mpboot,"mpboot -s ",output_dir,patient_ID,"_",mut_id,"_for_MPBoot.fa -bb 1000"),ignore.stdout = T)

#----------------------------------
# Map Mutations on Tree using treemut
#----------------------------------

print("Mapping mutations...")
print("Assigning mutation without an ancestral branch")

tree=read.tree(paste0(output_dir,patient_ID,"_",mut_id,"_for_MPBoot.fa.treefile"))
tree=drop.tip(tree,"Ancestral")
if(!keep_ancestral){
  tree$edge.length=rep(1,nrow(tree$edge))
  NR_tree=NR_filtered[present_vars_full,]
  NV_tree=NV_filtered[present_vars_full,]
  res=assign_to_tree(tree,
                     mtr=as.matrix(NV_tree),
                     dep=as.matrix(NR_tree))
}else{
  tree <- add_ancestral_outgroup(tree) #Re add the ancestral outgroup after making tree dichotomous - avoids the random way that baseline polytomy is resolved
  tree$edge.length = rep(1, nrow(tree$edge)) 
  
  NR_tree=NR_filtered[present_vars_full,]
  NR_tree$Ancestral=30
  NV_tree=NV_filtered[present_vars_full,]
  NV_tree$Ancestral=0
  
  p.error = rep(0.01,ncol(NV_tree))
  p.error[colnames(NV_tree)=="Ancestral"]=1e-6
  res=assign_to_tree(tree,
                     mtr=as.matrix(NV_tree),
                     dep=as.matrix(NR_tree),
                     error_rate = p.error)
}

edge_length_nonzero = table(res$summary$edge_ml[res$summary$p_else_where<tree_mut_pval])
edge_length = rep(0,nrow(tree$edge))
names(edge_length)=1:nrow(tree$edge)
edge_length[names(edge_length_nonzero)]=edge_length_nonzero
tree$edge.length=as.numeric(edge_length)

if(create_multi_tree){
  print("Converting to a multi-furcating tree structure")
  if(keep_ancestral) {
    #Maintain the dichotomy with the ancestral branch
    ROOT=tree$edge[1,1]
    current_length<-tree$edge.length[tree$edge[,1]==ROOT & tree$edge[,2]!=which(tree$tip.label=="Ancestral")]
    new_length<-ifelse(current_length==0,1,current_length)
    tree$edge.length[tree$edge[,1]==ROOT & tree$edge[,2]!=which(tree$tip.label=="Ancestral")]<-new_length
  }
  tree<-di2multi(tree) #Now make tree multifurcating
  #Re-run the mutation assignment algorithm from the new tree
  res=assign_to_tree(tree,
                     mtr=as.matrix(NV_tree),
                     dep=as.matrix(NR_tree))
  edge_length_nonzero = table(res$summary$edge_ml[res$summary$p_else_where<tree_mut_pval])
  edge_length = rep(0,nrow(tree$edge))
  names(edge_length)=1:nrow(tree$edge)
  edge_length[names(edge_length_nonzero)]=edge_length_nonzero
  tree$edge.length=as.numeric(edge_length)
}

saveRDS(res,paste0(output_dir,patient_ID,"_",mut_id,"_assigned_to_tree.Rdata"))
write.tree(tree, paste0(output_dir,patient_ID,"_",mut_id,"_tree_with_branch_length.tree"))

if(split_trees&mut_id=="both"){
  Muts_coord=matrix(ncol=4,unlist(strsplit(rownames(NV_filtered)[present_vars_full],split="_")),byrow = T)
  is.indel=nchar(Muts_coord[,3])>1|nchar(Muts_coord[,4])>1
  
  edge_length_nonzero = table(res$summary$edge_ml[res$summary$p_else_where<tree_mut_pval&!is.indel])
  edge_length = rep(0,nrow(tree$edge))
  names(edge_length)=1:nrow(tree$edge)
  edge_length[names(edge_length_nonzero)]=edge_length_nonzero
  tree$edge.length=as.numeric(edge_length)
  pdf(paste0(output_dir,patient_ID,"_snv_tree_with_branch_length.pdf"))
  plot(tree)
  axisPhylo(side = 1,backward=F)
  dev.off()
  write.tree(tree, paste0(output_dir,patient_ID,"_snv_tree_with_branch_length.tree"))
  
  edge_length_nonzero = table(res$summary$edge_ml[res$summary$p_else_where<tree_mut_pval&is.indel])
  edge_length = rep(0,nrow(tree$edge))
  names(edge_length)=1:nrow(tree$edge)
  edge_length[names(edge_length_nonzero)]=edge_length_nonzero
  tree$edge.length=as.numeric(edge_length)
  pdf(paste0(output_dir,patient_ID,"_indel_tree_with_branch_length.pdf"))
  plot(tree)
  axisPhylo(side = 1,backward=F)
  dev.off()
  write.tree(tree, paste0(output_dir,patient_ID,"_indel_tree_with_branch_length.tree"))
  
}else{
  pdf(paste0(output_dir,patient_ID,"_",mut_id,"_tree_with_branch_length.pdf"))
  plot(tree)
  axisPhylo(side = 1,backward=F)
  dev.off()
  
  tree_collapsed=tree
  tree_collapsed$edge.length=rep(1,nrow(tree_collapsed$edge))
  pdf(paste0(output_dir,patient_ID,"_",mut_id,"_tree_with_equal_branch_length.pdf"))
  plot(tree_collapsed)
  dev.off()
}

Mutations_per_branch=as.data.frame(matrix(ncol=4,unlist(strsplit(rownames(NR_tree),split="_")),byrow = T))
colnames(Mutations_per_branch)=c("Chr","Pos","Ref","Alt")
Mutations_per_branch$Branch = tree$edge[res$summary$edge_ml,2]
Mutations_per_branch=Mutations_per_branch[res$summary$p_else_where<tree_mut_pval,]
Mutations_per_branch$Patient = patient_ID
Mutations_per_branch$SampleID = paste(patient_ID,Mutations_per_branch$Branch,sep="_")
write.table(Mutations_per_branch,paste0(output_dir,patient_ID,"_",mut_id,"_assigned_to_branches.txt"),quote=F,row.names=F,sep="\t")



