#' hmmratac
#'
#' Dedicated peak calling based on Hidden Markov Model for ATAC-seq
#' data.
#' @param bam Sorted BAM files containing the ATAC-seq reads. If
#'     multiple files are given as '-t A B C', then they will all be
#'     read and pooled together. REQUIRED.
#' @param outdir If specified all output files will be written to that
#'     directory. Default: the current working directory
#' @param name Name for this experiment, which will be used as a
#'     prefix to generate output file names. DEFAULT: "NA"
#' @param verbose Set verbose level of runtime message. 0: only show
#'     critical message, 1: show additional warning message, 2: show
#'     process information, 3: show debug messages. DEFAULT:2
#' @param log Whether to capture logs.
#' @param cutoff_analysis_only Only run the cutoff analysis and output
#'     a report.  After generating the report, the process will stop.
#'     The report will help user decide the three crucial parameters
#'     for `-l`, `-u`, and `-c`. So it's highly recommanded to run
#'     this first! Please read the report and instructions in `Choices
#'     of cutoff values` on how to decide the three crucial
#'     parameters.
#' @param em_skip Do not perform EM training on the fragment
#'     distribution. If set, EM_MEANS and EM.STDDEVS will be used
#'     instead. Default: False
#' @param em_means Comma separated list of initial mean values for the
#'     fragment distribution for short fragments, mono-, di-, and
#'     tri-nucleosomal fragments. Default: 50 200 400 600
#' @param em_stddevs Comma separated list of initial standard
#'     deviation values for fragment distribution for short fragments,
#'     mono-, di-, and tri-nucleosomal fragments. Default: 20 20 20 20
#' @param min_frag_p We will exclude the abnormal fragments that can't
#'     be assigned to any of the four signal tracks. After we use EM
#'     to find the means and stddevs of the four distributions, we
#'     will calculate the likelihood that a given fragment length fit
#'     any of the four using normal distribution. The criteria we will
#'     use is that if a fragment length has less than MIN_FRAG_P
#'     probability to be like either of short, mono, di, or tri-nuc
#'     fragment, we will exclude it while generating the four signal
#'     tracks for later HMM training and prediction.  The value should
#'     be between 0 and 1. Larger the value, more abnormal fragments
#'     will be allowed. So if you want to include more 'ideal'
#'     fragments, make this value smaller. Default = 0.001
#' @param hmm_binsize Size of the bins to split the pileup signals for
#'     training and decoding with Hidden Markov Model. Must >=
#'     1. Smaller the binsize, higher the resolution of the results,
#'     slower the process. Default = 10
#' @param hmm_lower Upper limit on fold change range for choosing
#'     training sites. Default: 20
#' @param hmm_upper Lower limit on fold change range for choosing
#'     training sites. Default: 10
#' @param hmm_maxTrain Maximum number of training regions to
#'     use. Default: 1000
#' @param hmm_training_flanking Training regions will be expanded to
#'     both side with this number of basepairs. The purpose is to
#'     include more background regions. Default: 1000
#' @param hmm_file A JSON file generated from previous HMMRATAC run to
#'     use instead of creating new one. When provided, HMM training
#'     will be skipped. Default: NA
#' @param hmm_training_regions Filename of training regions
#'     (previously was BED_file) to use for training HMM, instead of
#'     using foldchange settings to select. Default: NA
#' @param hmm_randomSeed Seed to set for random sampling of training
#'     regions. Default: 10151
#' @param hmm_modelonly Stop the program after generating model. Use
#'     this option to generate HMM model ONLY, which can be later
#'     applied with `--model`. Default: False
#' @param prescan_cutoff The fold change cutoff for prescanning
#'     candidate regions in the whole dataset. Then we will use HMM to
#'     predict states on these candidate regions. Higher the prescan
#'     cutoff, fewer regions will be considered. Must > 1. Default:
#'     1.2
#' @param openregion_minlen Minimum length of open region to call
#'     accessible regions. Must be larger than 0. If it is set as 0,
#'     it means no filtering on the length of the open regions
#'     called. Please note that, when bin size is small, setting a too
#'     small OPENREGION_MINLEN will bring a lot of false
#'     positives. Default: 100
#' @param pileup_short By default, HMMRATAC will pileup all fragments
#'     in order to identify regions for training and candidate regions
#'     for decoding. When this option is on, it will pileup only the
#'     short fragments to do so. Although it sounds a good idea since
#'     we assume that open region should have a lot of short
#'     fragments, it may be possible that the overall short fragments
#'     are too few to be useful. Default: False
#' @param keepduplicates Keep duplicate reads from analysis. By
#'     default, duplicate reads will be removed. Default: False
#' @param blacklist Filename of blacklisted regions to exclude
#'     (previously was BED_file). Examples are those from
#'     ENCODE. Default: NA
#' @param save_digested Save the digested ATAC signals of short-,
#'     mono-, di-, and tri- signals in three BedGraph files with the
#'     names NAME_short.bdg, NAME_mono.bdg, NAME_di.bdg, and
#'     NAME_tri.bdg. DEFAULT: False
#' @param save_likelihoods Save the likelihoods to each state
#'     annotation in three BedGraph files, named with NAME_open.bdg
#'     for open states, NAME_nuc.bdg for nucleosomal states, and
#'     NAME_bg.bdg for the background states. DEFAULT: False
#' @param save_states Save all open and nucleosomal state annotations
#'     into a BED file with the name NAME_states.bed. DEFAULT: False
#' @param save_train Save the training regions and training data into
#'     NAME_training_regions.bed and NAME_training_data.txt. Default:
#'     False
#' @param decoding_steps Number of candidate regions to be decoded at
#'     a time. The HMM model will be applied with Viterbi to find the
#'     optimal state path in each region. bigger the number,
#'     'possibly' faster the decoding process, 'definitely' larger the
#'     memory usage. Default: 1000.
#' @param buffer_size Buffer size for incrementally increasing
#'     internal array size to store reads alignment information. In
#'     most cases, you don't have to change this parameter. However,
#'     if there are large number of chromosomes/contigs/scaffolds in
#'     your alignment, it's recommended to specify a smaller buffer
#'     size in order to decrease memory usage (but it will take longer
#'     time to read alignment files). Minimum memory requested for
#'     reading an alignment file is about # of CHROMOSOME *
#'     BUFFER_SIZE * 8 Bytes. DEFAULT: 100000
#' @param ... More options for macs2.
#' @export
hmmratac <- function(bam, outdir = ".", name = "NA", verbose = 2L, log = TRUE,
                     cutoff_analysis_only = FALSE,
                     em_skip = FALSE,
                     em_means = list(50, 200, 400, 600),
                     em_stddevs = list(20, 20, 20, 20),
                     min_frag_p = 0.001,
                     hmm_binsize = 10L,
                     hmm_lower = 10L,
                     hmm_upper = 20L,
                     hmm_maxTrain = 1000,
                     hmm_training_flanking = 1000,
                     hmm_file = NULL,
                     hmm_training_regions = NULL,
                     hmm_randomSeed = 10151,
                     hmm_modelonly = FALSE,
                     prescan_cutoff = 1.2,
                     openregion_minlen = 100,
                     pileup_short = FALSE,
                     keepduplicates = FALSE,
                     blacklist = NULL,
                     save_digested = FALSE,
                     save_likelihoods = FALSE,
                     save_states = FALSE,
                     save_train = FALSE,
                     decoding_steps = 1000,
                     buffer_size = 100000,
                     ...){
    if(is.character(bam)){
        bam_file <- as.list(normalizePath(bam))
    }
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    cl <- basiliskStart(env_macs)
    on.exit(basiliskStop(cl))
    res <- basiliskRun(cl, function(.namespace, outdir,
                                    ...){
        opts <- .namespace()$Namespace(bam_file = bam_file,
                                       name = name,
                                       outdir = outdir,
                                       cutoff_analysis_only = cutoff_analysis_only,
                                       em_skip = em_skip,
                                       em_means = em_means,
                                       em_stddevs = em_stddevs,
                                       min_frag_p = min_frag_p,
                                       hmm_binsize = hmm_binsize,
                                       hmm_lower = hmm_lower,
                                       hmm_upper = hmm_upper,
                                       hmm_maxTrain = hmm_maxTrain,
                                       hmm_randomSeed = hmm_randomSeed,
                                       hmm_training_flanking = hmm_training_flanking,
                                       hmm_file = hmm_file,
                                       hmm_training_regions = hmm_training_regions,
                                       hmm_modelonly = hmm_modelonly,
                                       prescan_cutoff = prescan_cutoff,
                                       openregion_minlen = openregion_minlen,
                                       pileup_short = pileup_short,
                                       misc_keep_duplicates = keepduplicates,
                                       blacklist = blacklist,
                                       save_likelihoods = save_likelihoods,
                                       save_digested = save_digested,
                                       save_states = save_states,
                                       save_train = save_train,
                                       decoding_steps = decoding_steps,
                                       buffer_size = buffer_size,
                                       verbose = verbose,
                                       ratio = NA)

        .hmmratac <- reticulate::import("MACS3.Commands.hmmratac_cmd")
        if(log){
            reticulate::py_capture_output(.hmmratac$run(opts))
        }else{
            .hmmratac$run(opts)
        }
    }, .namespace = .namespace, outdir = outdir)

    if(log){
        message(res)
    }
    outputs <- list.files(path = outdir, pattern = paste0(name, "_.*"), full.names = TRUE)
    args <- as.list(match.call())
    macsList(arguments = args, outputs = outputs, log = res)
}
