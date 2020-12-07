
r1 <- bdgpeakcall("/Users/qi28068/Workspace/MACS/temp/0_run_bdgcmp/run_bdgcmp_FE.bdg", cutoff = 2, outdir = "/tmp/", outputfile = "bdg")

r2 <- bdgbroadcall("/Users/qi28068/Workspace/MACS/temp/0_run_bdgcmp/run_bdgcmp_FE.bdg", cutoffpeak = 2, cutofflink = 1.5, outdir="/tmp", outputfile = "bdgbroad")

r3 <- bdgcmp("/Users/qi28068/Workspace/MACS/temp/0_run_pileup/run_pileup_ChIP.bed.bdg", "/Users/qi28068/Workspace/MACS/temp/0_run_pileup/run_pileup_CTRL.bed.bdg", method = c("ppois", "FE"), pseudocount = 1, outdir = "/tmp", oprefix = "bdgcmp")

r4 <- bdgopt("/Users/qi28068/Workspace/MACS/temp/0_run_callpeak_narrow/run_callpeak_narrow0_treat_pileup.bdg",
             method = "min", extraparam = 10, outdir = "/tmp", outputfile = "bdgopt_min.bdg")

r5 <- cmbreps(list("/Users/qi28068/Workspace/MACS/temp/0_run_callpeak_narrow/run_callpeak_narrow0_treat_pileup.bdg",
                   "/Users/qi28068/Workspace/MACS/temp/0_run_callpeak_narrow/run_callpeak_narrow0_control_lambda.bdg",
                   "/Users/qi28068/Workspace/MACS/temp/0_run_bdgcmp/run_bdgcmp_ppois.bdg"),
              method = "max", outdir = "/tmp", outputfile = "cmbreps_max.bdg")

r6 <- bdgdiff(t1bdg = "/Users/qi28068/Workspace/MACS/temp/0_run_callpeak_narrow/run_callpeak_narrow0_treat_pileup.bdg",
              c1bdg = "/Users/qi28068/Workspace/MACS/temp/0_run_callpeak_narrow/run_callpeak_narrow0_control_lambda.bdg",
              t2bdg = "/Users/qi28068/Workspace/MACS/temp/0_run_callpeak_narrow_revert/run_callpeak_narrow_revert_treat_pileup.bdg",
              c2bdg = "/Users/qi28068/Workspace/MACS/temp/0_run_callpeak_narrow_revert/run_callpeak_narrow_revert_control_lambda.bdg",
              oprefix = "bdgdiff",
              outdir = "/tmp")

r7 <- pileup("/Users/qi28068/Workspace/MACS/test/CTCF_SE_ChIP_chr22_50k.bed.gz",
             format = "BED", outdir = "/tmp", outputfile = "pileup_bed.bdg")

r8 <- randsample("/Users/qi28068/Workspace/MACS/test/CTCF_SE_ChIP_chr22_50k.bed.gz",
                 number = 10000L, seed = 31415926,
                 outdir = "/tmp", outputfile = "randsample.bed")

r9 <- refinepeak(bedfile = "/Users/qi28068/Workspace/MACS/temp/0_run_callpeak_narrow/run_callpeak_narrow0_peaks.narrowPeak",
                 ifile = "/Users/qi28068/Workspace/MACS/test/CTCF_SE_ChIP_chr22_50k.bed.gz",
                 outdir = "/tmp", outputfile = "refinepeak.bed")

macsdir <- "/Users/qi28068/Workspace/MACS/test"
CALLVAR_SE_BED= file.path(macsdir, "callvar_examples/PE_demo/PEsample_peaks_sorted.bed")
CALLVAR_SE_TBAM= file.path(macsdir, "callvar_examples/PE_demo/PEsample_peaks_sorted.bam")
CALLVAR_SE_CBAM= file.path(macsdir, "callvar_examples/PE_demo/PEcontrol_peaks_sorted.bam")

r10 <- callvar(peakbed = CALLVAR_SE_BED, tfile = CALLVAR_SE_TBAM, cfile = CALLVAR_SE_CBAM,
               outdir = "/tmp", outputfile = "var.vcf")
