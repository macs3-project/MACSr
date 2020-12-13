import MACS3
import argparse
from MACS3.Utilities import OptValidator
import MACS3.Commands.callvar_cmd as callvar

class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

options = Namespace(peakbed = '../../../MACS/test/callvar_examples/PE_demo/PEsample_peaks_sorted.bed',
                    tfile = '../../../MACS/test/callvar_examples/PE_demo/PEsample_peaks_sorted.bam',
                    cfile = '../../../MACS/test/callvar_examples/PE_demo/PEcontrol_peaks_sorted.bam',
                    ofile = '/tmp/out.vcf',
                    verbose = 2,
                    GQCutoffHetero = 0,
                    GQCutoffHomo = 0,
                    Q = 20,
                    maxDuplicate = 1,
                    fermi = 'auto',
                    fermiMinOverlap = 30,
                    top2allelesMinRatio = 0.8,
                    altalleleMinCount = 2,
                    maxAR = 0.95,
                    np = 1)

callvar.run(options)
