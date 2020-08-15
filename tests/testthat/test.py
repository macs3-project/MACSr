import MACS2
import argparse
from MACS2 import filterdup_cmd
from MACS2 import OptValidator

class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

options = Namespace(gsize = 2.7e9,
                    format = "BEDPE", 
                    keepduplicates = "1",
                    verbose = 1,
                    outputfile = "CTCF_ChIP_200K_filterdup.bed",
                    outdir = "/tmp",
                    ifile = ["/Users/qi28068/Workspace/MACSr/tests/testthat/CTCF_PE_CTRL_chr22_50k.bedpe.gz"],
                    buffer_size = 10000,
                    dryrun = False)

opts = OptValidator.opt_validate_filterdup(options)
filterdup_cmd.run(options)
