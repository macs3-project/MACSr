import logging
import sys

def testFun():
    logging.basicConfig(level=20, stream=sys.stderr)
    logging.info("test")