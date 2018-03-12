# ******************************************************************************
# Name:    BLANK
# Author:  Vadim Korotkikh
# Date:    2018
# Description:  BLANK
#
# ******************************************************************************


import sys

pyver = sys.version_info[0:3]
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

#>******************************************************************************
def main(maxltsize=31500, nslicenum=10 ):

    genlt = list(range(0, int(maxltsize)))
    repklt = do_ltslice_repack(genlt, int(nslicenum))



#>******************************************************************************
def do_ltslice_repack(lt_arg, numpaks):

    logger.info("Repacking list into %i lists" % numpaks)
    ltlen = len(lt_arg)
    paklen = int(ltlen/numpaks) # Calculates the size of each list slice
    ltrema = ltlen%numpaks # Remainder
    if ltrema > 0:
        logger.info(" Remainder %i " % ltrema)
    nltinds = list(range(0, ltlen, paklen))
    print(nltinds)
    lt_repack = [lt_arg[n:n+paklen] for n in nltinds if (n+paklen) < ltlen]
    if ltrema > 0 and (nltinds[-1]+int(ltlen/2)) > ltlen:
        logger.info(" Last pack too small - Range - %i : %i " % (nltinds[-1], ltlen))


    return lt_repack


#>******************************************************************************
if __name__ == "__main__":

    try:
        main(sys.argv[1], sys.argv[2])
    except IndexError:
        main()
