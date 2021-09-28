import numpy as np
import sys
from astropy.time import Time

def gps2mjd(gps):
    return Time(gps, format='gps').mjd

if __name__ == "__main__":
    times = gps2mjd(sys.argv[1:])
    for t in times:
        print(t)
