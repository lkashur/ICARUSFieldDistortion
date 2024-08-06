# imports
import numpy as np
from scipy.stats import bootstrap
import argparse

def main(args):
    
    # output
    with open(f'aggoffsets_{args.tpc}.txt', 'w') as outf:
        outf.write('gbn,med_low,med,med_high\n')
        
        # load all offsets
        with open(args.in_file, 'r') as f:
            lines = f.readlines()
            for i,l in enumerate(lines):
                # skip header
                if i < 1:
                    continue

                vals = l.split(',')
                vals = [float(v) for v in vals]
                gbn = vals[0]
                offsets = vals[1:]
                offsets = (offsets,)

                # minimum of two values for bootstrapping
                if len(offsets[0]) < 2:
                    med = -99999.0
                    bootstrap_ci = (-99999.0,-99999.0)
                else:
                    # median
                    med = np.median(offsets)
                    # set up confidence interval
                    bootstrap_ci = bootstrap(offsets, np.median, n_resamples=1000, method='basic', batch=100, confidence_level=0.6827).confidence_interval
                outstring = str(gbn) + ',' + str(bootstrap_ci[0]) + ',' + str(med) + ',' + str(bootstrap_ci[1]) + '\n'
                outf.write(outstring)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_file', type=str, required=True)
    parser.add_argument('--tpc', type=str, required=True)
    args = parser.parse_args()
    main(args)
