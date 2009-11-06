#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Brant Faircloth on 2009-11-04.
Copyright (c) 2009 Brant Faircloth. All rights reserved.
"""

import pdb
from Bio import SeqIO
from p3wrapr import primer


def prettyPrint(output, start, primer3):
    pl_start = start + int(primer3.primers[0]['PRIMER_LEFT'].split(',')[0])
    pl_length = primer3.primers[0]['PRIMER_LEFT'].split(',')[1]
    pr_start = start + int(primer3.primers[0]['PRIMER_RIGHT'].split(',')[0])
    pr_length = primer3.primers[0]['PRIMER_RIGHT'].split(',')[1]
    output.write(('%s,%s\t%s\t%s\t%s\t%s,%s\t%s\t%s\t%s\t%s\n') % (\
    pl_start, \
    pl_length, \
    primer3.primers[0]['PRIMER_LEFT_SEQUENCE'], \
    primer3.primers[0]['PRIMER_LEFT_TM'], \
    primer3.primers[0]['PRIMER_LEFT_GC_PERCENT'], \
    pr_start, \
    pr_length, \
    primer3.primers[0]['PRIMER_RIGHT_SEQUENCE'], \
    primer3.primers[0]['PRIMER_RIGHT_TM'], \
    primer3.primers[0]['PRIMER_RIGHT_GC_PERCENT'], \
    primer3.primers[0]['PRIMER_PAIR_PRODUCT_SIZE']\
    ))

def design(record, output, start = 0, stop = 150, target = '50,30'):
    output.write('Left Position\tLeft Primer\tLeft TM\tLeft GC\tRight Position\tRight Primer\tRight TM\tRight GC\tProduct Size\t\n')
    #pdb.set_trace()
    while stop < len(record.seq):
        seq_slice = str(record.seq[start:stop])
        primer3_settings = primer.Settings()
        primer3_settings.basic()
        # override some of the basic parameters
        primer3_settings.params['PRIMER_PRODUCT_SIZE_RANGE']        = '80-130'
        primer3_settings.params['PRIMER_MIN_GC'] = 20
        primer3_settings.params['PRIMER_MAX_GC'] = 80
        primer3_settings.params['PRIMER_GC_CLAMP'] = 0
        primer3_settings.params['PRIMER_MIN_TM']                    = 52.   # deg C
        primer3_settings.params['PRIMER_MAX_TM']                    = 60.   # deg C
        primer3_settings.params['PRIMER_OPT_TM']                    = 56.   # deg C
        primer3_settings.params['PRIMER_MAX_END_STABILITY']         = 10.0  # delta G
        primer3_settings.params['PRIMER_NUM_RETURN']                = 2     # count
        primer3_settings.params['PRIMER_MAX_POLY_X']                = 4     # nt
        primer3 = primer.Primers()
        primer3.pick(primer3_settings, sequence=seq_slice, target=target, name = 'dunno')
        #pdb.set_trace()
        # move through sequence iteratively designing primers
        #primer3.design()
        if 0 in primer3.primers.keys():
            #pdb.set_trace()
            prettyPrint(output, start, primer3)
            # do something with first set
            rp_end = start + int(primer3.primers[0]['PRIMER_RIGHT'].split(',')[0])
            rp_len = primer3.primers[0]['PRIMER_RIGHT'].split(',')[1]
        else:
            output.write("No primers found !!!!\n")
            rp_end = start + 100
            rp_len = '20'
        rp_start = int(rp_end) - int(rp_len)
        #if int(rp_len) >= 25:
        #    pdb.set_trace()
        rp_plus_buffer = rp_start - (25-int(rp_len))
        # need to maintain start and stop for slicing
        start = rp_plus_buffer - 30
        stop = start + 150
        # subtract start to 0-index
        target = ('%s,%s') % (rp_plus_buffer-start, 80)
        #pdb.set_trace()

def main():
    output = open('primersDesigned.txt','w')
    count = 0
    # read the sequence
    for record in SeqIO.parse(open('dog_mtdna_entire.fasta', 'r'), 'fasta'):
        design(record, output)
    output.close()

if __name__ == '__main__':
    main()

