#! /usr/bin/env python
import sys
import argparse
import screed
import khmer
import pickle
import shelve
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('transcriptomes', nargs='+')
    parser.add_argument('-k', '--ksize', type=int, default=31)
    parser.add_argument('-o', '--output')
    args = parser.parse_args()

    assert args.output

    kh = khmer.Nodetable(args.ksize, 1, 1)

    family_ids = {}
    family_counter = 0

    kmer_to_family_id = defaultdict(set)

    n = 0
    for tr_filename in args.transcriptomes:
        for record in screed.open(tr_filename):
            n += 1
            if n % 1000 == 0:
                print('...', n)

            family_name = record.name.split('|')[1]

            family_id = family_ids.get(family_name)
            if family_id is None:
                family_id = family_counter
                family_counter += 1
                family_ids[family_name] = family_id
            
            hashvals = kh.get_kmer_hashes(record.sequence)

            for hashval in hashvals:
                kmer_to_family_id[hashval].add(family_id)

    with open(args.output, 'wb') as fp:
        pickle.dump((family_ids, kmer_to_family_id), fp)


if __name__ == '__main__':
    main()
