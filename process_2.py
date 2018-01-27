#! /usr/bin/env python
import sys
import argparse
import screed
import khmer
import pickle
import numpy
import bbhash


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('reads')
    parser.add_argument('-k', '--ksize', type=int, default=31)
    args = parser.parse_args()

    if args.reads == '-':
        args.reads = sys.stdin

    kh = khmer.Nodetable(args.ksize, 1, 1)

    mphf_filename = args.database + '.mphf'
    array_filename = args.database + '.arr'
    print('loading database {}'.format(args.database))
    
    with open(array_filename, 'rb') as fp:
        mphf_to_kmer, mphf_to_cdbg, family_ids, cdbg_to_family_id = pickle.load(fp)
    mphf = bbhash.load_mphf(mphf_filename)

    print('done!')

    def get_kmer_to_family_ids(hashval):
        mphf_hash = mphf.lookup(hashval)
        if mphf_hash is None:
            return set()
        
        kmer_hash = mphf_to_kmer[mphf_hash]
        if kmer_hash != hashval:
            return set()

        cdbg_id = mphf_to_cdbg[mphf_hash]
        id_list = cdbg_to_family_id[cdbg_id]
        return id_list

    n_same = 0
    n_different = 0

    n = 0
    for record in screed.open(args.reads):
        n += 1
        if n % 1000 == 0:
            print('...', n)
            if n > 5000:
                break

        hashvals = kh.get_kmer_hashes(record.sequence)
        if len(hashvals) <= 1:
            continue

        first = hashvals[0]
        last = hashvals[-1]

        # find the first unambiguously assigned k-mer
        first_ids = get_kmer_to_family_ids(first)
        idx = 1
        while idx < len(hashvals)/2 and len(first_ids) != 1:
            first = hashvals[idx]
            idx += 1

        # find the last unambiguously assigned k-mer
        last_ids = get_kmer_to_family_ids(last)
        idx = len(hashvals) - 2
        while idx > len(hashvals) / 2 and len(last_ids) != 1:
            last = hashvals[idx]
            idx -= 1

        if len(first_ids) == 1 and len(last_ids) == 1 and \
           first_ids == last_ids:
            n_same += 1
        else:
            print('different {} {}'.format(first_ids, last_ids))
            n_different += 1

    print('same:', n_same)
    print('different:', n_different)


if __name__ == '__main__':
    main()
