#! /usr/bin/env python
import sys
import argparse
import screed
import khmer
import pickle


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('reads')
    parser.add_argument('-k', '--ksize', type=int, default=31)
    args = parser.parse_args()

    if args.reads == '-':
        args.reads = sys.stdin

    kh = khmer.Nodetable(args.ksize, 1, 1)

    print('loading database {}'.format(args.database))
    with open(args.database, 'rb') as fp:
        family_ids, kmer_to_family_id = pickle.load(fp)

    print('done!')

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
        first_ids = kmer_to_family_id.get(first, set())
        idx = 1
        while idx < len(hashvals)/2 and len(first_ids) != 1:
            first = hashvals[idx]
            idx += 1

        # find the last unambiguously assigned k-mer
        last_ids = kmer_to_family_id.get(last, set())
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
