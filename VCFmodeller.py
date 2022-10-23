#!/usr/bin/env python3

import re
import argparse
import random
import atexit
import gzip


class VirtualVCF:
    def __init__(self, file_name, n=1000):
        self._num = n
        self._title = '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] +
                                ['S' + str(n) for n in range(1, self._num + 1)])
        self._name = file_name
        self._file = None
        atexit.register(self.cleanup)

    def set_header(self, header='##fileformat=VCFv4.2'):
        header = header.strip()
        try:
            self._file = open(self._name, 'w')
        except IOError:
            print("Can't open file: {}".format(self._name))

        if not re.match(r'#CHROM', header):
            header += '\n' + self._title
        print(header, file=self._file)

        return header

    def add_variation(self, chrom, pos, var_id, ref, alt, quality, flt, info, af=0, homo=1):
        haplo_bank = []  # a set of sample haplotypes
        if not self._file:
            raise IOError("ERROR: VCF header was not set for file {}".format(self._name))

        if af < 0 or af > 1:
            raise ValueError("ERROR: Allele frequency (AF) must be in the range: 0 <= AF <= 1")

        if homo == 1:
            n_homo_1 = round(self._num * af)  # number of homozygous samples with alternative allele
            n_homo_0 = self._num - n_homo_1  # number of homozygous samples with reference allele

            for n in range(n_homo_1):
                haplo_bank.append(Haplotype(1, 1))
            for n in range(n_homo_0):
                haplo_bank.append(Haplotype(0, 0))
        elif 0 <= homo < 1:
            n_1 = round(2 * self._num * af)
            n_0 = 2 * self._num - n_1

            n_hetero = self._num - round(homo * self._num)  # number of heterozygous samples
            n_homo_1 = round((n_1 - n_hetero) / 2)
            n_homo_0 = round((n_0 - n_hetero) / 2)
            if n_homo_1 < 0 or n_homo_0 < 0:
                raise ValueError("ERROR: wrong homozygosity rate (HR): {}".format(homo))
            n_hetero = self._num - n_homo_1 - n_homo_0

            for n in range(n_hetero):
                haplo_bank.append(Haplotype(0, 1))
            for n in range(n_homo_1):
                haplo_bank.append(Haplotype(1, 1))
            for n in range(n_homo_0):
                haplo_bank.append(Haplotype(0, 0))
        else:
            raise ValueError("ERROR: homozygosity rate (HR) must be in the range: 0 <= HR <= 1")

        random.shuffle(haplo_bank)
        record = '\t'.join([
            str(chrom),
            str(pos),
            str(var_id),
            ref,
            alt,
            str(quality),
            str(flt),
            str(info),
            Haplotype.format,
            '\t'.join(hp.get() for hp in haplo_bank)
        ])
        print(record, file=self._file)

        return self.check(record)

    def check(self, rec):  # self check-up function
        n = 0
        n_alt = 0
        n_homo = 0
        for match in re.finditer(r'\t([01])/([01])', rec):
            h1, h2 = int(match.group(1)), int(match.group(2))
            n_alt += h1 + h2
            if h1 == h2:
                n_homo += 1
            n += 2

        af = round(n_alt / n, 4)
        hr = round(2 * n_homo / n, 4)

        if self._num != n / 2:
            raise ValueError("ERROR: wrong number of samples in the output!")

        return af, hr

    def cleanup(self):
        try:
            self._file.close()
        except IOError:
            print("Can't close file {}".format(self._name))

# end of class VirtualVCF


class Haplotype:
    format = 'GT'

    def __init__(self, hap_1=0, hap_2=0):
        self.gt = str(hap_1) + '/' + str(hap_2)

    def set(self, hap_1=0, hap_2=0):
        self.gt = str(hap_1) + '/' + str(hap_2)
        return self.gt

    def get(self):
        return self.gt

# end of class Haplotype


class HR:
    def __init__(self, af):
        if af < 0 or af > 1:
            raise ValueError("ERROR: wrong AF value: {}".format(af))
        self.af = af
        self.max_hr = self.get_max()
        self.optimal_hr = self.get_optimal()

    def get_max(self):
        if 0 <= self.af <= 0.5:
            return 2 * (1 - self.af)

        return 2 * self.af - 1

    def get_optimal(self):
        # hetero = 1 - self.af ** 2 - (self.af - 1) ** 2
        homo = self.af ** 2 + (self.af - 1) ** 2
        return homo

    def get(self, hr=None):
        if hr is None:
            return self.optimal_hr
        elif hr > self.max_hr:
            return self.max_hr
        elif hr < 0:
            raise ValueError("ERROR: wrong HR value: {}".format(af))

        return hr

# end of class HR


def main():
    parser = argparse.ArgumentParser(description=
                                     "A program to make a VCF file with virtual haplotypes \
                                     based on Allele Frequency (AF) and Homozygosity Rate (HR)")

    parser.add_argument('-i', help='input gzipped decomposed VCF file (INFO must have `AF` and `HR` values)')
    parser.add_argument('-o', default='output.vcf.gz', help='output gzipped VCF file with virtual haplotypes')
    parser.add_argument('-n', default='1000', help='number of virtual samples')

    args = parser.parse_args()

    in_file = args.i
    out_file = args.o
    n_samples = args.n

    virtual_vcf = VirtualVCF(out_file, n_samples)
    with gzip.open(in_file, 'wt') as f_in:
        header = ''
        for line in f_in:
            if re.match(r'#', line):
                header += line
            else:
                virtual_vcf.set_header(header)
                break

        for line in f_in:
            fields = re.split(r'\t', line)
            if len(fields) > 7:
                chrom, pos, var_id, ref, alt, quality, flt, info = fields[0:8]

                found = re.search(r'[\t; ]AF=([0-9.]+)', info)
                af = float(found.group(1)) if found else 0
                homo = HR(af)

                found = re.search(r'[\t; ]HR=([0-9.]+)', info)
                hr = float(found.group(1)) if found else None
                hr = homo.get(hr)

                af_chk, hr_chk = virtual_vcf.add_variation(chrom, pos, var_id, ref, alt, quality, flt, info, af, hr)
                print('{}:{} => AF={} ({}), HR={} ({})'.format(chrom, pos, af, af_chk, hr, hr_chk))

    print('...done')

# end of main()


if __name__ == '__main__':
    main()

    af = 0.4
    hr = 0.80
    chrom = 'chr1'
    pos = 123456

    var = VirtualVCF('virtual.vcf', 100)
    var.set_header()
    af_chk, hr_chk = var.add_variation(chrom, pos, 'rs123456', 'T', 'C', '.', 'PASS', 'info=HaHa;', af, hr)

    print('{}:{} => AF={} ({}), HR={} ({})'.format(chrom, pos, af, af_chk, hr, hr_chk))

