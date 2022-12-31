#!/usr/bin/env python3

import re
import sys
import argparse
import random
import atexit
import gzip


class VirtualVCF:
    def __init__(self, file_name, n=1000):
        self._num = n
        self._title = '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] +
                                ['V' + str(n) for n in range(1, self._num + 1)])
        self._name = file_name
        self._file = None
        atexit.register(self.cleanup)

    def set_header(self, header='##fileformat=VCFv4.2'):
        header = header.strip()
        try:
            self._file = open(self._name, 'w')
        except IOError:
            print("Can't open file: {}".format(self._name))

        if not re.search(r'#CHROM', header):
            header += '\n' + self._title
        print(header, file=self._file)

        return header

    def add_variation(self, chrom, pos, var_id, ref, alt, quality, flt, haplo_bank):
        if not self._file:
            raise IOError("ERROR: VCF header was not set for file {}".format(self._name))

        an = len(haplo_bank) * 2
        ac_set = {}
        for k, hap in enumerate(haplo_bank):
            allele1, allele2 = [int(x) for x in re.split(r'[/|]', hap)]

            if allele1 in ac_set.keys():
                ac_set[allele1] += 1
            else:
                ac_set[allele1] = 1

            if allele2 in ac_set.keys():
                ac_set[allele2] += 1
            else:
                ac_set[allele2] = 1

        ac = []
        af = []
        hr = []
        num_samples = len(haplo_bank)
        accuracy = len(str(num_samples))
        hr_set = [0] * (max(ac_set.keys()) + 1)
        maf = 0
        for i in range(1, len(hr_set)):
            if i not in ac_set.keys():
                raise ValueError("ERROR - wrong number of allele counts.")
            af.append(str(round(ac_set[i] / an, accuracy)))
            ac.append(str(ac_set[i]))
            if i == 1:
                maf = af[0]

            for h in haplo_bank:
                a1, a2 = re.split(r'[|/]', h)
                if int(a1) == i and int(a1) == int(a2):
                    hr_set[i] += 1
        for h in hr_set[1:]:
            hr.append(str(round(h / num_samples, accuracy)))

        info = 'AC={};AN={};AF={};HR={};'.format(','.join(ac), str(an), ','.join(af), ','.join(hr))

        record = '\t'.join([
            str(chrom),
            str(pos),
            str(var_id),
            ref,
            alt,
            str(quality),
            str(flt),
            info,
            'GT',
            '\t'.join(haplo_bank)
        ])
        print(record, file=self._file)

        return maf

    def cleanup(self):
        try:
            self._file.close()
        except IOError:
            print("Can't close file {}".format(self._name))

# end of class VirtualVCF


class HaploRevolver:
    def __init__(self, num_of_samples, hr_type=True):
        self._n_samples = num_of_samples
        self._n_alleles_total = num_of_samples * 2
        self._n_alleles_ref = self._n_alleles_total
        self._revolver = []
        self._haplotypes = []
        self._hr_type = hr_type
        self._num_alt = 0

    def add_alt_allele(self, af, hr=None):
        num = self._num_alt + 1
        n_alleles = round(self._n_alleles_total * af)
        if not n_alleles:
            return False

        self._n_alleles_ref -= n_alleles
        if self._n_alleles_ref < 0:
            self._n_alleles_ref = 0

        if self._hr_type:
            if hr is None:
                raise ValueError("ERROR - no HR value, set 'hr_type=False' before.")

            n_homo = round(self._n_samples * hr)
            if n_alleles < n_homo * 2:
                n_homo = int(n_alleles / 2)

            if n_alleles - n_homo * 2 < 0:
                raise ValueError("ERROR - wrong ALT allele homozygosity rate!")
            self._haplotypes += [str(num) + '/' + str(num)] * n_homo
            if n_alleles - n_homo * 2 > 0:  # important!
                self._revolver.append([num] * (n_alleles - n_homo * 2))
        else:
            if n_alleles > 0:  # important!
                self._revolver.append([num] * n_alleles)

        self._num_alt += 1

        return True

    def get(self):
        if self._n_alleles_ref > 0:  # important!
            self._revolver.append([0] * self._n_alleles_ref)
        empty_bullets = set()
        revolver_size = len(self._revolver)
        while len(empty_bullets) < revolver_size:
            allele = []
            exclude = set()
            for n in range(2):
                if len(list(set(range(revolver_size)) - empty_bullets - exclude)) < 1:
                    break
                k = random.choice(list(set(range(revolver_size)) - empty_bullets - exclude))
                if len(self._revolver[k]):
                    allele.append(str(self._revolver[k][0]))
                    self._revolver[k].pop()
                    if self._hr_type and len(empty_bullets) + 1 < revolver_size:
                        # exclude if only one allele in the revolver or no HR
                        exclude.add(k)
                else:
                    raise ValueError("ERROR - no bullet in the revolver!")

                if not len(self._revolver[k]):
                    empty_bullets.add(k)

            if len(allele) == 2:
                if len(self._haplotypes) < self._n_samples:
                    self._haplotypes.append('/'.join(allele))

        random.shuffle(self._haplotypes)
        if len(self._haplotypes) != self._n_samples:
            raise ValueError(
                "ERROR - wrong number of haplotipes: {} (mast be {})".format(len(self._haplotypes), self._n_samples))

        return self._haplotypes

# end of class HaploRevolver


def main():
    parser = argparse.ArgumentParser(description=
                                     "A program to make a VCF file with virtual haplotypes \
                                     based on Allele Frequency (AF) and Homozygosity Rate (HR)")

    parser.add_argument('-i', help='input gzipped decomposed VCF file (INFO must have `AF` and `HR` values)')
    parser.add_argument('-o', default='output.vcf', help='output VCF file with virtual haplotypes')
    parser.add_argument('-n', default='1000', help='number of virtual samples')

    args = parser.parse_args()

    in_file = args.i
    out_file = args.o
    n_samples = int(args.n)

    virtual_vcf = VirtualVCF(out_file, n_samples)
    with gzip.open(in_file, 'rt') as f_in:
        header = ''
        for line in f_in:
            if re.match(r'#CHROM', line):
                if not re.search(r'^<ID=GT,', header):
                    header += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                if not re.search(r'<ID=HR', header):
                    header += '##INFO=<ID=HR,Number=A,Type=Float,Description="Homozigous Rate">\n'
                virtual_vcf.set_header(header)
                break
            elif re.match(r'#', line):
                header += line

        for line in f_in:
            fields = re.split(r'\t', line)
            if len(fields) > 7:
                chrom, pos, var_id, ref, alt, quality, flt, info = fields[0:8]

                if int(pos) == 15388335:
                    pass

                found = re.search(r'[\t; ]AF=([^;]+);', info)
                af_str = found.group(1) if found else ''
                af_set = re.split(r',', af_str)
                alt_set = re.split(r',', alt)

                if len(alt_set) != len(af_set):
                    print("ERROR - AF/Allele number mismatch:\n{}\n".format(line))
                    continue

                found = re.search(r'[\t; ]HR=([^;]+);', info)
                hr_str = found.group(1) if found else ''
                hr_set = re.split(r',', hr_str)

                hr_type = True
                if len(hr_set) > 0 and len(hr_set) != len(af_set):
                    print("ERROR - AF/HR number mismatch:\n{}\n".format(line))
                    hr_type = False
                    continue

                try:
                    revolver = HaploRevolver(n_samples, hr_type)
                    alt_alleles = []
                    for n in range(len(af_set)):
                        if hr_type:
                            if revolver.add_alt_allele(float(af_set[n]), float(hr_set[n])):
                                alt_alleles.append(alt_set[n])
                        else:
                            if revolver.add_alt_allele(float(af_set[n])):
                                alt_alleles.append(alt_set[n])

                    hamplo_set = revolver.get()
                    ptrn = re.search(r'AF=([^;,]+)[,;]', info)
                    real_maf = round(float(ptrn.group(1)), len(str(len(hamplo_set)))) if ptrn else 0
                    if len(alt_alleles):
                        maf = virtual_vcf.add_variation(chrom, pos, var_id, ref, ','.join(alt_alleles), quality, flt, hamplo_set)
                        print('{}\t{}\treal_AF= {}\tAF= {}'.format(chrom, pos, real_maf, maf))

                except ValueError:
                    print("ERROR in the record: {} {} {} {}".format(chrom, pos, var_id, ref), file=sys.stderr)

    print('...done')

# end of main()


if __name__ == '__main__':
    main()

