import re

import gzip
import pandas as pd

from vcf_to_dataframe.helpers import make_chromosome_series_categorical


GENO_REGEX = re.compile(r'([\d|\.](?:[/|\|][\d|\.])?)')


def vcf_to_dataframe(vcf_path, keep_samples=None, keep_format_data=False):
    """
    Generates a pandas.DataFrame of the variants present in a VCF. Removes
    the duplicated rows.

    To avoid high consumption of RAM and cycles, the default behavior is not to
    read any of the genotypes. In case you want to keep the genotypes, you have
    to set keep_samples explicitely with a sample ID or a list of sample IDs.

    If you set keep_samples and keep_format_data, it will keep the metadata for
    each genotype call, e.g. AD, DP, GQ, etc. If not, it will only keep the
    genotypes (GT).
    """
    header = _header_from_vcf(vcf_path)

    if keep_samples:
        if isinstance(keep_samples, str):
            keep_samples = [keep_samples]

        for sample in keep_samples:
            if sample not in header:
                raise ValueError('"{}" not found in this VCF'.format(sample))
    else:
        keep_samples = []

    usecols = header[:9] + keep_samples

    df = pd.read_table(vcf_path, comment='#', names=header, low_memory=False,
                       usecols=usecols)
    df = df.rename(columns={col: col.lower() for col in df.columns[:9]})

    if keep_samples and not keep_format_data:
        _extract_genos_and_make_them_categorical(df)

    for col in 'ref alt filter'.split():
        df[col] = df[col].astype('category')

    df['chrom'] = make_chromosome_series_categorical(df['chrom'])

    return df.drop_duplicates()


def _header_from_vcf(vcf_path):
    """Read the header from a VCF file and return the found field names. It
    will gunzip on the fly if the path ends with '.gz' or if gzipped is set."""
    fn_open = gzip.open if vcf_path.endswith('.gz') else open

    with fn_open(vcf_path, 'rb') as vcf_file:
        for line in vcf_file:
            if isinstance(line, bytes):
                line = line.decode('utf-8')
            if line.startswith('#CHROM'):
                return line.strip().replace('#CHROM', 'CHROM').split('\t')


def available_samples(vcf_path):
    """Return a list of the sample IDs present in a VCF header."""
    return _header_from_vcf(vcf_path)[9:]


def _extract_genos_and_make_them_categorical(df):
    """Extract the genotypes from the format fields of the VCF and make them
    categorical. It changes the input DataFrame!"""
    # Genotype columns range from the 10th to the last one
    df.iloc[:, 9:] = df.iloc[:, 9:].applymap(_extract_genotype)

    # Cast genotypes as category since it's MUCH more memory efficient
    df.iloc[:, 9:] = df.iloc[:, 9:].apply(
            lambda series: series.astype('category'))


def _extract_genotype(geno_field):
    """Extract the genotype from a format field."""
    # Assume the genotype is the first format field and raise if it's not
    geno = geno_field.split(':')[0]
    if not GENO_REGEX.search(geno):
        raise ValueError('"{}" does not look like a genotype'.format(geno))
    return geno

