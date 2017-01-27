import re
from functools import partial
import warnings

import gzip
import pandas as pd

from .helpers import (
    make_chromosome_series_categorical,
    nan_to_None,
    dot_to_None,
)


GENO_REGEX = re.compile(r'([\d|\.](?:[/|\|][\d|\.])?)')


def vcf_to_dataframe(vcf_path, keep_samples=None, keep_format_data=False):
    """
    Generates a pandas.DataFrame of the variants present in a VCF. Removes
    the duplicated rows.

    To avoid high consumption of RAM and cycles, the default behavior is not to
    read any of the genotypes. In case you want to keep the genotypes, you have
    to set keep_samples explicitely with a sample ID or a list of sample IDs.

    If keep_format_data=False, it will only keep the genotypes (GT), with one
    column per sample. This option only makes sense if keep_samples is set.

    If keep_format_data=True, it will keep the metadata for each
    genotype call, e.g. AD, DP, GQ, etc. and for each variant. This means that
    for each variant there will now be many rows, one per sample, with all the
    genotype metadata (GQ, AD, DP) in new columns.
    """
    vcf_header = _header_from_vcf(vcf_path)

    if keep_samples:
        if isinstance(keep_samples, str) or isinstance(keep_samples, int):
            keep_samples = [str(keep_samples)]

        for sample in keep_samples:
            if sample not in vcf_header:
                raise ValueError('"{}" not found in this VCF'.format(sample))
    else:
        keep_samples = []

        if keep_format_data:
            warnings.warn('the option to "keep_format_data" only makes sense '
                          'if "keep_samples" is set. If not, the option is '
                          'just ignored, since there are no genotypes with '
                          'FORMAT data to display.')

    usecols = vcf_header[:9] + keep_samples  # 9 = number of standard VCF fields
    df = pd.read_table(vcf_path, comment='#', names=vcf_header,
                       low_memory=False, usecols=usecols)
    df = df.rename(columns={col: col.lower() for col in df.columns[:9]})

    if keep_samples:
        if keep_format_data:
            df = _unfold_genotype_data(df)
            # Cast genotypes as category since it's MUCH more memory efficient
            # In this case, genotypes will all be in the column 'GT':
            df['GT'] = df['GT'].astype('category')
        else:
            _extract_genos(df)
            # Cast genotypes as category since it's MUCH more memory efficient
            # In this case, there will be many genotype columns:
            df.iloc[:, 9:] = df.iloc[:, 9:].apply(
                lambda series: series.astype('category'))

    for col in 'ref filter'.split():
        df[col] = df[col].astype('category')

    df['chrom'] = make_chromosome_series_categorical(df['chrom'])

    # This drop_duplicates() must happen before the following parsings, or else
    # it will raise an Exception:
    df = df.drop_duplicates()

    df['info'] = df['info'].map(_info_line_as_dict)

    # Make the ALT alleles field always a list:
    df['alt'] = df['alt'].map(lambda alleles: alleles.split(','))

    integer_fields = [
        'DP',
        'GQ',
        'RGQ',
    ]
    for field in integer_fields:
        if field in df:
            df[field] = (df[field]
                         .map(nan_to_None).map(dot_to_None)
                         .map(int, na_action='ignore'))

    list_fields = [
        'AD',
        'PL',
    ]
    for field in list_fields:
        if field in df:
            df[field] = (df[field]
                         .map(nan_to_None).map(dot_to_None)
                         .map(lambda value: tuple(value.split(',')),
                              na_action='ignore'))

    return df


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


def _extract_genos(df):
    """Extract the genotypes from the format fields of the VCF. It changes the
    input DataFrame!"""

    def _extract_genotype(geno_field):
        """Extract the genotype from a format field."""
        # Assume the genotype is the first format field and raise if it's not
        geno = geno_field.split(':')[0]
        if not GENO_REGEX.search(geno):
            raise ValueError('"{}" does not look like a genotype'.format(geno))
        return geno

    # Genotype columns range from the 10th to the last one
    df.iloc[:, 9:] = df.iloc[:, 9:].applymap(_extract_genotype)


def _info_line_as_dict(info_line):
    """
    Convert one info line from a VCF to a dictionary. Example:
    > "DP=20;AN=16;EX_TARGET" => {'DP': 20, 'AN': 16, 'EX_TARGET': None}
    """
    info_dict = {}

    for element in info_line.split(';'):
        try:
            key, value = element.split('=')
        except ValueError:
            key = element
            value = None

        info_dict[key] = value

    return info_dict


def _unfold_genotype_data(df):
    """
    Given a dataframe with the standard VCF fields and genotypes columns,
    unfold each genotype column's data (GT, DP, GQ, etc.) adding a new column
    for each metadata category.

    The result will be a tidy dataframe with one row per combination of sample
    and variant.
    """
    sample_ids = list(df.columns[9:])
    # ^ The sample IDs come after the 9 VCF standard fields

    def parse_genotypes(variant, sample_id):
        """
        Make a dict out of the genotype data for one sample & variant.
        Example: ['GT:GQ', '1/1,99'] => {'GT': '1/1', 'GQ': 99}
        """
        genotype_fields = variant['format'].split(':')
        genotype_data = variant[sample_id].split(':')
        return dict(zip(genotype_fields, genotype_data))

    sample_frames = []
    for sample_id in sample_ids:
        genotypes = df.apply(partial(parse_genotypes, sample_id=sample_id), axis=1)
        genotypes = pd.DataFrame(list(genotypes))
        genotypes['sample_id'] = sample_id
        sample_df = pd.concat([df.drop(['format'] + sample_ids, axis=1), genotypes],
                              axis=1)
        sample_frames.append(sample_df)

    unfolded_df = pd.concat(sample_frames, ignore_index=True)
    return unfolded_df.sort_values(['chrom', 'pos']).reset_index(drop=True)

