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

    keep_samples='All' (case insensitive) is treated specially: all samples
    will be kept. Don't use this option with big VCFs, like 1000 Genomes that
    has 2,504 samples.

    If keep_format_data=False, it will only keep the genotypes (GT), with one
    column per sample. This option only makes sense if keep_samples is set.

    If keep_format_data=True, it will keep the metadata for each
    genotype call, e.g. AD, DP, GQ, etc. and for each variant. This means that
    for each variant there will now be many rows, one per sample, with all the
    genotype metadata (GQ, AD, DP) in new columns.
    """
    vcf_header = _header_from_vcf(vcf_path)

    if isinstance(keep_samples, str) and keep_samples.lower() == 'all':
        keep_samples = available_samples(vcf_path)

    keep_samples = _parse_samples(keep_samples, vcf_header)

    if not keep_samples and keep_format_data:
        warnings.warn('the option to "keep_format_data" only makes sense '
                      'if "keep_samples" is set. If not, the option is '
                      'just ignored, since there are no genotypes with '
                      'FORMAT data to display.')

    usecols = vcf_header[:9] + keep_samples  # 9 = number of standard VCF fields
    skiprows = _count_comment_rows(vcf_path)

    df = pd.read_table(vcf_path, names=vcf_header, low_memory=False,
                       usecols=usecols, skiprows=skiprows)
    df = df.rename(columns={col: col.lower() for col in df.columns[:9]})

    if keep_samples:
        if keep_format_data:
            df = _unfold_genotype_data(df)
            if 'GT' in df:
                # Cast genotypes as category since it's MUCH more memory
                # efficient:
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

    integer_list_fields = [
        'AD',
        'PL',
    ]
    for field in integer_list_fields:
        if field in df:
            df[field] = \
                (df[field].map(nan_to_None).map(dot_to_None)
                 .map(lambda value: tuple(int(n) for n in value.split(',')),
                      na_action='ignore'))

    return df


def _header_from_vcf(vcf_path):
    """Read the header from a VCF file and return the found field names. It
    will gunzip on the fly if the path ends with '.gz' or if gzipped is set."""
    vcf_lines_generator = lines_from_vcf(vcf_path)

    for line in vcf_lines_generator:
        if line.startswith('#CHROM'):
            vcf_lines_generator.close()  # Don't leave the file handle opened
            return line.strip().replace('#CHROM', 'CHROM').split('\t')


def lines_from_vcf(vcf_path):
    """Generator. Yield line by line of a given (possibly gzipped) VCF."""
    fn_open = gzip.open if vcf_path.endswith('.gz') else open

    with fn_open(vcf_path) as vcf_file:
        for line in vcf_file:
            if isinstance(line, bytes):
                line = line.decode('utf-8')

            yield line


def _parse_samples(sample_list, vcf_header):
    """
    Given the fields in *vcf_header*, check the samples in *sample_list*
    and raise if any doesn't exist in the header. Also, stringify the
    sample IDs.
    """
    if sample_list:
        if isinstance(sample_list, str) or isinstance(sample_list, int):
            sample_list = [str(sample_list)]

        for sample in sample_list:
            if sample not in vcf_header:
                raise ValueError('"{}" not found in this VCF'.format(sample))
    else:
        sample_list = []

    return sample_list


def _count_comment_rows(vcf_path):
    """Count how many comment rows sit on top of the passed VCF file."""
    vcf_lines_generator = lines_from_vcf(vcf_path)

    comment_lines_count = 0
    for line in vcf_lines_generator:
        if line.startswith('#'):
            comment_lines_count += 1
        else:
            vcf_lines_generator.close()  # Don't leave the file handle opened
            # Don't continue reading the VCF once the comments section ended
            break

    return comment_lines_count


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

