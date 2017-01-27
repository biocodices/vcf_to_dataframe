from os.path import dirname, realpath, join
import re
import gzip

import pytest
import pandas as pd

from vcf_to_dataframe import vcf_to_dataframe, available_samples
from vcf_to_dataframe.helpers import (
    nan_to_None,
    dot_to_None,
)


TEST_DIR = dirname(realpath(__file__))


SAMPLE_IDS = ['HG00096', 'HG00097']
TEST_PARAMS = [
    ('sample.vcf.gz', True),
    ('sample.vcf', False)
]


@pytest.mark.parametrize('filename,gzipped', TEST_PARAMS)
def test_basic_functionality(filename, gzipped):
    df = vcf_to_dataframe(_test_file(filename))
    assert len(df['chrom'].cat.categories) == 5

    # Genotypes *not* read by default:
    assert all(sample_id not in df.columns for sample_id in SAMPLE_IDS)

    # Read the file "manually" to compare the data
    rows = read_file(_test_file(filename), gzipped)
    records = [re.split(r'\s+', row) for row in rows]
    seen_ids = set(record[2] for record in records)
    assert seen_ids == set(df['id'].unique())

    # Check the INFO field is read as a dictionary
    first_info_dict = df['info'].iloc[0]
    assert isinstance(first_info_dict, dict)
    assert first_info_dict['AN'] == '5008'


@pytest.mark.parametrize('filename,gzipped', TEST_PARAMS)
def test_genotypes_are_read_ok(filename, gzipped):
    sample_id = SAMPLE_IDS[0]

    df = vcf_to_dataframe(_test_file(filename), keep_samples=sample_id)
    assert sample_id in df.columns

    # Check category dtypes
    for col in ['chrom', 'ref', 'filter', sample_id]:
        assert df[sample_id].dtype == 'category'


@pytest.mark.parametrize('filename,gzipped', TEST_PARAMS)
def test_error_if_sample_not_present(filename, gzipped):
    with pytest.raises(ValueError):
        vcf_to_dataframe(_test_file(filename), keep_samples='non_existent')


@pytest.mark.parametrize('filename,gzipped', TEST_PARAMS)
def test_genotype_metadata_are_read_ok(filename, gzipped):
    df = vcf_to_dataframe(_test_file(filename), keep_samples=SAMPLE_IDS,
                          keep_format_data=True)
    assert 'GT' in df.columns
    assert all(sample_id not in df.columns for sample_id in SAMPLE_IDS)

    # Check there's one row per variant & sample
    variant_id = 'rs870124'
    variant_rows = df[df['id'] == variant_id]
    assert len(variant_rows) == len(SAMPLE_IDS)
    assert set(variant_rows['sample_id']) == set(SAMPLE_IDS)

    # Check one random genotype to make sure we're reading it correctly
    heterozygous_variant = df.loc[(df['id'] == variant_id) &
                                  (df['sample_id'] == 'HG00097')]
    assert heterozygous_variant.iloc[0]['GT'] == '0|1'


@pytest.mark.parametrize('filename,gzipped', TEST_PARAMS)
def test_available_samples(filename, gzipped):
    vcf_path = _test_file(filename)

    found_samples = available_samples(vcf_path)

    # Read the file "manually" to compare the data
    header = [line for line in read_file(vcf_path, gzipped, keep_header=True)
              if line.startswith('#CHROM')][0]
    expected_samples = header.split('\tFORMAT\t')[-1].split('\t')

    assert found_samples == expected_samples


def read_file(path, gzipped, keep_header=False):
    if gzipped:
        with gzip.open(path) as f:
            lines = [line.decode('utf-8').strip() for line in f.readlines()]
    else:
        with open(path) as f:
            lines = [line.strip() for line in f.readlines()]

    if not keep_header:
        lines = [line for line in lines if not line.startswith('#')]

    return lines


def test_dot_to_None():
    series = pd.Series(['foo', 'bar', '.', 'baz', '.'])
    result = series.map(dot_to_None)
    assert list(result) == ['foo', 'bar', None, 'baz', None]


def test_nan_to_None():
    series = pd.Series([1.0, 2.0, None])  # Nones converted to NaN by pandas
    series.loc[0] = '1'  # Now it's an 'object' series
    result = series.map(nan_to_None)
    assert list(result) == ['1', 2.0, None]


def _test_file(filename):
    return join(TEST_DIR, 'files', filename)

