from os.path import dirname, realpath, join
import re
import gzip

import pytest

from vcf_to_dataframe import vcf_to_dataframe, available_samples


TEST_DIR = dirname(realpath(__file__))


TEST_PARAMS = [
    ('sample.vcf.gz', True),
    ('sample.vcf', False)
]


@pytest.mark.parametrize('filename,gzipped', TEST_PARAMS)
def test_vcf_to_dataframe(filename, gzipped):
    vcf_path = _test_file(filename)

    df = vcf_to_dataframe(vcf_path)
    sample_id = 'HG00096'
    assert len(df['chrom'].cat.categories) == 5
    assert sample_id not in df.columns

    # Read the file "manually" to compare the data
    rows = read_file(vcf_path, gzipped)

    records = [re.split(r'\s+', row) for row in rows]
    seen_ids = set(record[2] for record in records)

    assert seen_ids == set(df['id'].unique())

    df = vcf_to_dataframe(vcf_path, keep_samples=sample_id)
    assert sample_id in df.columns

    with pytest.raises(ValueError):
        vcf_to_dataframe(vcf_path, keep_samples='non existent')


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


def _test_file(filename):
    return join(TEST_DIR, 'files', filename)

