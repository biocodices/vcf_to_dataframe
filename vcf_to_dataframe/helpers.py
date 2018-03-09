from pandas.api.types import CategoricalDtype
import numpy as np


def make_chromosome_series_categorical(series):
    """
    Receives a pandas Series with chromosomes and returns a new categorical
    Series with the same values.
    """
    new_series = series.astype(str)  # Creates a new copy of the series

    if all(new_series.str.startswith('chr')):
        new_series = new_series.str.replace('chr', '')

    sorted_chromosomes = [str(chrom)
                          for chrom in list(range(1, 23)) + ['X', 'Y']
                          if str(chrom) in new_series.values]
    chromosome_dtype = CategoricalDtype(categories=sorted_chromosomes,
                                        ordered=True)

    return (new_series
            .replace('23', 'X')
            .replace('24', 'Y')
            .astype(chromosome_dtype))


def nan_to_None(value):
    """Given an arbitrary value, replace NaN with None or return unchanged."""
    try:
        is_nan = np.isnan(value)
    except TypeError:
        is_nan = False

    return None if is_nan else value


def dot_to_None(value):
    """Given an arbitrary value, replace '.' with None or return unchanged."""
    return None if value == '.' else value

