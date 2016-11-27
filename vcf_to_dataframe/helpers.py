def make_chromosome_series_categorical(series):
    """
    Receives a pandas Series with chromosomes and returns a new categorical
    Series with the same values.
    """
    new_series = series.astype(str)  # Creates a new copy of the series
    if all(new_series.str.startswith('chr')):
        new_series = new_series.str.replace('chr', '')

    new_series = new_series.replace('23', 'X').replace('24', 'Y').astype(str)
    chromosomes = [str(chrom) for chrom in list(range(1, 23)) + ['X', 'Y']
                   if str(chrom) in new_series.values]
    new_series = new_series.astype('category', categories=chromosomes,
                                   ordered=True)
    return new_series

