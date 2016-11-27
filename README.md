# VCF to DataFrame

Get a pandas DataFrame from a possibly gzipped VCF file. Keep only the variants
or also the genotypes from the selected samples.

Usage:

```python
from vcf_to_dataframe import vcf_to_dataframe

df = vcf_to_dataframe('data.vcf')
# => DataFrame with the variants in the data.vcf file

df = vcf_to_dataframe('data.vcf.gz')
# => DataFrame with the variants in the *gzipped* data.vcf.gz file

df = vcf_to_dataframe('data.vcf.gz', keep_samples='HG00096')
# => Same, but now it keeps the genotypes of the selected sample(s)

df = vcf_to_dataframe('data.vcf.gz', keep_samples=['HG00096', 'HG00097'],
                      keep_format_data=True)
# => Same, but now you have the genotypes and the each call metadata,
#    like AD, DP, GQ, and whatever there is in the VCF
```

