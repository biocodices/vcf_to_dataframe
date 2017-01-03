from setuptools import setup, find_packages


with open('requirements.txt') as f:
    requirements = f.read().split('\n')

with open('README.md') as f:
    long_description = f.read()

setup(
    name='vcf_to_dataframe',
    version='1.1.0',
    description='Get a pandas DataFrame from a possibly gzipped VCF file.',
    author='Juan Manuel Berros',
    author_email='juanma.berros@gmail.com',
    url='https://github.com/biocodices/vcf_to_dataframe',
    license='MIT',
    install_requires=requirements,
    packages=find_packages(),
    long_description=long_description,
)
