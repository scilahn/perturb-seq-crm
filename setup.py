from setuptools import setup, find_packages

setup(
    name='crm',
    version='0.1.0',
    description='Cardiovascular/Renal/Metabolic Perturb-seq analysis library',
    author='Richard Ahn',
    author_email='sungho.richard@gmail.com',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    python_requires='>=3.11',
    install_requires=[
        'pertpy>=0.7.0',
        'scvi-tools>=1.1.0',
        'scanpy>=1.10.0',
        'anndata>=0.10.0',
        'decoupler>=1.6.0',
        'pydeseq2>=0.4.0',
        'boto3>=1.34.0',
    ],
)
