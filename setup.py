import setuptools

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

setuptools.setup(
    version='0.1.8',
    name='cpiextract',
    author='Michael Sebek, Andrea Piras, Lucía Prieto Santamaría, Shi Chenghao, Gordana Ispirova, Giulia Menichetti',
    author_email='giulia.menichetti@channing.harvard.edu',
    description='CPIExtract is a software package to collect and harmonize small molecule and protein interactions.',
    keywords='data-science, bioinformatics, cheminformatics, proteins, network-science, data-harmonization, binding-affinity, chemical-compounds',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/menicgiulia/CPIExtract',
    packages=setuptools.find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.8',
    install_requires=[
        'pandas>=1.1.0; python_version<"3.8"',
        'pandas>=1.3.0; python_version>="3.8"',
        'numpy>=1.21.0,<2.0.0',
        'mysql-connector-python>=8.3.0',
        'biomart>=0.9.2',
        'chembl_webresource_client>=0.10.8',
        'pubchempy>=1.0.5'
    ],
    extras_require={
        'interactive': ['jupyter'],
        'performance': ["duckdb>=0.9.0"]
    }
)
