"""
"""
import setuptools

setuptools.setup(
    name='gaston',
    version='v0.1.1',
    python_requires='>=3.8',
    packages=['gaston'],
    package_dir={'': 'src'},
    author='Uthsav Chitra',
    author_email='uchitra@princeton.edu',
    description='GASTON: interpretable deep learning model of gene expression topography',
    url='https://github.com/raphael-group/GASTON',
    install_requires=[
        'torch',
        'matplotlib',
        'numpy',
        'pandas',
        'seaborn',
        'scikit-learn',
        'tqdm',
        'scipy',
        'jupyterlab',
        'glmpca',
        'scanpy'
    ],
    include_package_data = True,
    package_data = {
        '' : ['*.txt']
        },
    license='BSD',
    platforms=["Linux", "MacOs", "Windows"],
    classifiers=[
        'Programming Language :: Python :: 3.8',
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords=[
        'spatial transcriptomics',
        'neural network',
        'tissue layer'],
    entry_points={'console_scripts': 'gaston=gaston.__main__:main'}
)
