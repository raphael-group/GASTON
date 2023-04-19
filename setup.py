"""
"""
import setuptools

setuptools.setup(
    name='spatialnn',
    version='v0.0.1',
    python_requires='>=3.8',
    packages=['spatialnn'],
    package_dir={'': 'src'},
    author='Uthsav Chitra',
    author_email='uchitra@princeton.edu',
    description='SpatialNN does some things that are useful',
    url='https://github.com/raphael-group/SpatialNN',
    install_requires=[
        'torch',
        'matplotlib',
        'numpy',
        'pandas',
        'seaborn',
        'scikit-learn'
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
    entry_points={'console_scripts': 'spatialnn=spatialnn.__main__:main'}
)
