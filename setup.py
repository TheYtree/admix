from setuptools import setup

setup(
    name='admix',

    version='1.3',

    description='An admixture analysis tool that supports raw data from 23andme, AncestryDNA, etc.',
    long_description=open('admix/README.rst','r').read(),

    url='https://github.com/TheYtree/admix',

    author='Steven Liu',
    author_email='me@yliu.io',

    license='GNU General Public License v3.0',

    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Natural Language :: Chinese (Simplified)',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],

    keywords='bio DNA SNP ancestry admixture position-based',

    install_requires=['numpy','scipy'],

    packages=['admix'],

    package_data={
        'admix':['data/*', 'data/position_models/*', 'README.rst']
    },

    entry_points={
        'console_scripts': ['admix=admix.admix:main'],
    },
)
