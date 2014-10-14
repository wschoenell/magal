from distutils.core import setup

setup(
    name='magal',
    version='0.1_beta',
    packages=['magal', 'magal.io', 'magal.fit', 'magal.core', 'magal.util', 'magal.plots', 'magal.library', 'magal.photometry'],
    requires=['pystarlight', 'numpy', 'astropy', 'atpy', 'h5py', 'matplotlib'],
    scripts = ['src/scripts/magal_filterset',
               'src/scripts/magal_library',
               'src/scripts/magal_input',
               'src/scripts/magal_fit'
                   ],
    package_dir={'': 'src'},
    url='http://www.iaa.es/~william/magal/',
    license='3-clause BSD',
    author='William Schoenell',
    author_email='wschoenell@gmail.com',
    description='Magal is a tool to retrieve galaxy information out of photometry.'
)
