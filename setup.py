from distutils.core import setup, Extension
import os, glob

__version__ = '0.0.1'

def indir(dir, files): return [dir+f for f in files]
def globdir(dir, files):
    rv = []
    for f in files: rv += glob.glob(dir+f)
    return rv

setup(name = 'omnical',
    version = __version__,
    description = __doc__,
    long_description = __doc__,
    license = 'GPL',
    author = 'Jeff Zheng, Adrian Liu, Aaron Parsons',
    author_email = '',
    url = '',
    package_dir = {'omnical':'src'},
    packages = ['omnical'],
    ext_modules = [
        Extension('omnical._omnical',
            globdir('src/_omnical/',
                ['*.cpp', '*.c']),
            include_dirs = ['src/_omnical/include'],
        )
    ],
    scripts = glob.glob('scripts/*'),
)
