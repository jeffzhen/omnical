from distutils.core import setup, Extension


#module = Extension('calibration_omni_extension',
                    #define_macros = [('MAJOR_VERSION', '0'),
                                     #('MINOR_VERSION', '1')],
                    #include_dirs = ['/usr/local/include','.'],
                    #libraries = [],
                    #library_dirs = ['/usr/local/lib','/usr/local/lib64'],
                    #sources = ['calibration_omni_extension.cc','calibration_omni.cc'])

#setup (name = 'calibration_omni_extension',
       #version = '0.1',
       #description = 'This is a python extension for calibration_omni package',
       #author = 'Jeff Haoxuan Zheng, Eric Yang, MITEoR group @ MIT',
       #author_email = 'jeff_z@mit.edu',
       #url = 'github.com/jeffzhen/omnical',
       #long_description = '''
#This is a really long description.
#''',
       #ext_modules = [module])

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

