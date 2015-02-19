from distutils.core import setup, Extension
import os, glob, numpy



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

__version__ = '3.4.1'

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
    author = 'Jeff Zheng, Eric Yang, Aaron Parsons, Shana Tribiano, Adrian Liu',
    author_email = '',
    url = '',
    package_dir = {'omnical':'src'},
    packages = ['omnical'],
    ext_modules = [
        Extension('omnical._omnical',
            globdir('src/_omnical/',
                ['*.cpp', '*.c', '*.cc']),
            include_dirs = ['src/_omnical/include', numpy.get_include()],
            extra_compile_args=['-Wno-write-strings', '-O3']
        )
    ],
    scripts = ['scripts/omnical_PSA128.py', 'scripts/omnical2npz.py', 'scripts/plot_omni.py', 'scripts/first_cal.py', 'scripts/bury_treasure.py', 'scripts/apply_omnigain.py'],
)

