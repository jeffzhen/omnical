from distutils.core import setup, Extension

module = Extension('calibration_omni_extension',
                    define_macros = [('MAJOR_VERSION', '0'),
                                     ('MINOR_VERSION', '1')],
                    include_dirs = ['/usr/local/include','.'],
                    libraries = [],
                    library_dirs = ['/usr/local/lib','/usr/local/lib64'],
                    sources = ['calibration_omni_extension.cc','calibration_omni.cc'])

setup (name = 'calibration_omni_extension',
       version = '0.1',
       description = 'This is a python extension for calibration_omni package',
       author = 'Jeff Haoxuan Zheng, Eric Yang, MITEoR group @ MIT',
       author_email = 'jeff_z@mit.edu',
       url = 'github.com/jeffzhen/omnical',
       long_description = '''
This is a really long description.
''',
       ext_modules = [module])
