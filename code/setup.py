import os
from os.path import join as pjoin
import glob
from distutils.core import setup, Extension
import numpy.distutils.misc_util
import subprocess
from subprocess import Popen

# Locate prefix based on Python path
ins_prefix = subprocess.check_output(['which python'], shell=True )
ins_prefix = ins_prefix.split('/')[:-2]

prefix = "/"
for p in ins_prefix:
    prefix = os.path.join( prefix, p )
    
    print 'installing in', prefix

extensions = []
# Add extension source files.
BEML_source_files = glob.glob( './extensions/src/*c' )
# Add external library path
BEML_lib_dirs = ['/opt/conda/anaconda/envs/bempp/lib/']
BEML_lib_dirs.append( os.path.join( prefix, 'lib' ) )
# Add include directories
BEML_inc_dirs = ['/opt/conda/anaconda/envs/bempp/include']
BEML_inc_dirs.append( os.path.join( prefix, 'include' ) )
BEML_inc_dirs.append( numpy.distutils.misc_util.get_numpy_include_dirs()[0] )
BEML_inc_dirs.append( './extensions/include' )

# Setup BEML Extension
_BEML =  Extension(   "_beml",
                            BEML_source_files ,
                            language = "c",
                            library_dirs = BEML_lib_dirs ,
                            include_dirs = BEML_inc_dirs ,
                            libraries=["gsl", "gslcblas"],
                            extra_compile_args=["-fopenmp", '-std=c99', '-lm', '-fPIC', '-Wall'] , )

extensions.append( _BEML )

# *****************************************************************************                               
# Disutils setup specs                                                                                        
# *****************************************************************************                               
setup(                                                                                                        
    name = 'beml',
    package_dir =  { 'beml'      : 'pythonsrc' },
    packages = ['beml',
                'beml.beml'],
    author  = 'Maria Ignacia Fierro Piccardo',                                                                        
    author_email = 'ignaciapiccardo@gmail.com',                                                                    
    version = '0.1',                                                                                          
    ext_modules = extensions,                                                                                 
    )                                                                                                         
# *****************************************************************************

