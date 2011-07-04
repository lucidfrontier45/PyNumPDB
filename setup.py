#!/usr/bin/python

#from distutils.core import setup
from numpy.distutils.core import setup, Extension

ext = Extension(name="pynumpdb._anm",sources=["ext_module/_anm.f90",])

setup(name='pynumpdb',
      version='0.2',
      description='Python Library for Treating Protein Data Bank format',
      author='Shiqiao',
      author_email='lucidfrontier.45@gmail.com',
      url='http://frontier45.web.fc2.com/',
      packages=['pynumpdb'],
      ext_modules = [ext,]
     )
