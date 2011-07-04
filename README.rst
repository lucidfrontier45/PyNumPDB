.. -*- mode: rst -*-

About
=====

PyNumPDB is a python module for treating Protein Data Bank (PDB) file.

The module is designed especially for the aids of numerical calculation 
such as pca or clustering.

The project was started in 2010 by Shiqiao Du for the sake of education.

It is currently maintained by Du.


Download
========

You can download source code from my homepage:

http://sourceforge.net/projects/scikit-learn/files/


Dependencies
============

The required dependencies to build the software are python >= 2.5,
NumPy >= 1.1, SciPy and scikits.learn >= 0.4.

Optional dependencies are MODELLER and fortran compiler.


Install
=======

This packages uses distutils, which is the default way of installing
python modules. The install command is::

  python setup.py install


Bug report
============

If you want to ask something, feel free to mail to:

lucidfrontier.45@gmail.com

Any bug reports are welcome.

You can also ask in the forum of my homepage:
http://frontier45.web.fc2.com/index.html

Testing
-------

There are some sample scripts in "samples" and test pdb fiels in "samples/test_pdbs".
For example, you can calculate RMSD of decoyes to native structure in that folder::

  cd samples
  python pdbrmsd.py test_pdbs/* test_pdbs/native.pdb


