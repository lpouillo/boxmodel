IsotopicBoxModel
================

A box model written in Python.
Numerical development is an adaptation of methods described in:
http://perso.ens-lyon.fr/francis.albarede/geochemodel.html

The current engine allow to compute the evolution of a multibox system,
with different box mass, isotopic ratio, partition coefficient and flux. 

It was initially developped for the Fe ratio modeling in a human body.
A detailed explanation of the model can be found in the PhD manuscript of
Klervia Jaouen (http://tel.archives-ouvertes.fr/tel-00781645).


Installation
------------
Packages requirements
- python 			2.7
- numpy 			1.7.1
- scipy				0.12
- matplotlib		1.2.1
- pydot				1.0.28
- execo				2.4.3


Installing packages on Debian/Ubuntu, with root privileges:

    apt-get install graphviz python-setuptools python-graphviz python-scipy \
    python-numpy python-matplotlib python-networkx
      
And then
    
    easy_install execo
    

Finally clone the repository

    git clone https://github.com/lpouillo/boxmodel.git
    
Usage
-----
The IsotopicBoxModel is a basic engine that can be derived to build specific element model. It does not work by itself.
You can run the demos located in src/ directory by typing:

    # Evolution of the iron ratio
    python FeSimple.py -ML
    # Evolution of the zinc ratio
    python ZnCalibration.py  -ML
    # Evolution of the zinc ratio for a given parameter range (Diet ratio variation)
    python ZnDietRatio.py  -ML



