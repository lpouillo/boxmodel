BoxModel
========

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
- python 	2.7
- numpy 	1.7.1
- scipy		0.12
- matplotlib	1.2.1
- pydot		1.0.28
- execo		2.1


Installing required packages
Debian

    sudo apt-get install git graphviz python-setuptools libatlas3-base libatlas-dev libblas3 libblas-dev liblapack3 liblapack-dev build-essential gfortran libfreetype6-dev libpng1.2-dev

Ubuntu

    sudo apt-get install git graphviz python-setuptools libatlas3gf-base libatlas-dev libblas3gf libblas-dev liblapack3gf liblapack-dev build-essential gfortran libfreetype6-dev libpng12-dev
      
And then
    
    sudo easy_install numpy==1.7.1 scipy==0.12 matplotlib==1.2.1 pydot
    wget http://execo.gforge.inria.fr/downloads/execo-2.1.tar.gz && tar xzf execo-2.1.tar.gz && cd execo-2.1/ && sudo make install

Finally clone the repository

    git clone https://github.com/lpouillo/boxmodel.git
    
Usage
-----
The BoxModel is a basic engine that can be derived to build specific element model. It does not work by itself.
You can run the demo by typing:

    execo-run FeSimple -ML 	        # Evolution of the iron ratio
    execo-run ZnCalibration -ML		# Evolution of the zinc ratio
    execo-run ZnDietRatio -ML		# Evolution of the zinc ratio pour un espace de paramètre donné 

Demo
----
The proposed model is a simple simulation of the Fe ratio evolution in a human body. 
It can be run using:
	execo-run FeModel -ML

Building customized model
-------------------------
