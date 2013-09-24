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


=Requirements=
- python 	2.7
- numpy 	1.7.1
- scipy		0.12
- matplotlib	1.2.x
- pydot		1.0.x
- execo		2.1

=Usage=
The BoxModel is a basic engine that can be derived to build specific element model. The only thing you need is to create the simulation parameters in a parameters method and run the engine.

==Demo==
The proposed model is a simple simulation of the Fe ratio evolution in a human body. 
It can be run using:
	execo-run FeModel -ML

==Building customized model==

