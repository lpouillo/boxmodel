boxmodel
========

A box model written in Python.

The current script allow to compute the evolution of a multibox system,
with different box mass, isotopic ratio, partition coefficient and flux. 

It was initially developped for the Fe ratio modeling in a human body.
A detailed explanation of the model can be found in the PhD manuscript of
Klervia Jaouen (http://tel.archives-ouvertes.fr/tel-00781645).

Numerical development is an adaptation of methods described in:
http://perso.ens-lyon.fr/francis.albarede/geochemodel.html


=Requirements=
- python 	2.7
- numpy 	1.7.1
- scipy		0.12
- execo		2.2

=Usage=
The proposed model is a simple simulation of the Fe ratio evolution. 
It can be run using:
execo-run BoxModel -ML


