# Code detecting Motifs with signed uni-directed links

 A motif is a pattern which observed higher/lower than expected in random networks with conserved the degree properties.
We calculate the z-score to compare the frequencies of motifs in networks consisting of different nodes and links.
The z-score indicates the number of standard deviations from the mean a data point is.

<p align="center">
  <img src="images/motifs.png" width="80%">
</p>

 Three-node motifs should have two or three links to connect all nodes. 
There are five types of motifs considering only the direction of each link because we regard the mirror and rotational symmetry. 
However, we can obtain 22 types of motifs by considering both the direction and sign. 

 Our codes are based on C language with gcc compiler.
