# Material_Cavity
This is a Python Script which aids in the development of a given crystallized structure with an XYZ file of a molecule inserted in a parameterized elliptical cavity within the crystal.

# Motivations
We are hoping to utilize this script in conjunction with a research project at the University of Arizona to use Frozen Density Embedded Theory (FDET) calculations for certain molecules within a crystallized structure.


# How to use this script
This python script will use the publicly-available ASE library to create a crystallized structure, with the default being diamond. It will then create an elliptical cavity within the structure itself, with a default setting but can be parameterized by the user themselves.

The first step is to clone the code yourself.

If you have an .xyz file (we use water in our example), please put it in the same directory as the cloned project.

I will show an example script, in bash, of how to run the code, assuming the code is called "test2.py", and with an example water.xyz code (assuming a water.xyz file is in your directory. This file is available on the web).

To activate the script, you will need to go to the folder in which it is available, and run the script (using Bash) as follows, in the console/terminal:

"python3 test2.py --xyz butynal.xyz --element C --crystalstructure diamond --a 3.57 --size 5 5 5 --ellipse 4.5 4.5 5"

(Note: Dependent on the Python version downloaded on your local machine/cluster, you can also use "python" instead of "python3").

The ASE documentation library is well-documented and open-source, and encouraged to be scrolled through here, in case you have an idea to optimize the code for your own reasons:

https://wiki.fysik.dtu.dk/ase/




