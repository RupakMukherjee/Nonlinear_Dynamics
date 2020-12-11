This is a PyClaw input file to solve acoustics in 1D inhomogeneous medium.
The user need to run only the file script.sh from terminal using a command ./script.sh
Make sure you have execute permission for the file script.sh
If you do not have the permission, you should change the permission by running this command within the local directory: chmod ugo+rwx script.sh
Once the run is completed, it will generate a file called "SlopePlot.vtr" which you can visualize using Paraview.
Also it will generate another file called "ParameterPlot.dat" which you can visualize via GNUplot or Python.
The raw data will be written in ASCII format in a file called "SlopePlot.dat" and one can use python or any other visualization tool to plot these data.
