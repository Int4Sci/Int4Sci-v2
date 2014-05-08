Int4Sci-v2
==========

Int4Sci is interval analysis toolbox for Scilab


Basic Int4Sci binary versioon installation

1 Unpack the archive with the required version

2 Copy the folder "int4sci-v2" to the folder $SCIHOME/contrib/

3 Run Scilab.

Installation only tested with :

Scilab 5.4.1 binary files downloaded ( http://www.scilab.org/download/index_download.php) and compiled on Linux (Ubuntu 12.04).


Basic Int4Sci installation
--------------------------

1. You must download the source Profil-2.0.8 and apply to them the patch ( https://github.com/Int4Sci/Profil-2.0.8-patch )
2. Copy the folder "lib" and "include" from Profil-2.0.8 to Int4Sci_Folder/Profil
3. You should have installed gcc and g++ compiler on Linux or you shoud have MinGW compiler for Windows.
4. Copy the folder "int4sci-v2" to the folder $SCIHOME/contrib/
5. Run Scilab
6. Install Mingw Compiler for Scilab
7. Run 
	--> exec builder.sce;
8. To run the tests for Int4sci
   --> I4Stest();



Basic use
To get started with Int4Sci, see the Tutorial.

You can access Int4Sci on line documentation by typing in your Scilab window :

--> apropos interval;

Finally, a set of tests for Int4sci is available by typing in your Scilab window :

--> I4Stest();

