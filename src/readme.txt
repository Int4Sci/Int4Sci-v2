Basic Use
--------

When receiving the Int4Sci archive file, extract it. You will have a Int4Sci directory.
Go to the Int4Sci directory and open scilab.

Under scilab, if you have a source distribution, type :
    --> exec builder.sce; 
    --> exec int4sciloader.sce;


After the first use, you will need only to load Int4Sci by :
    --> exec int4sciloader.sce;

To run the tests for Int4sci
   --> I4Stest();

For automatic use from any directory, see Install.htm file.

Look at the Int4Sci_Tutorial.htm file for detailed use.

Scilab on line help :
-------------------

Scilab on line documentation for Int4Sci is available.

Just type :
	--> apropos interval
to have an instant access to the whole documentation.

There are seven help files than you can consult from your scilab
window :

	--> help interval
for documentation on the basic functions of Int4Sci.

	--> help I4S_Basic_Func
for documentation on classical functions dealing with intervals
(intersection, inclusion, radius, ...) .

	--> help I4S_Functions
for documentation on Scilab functions applied to intervals
(norms, matrix reshaping, size, ...).

	--> help I4S_Arith_Op
for documentation on interval arithmetic operations available with Int4Sci.

	--> help I4S_Arith_Func
for documentation on elementary functions applied to intervals (trigonometric,
hyperbolic, ...).

	--> help I4Slinearsolve
for documentation about interval linear systems solving.

	--> help I4S_Error
for documentation about Int4Sci errors.

Installation only tested with :
-------------------------------

Scilab 5.1.1 binary files downloaded ( http://www.scilab.org/download/index_download.php)
and compiled on Linux (fedora core 5 and 7).


More installation tests still to come !

License :
---------

This toolbox is under the GPL license (http://www.gnu.org/copyleft/gpl.html) .

See the LICENSE file for more details.

Bugs :
------
Watch the RELEASE_NOTE file for more details.
 
Authors
-------
Raphael PEREIRA
David DANEY
Jean Pierre MERLET
 INRIA Sophia Antipolis
 COPRIN project-team

Contacts and technical support
------------------------------

int4sci@lists-sop.inria.fr



