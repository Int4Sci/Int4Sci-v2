//C Builder file
//------------

//    <Int4Sci, a scilab interface for Interval Analysis.>
//    Copyright (C) <2010>  <D. DANEY, B. NEVEU>

//    This program is free software: you can redistribute it and/or
//    modify
//    it under the terms of the GNU General Public License as published
//    by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see
//   <http://www.gnu.org/licenses/>.


//   Contact : D. Daney <int4sci@lists-sop.inria.fr>
//   COPRIN team
//   INRIA Sophia Antipolis
//   2004 route des lucioles - BP 93
//   06902 Sophia Antipolis Cedex
//   FRANCE

// Developer only executions



src_cpp_path = get_absolute_file_path('builder_cpp.sce');

	if getos() == "Windows" then
    //** ----------- Windows section  -----------------
		src_cpp_path = strsubst(src_cpp_path,'\','/');
		CURRENT_PATH = strsubst(CURRENT_PATH,'\','/');

		CFLAGS = "-I""" + src_cpp_path + ..
			""" -I"""+CURRENT_PATH+"/includes""" + ..
			" -I"""+CURRENT_PATH+"/Profil/include""" + ..
			" -I"""+CURRENT_PATH+"/Profil/include/lr""";

		LDFLAGS = "-g "+ ..                  
			" -L"""+CURRENT_PATH+"/Profil/lib"""+...
			" -lProfilPackages -lProfil -llr -lBias -lm";

	else
    //** ---------- Linux/MacOS/Unix section ---------------------
		CFLAGS = "-I" + src_cpp_path + ..
			" -I"+CURRENT_PATH+"/includes " + ..
			" -I"+CURRENT_PATH+"/Profil/include " + ..
			" -I"+CURRENT_PATH+"/Profil/include/lr ";

		LDFLAGS = "-g "+ ..                  
			" -L"+CURRENT_PATH+"/Profil/lib "+...
			" -lProfilPackages -lProfil -llr -lBias -lm";
	end
	
tbx_build_src([ 'sciILSR', 'sciI4Svarsend', 'sciI4Svarget'], [ 'ProfilILSR.cpp'], 'cpp', src_cpp_path, '', LDFLAGS, CFLAGS);

clear tbx_build_src;
 clear src_c_path;
clear LDFLAGS
clear CFLAGS;
