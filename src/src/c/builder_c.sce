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



src_c_path = get_absolute_file_path('builder_c.sce');                           
                                                                                
CFLAGS = "-I" + src_c_path + " -I"+CURRENT_PATH+"/includes "+ " -I"+CURRENT_PATH+"/Profil/include " ;                                                          
                                                                                
LDFLAGS = " -L"+CURRENT_PATH+"/Profil/lib "+...   
    " -lstdc++ -lProfilPackages -lProfil -llr -lBias";                           
LDFLAGS = "-g"+ LDFLAGS +" -lm";                                           
                                                                                
                                                                                
tbx_build_src(['scioperator','scipower','sciunary'], ['Biasoperator.c'], 'c', ...
src_c_path, '', LDFLAGS, CFLAGS);                                 
                                                                                
clear tbx_build_src;                                                            
clear src_c_path;                                                               
clear LDFLAGS                                                                   
clear CFLAGS;   
