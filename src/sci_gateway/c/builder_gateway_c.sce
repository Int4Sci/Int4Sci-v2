

libname = "libI4S";
makename = "Makelib";

table = ["scioperator","intoperator";"scipower","intpower" ; "sciunary", ..
         "intunary";"sciILSR","intILSR"; ..
         "sciI4Svarsend","intI4Svarsend"; ..
	 "sciI4Svarget","intI4Svarget"];

objs = ["intoperator.c","intILSR.c"];

CFLAGS = " -I"+CURRENT_PATH+"/includes "+ ..
" -I"+CURRENT_PATH+"/Profil/include "+ ..
" -I"+CURRENT_PATH+"/Profil/include/BIAS ";

LDFLAGS = " -L"+CURRENT_PATH+"/Profil/lib "+..
	"-L"+CURRENT_PATH+"/src/c " + ..
        "-L"+CURRENT_PATH+"/src/cpp "+ ..
	" -lstdc++ -lProfilPackages"+..
	" -lProfil -lBias -lm ";

tbx_build_gateway (TOOLBOX_NAME,table,objs,get_absolute_file_path('builder_gateway_c.sce'),['../../src/c/libscioperator','../../src/cpp/libsciILSR'],LDFLAGS,CFLAGS);

clear tbx_build_gateway;   