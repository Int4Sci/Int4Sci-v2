

libname = "libI4S";
makename = "Makelib";

table = ["scioperator","intoperator";"scipower","intpower" ; "sciunary", ..
         "intunary";"sciroundmode","introundmode";"sciILSR","intILSR"; ..
         "sciI4Svarsend","intI4Svarsend"; ..
	 "sciI4Svarget","intI4Svarget"];

objs = ["intoperator.c","intILSR.c"];

	if getos() == "Windows" then
    //** ----------- Windows section  -----------------
		CURRENT_PATH = strsubst(CURRENT_PATH,'\','/');
		CFLAGS = " -I"""+CURRENT_PATH+"/includes"""+ ..
			" -I"""+CURRENT_PATH+"/Profil/include"""+ ..
			" -I"""+CURRENT_PATH+"/Profil/include/BIAS""";

		LDFLAGS = " -L"""+CURRENT_PATH+"/Profil/lib"""+..
			" -L"""+CURRENT_PATH+"/src/c""" + ..
			" -L"""+CURRENT_PATH+"/src/cpp"""+ ..
			" -lstdc++ -lProfilPackages"+..
			" -lProfil -lBias -lm ";
	else
    //** ---------- Linux/MacOS/Unix section ---------------------
		CFLAGS = " -I"+CURRENT_PATH+"/includes "+ ..
			" -I"+CURRENT_PATH+"/Profil/include "+ ..
			" -I"+CURRENT_PATH+"/Profil/include/BIAS "+ ..
			"-D __USE_DEPRECATED_STACK_FUNCTIONS__";

		LDFLAGS = " -L"+CURRENT_PATH+"/Profil/lib "+..
			"-L"+CURRENT_PATH+"/src/c " + ..
			"-L"+CURRENT_PATH+"/src/cpp "+ ..
			" -lstdc++ -lProfilPackages"+..
			" -lProfil -lBias -lm ";
	end
	
tbx_build_gateway (TOOLBOX_NAME,table,objs,get_absolute_file_path('builder_gateway_c.sce'),['../../src/c/libscioperator','../../src/cpp/libsciILSR'],LDFLAGS,CFLAGS);

clear tbx_build_gateway;   

