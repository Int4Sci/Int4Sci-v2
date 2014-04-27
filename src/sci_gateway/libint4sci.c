#include <mex.h> 
#include <sci_gateway.h>
#include <api_scilab.h>
static int direct_gateway(char *fname,void F(void)) { F();return 0;};
extern Gatefunc intoperator;
extern Gatefunc intpower;
extern Gatefunc intunary;
extern Gatefunc intI4Svarsend;
extern Gatefunc intI4Svarget;
static GenericTable Tab[]={
  {(Myinterfun)sci_gateway,intoperator,"scioperator"},
  {(Myinterfun)sci_gateway,intpower,"scipower"},
  {(Myinterfun)sci_gateway,intunary,"sciunary"},
  {(Myinterfun)sci_gateway,intI4Svarsend,"sciI4Svarsend"},
  {(Myinterfun)sci_gateway,intI4Svarget,"sciI4Svarget"},
};
 
int C2F(libint4sci)()
{
  Rhs = Max(0, Rhs);
  if (*(Tab[Fin-1].f) != NULL) 
  {
     pvApiCtx->pstName = (char*)Tab[Fin-1].name;
    (*(Tab[Fin-1].f))(Tab[Fin-1].name,Tab[Fin-1].F);
  }
  return 0;
}
