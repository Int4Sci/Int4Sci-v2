lines(0);

function mygrid(xa,ya)
  a=get("current_axes")//get the handle of the newly created axes
  a.axes_visible="on"; // makes the axes visible
  a.font_size=3; //set the tics label font size
  a.x_location="top"; //set the x axis position
  a.data_bounds=[xa,ya,[0,0]];
  a.isoview= 'on'
  a.grid=[2,2];
  a.box="off";
endfunction
function mycircles(c)
  xset('color',c);
  xarc(-7,7,14,14,0,360*64)
  xarc(-6,6,12,12,0,360*64)
  xarc(-7+5,7+5,14,14,0,360*64)
  xarc(-6+5,6+5,12,12,0,360*64)
endfunction
function my2dbox(v,c)
xset('color',c);
xfrect(inf(v(1)),sup(v(2)),width(v(1)),width(v(2)));
endfunction

// Newton Schem
function f=F(z)
f=[(z(1))^2+(z(2))^2-#(6,7)^2 ;
   (z(1)-5)^2+(z(2)-5)^2-#(6,7)^2];
endfunction
function M=J(z)
   M=[z(1)+mid(z(1))  z(2)+mid(z(2)) ;
      z(1)+mid(z(1))-10 z(2)+mid(z(2))-10];
endfunction


function res=test(z)                                                 
  r=#(6,7)^2;                                                        
  x1=(z(1))^2;                                                       
  y1=(z(2))^2;                                                       
  x2=(z(1)-5)^2;                                                     
  y2=(z(2)-5)^2;                                                     
  c1=x1+y1;                                                          
  c2=x2+y2;                                                          
  if (in(c1,r) & in(c2,r)) | ~in(0,c1-r) | ~in(0,c2-r) then          
    res=[];                                                           
    return;                                                           
  end                                                                 
  res=z;  
endfunction    

// FILTER
//Newton

function t=Newton(z)
   t=  I4Slinearsolve(J(z),-F(mid(z)),"PGE")+mid(z);
endfunction
// Basic projections
function res=twoBfilter(z)
  r=#(6,7)^2;
  x1=(z(1))^2;
  y1=(z(2))^2;
  x2=(z(1)-5)^2;
  y2=(z(2)-5)^2;
  c1=x1+y1;
  c2=x2+y2;
  tmp=sqrt(r-z(2)^2);
  z(1)=intersection(z(1),#(-sup(tmp),sup(tmp)));
  tmp=sqrt(r-(z(1)-5)^2);
  z(2)=intersection(z(2),#(-sup(tmp),sup(tmp))+5);
  tmp=sqrt(r-(z(2)-5)^2);
  z(1)=intersection(z(1),#(-sup(tmp),sup(tmp))+5);
  tmp=sqrt(r-z(1)^2);
  z(2)=intersection(z(2),#(-sup(tmp),sup(tmp)));
  res=z;
endfunction    
    
    
    
    