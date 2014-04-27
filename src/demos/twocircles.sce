
exec ../loader.sce;


function mygrid(xa,ya)
a=get("current_axes")//get the handle of the newly created axes
a.axes_visible="on"; // makes the axes visible
a.font_size=3; //set the tics label font size
a.x_location="top"; //set the x axis position
a.data_bounds=[xa,ya,[0,0]];
 //set the boundary values for the x, y and z coordinates.
//a.sub_tics=[tics,tics];
//a.labels_font_color=5;
a.isoview= 'on'
a.grid=[2,2];
a.box="off";
endfunction


function my2dbox(v,c)
xset('color',c);
xfrect(inf(v(1)),sup(v(2)),width(v(1)),width(v(2)));
endfunction




function f=F(z)
f=[(z(1))^2+(z(2))^2-#(6,7)^2 ;
   (z(1)-5)^2+(z(2)-5)^2-#(6,7)^2];
endfunction

function M=J(z)
M=[z(1)+mid(z(1))-2*0  z(2)+mid(z(2))-2*0 ; 
   z(1)+mid(z(1))-2*5 z(2)+mid(z(2))-2*5];
endfunction

function r=innertest(z) 
 r=prod(in([(z(1)-0)^2+(z(2)-0)^2 ; (z(1)-5)^2+(z(2)-5)^2],[#(5,9)^2;#(5,9)^2]));
endfunction

function r=outertest(z) 
 r=~prod(in([0;0],F(z)));
endfunction

//function t=testnewton(z)
//	 t=  I4Slinearsolve(J(z),-F(mid(z)),z-mid(z),"PGS")+mid(z);
//endfunction


function t=testnewton(z)
	 t=  I4Slinearsolve(J(z),-F(mid(z)))+mid(z);

	 if t==[] then
	    t=z;
	 else
	    t=intersection(t+mid(z),z);
	    //if strictin(t,z) then
	    //   t=testnewton(z);
	    //end	    
	 end
endfunction


function z=likea2B(z)
	 z(1)=intersection(z(1),hull(-sqrt(#(5,9)^2-z(2)^2),sqrt(#(5,9)^2-z(2)^2)));
	 if isnan(z(1)) then
	    r=[];
	    return;
	 end
	 z(2)=intersection(z(2),hull(-sqrt(#(5,9)^2-z(1)^2),sqrt(#(5,9)^2-z(1)^2)));
	 if isnan(z(2)) then
	    r=[];
	    return;
	 end
	 z(1)=intersection(z(1),hull(-sqrt(#(5,9)^2-(z(2)-5)^2)+5,sqrt(#(5,9)^2-(z(2)-5)^2)+5));
	 if isnan(z(1)) then
	    r=[];
	    return;
	 end
	 z(2)=intersection(z(2),hull(-sqrt(#(5,9)^2-(z(1)-5)^2)+5,sqrt(#(5,9)^2-(z(1)-5)^2)+5));
	 if isnan(z(2)) then
	    r=[];
	    return;
	 end
endfunction




clf;
lines(0)

mygrid([-5,3],[2,10])

xarc(-9,9,18,18,0,360*64)
xarc(-5,5,10,10,0,360*64) 

xarc(-9+5,9+5,18,18,0,360*64)
xarc(-5+5,5+5,10,10,0,360*64) 


Z=[#(-5,3) #(2,10)]';


Zlst=[Z];

while length(Zlst) > 0
  Ztmp = Zlst(:,1);
  Zlst= Zlst(:,2:$);

  c=0;

  if outertest(Ztmp) then
     my2dbox(Ztmp,2);
  elseif innertest(Ztmp) then
     my2dbox(Ztmp,5);
  elseif prod(width(Ztmp) < [0.5;0.5]) then
     my2dbox(Ztmp,7);
  else
     //Ztt=testnewton(Ztmp);
     Ztt=likea2B(Ztmp);  
     if Ztt == Ztmp then
       t=Ztmp:2;
       Zlst=[Zlst t(1) t(2)];
     else
       Zlst=[Zlst Ztt];	
     end
  end

  
end





