
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

mygrid([-30,30],[-30,30])

function my2dbox(v)
xfrect(inf(v(1)),sup(v(2)),width(v(1)),width(v(2)));
endfunction
