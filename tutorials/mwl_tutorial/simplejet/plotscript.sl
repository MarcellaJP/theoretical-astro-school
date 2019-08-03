clear;
xlog;
ylog;

fancy_plot_unit("hz","mjy");

variable x,y;
variable poptf=@popt;
poptf.power=2; %kev * photon flux
poptf.xrange={5.e8,1.e20};
poptf.yrange={1.e-3,60.,-10,10};

point_size(1.4);
charsize(1.12);
variable wd=4;
set_plot_widths(;m_width=5, d_width=5, r_width=5,de_width=3, re_width=3);
set_frame_line_width(wd); set_line_width(3);

multiplot(3,1);                % Traditional 2 panel plot, 3:1 ratio
mpane(1);                      % First panel
poptf.res=0;                   % popt.res=0, don't plot residuals

poptf.dsym={5,-4,4};
poptf.dcol={12,4,8}; 
poptf.decol={3,5,7};
poptf.rcol={12,4,8};
poptf.recol={3,5,7};
poptf.mcol={2,2,2}; 

plot_unfold({rid,pid,hid},poptf;no_reset=1);  % Same data plot as before

(x,y)=readcol("./presyn.dat",1,2);
connect_points(1); oplot(10^x,10^y,5);  % Overplot the second model component (5=cyan)

linestyle(2);
(x,y)=readcol("./postsyn.dat",1,2);
connect_points(1); oplot(10^x,10^y,10);  % 10=green/cyan

linestyle(1); 
(x,y)=readcol("./total.dat",1,2);
connect_points(1); oplot(10^x,10^y,15);   % 15=light grey
set_line_width(1);

mpane(2);                  % Second panel
poptf.res=1;
poptf.rsym=[5,-4,4];
poptf.rcol={12,4,8};
poptf.recol={3,5,7};


plot_residuals({rid,pid,hid},poptf; no_reset=0);




