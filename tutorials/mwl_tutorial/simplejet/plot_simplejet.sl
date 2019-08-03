xlog;ylog;
set_plot_widths(;m_width=5, r_width=5, re_width=5, d_width=5);

%plot model 
fancy_plot_unit("kev","mjy");
plot_unfold({pid,hid};dsym={4,4},dcol={4,8},decol={4,8},rcol={4,8},recol={4,8},xrange={2,200},res=1);

(x,y)=readcol(“Output/Presyn.dat",1,2);
connect_points(1); oplot(10^x,10^y,5);  % Overplot the second model component (5=cyan)

linestyle(2);
(x,y)=readcol(“Output/Postsyn.dat",1,2);
connect_points(1); oplot(10^x,10^y,10);  % 10=green/cyan

linestyle(1); 
(x,y)=readcol(“./Total.dat",1,2);
connect_points(1); oplot(10^x,10^y,15);   % 15=light grey
