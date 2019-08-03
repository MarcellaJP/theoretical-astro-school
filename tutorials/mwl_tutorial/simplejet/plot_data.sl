xlog;ylog;
set_plot_widths(;m_width=5, r_width=5, re_width=5, d_width=5);

% plot just Xrays, folded with detector, with and without background
%(UNCOMMENT BOTH AT ONCE!)
%plot_data({pid,hid};dsym=4,dcol=8,decol=8, bkg=1, yrange={0.005,30.}); % plot data+bkg
%plot_data({pid,hid};dsym=4,dcol=4,decol=4, bkg=0, oplt=1, yrange={0.005,30.}); % plot data-bkg

%plot just Xrays, unfolded
plot_unfold({pid,hid};dsym={4,4},dcol={4,8}, decol={4,8},res=0);


%fancy_plot_unit(“kev”,”ergs”);
%plot_unfold({pid,hid};dsym={4,4},dcol={4,8}, decol={4,8},res=0);

%plot with radio/IR
%fancy_plot_unit("hz","mjy");
%plot_unfold({rid,pid,hid};dsym={4,4,4},dcol={6,4,8}, decol={6,4,8},xrange={3e8,1e20},res=0);

