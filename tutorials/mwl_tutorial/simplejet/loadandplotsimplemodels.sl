%add a simple absorbed broken powerlaw model, grid data, and plot just Xray
%exclude(rid);
%usr_grid([rid,pid],-9,3,0.001,1);
%usr_grid([hid],0,3,0.001);
%fit_fun("constant(Isis_Active_Dataset)*phabs(1)*bknpower(1)");
%set_par("phabs(1).nH",0.6,1); % freeze column density to well known value
%set_par("constant(1).factor",1,1); % freeze detector constant of radio and soft X-rays to 1, hard X_ray dataset free to vary
%fit_counts; % should give chi^2 ~ 830 at 85 degrees of freedom (#bins - #free paramaters)
%fancy_plot_unit("kev","mjy");
%plot_unfold({pid,hid};dsym={4,4},dcol={4,8},decol={4,8},rcol={4,8},recol={4,8},xrange={2,200},res=1);

%%  add some other phenomenological components
%fit_fun("constant(Isis_Active_Dataset)*phabs(1)*highecut(1)*(bknpower(1)+diskbb(1)+ gaussian(1))");
%%list_par;
%fit_counts;
%plot_unfold({pid,hid};dsym={4,4},dcol={4,8},decol={4,8},rcol={4,8},recol={4,8},xrange={2,200},res=1);

%% downward curvature of residuals above ~100keV: cutoff naturally explained by a physical continuum due to
%% Compton upscattering by hot electron plasma. "nthcomp" good model: https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node203.html
%fit_fun("constant(Isis_Active_Dataset)*phabs(1)*(gaussian(1)+nthcomp(1))");
%%list_par;  
%fit_counts;
%plot_unfold({pid,hid};dsym={4,4},dcol={4,8},decol={4,8},rcol={4,8},recol={4,8},xrange={2,200},res=1);


%%add a reflection convolution model
%fit_fun("constant(Isis_Active_Dataset)*phabs(1)*reflect(1,(nthcomp(1)+gaussian(1)))");
%fit_counts; 
%plot_unfold({pid,hid};dsym={4,4},dcol={4,8},decol={4,8},rcol={4,8},recol={4,8},xrange={2,200},res=1);


%now add back in radio and fit with simple bknpowerlaw again
%include(rid);
%set_par("constant(3).factor",1,1); % radio and IR data shall have same normalization as dataset #1
%tie(18,1);
%fancy_plot_unit("hz","mjy");

% just evaluate model that we found when fitting X-rays first
%eval_counts;
%plot_unfold({rid,pid,hid};dsym={4,4,4},dcol={6,4,8}, decol={6,4,8},xrange={3e8,1e20},res=1);

%%go back to simple bknpowerlaw
fit_fun("constant(Isis_Active_Dataset)*phabs(1)*highecut(1)*(gaussian(1)+ bknpower(1))");
%list_par;
set_par(2,.6,1);  %fix absorption to measured value
set_par(3,75); %set cutoff energy
set_par(8,1.e4); %give starting normalisation 
set_par(9,0.75); %inverted PL to fit radio spectrum
set_par(10,0.001,0,1.e-5,1.e5); %give starting value/reset break range
set_par(11,1.6);
renorm_counts;
fit_counts;
plot_unfold({rid,pid,hid};dsym={4,4,4},dcol={6,4,8},decol={6,4,8},xrange={3e8,1e20},res=1);
