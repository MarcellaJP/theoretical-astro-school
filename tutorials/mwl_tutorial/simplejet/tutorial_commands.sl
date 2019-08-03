%% before starting: get to know ISIS, how to find help, syntax, etc.
%% -> Remeis tutorial https://www.sternwarte.uni-erlangen.de/wiki/index.php/Isis_tutorial
%%    "A walk through ISIS"



%% load PCA and HEXTE (small recab from previous tutorials) %%

%%%%%%%%%%%%%%
%% load PCA %%
%%%%%%%%%%%%%%

variable datapath = "/data/mwl_tutorial/";

variable pid = load_data(datapath+"pca.pha");
set_sys_err_frac(pid,0.005);
group(pid;min_sn=4.5,bounds=3.,unit="kev");
% why min_sn=4.5? Chi2 stat per bin -> need min 20-25 counts/bin
% -> Poisson -> Chi2 stat, while Poisson uncertainty per bin 
% S/N = S/sqrt(S) = sqrt(S) -> sqrt(20)=4.5
xnotice_en(pid,3,22);

%% plot, learn about units, folding/unfolding
%% (-> see also Remeis Tutorials and Mike's famous talks I,II
%% http://www.black-hole.eu/media/summerschool2/X-ray_Spectra_Part_I.pdf
%% http://www.black-hole.eu/media/summerschool2/X-ray_Spectra_Part_II.pdf)

xlog;ylog;
set_plot_widths(;m_width=5, r_width=5, re_width=5, d_width=5);

% different ways to plot data
plot_counts(pid;dsym=4,dcol=8,decol=8); % Counts / bin
plot_data(pid;dsym=4,dcol=8,decol=8); % Counts / s / keV
% spectrum folded with detector features, detector space
%  preferred;
% source spectrum as seen by the detector / intrinsic spectrum folded by detector response, i.e.,
% effective area [how sensitive at what energy due to mirror effects etc]
% redistribution function [what channel belongs to which energy]
% Advantage: see detector features (bumps around absorption edges, less counts at edges of energy range (see shape effective area)

%% PLT: plot_tutorial_pca_plot-data.eps %%

plot_unfold(pid;dsym=4,dcol=8,decol=8); % Photons / cm2 / s / keV
% spectrum unfolded with detector response -> source-intrinsic spectrum
% Guess of intrinsic spectrum just as good as we know detector, unfolding not always 100% successful, sometimes still features at absorption edges of effective area
% good to get idea what model we want to use (spectral model defined in source space of course)

%% PLT: plot_tutorial_pca_plot-unfold.eps %%

% check out background contribution
plot_data(pid;dsym=4,dcol=8,decol=8,bkg=-1); % background only
%% PLT: plot_tutorial_pca_plot-data_bkg.eps 
plot_data(pid;dsym=4,dcol=8,decol=8,bkg=0); % background subtracted, default
%% PLT: plot_tutorial_pca_plot-data_nobkg.eps
plot_data(pid;dsym=4,dcol=4,decol=4,bkg=1,oplt=1); % overplot data without bkg + data with bkg 
%% PLT: plot_tutorial_pca_plot-data_data+bkg.eps


%%%%%%%%%%%%%%%%%
%% load HEXTE %%%
%%%%%%%%%%%%%%%%%

variable hid = load_data(datapath+"hxt.pha");
group(hid;min_sn=4.5,bounds=20.,unit="kev");
xnotice_en(hid,20,200);
plot_data(hid;dsym=4,dcol=8,decol=8, bkg=1, yrange={0.005,2.2}); % plot data+bkg
plot_data(hid;dsym=4,dcol=4,decol=4, bkg=0, oplt=1, yrange={0.005,2.2}); % plot data-bkg
%% PLT: plot_tutorial_hexte_plot-data_data-bkg_data+bkg.eps


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% work on both datasets %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_data({pid,hid};dsym={4,4},dcol={4,8}, decol={4,8}); % plot both in detector space / folded spectrum
%% PLT: plot_tutorial_pca_hexte_plot-data.eps
plot_unfold({pid,hid};dsym={4,4},dcol={4,8}, decol={4,8}); % source-intrinsic spectrum / unfolded spectrum
%% PLT: plot_tutorial_pca_hexte_plot-unfold.eps

% change units to energy flux:
fancy_plot_unit("kev","ergs");
plot_unfold({pid,hid};dsym={4,4},dcol={4,8}, decol={4,8});
%% PLT: plot_tutorial_pca_hexte_plot-unfold_nuFnu.eps


%%%%%%%%%%%%%%%%%
%% load radio %%%
%%%%%%%%%%%%%%%%%

%% add radio + IR
variable rid=load_radio2(datapath+"53082_radio-OIR.dat");
fancy_plot_unit("hz","mjy");
plot_unfold({rid,pid,hid};dsym={4,4,4},dcol={6,4,8}, decol={6,4,8},xrange={3e8,1e20});

%% define custom grid at which model evaluated. If not given, model always evaluated on ARF grid,
%% which is usually finer than the data grid
%% Even if data is grouped, will be evaluated on ARF grid and model then rebinned to match grouped data grid
%% This binned model is then folded with the response and can be directly fitted to the observed data
usr_grid([rid,pid],-9,3,0.001,1);
usr_grid([hid],0,3,0.001);


%%%%%%%%%%%%%%%%%%%
%% Fit the data %%%
%%%%%%%%%%%%%%%%%%%


%% define model, refine model step by step by looking at chi-residuals



%% fit X-ray only %%%%%%%%%%%%%%%%%%%%%%

exclude(rid);
fit_fun("constant(Isis_Active_Dataset)*phabs(1)*bknpower(1)");
set_par("phabs(1).nH",0.6,1); % freeze column density to well known value
set_par("constant(1).factor",1,1); % freeze detector constant of radio and soft X-rays to 1, hard X_ray dataset free to vary
fit_counts; % should give chi^2 ~ 830 at 85 degrees of freedom (#bins - #free paramaters)
fancy_plot_unit("kev","mjy");
plot_unfold({pid,hid};dsym={4,4},dcol={4,8},decol={4,8},rcol={4,8},recol={4,8},xrange={2,200},res=1);
%% PLT: plot_tutorial_pca_hexte_plot-unfold_bknpower.eps


%%  add some other phenomenological components
fit_fun("constant(Isis_Active_Dataset)*phabs(1)*highecut(1)*(bknpower(1)+diskbb(1)+ gaussian(1))");
plot_unfold({pid,hid};dsym={4,4},dcol={4,8},decol={4,8},rcol={4,8},recol={4,8},xrange={2,200},res=1);

%% downward curvature of residuals above ~100keV: cutoff naturally explained by a physical continuum due to
%% Compton upscattering by hot electron plasma. "nthcomp" good model: https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node203.html
fit_fun("constant(Isis_Active_Dataset)*phabs(1)*(gaussian(1)+nthcomp(1))");
list_par;  
fit_counts;
plot_unfold({pid,hid};dsym={4,4},dcol={4,8},decol={4,8},rcol={4,8},recol={4,8},xrange={2,200},res=1);
%% PLT: plot_tutorial_pca_hexte_plot-unfold_nthcomp.eps

%% slope in residuals >10keV, looks like Compton reprocessing, for example in accretion disk
%% - Broad Compton hump peaking around 20-30keV
%% - line feature around 6 keV (reprocessing in disk, flourescence line)
%%
%% Possible physical explanations: Reflection off
%% 1) inner disk: relativistic effects, use "relxill" family of reprocessing models linking the emissivities at each point of the disk
%%                with correct relativistic kernel (light bending, grav. redshift, rel. Doppler, naturally explains also rel. blurred
%%                K alpha line around 6-7 keV)
%%                http://www.sternwarte.uni-erlangen.de/~dauser/research/relxill/
%% 2) outer, cool disk: easier, relativistic effects not important, "reflect" good and simple convolution model without
%%                line features

fit_fun(“constant(Isis_Active_Dataset)*phabs(1)*reflect(1,(nthcomp(1)+gaussian(1)))”);
fit_counts; % well, flattens Compton hump a little bit, but not everything (would need to run more complicated fitting)
plot_unfold({pid,hid};dsym={4,4},dcol={4,8},decol={4,8},rcol={4,8},recol={4,8},xrange={2,200},res=1);

%% logical next step, if model "relxill" was installed:
%% account for remaining Compton hump emission and also broad Fe line feature around 6keV.
%% Both natural features of relativistic reflection, use state-of-the-art model relxill that also accounts for ionization of disk,
%% where insident radiation is getting reprocessed
%% http://www.sternwarte.uni-erlangen.de/~dauser/research/relxill/)
fit_fun("constant(Isis_Active_Dataset)*phabs(1)*relxill(1))");  % relxill includes nthcomp


%%%fitting algorithm
% the default fitting algorithm is lmdif, which is not the most robust but it's fast.  You can change to another one like subplex 
% which will take much, much longer but will explore parameter space more thoroughly, but do that at HOME! 


%% fit radio + IR + X-ray %%%%%%%%%%%%%%%%%%%%%%
include(rid);
set_par("constant(3).factor",1,1); % radio and IR data shall have same normalization as dataset #1
%this is one way to link parameters for a given fit:
%set_par_fun("constant(3).factor","constant(1).factor");
%this is nice because you can also define functional forms!  
%But you can also just tie the parameters, which we'll do 
tie(18,1);
fancy_plot_unit("hz","mjy");
% just evaluate model that we found when fitting X-rays first
eval_counts;
plot_unfold({rid,pid,hid};dsym={4,4,4},dcol={6,4,8}, decol={6,4,8},xrange={3e8,1e20},res=1);
%% PLT: plot_tutorial_radio_pca_hexte_plot-unfold_nthcomp.eps


%%will skip below in class:  
% doesn't really fit lower energies;
% right now we only see model plotted on data grid.
% overplot model evaluated on arbitrary, much finer grid independent of data binning
fancy_plot_unit("keV"); % switch back to keV on x-axis and photon flux on y axis
plot_unfold({rid,pid,hid};dsym={4,4,4},dcol={6,4,8}, decol={6,4,8},xrange={1e-9,1000}); % plot data and model
(lo,hi)=log_grid(1e-9,1e3,1e4); % define logarithmic grid
f=eval_fun_keV(lo,hi); % and evaluate model on it, see 'help eval_fun_keV' in which units function returns the model
ohplot(lo,hi,f/(hi-lo)); % overplot evaluated model over binned data in correct units 
% can see that photoabsorption model only defined down to UV,
% in optical would need to account for extinction / reddening (don't do this here)
% and in radio for Synchrotron Self Absorption (also not done here)

%% PLT: plot_tutorial_radio_pca_hexte_plot-unfold_nthcomp_finemodelgrid.eps


% model with inverse-Comptonized continuum (reflected off cold outer disk)
% fits X_rays perfectly (besides missing relativistic reflection component)
% but radio and IR not well modelled -> this is where jet kicks in!
% Can model broad-band spectrum in different ways:
% 1) phenomenological model (remember first fit of broken powerlaw)
% 2) combination of the model just found for the X-rays and
%    extra component for the jet-synchrotron and disk-thermal emission + companion star (IR)
% 3) physical jet model predicting all that self-consistently from radio - gamma rays (developed at API, but too complicated for now)
%    but: lacking relativistic effects that lead to broad Fe line at ~6keV and parts of Compton hump peaking at ~20-30 keV

% for the given time, try 1) and recall first fit, give some values to start off and save time, from 'eyeballing'
fit_fun(“constant(Isis_Active_Dataset)*phabs(1)*highecut(1)*(gaussian(1) + bknpower(1))”);
list_par;
set_par(2,.6,1);  %fix absorption to measured value
set_par(3,75); %set cutoff energy
set_par(8,1.e4); %give starting normalisation 
set_par(9,0.75); %inverted PL to fit radio spectrum
set_par(10,0.001,0,1.e-5,1.e5); %give starting value/reset break range
set_par(11,1.6);
eval_counts;  % just evaluate and see how it fits in the radio and IR (we made it fit the X-rays well before; ISIS remembers these parameters)
renorm_counts; fit_counts;
plot_unfold({rid,pid,hid};dsym={4,4,4},dcol={6,4,8},decol={6,4,8},xrange={3e8,1e20},res=1);
%% PLT: plot_tutorial_radio_pca_hexte_plot-unfold_bknpower.eps


% General guidelines:  tweak parameters to make entire continuum fit
% evaluate (eval_counts)
% if model seems to come close to data, run fit_counts

%%%save and exit Isis to slirp in a custom model
save_input(“mwfitting.sl”); %history!!
save_par(“mwfitting.par”);
quit;

%%%%Follow steps to "slirp" model in /data/mwl_tutorial directory
unix> slirp -make simplejet.f
unix> make
unix> make test

%%%% Make sure Isis paths set to find module
% for this tutorial add lines to .isisrc:

prepend_to_isis_load_path(“/home/ataschool/data/mwl_tutorial“); 
                             
prepend_to_isis_module_path(“/home/ataschool/data/mwl_tutorial“);                            

%%%Go back into isis and load the model, cache it to grid, evaluate and plot:
%%% have to copy in all the load/setup commands then:

()=evalfile(“simplejet.sl”);
loval = _A(10^[-9:3:0.001]);
hival = make_hi_grid (loval);
cache_fun("simplejet",loval,hival);
fit_fun(“simplejet_cache”);

eval_counts;

%%%use precooked plotting script that also overplots components
.load plotscript

%%%can also get started with some better parameters to get you started
load_par("sj_cache_start.par");
eval_counts;
.load plotscript





