
% Create a S-lang library of the function using SLIRP:

% unix%> slirp -make simplejet.f
% unix%> make 

% Move mffjet-module.so to your S-lang modules directory (for me,
% /usr/local/isis/modules).  In your ISIS .isisrc, or other S-lang
% input file, run all the code below.

import("simplejet");

define simplejet_fit(lo,hi,par)
{
   variable ear = [_A(hi),_A(lo[0])];
   variable nein = length(lo), photar=Float_Type[nein], photer=@photar;

   variable pars = par[[1:15]];

   % !!!  It is important to match the type for arrays that are input, !!!
   % !!!  and recieve return values.  I.e., phot & energ go in and     !!!
   % !!!  out as real*8 == Double_Type, not real*4 == Float_Type       !!!

   variable ne = 600, ifl=1, phot = Double_Type[ne], energ = @phot;
   variable ebins = [-10:7.62:17.617/ne];
   
   ebins = 10.^ebins;

   % Model called on a fixed grid, then interpolated

   xrbjet(ebins,ne,pars,ifl,phot,energ);

   % Convert to keV and ph/cm^2/s/keV

   energ = (10.^energ)*(6.6260755e-18/1.60217733);
   phot = (10.^phot)*(10./6.6260755/energ); 

   % Interpolate back to original grid size
   xrbinterp(ear, energ, phot, photar, ne, nein);

   % Multiply by normalization, and return in wavelength ascending order

   return reverse(par[0]*photar);
}

static variable jet_input = { 
"norm",              1.,     0.,   1.e6,
"mbh [msun]",       10.,     1.,   1.e9,
"jetrat [L_edd]", 1.e-3,  1.e-7,     1.,
"pspec",             2.,    1.7,     4.,
"zsh [r_g]",        10.,     5.,   1.e4,
"r0 [r_g]",          5.,     2.,   100.,
"hratio",           1.5,    0.1,   100.,
"incl [deg]",       40.,     0.,    90.,
"eltemp [K]",     2.e10,   2.e9,  5.e11,
"plfrac",          0.75,  1.e-4,     1.,
"dkpc [kpc]",        1.,    0.1,   1.e4,
"comsw",             0.,    0.,      1.,
"plotsw",            1.,     0.,     1.,
"fsc",           3.6e-3,  1.e-7, 3.6e-2,
"zmax [lg cm]",     14.,    13.,    16.,
"equip",             2.,    0.1,   100.,
};

variable npar=length(jet_input)/4;
variable jet_pars=String_Type[npar];
variable jet_def=Float_Type[npar];
variable jet_min=Float_Type[npar];
variable jet_max=Float_Type[npar];

variable i;
for(i=0;i<=npar-1;i++)
{
    jet_pars[i]=jet_input[4*i];
    jet_def[i]=jet_input[4*i+1];
    jet_min[i]=jet_input[4*i+2];
    jet_max[i]=jet_input[4*i+3];
}

static variable jet_frz = Integer_Type[npar];

% Default above is no frozen parameters, now set the frozen ones

jet_frz[[0,1,6,7,10,11,12,14]]=1;

add_slang_function("simplejet",jet_pars);

define simplejet_defaults(i)
{
   return (jet_def[i],jet_frz[i],jet_min[i],jet_max[i]);
}

set_param_default_hook("simplejet","simplejet_defaults");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Just making some wavelength grids and parameters to test the output
% of the above

variable par = @jet_def;
variable pars = par[[1:15]];
variable lo = _A([0.01:10:0.01]);
variable hi = make_hi_grid(lo);
