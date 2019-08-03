variable datapath = "/ata-school/mwl_tutorial/";

%%PCA
variable pid = load_data(datapath+"pca.pha");
set_sys_err_frac(pid,0.005);
group(pid;min_sn=4.5,bounds=3.,unit="kev");
xnotice_en(pid,3,22);

%%HEXTE
variable hid = load_data(datapath+"hxt.pha");
group(hid;min_sn=4.5,bounds=20.,unit="kev");
xnotice_en(hid,20,200);

%%Radio/IR
variable rid=load_radio2(datapath+"53082_radio-OIR.dat");


