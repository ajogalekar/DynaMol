$RunGauss
%Chk=1OKL_small_opt.chk
%Mem=3000MB
%NProcShared=2
#N B3LYP/6-31G* Freq=NoRaman Geom=AllCheckpoint Guess=Read
Integral=(Grid=UltraFine) SCF=XQC IOp(7/33=1)
 
 
