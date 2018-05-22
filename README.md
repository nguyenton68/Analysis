# Analysis
small angle GDH analysis

Detector analysis:

1. VDC efficiency:
vdc_eff = one_track/total_good_track
Modifications:
- Change fabs --> abs (04/03/2018): the code already has #include <TMath.h>, so abs(double)=double
- Subroutine might missing argument: need to make sure it match

Command to run:
 analyzer -q -d best_vdc_multi_track_eff_list_v4_final.C+


2. N2 dilution and pressure curves:
- The study is for 2.2 GeV elastic kinematic.
- Purpose: 
  a, Get n2 density inside polarized he3 cell.
  b, Get n2 dilution in asymmetry analysis.
- Details:
  + Simulation: n2 cross section is used in this simulation. You need to be careful with inputs:
   a, Radiation length: is it He3 or reference cell and at what density?
   b, Before apply W cut, make sure that all W peaks are aligned
  + Extract data code:
   a, Get absolute n2 cross section: /n2_absolute_xs/get_xs.C
   b, Get n2 yield from reference cell: /pressure_curve/get_yield.C
   c, Plot all n2 reference cells in order to determine how many MeV need to be shifted before apply cut on W: /pressure_curve/plot_all_refcell.C
   d, The output dat file for pressure curve, the yield is already corrected for different in radiation length: /pressure_curve/yield_atm.txt
   e, Analyze the raw root file and save it into a new root file to analyze further: /pressure_curve/temp.C
   f, Normalization file for all runs involved in this study: 2080 (glass), 2081-2085 (n2 ref cell), 2087 (he3).
   Run, p0, ps1, livetime, charge, preshower lower limit, shower lower limit, cherenkov lower limit, scintillator eff, pid eff, vdc eff
