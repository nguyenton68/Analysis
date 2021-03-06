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
    - Get n2 density inside polarized he3 cell.
    - Get n2 dilution in asymmetry analysis.
- Details:
   - Simulation: n2 cross section is used in this simulation. You need to be careful with inputs:
     - Radiation length: is it He3 or reference cell and at what density?
     - Before apply W cut, make sure that all W peaks are aligned
   - Extract data code:
     - Get absolute n2 cross section: /n2_absolute_xs/get_xs.C
     - Get n2 yield from reference cell: /pressure_curve/get_yield.C
     - Plot all n2 reference cells in order to determine how many MeV need to be shifted before apply cut on W: /pressure_curve/plot_all_refcell.C
     - The output dat file for pressure curve, the yield is already corrected for different in radiation length: /pressure_curve/yield_atm.txt
     - Analyze the raw root file and save it into a new root file to analyze further: /pressure_curve/temp.C
     - Normalization file for all runs involved in this study: 2080 (glass), 2081-2085 (n2 ref cell), 2087 (he3).
     Run, p0, ps1, livetime, charge, preshower lower limit, shower lower limit, cherenkov lower limit, scintillator eff, pid eff, vdc eff\
     
3. He3 elastic cross section:
- get_xs.C: run by type analyzer get_xs.C
  - Need input root file for he3, n2, glass runs. These root files have been applied analysis cut. Analysis cut is done by temp.C.
  - Output is cross section in nb.
  - Need to check: density, target length, solid angle,...
- temp.C: apply analysis cut: 2D cut on theta, phi reconstructed angles. PID cut, W cut. Output is root files used for get cross section later
- MC: typical input file for he3 simulation. And main code for simulate he3 elastic
