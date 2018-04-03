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
