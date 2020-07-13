msconvert Z:\1_QEdata\LTC\DATA\HILIC\190718_DepthProfiles_FK180310\*.raw --mzML --filter "peakPicking true 1-" --filter "polarity positive" -z -o mzMLs/pos/
msconvert Z:\1_QEdata\LTC\DATA\HILIC\190718_DepthProfiles_FK180310\MSMS\*DDApos*.raw --mzML --filter "peakPicking true 1-" --filter "polarity positive" -z -o mzMLs/pos/MSMS/
msconvert Z:\1_QEdata\LTC\DATA\HILIC\190718_DepthProfiles_FK180310\*.raw --mzML --filter "peakPicking true 1-" --filter "polarity negative" -z -o mzMLs/neg/
msconvert Z:\1_QEdata\LTC\DATA\HILIC\190718_DepthProfiles_FK180310\MSMS\*DDAneg*.raw --mzML --filter "peakPicking true 1-" --filter "polarity negative" -z -o mzMLs/neg/MSMS/

msconvert Z:\1_QEdata\LTC\DATA\HILIC\HILIC_2018\180205_Wei_MESO-SCOPE_HRM_DepthProfile\*.raw --mzML --filter "peakPicking true 1-" --filter "polarity positive" -z -o mzMLs/pos/
msconvert Z:\1_QEdata\LTC\DATA\HILIC\HILIC_2018\180205_Wei_MESO-SCOPE_HRM_DepthProfile\*Pos_dda*.raw --mzML --filter "peakPicking true 1-" --filter "polarity positive" -z -o mzMLs/pos/MSMS/
msconvert Z:\1_QEdata\LTC\DATA\HILIC\HILIC_2018\180205_Wei_MESO-SCOPE_HRM_DepthProfile\*.raw --mzML --filter "peakPicking true 1-" --filter "polarity negative" -z -o mzMLs/neg/
msconvert Z:\1_QEdata\LTC\DATA\HILIC\HILIC_2018\180205_Wei_MESO-SCOPE_HRM_DepthProfile\*Neg_dda*.raw --mzML --filter "peakPicking true 1-" --filter "polarity negative" -z -o mzMLs/neg/MSMS/

msconvert Z:\1_QEdata\LTC\DATA\HILIC\HILIC_2017\170706_Gradients2.0_CruiseFilters_2\170706_Blk_Blk0p2_?.raw --mzML --filter "peakPicking true 1-" --filter "polarity negative" -z -o mzMLs/neg/
msconvert Z:\1_QEdata\LTC\DATA\HILIC\HILIC_2017\170706_Gradients2.0_CruiseFilters_2\170706_Blk_Blk0p2_?.raw --mzML --filter "peakPicking true 1-" --filter "polarity positive" -z -o mzMLs/pos/

msconvert Z:\1_QEdata\Will\falkor_MSMS_inclusion\*DDApos*.raw --mzML --filter "peakPicking true 1-" --filter "polarity positive" -z -o mzMLs/pos/MSMS/