msconvert Z:\1_QEdata\LTC\DATA\HILIC\190718_DepthProfiles_FK180310\*.raw --mzML --filter "peakPicking true 1-" --filter "polarity positive" -z -o mzMLs_pos/
msconvert Z:\1_QEdata\LTC\DATA\HILIC\HILIC_2017\170706_Gradients2.0_CruiseFilters_2\170706_Blk_Blk0p2_?.raw --mzML --filter "peakPicking true 1-" --filter "polarity positive" -z -o mzMLs_pos/
msconvert Z:\1_QEdata\LTC\DATA\HILIC\190718_DepthProfiles_FK180310\MSMS\*.raw --mzML --filter "peakPicking true 1-" --filter "polarity positive" -z -o mzMLs_pos\MSMS\
msconvert Z:\1_QEdata\LTC\DATA\HILIC\190718_DepthProfiles_FK180310\*.raw --mzML --filter "peakPicking true 1-" --filter "polarity negative" -z -o mzMLs_neg/
msconvert Z:\1_QEdata\LTC\DATA\HILIC\HILIC_2017\170706_Gradients2.0_CruiseFilters_2\170706_Blk_Blk0p2_?.raw --mzML --filter "peakPicking true 1-" --filter "polarity negative" -z -o mzMLs_neg/
msconvert Z:\1_QEdata\LTC\DATA\HILIC\190718_DepthProfiles_FK180310\MSMS\*.raw --mzML --filter "peakPicking true 1-" --filter "polarity negative" -z -o mzMLs_neg\MSMS\