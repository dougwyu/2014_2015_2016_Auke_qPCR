2014_2015_2016_Auke_qPCR_supersimple_revision2.R contains the code used for the post-review version of the manuscript. The only real change is the use of log(Qcorr_qPCR) instead of Qcorr_qPCR as the predictor. This improved fit to zero count values.

2014_2015_2016_Auke_qPCR_supersimple.R contains the final code used for the analyses and figures included in "Environmental DNA for the enumeration and management of Pacific salmon" by Levi et al.

2014_2015_2016_Auke_qPCR.Rmd contains much additional code used during exploration of the two datasets, which are also included in this repo.

The manuscript is now published:

Levi, T., Allen, J.M., Bell, D., Joyce, J., Russell, J.R., Tallmon, D.A., Vulstek, S.C., Yang, C.Y., Yu, D.W. 2019. Environmental DNA for the enumeration and management of Pacific salmon. *Molecular Ecology Resources*. 19:597-608. [https://doi.org/10.1111/1755-0998.12987]


Datafile column definitions (Sockeyelong_all.tsv and Coholong_all.tsv)

Date:  Date of sample

Gage_Height:  Height of river in inches

Depth_in:  Depth of river in inches

Q_cfs:  Water flux in cubic feet per second

Temp_C:  water temperature, Celsius

CT_mean:  mean of the CT (qPCR 'threshold cycle')  A relative measure of the concentration of the target molecule in the sample (in our case, the PCR product)

CT_sd:  standard deviation of the CT

QUANT_mean:  mean of the estimated concentration of the target molecule (ng/ul), estimated using the CT_mean and a standardisation curve

QUANT_sd:  standard deviation of the estimated concentration of the target molecule

n_qpcrs:  number of qPCRs carried out per sample

Qcorr_qPCR:  QUANT_mean * Q_cfs (flow-corrected eDNA rate)

date.n:  index number for date

Qcorr_qpcr.lag:  Qcorr_qPCR shifted forward by one day

Qcorr_qpcr.lead:  Qcorr_qPCR shifted backward by one day

Qcorr_qpcr.lead2:  Qcorr_qPCR shifted forward by two days

Qcorr_qpcr.lead3:  Qcorr_qPCR shifted backward by two days

Sockeyetype/Cohotype:  life-history stage

Count:  Number of individuals (the response variable)

Year:  year of sample
