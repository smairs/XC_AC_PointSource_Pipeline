import numpy as np
import matplotlib.pyplot as plt

##################################################
# Set up plot style for consistency for paper!
from matplotlib.pyplot import rc
rc('text', usetex=True)
rc('font',family = 'serif')
rc('xtick', labelsize=14)
rc('xtick', top=True)
rc('xtick', direction='in')
rc('ytick', labelsize=14)
rc('ytick', right=True)
rc('ytick', direction='in')
rc('legend', markerscale=1)
rc('legend', fontsize='medium')
rc('axes',labelsize=16)
##################################################


def weighted_cal(cut_off,variables,filename,outdir,date_cutoff=20210410,numiter=5,fid_450_perc=0.05,goodmaps=False):
    '''

    :param cut_off:  how deep into the file of sources do we want to look? 0.1Jy @ 850, 1 Jy @450 -- no deeper!
    :param variables: array of soureceinfo indices -- sources you don't want to include (can be empty)!
    :param filename: sourceinfo file
    :param numiter: The number of iterations - we iterate 5 times by default just because it converges there
    :param fid_450_perc: The fiducial threshold with which to multiply the weighted mean of the source flux (0.05 = 5%)
    :return: lots
    '''

    ###########
    # Figure out what wavelength and region we are dealing with
    ###########
    namecomponents = np.array(filename.split('/')[-1].split('_'))
    for eachcomp in namecomponents:
        if eachcomp[1:]=='50':
            wave = eachcomp

    region = filename.split('/')[-1].split('_')[0]


    ############
    # Read in sources and remove those that are known to vary
    # Get number of epochs and sources, too
    ############

    data    = np.genfromtxt(filename,dtype=None,names=True)
    SR_all  = len(data['Index'])
    Nvars   = len(variables)
    SR      = SR_all-Nvars  	           # Number of sources (minus variables)
    EPtot      = int((len(data.dtype.names)-15)/2) # Number of epochs in total
    caldates     = []
    alldatescans = []
    for i in list(data.dtype.names)[15:15+EPtot]:
        thisdate = i.split('_')[1]
        alldatescans.append(i.split('_')[1]+'_'+i.split('_')[2])
        if int(thisdate)<=date_cutoff:
            caldates.append(thisdate)
    calEP = len(caldates) # Number of dates to consider during calibration (up to date_cutoff) 
    allEP = len(alldatescans) # Number of epochs in total
    srcn    = [i for i in data['Index'] if i not in variables] # Source Indices that are non-variable

    #################
    # Fix up the units and construct array of fluxes in each epoch indexed by source
    #################

    #Units for 850 should be in Jys while 450 and High mass will be in mJy( for now)
    if np.logical_or(wave == '450',region in ['DR21C','DR21N','DR21S','S255','M17','M17SWex']):
        unit_conv = 1.0 #1000.0
    else:
        unit_conv = 1.0

    calflux = np.array([list(data[i])[15:15 + calEP] for i in srcn])/unit_conv # Only include dates <= Date Cutoff
    allflux = np.array([list(data[i])[15:15 + allEP] for i in srcn])/unit_conv # Include All Dates

    #####################
    # Define arrays for the 6 Variables to keep track of:
    #####################

    # These first three should have a value for every source and only rely on the dates we wish to use for the calibration
    means = np.zeros(SR)  #(weighted) mean flux of a given source
    sdws  = np.zeros(SR) #(weighted) flux stand dev of given source
    sdfs  = np.zeros(SR) #fiducial flux stand dev of given source -- this differs at 450 and 850

    # These next three should have a value for every Epoch and should include ALL eopchs
    mult  = np.zeros(allEP)+1    # calibration multiplier -- DOUG CALIBRATES WITH MULTIPLIER -- BUT I CHANGE THIS TO A DIVISOR LATER IN THE OUTPUT FILES FOR CONSISTENCY
    sde   = np.zeros(allEP)+0.02 # formal cal unc (comes out mult calc) - not correct to use in further calculations since it assumes all epochs are the same!
    sdwe  = np.zeros(allEP)+0.02 # (weighted) calibration uncertainty -- for the epoch


    ##################
    ##################
    ##################
    # Begin Iterations!
    ##################
    ##################
    ##################

    for i in range(numiter):

        #####################
        #####################
        # Measurement estimates for all ****Sources**** -- HERE WE ARE FINDING WHAT SOURCES TO USE FOR EPOCH CALIBRATION, SO ONLY USE CAL DATES
        #####################
        #####################

        for s in range(SR):

            ############
            # Weighted mean calculation
            ############

            dw = sum(1/(sdwe[0:calEP]**2)) # From Chris Beaumont's wmean (2009)
            numerator = sum(mult[0:calEP]*np.array(calflux[s])/(sdwe[0:calEP]**2)) # From Chris Beaumont's wmean (2009) - no square on top
            #print(i,s,numerator,dw,sdwe[0:calEP])
            weighted_mean = numerator/dw # From Chris Beaumont's wmean (2009)
            means[s] = weighted_mean # w. mean of the best multiplier we have right now -- sdwe, (weight of epoch)

            ############
            # Delta between mean value and the calibrated flux
            ############

            delt     = (means[s] - mult[0:calEP]*np.array(calflux[s]))

            ############
            # Now do weighting for what we think uncertainty in source is!
            ############

            numerator2 = sum(delt**2/sdwe[0:calEP]**2) # The weighted mean here is over delt**2 which is why square on top
            weighted_mean2 = numerator2/dw
            sdws[s] = np.sqrt(weighted_mean2*(calEP/(calEP-1.))) # weighted measure of what we think is the src uncertainty


            #############
            # Fiducial Threshold to prevent Runaways
            #############

            # Determination fiducial to compare against. The only reason it comes in here is that
            # We compare the standard deviation to the fiducial and say it can't get bettter than dfiducial
            # (to stop runaways where 1 source is doing ALL the work)
            # For 450 - we just set it to some percentage (default 5) of the mean flux as the cutoff.
            if wave == '850':
                sdfs[s]  = np.sqrt(14*14 + 0.02*0.02*means[s]**2) # should be in millijanskys already, so 14 millijanskyes instead of 0.014
            else:
                sdfs[s] = fid_450_perc*means[s]
            #print("Src: ", i, s, means[s], sdws[s], sdfs[s], sdws[s]/sdfs[s])
            if sdws[s]/sdfs[s] < 1:
                sdws[s] = sdfs[s]

        #######################
        # Determination of which sources to use for Epoch Calibration
        #######################

        #print('means',means)
        #print('cut_off',cut_off)
        slist = np.where(means > cut_off) # which sources gt than cutoff
        SRC = len(means[slist]) # number of sources gt than cutoff
        tw   = sum((means[slist]*means[slist])/(sdws[slist]*sdws[slist]))	#total weighting - for fractional weights


        #######################
        #######################
        # Get measurement estimates for all **Epochs**
        #######################
        #######################
        for e in range(allEP):

            #############
            # least squares fit with b=0 -- sxy flux/mean -- sxx is flux^2/mean
            #############

            #print('e',e)
            #print('slist',slist)
            #print('sdws[slist]',sdws[slist])
            sxy = sum((allflux[slist].T[e]*means[slist])/(sdws[slist]*sdws[slist]))
            sxx = sum((allflux[slist].T[e]*allflux[slist].T[e])/(sdws[slist]*sdws[slist]))

            #############
            # Extract Values
            #############

            mult[e] = sxy/sxx #THE MULTIPLIER!
            sde[e]  = np.sqrt(1./sxx)
            delt    = (means[slist]/allflux[slist].T[e] - mult[e]) # means/flux should be same as the multiplier, ideally!

            #############
            # Epoch Uncertainty
            #############

            dw = sum(1.0/(sdws[slist]/means[slist])**2)
            weighted_mean = sum(delt**2/((sdws[slist]/means[slist])**2))/dw
            sdwe[e] = np.sqrt(weighted_mean/(1.*SRC)) # over all sources - divide by sqrt of sources - Epoch unc

            # Don't let epoch uncertainty get 70% better than formal error -- 450 runaway prevention.
            # Any epoch can have only about twice the weighting it would naturally have,
            # Since this 0.7 comes in as a square -- it can get much less, but not much better.
            # This helps to have a measure for when 450 is lousy
            if sdwe[e]/sde[e] < 0.70:
                sdwe[e] = 0.70*sde[e]
            #print("Epo: ", i, e, mult[e], sde[e], sdwe[e], SRC)


        ################
        # Renormalisation
        ################

        # Renormalize to keep the over-all calibration from drifting
        # treated all epochs independently -- but you need to know where you are fixing the flux to
        # 2021-05-14 -- The first epoch was the original choice here, but I have changed this to be the average between all calibration dates for consistency with the PScal
        div    = np.nanmedian(mult[0:calEP+1]) # calEP is the number of calibration dates. The +1 is for python indexing, to be inclusive of all cal dates. Used to be: mult[0]
        mult   = mult/div
        sdwe   = sdwe/div
        sde    = sde/div

        #################
        #Printing Iteration Information To Screen
        #################

        weightpersource = np.sqrt(sum(sdws[slist]*sdws[slist]/(means[slist]*means[slist]))/(SRC)) # Will not change as we add epochs, based only on date_cutoff data
        # Iteration number, Number of Epochs, Number of Sources (sans vars),
        # Number of Calibrators, Mean mult, Error mult, Total Weight, Average Weight per source
        #print('Iteration, NumEpochs, NumSources_novar, NumCals, avgMult, stdMults, TotalWeight, AvgWeightperSource:')
        #print(i, allEP, SR, SRC, np.mean(mult), np.std(mult,ddof=1), tw, weightpersource)

    ####################
    ####################
    ####################
    # Creating Info Files
    ####################
    ####################
    ####################

    ########################
    # **All** Source Information - Using ***Uniform*** Weighting per Epoch
    ########################

    if goodmaps:
        fsource = open(outdir+filename.split('/')[-1].split('.txt')[0]+'_allsource_uniformweighting_GoodMaps.txt','w')
    else:
        fsource = open(outdir+filename.split('/')[-1].split('.txt')[0]+'_allsource_uniformweighting.txt','w')
    fsource.write('#SourceInd\tWMeanFlux\tWFluxSD\tWFidSD\tWFluxSD_div_WFidSD\n')
    for s in range(SR):
        means[s] = np.mean(mult*np.array(allflux[s]))
        sdws[s]  = np.std(mult*np.array(allflux[s]),ddof=1) # NOT REALLY WEIGHTED - JUST FOR COMPARISON WITH MY VALUES
        if wave == '850':
            sdfs[s] = np.sqrt(14 * 14 + 0.02 * 0.02 * means[s] ** 2) # should be in millijanskys already, so 14 millijanskyes instead of 0.014
        else:
            sdfs[s] = fid_450_perc * means[s]
        fsource.write(str(srcn[s])+'\t'+str(means[s])+'\t'+str(sdws[s])+'\t'+str(sdfs[s])+'\t'+
                      str(sdws[s]/sdfs[s])+'\n')

    fsource.close()


    ########################
    # **THE** Important file -- Mean fluxes, weightings, which sources to use -- essential for future epochs!
    ########################


    # Calibrator Source Information - Using Weighted Epochs
    if goodmaps:
        fcal = open(outdir+filename.split('/')[-1].split('.txt')[0]+'_calsource_weightedmean_GoodMaps.txt','w')
    else:
        fcal = open(outdir+filename.split('/')[-1].split('.txt')[0]+'_calsource_weightedmean.txt','w')
    fcal.write('#SourceInd\tWMeanFlux\tWFluxSD\tWFidSD\tWFluxSD_div_WFidSD\twt\n')

    for s in range(SR):

        ############
        # Weighted mean calculation
        ############
        dw = sum(1. / (sdwe[0:calEP] ** 2))  # From Chris Beaumont's wmean (2009)
        numerator = sum(mult[0:calEP] * np.array(calflux[s]) / (sdwe[0:calEP] ** 2)) # From Chris Beaumont's wmean (2009) - no square on top
        weighted_mean = numerator/dw  # From Chris Beaumont's wmean (2009)
        means[s] = weighted_mean  # w. mean of the best multiplier we have right now -- sdwe, (weight of epoch)

        ############
        # Delta between mean value and the calibrated flux
        ############

        delt = (means[s] - mult[0:calEP] * np.array(calflux[s]))

        ############
        # Now do weighting for what we think uncertainty in source is!
        ############

        numerator2 = sum(delt ** 2 / sdwe[0:calEP] ** 2)  # The weighted mean here is over delt**2 which is why square on top
        weighted_mean2 = numerator2 / dw
        sdws[s] = np.sqrt(weighted_mean2 * (calEP / (calEP - 1.)))  # weighted measure of what we think is the src uncertainty

        #############
        # Fiducial Threshold to prevent Runaways
        #############

        # Determination fiducial to compare against. The only reason it comes in here is that
        # We compare the standard deviation to the fiducial and say it can't get bettter than dfiducial
        # (to stop runaways where 1 source is doing ALL the work)
        # For 450 - we just set it to some percentage (default 5) of the mean flux as the cutoff.
        if wave == '850':
            sdfs[s] = np.sqrt(14 * 14 + 0.02 * 0.02 * means[s] ** 2) # should be in millijanskys already, so 14 millijanskyes instead of 0.014
        else:
            sdfs[s] = fid_450_perc * means[s]
        #print("Src: ", i, s, means[s], sdws[s], sdfs[s], sdws[s] / sdfs[s])
        if sdws[s] / sdfs[s] < 1:
            sdws[s] = sdfs[s]

    tw = sum((means[slist] * means[slist])/(sdws[slist] * sdws[slist]))  # total weighting - for fractional weights
    for s in slist[0]:
        wt = ((means[s]*means[s])/(sdws[s]*sdws[s]))/tw
        fcal.write(str(srcn[s])+'\t'+str(means[s])+'\t'+str(sdws[s])+'\t'+str(sdfs[s])+'\t'+
                   str(sdws[s]/sdfs[s])+'\t'+str(wt)+'\n')
    fcal.close()

    ##############
    #Printing out a file with the epoch information!
    ##############
    if goodmaps:
        fepoch = open(outdir+filename.split('/')[-1].split('.txt')[0]+'_calepoch_weightedmean_GoodMaps.txt','w')
    else:
        fepoch = open(outdir+filename.split('/')[-1].split('.txt')[0]+'_calepoch_weightedmean.txt','w')
    fepoch.write('#DateScan\tDivisor\tFormalUnc\tWeightedCalUnc\tNumCal\n')
    for e in range(allEP):
        # CHANGING TO DIVISOR FOR OUR RECORDS - NEED TO CHANGE UNCERTAINTY VIA RELATIVE ERROR CALCULATIONS
        # MULT +/- MULTErr ---> DIV +/- (MULTErr/MULT)*DIV
        fepoch.write(alldatescans[e]+'\t'+str(1.0/mult[e])+'\t'+str((sde[e]/mult[e])*(1.0/mult[e]))+
                     '\t'+str((sdwe[e]/mult[e])*(1.0/mult[e]))+'\t'+str(SRC)+'\n')
    fepoch.close()

    ##################
    ##################
    ##################
    # Printing Information to Screen!
    ##################
    ##################
    ##################

    #print(sorted(list(sdws[slist]/sdfs[slist]))[::-1])
    #order = sorted(list(sdws[slist]/sdfs[slist]))[::-1]
    #for j in range(SRC):
    #    s = slist[order[j]]
    #    wt = ((means[s]*means[s])/(sdws[s]*sdws[s]))/tw
    #    print('\n\nORDER OF SourceUnc/FidUnc:\n')
    #    print('SourceInd, WMeanFlux, WSourceUnc, Flux/Unc, Unc/Fiducial, FractionalWeight')
    #    print(srcn[s], means[s], sdws[s], means[s]/sdws[s], sdws[s]/sdfs[s], wt)


    #order = reverse(sort(((means[slist]*means[slist])/(sdws[slist]*sdws[slist]))))
    #for j in range(SRC):
    #    s = slist[order[j]]
    #    wt = ((means[s]*means[s])/(sdws[s]*sdws[s]))/tw
    #    print('\n\nORDER OF MeanFlux^2/SourceUnc^2:\n')
    #    print('SourceInd, WMeanFlux, WSourceUnc, Flux/Unc, Unc/Fiducial, FractionalWeight')
    #    print(srcn[s], means[s], sdws[s], means[s]/sdws[s], sdws[s]/sdfs[s], wt)


    #for e in range(EP):
    #    # CHANGING TO DIVISOR FOR OUR RECORDS - NEED TO CHANGE UNCERTAINTY VIA RELATIVE ERROR CALCULATIONS
    #    # MULT +/- MULTErr ---> DIV +/- (MULTErr/MULT)*DIV
    #    print('\n\nEpoch Ind, Divisor, FormalUnc, DivisorUnc, DivisorUnc/FormalUnc, NumCal')
    #    print(e, 1.0/mult[e], (sde[e]/mult[e])*(1.0/mult[e]), (sdwe[e]/mult[e])*(1.0/mult[e]), sdwe[e]/sde[e], SRC)

    #####################
    #####################
    #####################
    # Plots!
    #####################
    #####################
    #####################

    ##########
    # Set up 2x2 plot structure
    ##########

    fig = plt.figure()
    ax11 = fig.add_subplot(221)
    ax12 = fig.add_subplot(222)
    ax21 = fig.add_subplot(223)
    ax22 = fig.add_subplot(224)

    ############################
    # Plot 0: Source Flux versus Fiducial Unc
    ############################
    ax11.scatter(means,means/sdws,marker='.',color='b',s=4)
    ax11.set_xscale('log')
    ax11.set_yscale('log')
    flx = (np.arange(0,10000.5,1)+1)*0.001
    if wave == '850':
        ax11.plot(flx,flx/np.sqrt(14*14 + 0.02*0.02*flx*flx),linewidth=2,linestyle='dashed',color='k') # should be in millijanskys already, so 14 millijanskyes instead of 0.014
    else:
        ax11.plot(flx,flx/(fid_450_perc * flx))
    ax11.set_ylabel('Source Flux/Source Unc')
    ax11.set_xlabel('Source Flux')

    ############################
    # Plot 1 - Source flux vs Calibration Weight
    ############################
    ax12.scatter(means[slist],((means[slist]*means[slist])/(sdws[slist]*sdws[slist]))/tw,marker='.',s=4,color='b')
    ax12.set_xscale('log')
    ax12.set_yscale('log')
    ax12.set_ylabel('Source Cal Weight')
    ax12.set_xlabel('Source Flux')

    ############################
    # Plot 2 - Epoch Calibration Divisor
    ############################
    #print(np.arange(0,allEP+0.1,1))
    #print(1.0/mult)
    ax21.scatter(np.arange(0,allEP,1),(1.0/mult), color='r',marker='o')
    ax21.set_xlim(xmin=-1, xmax=allEP)
    ax21.set_ylim(ymin = np.mean(1.0 / mult) - 5 * np.std((1.0 / mult), ddof=1),
                  ymax = np.mean(1.0 / mult) + 5 * np.std((1.0 / mult), ddof=1))

    ax21.axhline(y=np.mean(1.0/mult),linestyle='solid',color='k')
    ax21.axhline(y=np.mean(1.0/mult)+np.std(1.0/mult,ddof=1), linestyle='dashed',color='k')
    ax21.axhline(y=np.mean(1.0/mult)-np.std(1.0/mult,ddof=1), linestyle='dashed',color='k')
    ax21.set_ylabel('Cal Divisor')
    xlabels = []
    for i in alldatescans:
        xlabels.append(i.split('_')[0])
    ax21.set_xticklabels(xlabels,rotation=20)

    ############################
    # Plot 3 - Epoch Calibration Uncertainty
    ############################
    # CHANGING TO DIVISOR FOR OUR RECORDS - NEED TO CHANGE UNCERTAINTY VIA RELATIVE ERROR CALCULATIONS
    # MULT +/- MULTErr ---> DIV +/- (MULTErr/MULT)*DIV
    sdwe_converted_to_div = (sdwe/mult)*(1.0/mult)
    ax22.scatter(np.arange(0,allEP,1),sdwe_converted_to_div,color='r')
    ax22.set_xlim(xmin=-1,xmax=allEP)
    ax22.set_ylim(ymin=np.mean(sdwe_converted_to_div)-5*(np.std(sdwe_converted_to_div,ddof=1)+0.01),
                  ymax = np.mean(sdwe_converted_to_div)+5*(np.std(sdwe_converted_to_div,ddof=1)+0.01))
    ax22.set_ylabel('Epoch Cal Unc')
    ax22.set_xticklabels(xlabels,rotation=20)

    if goodmaps:
        plt.savefig(outdir+filename.split('/')[-1].split('.txt')[0]+'_weightedcal_plots_GoodMaps.pdf',format='pdf')
    else:
        plt.savefig(outdir+filename.split('/')[-1].split('.txt')[0]+'_weightedcal_plots.pdf',format='pdf')
    plt.clf()

###############################################################
###############################################################
###############################################################
# RUN PROGRAM
###############################################################
###############################################################

##OMC 2/3 450 micron test:
#weighted_cal(0.1,[9],'OPHCORE_20200901_00015_450_sourceinfo.txt',numiter=5,fid_450_perc=0.05)
#
##OMC 2/3 850 micron test:
##weighted_cal(1.0,[1,40,57,97],'OMC23_20200222_00020_850_sourceinfo_CR3.txt',numiter=5,fid_450_perc=0.05)


# Variables (2021-04-01):

 #############
 #############
 #450 Microns:
 #############
 #############

 #       if region == 'IC348':
 #           known_vars = [1]
 #           var_names  = ['IC 348 1']
 #       if region == 'NGC1333':
 #           known_vars = [0,10]
 #           var_names  = ['IRAS 4A','BOLO 40']
 #       if region == 'NGC2024':
 #           known_vars = []
 #           var_names  = []
 #       if region == 'NGC2071':
 #           known_vars = [1,2]
 #           var_names  = ['HOPS 358','HOPS 373']
 #       if region == 'OMC23':
 #           known_vars = [1,12,29,100]
 #           var_names  = ['HOPS 88','HOPS 370','HOPS 383','V1017 ORI']
 #       if region == 'OPHCORE':
 #           known_vars = [9]
 #           var_names  = ['OPH 162624']
 #       if region == 'SERPM':
 #           known_vars = [0,10,11]
 #           var_names  = ['SMM 1','EC 53','SMM 10']
 #       if region == 'SERPS':
 #           known_vars = [12]
 #           var_names  = ['IRAS 18270']
 #       if region == 'DR21C':
 #           known_vars = []
 #           var_names = []
 #       if region == 'DR21N':
 #           known_vars = []
 #           var_names = []
 #       if region == 'DR21S':
 #           known_vars = []
 #           var_names = []
 #       if region == 'M17':
 #           known_vars = []
 #           var_names = []
 #       if region == 'M17SWex':
 #           known_vars = []
 #           var_names = []
 #       if region == 'S255':
 #           known_vars = []
 #           var_names = []

 ###########
 ###########
 # 850 microns!
 ###########
 ###########

 #       if region == 'IC348':
 #           known_vars = [1]
 #           var_names  = ['IC 348 1']
 #       if region == 'NGC1333':
 #           known_vars = [0,19,24]
 #           var_names  = ['IRAS 4A','BOLO 40']
 #       if region == 'NGC2024':
 #           known_vars = []
 #           var_names  = []
 #       if region == 'NGC2071':
 #           known_vars = [1,12]
 #           var_names  = ['HOPS 358','HOPS 373']
 #       if region == 'OMC23':
 #           known_vars = [40,1,57,97]
 #           var_names  = ['HOPS 88','HOPS 370','HOPS 383','V1017 ORI']
 #       if region == 'OPHCORE':
 #           known_vars = [28]
 #           var_names  = ['OPH 162624']
 #       if region == 'SERPM':
 #           known_vars = [0,2,3]
 #           var_names  = ['SMM 1','EC 53','SMM 10']
 #       if region == 'SERPS':
 #           known_vars = [46]
 #           var_names  = ['IRAS 18270']
 #       if region == 'DR21C':
 #           known_vars = []
 #           var_names = []
 #       if region == 'DR21N':
 #           known_vars = []
 #           var_names = []
 #       if region == 'DR21S':
 #           known_vars = []
 #           var_names = []
 #       if region == 'M17':
 #           known_vars = []
 #           var_names = []
 #       if region == 'M17SWex':
 #           known_vars = []
 #           var_names = []
 #       if region == 'S255':
 #           known_vars = []
 #           var_names = []
