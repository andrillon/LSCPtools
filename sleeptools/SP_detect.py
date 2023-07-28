#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed 28 2023
@author: Arthur
From @author: Thomas

SS_detection_vALC.py
"""
# %% 

import multiprocessing

def SS_detection(file_path):
    
        ### import libraries
    import numpy as np, scipy as sc
    import mne, localdef
    from datetime import date
    todaydate = date.today().strftime("%d%m%y")
    
    # init variables
    sleepSpindles=[]
    spindlesOutEnergy=[]
    spindlesOutDuration=[]
    channels = ["F3", "C3", "O1"]
    
    root_dir = localdef.DDE_path_PSG
    preproc_dir = f"{root_dir}/Preproc"
    fig_dir = f"{root_dir}/Figs"
    
    ### Selection criteria for Sleep Spindles
    minMax_SD_threshold = [2, 10, 1] # SD for threshold, SD for artefact rejection
    minMaxDurations = [0.5, 2.5]
    minFreqValue = 11
    maxFreqValue = 16
    NREMonlyForThr = 0
    filt_range = [11, 16] # in Hz
    
    sub_id = file_path[len(preproc_dir):]
    
    eeg_raw_data = mne.io.read_raw_fif(file_path, preload = True)
    sfreq = int(eeg_raw_data.info['sfreq'])
    scoring_ts = ""
    
    # warning scoring should be in sample
    nchan, nsamples = eeg_raw_data._data.shape
    #if np.size(scoring_ts)!=nsamples:
    #    return sleepSpindles, spindlesOutEnergy, spindlesOutDuration
    
    sfreq=int(sfreq)
    
    ### Prepare EEG for pre-processing
    eeg_ss_bp = eeg_raw_data.copy()
    eeg_broad_bp = eeg_raw_data.copy()
    eeg_high_bp = eeg_raw_data.copy()
    eeg_alpna_bp = eeg_raw_data.copy()

    ### Pre-processing of EEG data: Filtering
    print('Filter spindle range: ')
    eeg_ss_bp.filter(filt_range[0], filt_range[1], l_trans_bandwidth='auto', h_trans_bandwidth='auto',
               filter_length='auto', phase='zero')
    nchan, nsamples=eeg_ss_bp._data.shape
    
    print('Filter broad range: ')
    eeg_broad_bp.filter(0.1, 30, l_trans_bandwidth='auto', h_trans_bandwidth='auto',
               filter_length='auto', phase='zero')
    print('Filter high range: ')
    eeg_high_bp.filter(60, 90, l_trans_bandwidth='auto', h_trans_bandwidth='auto',
               filter_length='auto', phase='zero')
    print('Filter alpha range: ')
    eeg_alpna_bp.filter(8, 10, l_trans_bandwidth='auto', h_trans_bandwidth='auto',
               filter_length='auto', phase='zero')


    ### Loop across channels
    nchan, nsamples=eeg_ss_bp._data.shape
    for ch in range(0,len(channels)):
        print('Processing: ')
        print(channels[ch])
        this_eeg_chan=eeg_ss_bp.copy().pick_channels([channels[ch]])
        this_eeg_chan=this_eeg_chan[0]
        this_eeg_chan=this_eeg_chan[0][0,:]
        
        this_high_eeg_chan=eeg_high_bp.copy().pick_channels([channels[ch]])
        this_high_eeg_chan=this_high_eeg_chan[0]
        this_high_eeg_chan=this_high_eeg_chan[0][0,:]
        
        this_bp_eeg_chan=eeg_broad_bp.copy().pick_channels([channels[ch]])
        this_bp_eeg_chan=this_bp_eeg_chan[0]
        this_bp_eeg_chan=this_bp_eeg_chan[0][0,:]
        
        this_alpha_eeg_chan=eeg_alpna_bp.copy().pick_channels([channels[ch]])
        this_alpha_eeg_chan=this_alpha_eeg_chan[0]
        this_alpha_eeg_chan=this_alpha_eeg_chan[0][0,:]
        
        if np.max(this_eeg_chan)<1: # probably in V and not uV
            this_eeg_chan=this_eeg_chan*1000000
            this_high_eeg_chan=this_high_eeg_chan*1000000
            this_bp_eeg_chan=this_bp_eeg_chan*1000000
            this_alpha_eeg_chan=this_alpha_eeg_chan*1000000
            print('converting to uV')
        
        
        # identify artefacted epochs
        #artefactsEpochs=np.intersect1d(np.where(np.abs(sc.signal.hilbert(this_high_eeg_chan))>30),np.where(np.abs(this_bp_eeg_chan)>500))
        artefactsEpochs=np.where(np.abs(this_bp_eeg_chan)>500)[0]
        
        # Compute envelopes
        envelope = (np.abs(sc.signal.hilbert(this_eeg_chan)));
        envelope_high = (np.abs(sc.signal.hilbert(this_high_eeg_chan)));
        envelope_bp = (np.abs(sc.signal.hilbert(this_bp_eeg_chan)));
        envelope_alpha = (np.abs(sc.signal.hilbert(this_alpha_eeg_chan)));
        envelope_EpochsOfInterest = envelope;
        
        for k in range(0,np.size(artefactsEpochs)):
            envelope_EpochsOfInterest[
                np.max([int(artefactsEpochs[k]-sfreq),0]
                       ):np.min(
                           [int(artefactsEpochs[k]+sfreq),np.size(
                               envelope_EpochsOfInterest)])]=np.nan
    
        # Compute mean and SD (NREM only, not artefacted epochs)
        if NREMonlyForThr==1:
            envelope_EpochsOfInterest[
                np.union1d(np.where(scoring_ts<1),
                           np.where(scoring_ts>3))]=np.nan
        detectionThresholds=[
            np.nanmedian(envelope_EpochsOfInterest),
            np.nanstd(envelope_EpochsOfInterest)
            ]
    
        DetectionThreshold  = detectionThresholds[0] + detectionThresholds[1]*minMax_SD_threshold[0]
        RejectThreshold     = detectionThresholds[0] + detectionThresholds[1]*minMax_SD_threshold[1]
        StartEndThreshold   = detectionThresholds[0] + detectionThresholds[1]*minMax_SD_threshold[2]
    
        # Initialise variables to save
        these_spindles = []
        these_spindlesOutDuration = []
        these_spindlesOutEnergy = []
        nspin=0
        
        # Detection
        # Transform envelope so that every point above the threshold is a 'zero crossing'
        envelopeMinusThreshold = envelope - DetectionThreshold
        pos_index=np.transpose([0]*np.size(envelopeMinusThreshold))
        #index of all positive points for EEG
        pos_index[envelopeMinusThreshold>0]=1
        # Find positive zero-crossing (start of potential REM)
        poscross = np.where(np.diff(pos_index)>0)[0]
        poscross = [x+1 for x in poscross]
        # Find negative zero-crossing (end of potential REM)
        negcross = np.where(np.diff(pos_index)<0)[0]
        negcross = [x+1 for x in negcross]

        peaks = np.where(np.diff(pos_index)==1)[0]
        peaks = [x+1 for x in peaks]
        peaks = [x for x in peaks if envelopeMinusThreshold[x] > 0]
        
        if np.size(negcross)==0 | np.size(poscross)==0:
            sleepSpindles.append(these_spindles)
            print('... no spindle found')
            continue
        
        # Spindles Between 2 segments
        # Correct anomalies in negcross and poscross
        if poscross[-1]>negcross[-1]:
            poscross=poscross[0:np.size(poscross)-1]
        elif negcross[0]<poscross[0]:
            negcross=negcross[1:np.size(negcross)]
        if np.size(negcross)==0 | np.size(poscross)==0:
            sleepSpindles.append(these_spindles)
            print('... no spindle found')
            continue
        if np.size(poscross) != np.size(negcross) :
            if negcross[0]<poscross[0]:
                negcross=negcross[1:np.size(negcross)]
            elif poscross[-1]>negcross[-1]:
                poscross=poscross[0:np.size(poscross)-1]

        # Transform envelope relative to StartEnd threshold to compute spindle durations
        envelopeMinusBottomThreshold = envelope - StartEndThreshold
        pos_index2=np.transpose([0]*np.size(envelopeMinusBottomThreshold))
        pos_index2[envelopeMinusBottomThreshold>0]=1
        poscrossMinimal = np.where(np.diff(pos_index2)>0)[0]
        poscrossMinimal = [x+1 for x in poscrossMinimal]
        negcrossMinimal = np.where(np.diff(pos_index2)<0)[0]
        negcrossMinimal = [x+1 for x in negcrossMinimal]
        
        # Compute spindles candidates properties
        spindleStartTimeArray = [] # we build these arrays as a tool to avoid detecting the same spindle twice
        spindleEndTimeArray   = []
        sndx=-1
        lastSpindleStart=0
        while sndx<np.size(poscross)-1:
            
            sndx=sndx+1
            if np.size(poscrossMinimal)==0 | np.size(negcrossMinimal)==0: # Signal crossed the detection threshold but did not corss twice the StartEnd thr
                continue
          
            # Compute Spindle Start and End time
            if np.size(np.where(poscrossMinimal<poscross[sndx]))==0:
                continue
            tp = [x for x in poscrossMinimal if x < poscross[sndx]]
            spindleStartTime = tp[-1]
            tp = [x for x in negcrossMinimal if x > negcross[sndx]]
            spindleEndTime = tp[0]

            if spindleStartTime-sfreq<0 | spindleEndTime+sfreq>np.size(this_eeg_chan):
                continue
            if np.size(spindleStartTime)==0 | np.size(spindleEndTime)==0: # failed to get start time or end time for some reason
                continue
            spindleDuration  = spindleEndTime - spindleStartTime + 1;
            
            # Check duration is over 0.5s
            if spindleDuration/sfreq<0.5:
                continue
            
            # Check that this spindle candidate is not too close from the next one
            # (<500ms). If it is the case, spindles candidates are concatenated
            if sndx<np.size(poscross)-1 & sndx<np.size(negcross)-1:
                while sndx+1<np.size(poscross)+1 & (negcross[sndx+1]-poscross[sndx])<0.5*sfreq:
                    tp = [x for x in poscrossMinimal if x < poscross[sndx]]
                    spindleStartTime = tp[-1]
                    tp = [x for x in negcrossMinimal if x > negcross[sndx+1]]
                    spindleEndTime = tp[0]
                    del poscross[sndx+1]
                    del negcross[sndx]
            spindleDuration  = spindleEndTime - spindleStartTime + 1


            if (spindleStartTime-lastSpindleStart)<0.5*sfreq: # do not forget to update lastSpindleStart
                continue
            
            #Make sure this is not just another peak for the same spindle (as judged by its start and end times) if so continue and ignore..
            if np.size(np.intersect1d(spindleStartTime, spindleStartTimeArray)) > 0:
                continue
            elif np.size(np.intersect1d(spindleEndTime, spindleEndTimeArray)) > 0:
                continue
            spindleStartTimeArray = np.append(spindleStartTimeArray,spindleStartTime)
            spindleEndTimeArray = np.append(spindleEndTimeArray,spindleEndTime)
            
            # Compute Peak in amplitude time and power
            relevantPeakIdx = np.intersect1d(np.where(peaks > spindleStartTime),np.where(peaks < spindleEndTime))
            relevantPeakIdx = relevantPeakIdx[0]
            peakTimes = peaks[relevantPeakIdx]
            peakTime = peakTimes[envelope[peakTimes] == np.max(envelope[peakTimes])]
            if np.size(peakTime)==0:
                continue
            peakTime=peakTime[0]
            peakEnergy = envelope[peakTime]
            peakEnergyNorm = (envelope[peakTime]-detectionThresholds[0])/detectionThresholds[1]
            
            # discard spindle if spindles occurs during artefacted epochs
            if np.max(np.abs(envelope_high[int(np.max([spindleStartTime-sfreq,0])):int(np.min([spindleEndTime+sfreq,np.size(envelope_high)]))]))>30:
                continue

            # Compute Spindle candidate frequency using the power spectrum
            # (precision: 1Hz)
            temp_sig=this_bp_eeg_chan[int(peakTime-0.5*sfreq):int(peakTime+0.5*sfreq)]
            myPSDsFilt=np.real(np.fft.rfft(temp_sig))
            myPSDsFilt=myPSDsFilt[1:np.size(myPSDsFilt)-1]
            faxisFilt=np.arange(1/(np.size(temp_sig)/sfreq), sfreq/2, 1/(np.size(temp_sig)/sfreq))
            SPFilt=myPSDsFilt[np.intersect1d(np.where(faxisFilt>=minFreqValue),np.where(faxisFilt<=maxFreqValue))]
            SPFreq=faxisFilt[np.intersect1d(np.where(faxisFilt>=minFreqValue),np.where(faxisFilt<=maxFreqValue))]
            freqSpindle=SPFreq[np.where(SPFilt==np.max(SPFilt))][0]
            
            # compute ratio with alpha power
            PowerSP=np.mean((envelope_bp[spindleStartTime:spindleEndTime]));
            PowerAlpha=np.mean((envelope_alpha[spindleStartTime:spindleEndTime]));
            
            # Stage
            thisStage=scoring_ts[spindleStartTime]

            # Sum-up
            currentSpindle = (spindleStartTime, spindleEndTime, peakTime, peakEnergy, peakEnergyNorm, freqSpindle, spindleDuration, PowerSP, PowerAlpha, thisStage)
            currentSpindle = np.transpose(currentSpindle)
            
            # Sort Spindles
            # Candidate discarded: Peak Energy too highg
            if (peakEnergy > RejectThreshold):
                these_spindlesOutEnergy.append(currentSpindle)
                continue
            
            # Candidate discarded: Out Because of  Duration
            if (spindleDuration/sfreq > minMaxDurations[0]) & (spindleDuration/sfreq < minMaxDurations[1]):
                these_spindles.append(currentSpindle)
                nspin=nspin+1
            else:
                these_spindlesOutDuration.append(currentSpindle)
                continue
            lastSpindleStart=spindleStartTime
            
        sleepSpindles.append(these_spindles)
        spindlesOutDuration.append(these_spindlesOutDuration)
        spindlesOutEnergy.append(these_spindlesOutEnergy)

        print('Number of sleep spindles found: ')
        print(np.size(sleepSpindles[ch],0))

    sleepSpindles_savename = (f"{fig_dir}/allwaves/{sub_id}_{todaydate}.npy")
    spindlesOutEnergy_savename = (f"{fig_dir}/allwaves/{sub_id}_{todaydate}.npy")
    spindlesOutDuration_savename = (f"{fig_dir}/allwaves/{sub_id}_{todaydate}.npy")
    np.save(sleepSpindles_savename, np.asarray(sleepSpindles, dtype = object))
    np.save(spindlesOutEnergy_savename, np.asarray(spindlesOutDuration, dtype = object))
    np.save(spindlesOutDuration_savename, np.asarray(spindlesOutEnergy, dtype = object))


    return "\n\n...All files were correctly processed!\n\n"
             
if __name__ == '__main__':
    import glob
    # Get the list of EEG files
    eeg_files = glob.glob("/Volumes/DDE_ALC/PhD/EPISSE/CGC_Pilots/Preproc/*_concat_raw.fif")
    
    # Set up a pool of worker processes
    num_processes = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=8)
    
    # Process the EEG files in parallel
    pool.map(SS_detection, eeg_files)
    
    # Clean up the pool of worker processes
    pool.close()
    pool.join()
    
   
