cp data_mwk/SWDebug_20110607_RSVPNicole10x6_A1to1_RefAA_002.mwk data_merged_tmp/SWDebug_20110607_RSVPNicole10x6_A1to1_RefAA_002.mwk && bin/merge.py data_merged_tmp/SWDebug_20110607_RSVPNicole10x6_A1to1_RefAA_002.mwk /mindhive/dicarlolab/u/hahong/teleport/array_tmp/blackrock_default//SWDebug_20110607_RSVPNicole10x6_A1to1_RefAA_002.nev nowav   && rm -f data_merged_tmp/SWDebug_20110607_RSVPNicole10x6_A1to1_RefAA_002.mwk/*.bak && bin/collect_PS_firing.py data_merged_tmp/SWDebug_20110607_RSVPNicole10x6_A1to1_RefAA_002.mwk data_postproc_tmp/SWDebug_20110607_RSVPNicole10x6_A1to1_RefAA_002.psf.pk 300 128    exclude_img=circ_mask
