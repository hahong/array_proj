cp data/d002_Tito/mwk/Tito_20110822_RSVPVar06_nopos_S110720A_001.mwk data/d002_Tito/mwk_merged/Tito_20110822_RSVPVar06_nopos_S110720A_001.mwk && bin/merge.py data/d002_Tito/mwk_merged/Tito_20110822_RSVPVar06_nopos_S110720A_001.mwk data/d002_Tito/neudat/Tito_20110822_RSVPVar06_nopos_S110720A_001.nev nowav   && rm -f data/d002_Tito/mwk_merged/Tito_20110822_RSVPVar06_nopos_S110720A_001.mwk/*.bak && bin/collect_PS_firing.py data/d002_Tito/mwk_merged/Tito_20110822_RSVPVar06_nopos_S110720A_001.mwk data/d002_Tito/data_postproc/Tito_20110822_RSVPVar06_nopos_S110720A_001.psf.pk 300 110    ch_shift=20110720A extinfo c_success=images_shown ign_unregistered
