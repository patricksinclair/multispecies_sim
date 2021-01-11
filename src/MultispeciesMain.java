public class MultispeciesMain {

    public static void main(String[] args) {

        //todo - make sure alpha is changed to correspond with the new parameter values

        int nCores = Integer.parseInt(args[0]); //no. of cores used in parallel runs
        //changed nBlocks to 10, we'll do 10x10 runs to get more cores available
        //todo - for 14% resistance, make sure no. of cores is 10 (100 runs at a time)
        //todo - for 15% resistance, make sure no. of cores is 25 (25 runs at a time)
        //todo - for 16% resistance, make sure no. of cores is 25 (25 runs at a time)
        //todo - this will require a different runID offset
        //no. of times a parallel run is performed.  total no. of runs = nCores * nBlocks.
        int nBlocks_14 = 10; //10 runs on each of 10 cores
        int nBlocks_15 = 1; //1 run on each of 25 cores
        int nBlocks_16 = 1; //1 run on each of 25 cores
        //runID_offset is used to adjust the run ID for successive runs, so that it starts at runID_offset instead of 0
        //todo - make sure the runID_offset is correct (add on 100 for 14, 25 for 15%, 25 for 16%) - check the end of the counters dataframe
        int runID_offset_14 = 850; //session 8 WE NOW HAVE 1100 RUNS FOR 14%, SO CAN STOP RUNNING THESE FOR NOW
        //todo - trying these on the 24 hour queue for now
        int runID_offset_15 = 425; //session 10
        int runID_offset_16 = 300; //session 10
        String date = "-11-Jan-2021"; //session 10

        //Depending on our choices of N* and r_det we will either be in phase 2 or 4 of the bftt phase diagram
        //Need to save our values in the corresponding results directory
        //weren't getting any results with the initial values we chose, so instead lets try the values that gave the thickest biofilm
        //phase2 params are now chosen for the lowest det rate
        //we'll use an object array to store [directory_ID, N*, r_det_ratio]
        Object[] phase2_params = new Object[]{"_phase2", 0.75, 0.22};
        Object[] phase4_params = new Object[]{"_phase4", 0.625, 0.7};

        //also need folder ID depending on our values of the geno distbs.
        //we'll use an object array to store [subDirectory_ID, scale, sigma]
        //third attempt to get growth to occur. now c_max is set to 5.
        //Object[] params_14_resistant = new Object[]{"14_resistant"+date, 2.703747953786337, 0.5690825284230452};
        Object[] params_15_resistant = new Object[]{"15_resistant"+date, 2.6133256846855746, 0.6260058161550592};
        Object[] params_16_resistant = new Object[]{"16_resistant"+date, 2.47772924764521, 0.7060073500033884};

        //BioSystem.getEventCountersAndRunPopulations(nCores, nBlocks_14, params_14_resistant, phase2_params, runID_offset_14);
        BioSystem.getEventCountersAndRunPopulations(nCores, nBlocks_15, params_15_resistant, phase2_params, runID_offset_15);
        //BioSystem.getEventCountersAndRunPopulations(nCores, nBlocks_16, params_16_resistant, phase2_params, runID_offset_16);
        //BioSystem.timeToFailure(nCores, nBlocks, scale_93, sigma_93, folderID93);

        //time to failure params
        //[fileID, scale, sigma, c_max]
//        double c_max = 5.;
//        int nBlocks = 20;
//        Object[] ttf_14_resistant_params = new Object[]{"timeToFailure-14_pc_res", 2.703747953786337, 0.5690825284230452, c_max};
//        Object[] ttf_15_resistant_params = new Object[]{"timeToFailure-15_pc_res", 2.6133256846855746, 0.6260058161550592, c_max};
//        Object[] ttf_16_resistant_params = new Object[]{"timeToFailure-16_pc_res", 2.47772924764521, 0.7060073500033884, c_max};
//        //todo - make sure the update biofilm size method has the failure limit check included. (actually changed the thickness limit arguments so this might not be necessary).
//        BioSystem.timeToFailure_vs_c_max(nCores, nBlocks, ttf_14_resistant_params, phase2_params);
    }

//    these are the outdated values for the distributions from when c_max was 10.
//    now that c_max is 6, and the percentage resistant are slightly different
//    Object[] params_99_suscep = new Object[]{"99_suscep"+date, 2.71760274, 0.56002833};
//    Object[] params_98_suscep = new Object[]{"98_suscep"+date, 2.54669037, 0.66599239};
//    Object[] params_97_suscep = new Object[]{"97_suscep"+date, 2.37276256, 0.76486794};
//    Object[] params_96_suscep = new Object[]{"96_suscep"+date, 2.17434104, 0.87159677};
//    Object[] params_95_suscep = new Object[]{"95_suscep"+date, 1.9246899,  1.00179994};
//    Object[] params_94_suscep = new Object[]{"94_suscep"+date, 1.54590048, 1.20080013};
//    Object[] params_93_suscep = new Object[]{"93_suscep"+date, 1.01073016, 1.51389233};

    //these are the new updated versions for c_max = 6.
    //these were also found to be insufficient for growth to occur, so we'll lower c_max again
//    Object[] params_8_resistant = new Object[]{"8_resistant"+date, 2.70825993, 0.56614558};
//    Object[] params_10_resistant = new Object[]{"10_resistant"+date, 2.53711338, 0.6716238};
//    Object[] params_12_resistant = new Object[]{"12_resistant"+date, 2.22826069, 0.84302476};
}

