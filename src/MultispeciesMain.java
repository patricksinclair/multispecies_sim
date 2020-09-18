public class MultispeciesMain {

    public static void main(String[] args) {

        //todo - make sure alpha is changed to correspond with the new parameter values

        int nCores = Integer.parseInt(args[0]); //no. of cores used in parallel runs
        //changed nBlocks to 10, we'll do 10x10 runs to get more cores available
        int nBlocks = 10; //no. of times a parallel run is performed.  total no. of runs = nCores * nBlocks.
        //runID_offset is used to adjust the run ID for successive runs, so that it starts at runID_offset instead of 0
        int runID_offset = 0;
        String date = "-18-Sep-2020";

        //Depending on our choices of N* and r_det we will either be in phase 2 or 4 of the bftt phase diagram
        //Need to save our values in the corresponding results directory
        //we'll use an object array to store [directory_ID, N*, r_det_ratio]
        //weren't getting any results with the initial values we chose, so instead lets try the values that gave the thickest biofilm
        //phase2 params are now chosen for the lowest det rate
        Object[] phase2_params = new Object[]{"_phase2", 0.75, 0.22};
        Object[] phase4_params = new Object[]{"_phase4", 0.625, 0.7};

        //also need folder ID depending on our values of the geno distbs.
        //we'll use an object array to store [subDirectory_ID, scale, sigma]
        //these are the new updated versions for c_max = 6.
        Object[] params_8_resistant = new Object[]{"8_resistant"+date, 2.70825993, 0.56614558};
        Object[] params_10_resistant = new Object[]{"10_resistant"+date, 2.53711338, 0.6716238};
        Object[] params_12_resistant = new Object[]{"12_resistant"+date, 2.22826069, 0.84302476};


        BioSystem.getEventCountersAndRunPopulations(nCores, nBlocks, params_8_resistant, phase2_params);
        //BioSystem.timeToFailure(nCores, nBlocks, scale_93, sigma_93, folderID93);
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
}
