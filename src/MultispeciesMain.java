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
        //todo - accidentally did 250 runs on the 25-Nov-2020 14% run, make sure the runID is updated correctly accordingly
        //todo - do 150 runs for session 7 of 14%, this should balance it out.
        int runID_offset_14 = 750; //session 7
        int runID_offset_15 = 350; //session 7
        int runID_offset_16 = 225; //session 7
        String date = "-03-Dec-2020"; //session 7

        //Depending on our choices of N* and r_det we will either be in phase 2 or 4 of the bftt phase diagram
        //Need to save our values in the corresponding results directory
        //we'll use an object array to store [directory_ID, N*, r_det_ratio]
        //weren't getting any results with the initial values we chose, so instead lets try the values that gave the thickest biofilm
        //phase2 params are now chosen for the lowest det rate
        Object[] phase2_params = new Object[]{"_phase2", 0.75, 0.22};
        Object[] phase4_params = new Object[]{"_phase4", 0.625, 0.7};

        //also need folder ID depending on our values of the geno distbs.
        //we'll use an object array to store [subDirectory_ID, scale, sigma]
        //third attempt to get growth to occur. now c_max is set to 5.
        Object[] params_14_resistant = new Object[]{"14_resistant"+date, 2.703747953786337, 0.5690825284230452};
        Object[] params_15_resistant = new Object[]{"15_resistant"+date, 2.6133256846855746, 0.6260058161550592};
        Object[] params_16_resistant = new Object[]{"16_resistant"+date, 2.47772924764521, 0.7060073500033884};

        //BioSystem.getEventCountersAndRunPopulations(nCores, nBlocks_14, params_14_resistant, phase2_params, runID_offset_14);
        //BioSystem.getEventCountersAndRunPopulations(nCores, nBlocks_15, params_15_resistant, phase2_params, runID_offset_15);
        BioSystem.getEventCountersAndRunPopulations(nCores, nBlocks_16, params_16_resistant, phase2_params, runID_offset_16);
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

    //these are the new updated versions for c_max = 6.
    //these were also found to be insufficient for growth to occur, so we'll lower c_max again
//    Object[] params_8_resistant = new Object[]{"8_resistant"+date, 2.70825993, 0.56614558};
//    Object[] params_10_resistant = new Object[]{"10_resistant"+date, 2.53711338, 0.6716238};
//    Object[] params_12_resistant = new Object[]{"12_resistant"+date, 2.22826069, 0.84302476};
}
//    static void timeToNthMicrohabPhaseDiagram(Object[] params, int nCores, int microhab_lim){
//        //this method is used to make a colour plot of the time taken to reach the Nth microhabitat as
//        //a function of N* and r_det
//        //save all the values for each of the runs in seperate dataframes, then do averaging etc later manually
//        //todo - make sure the failure limit thing is set up correctly
//
//        String[] headers = new String[]{"n_thresh", "det_ratio", "time_to_n", "time_elapsed"};
//        String results_directory = "/Disk/ds-sopa-personal/s1212500/multispecies-sims/biofilm_threshold_theory/"+params[0];
//
//        double duration = 1000.;
//        int nMeasurements = 20;
//        double N_thresh_min = 0., N_thresh_max = 1.5;
//        double N_thresh_increment = (N_thresh_max - N_thresh_min)/nMeasurements;
//        double r_det_ratio_min = 0., r_det_ratio_max = 1.5;
//        double r_det_ratio_increment = (r_det_ratio_max - r_det_ratio_min)/nMeasurements;
//        //ArrayList<DataBox[]> dataBoxes = new ArrayList<>();
//
//        for(int nt = 0; nt <= nMeasurements; nt++){
//            for(int dr = 0; dr <= nMeasurements; dr++){
//                double n_thresh = N_thresh_min + (nt*N_thresh_increment);
//                double det_ratio = r_det_ratio_min + (dr*r_det_ratio_increment);
//
//                DataBox[] dataBoxes = timeToNthMicrohabPhaseDiagram_subroutine(nCores, duration, params, n_thresh, det_ratio, microhab_lim);
//                String filename = String.format("mhLim-%d_N^-%.3f_rDet-%.3f", microhab_lim, n_thresh, det_ratio);
//                Toolbox.writeTimeToNthMicrohabDataToFile(results_directory, filename, headers, dataBoxes);
//            }
//        }
//    }
