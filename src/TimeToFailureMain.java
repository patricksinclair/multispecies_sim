public class TimeToFailureMain {
    public static void main(String[] args) {
        //This class is used as the main for the time to failure simulations.
        //Should make it easier to manage things between the big geno runs and the time to failure code.
        int nCores = Integer.parseInt(args[0]); //no. of cores used in parallel runs
        // want to do 1000 total reps where n_reps = nCores*nBlocks
        double c_max = 5.;
        int nBlocks = 20;

        //time to failure params. these are varied to get ttf vs param plots
        double r_imm = 20; //default value is 20. lets do 18 -> 22 in steps of 1. -done
        double g_max = 0.083; //default value is 0.083.



        //format of the Object[] params used to be [fileID, scale, sigma, c_max, r_imm]
        // We're now varying other params, easiest (but not neatest) way is likely to just have a timeToFailure_vs_param method for each param of interest
        Object[] ttf_14_resistant_params = new Object[]{"timeToFailure-14_pc_res", 2.703747953786337, 0.5690825284230452, c_max, r_imm};
        Object[] ttf_15_resistant_params = new Object[]{"timeToFailure-15_pc_res", 2.6133256846855746, 0.6260058161550592, c_max, r_imm};
        Object[] ttf_16_resistant_params = new Object[]{"timeToFailure-16_pc_res", 2.47772924764521, 0.7060073500033884, c_max, r_imm};
        //we'll use an object array to store [directory_ID, N*, r_det_ratio]
        Object[] phase2_params = new Object[]{"_phase2", 0.75, 0.22};

        BioSystem.timeToFailure_vs_r_imm(nCores, nBlocks, ttf_14_resistant_params, phase2_params);
        //todo - make sure the update biofilm size method has the failure limit check included. (actually changed the thickness limit arguments so this might not be necessary).
        //BioSystem.timeToFailure_vs_c_max(nCores, nBlocks, ttf_14_resistant_params, phase2_params);
    }
}
