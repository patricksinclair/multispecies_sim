import java.util.HashMap;
import java.util.Map;

public class TimeToFailureMain {
    public static void main(String[] args) {
        //This class is used as the main for the time to failure simulations.
        //Should make it easier to manage things between the big geno runs and the time to failure code.
        int nCores = Integer.parseInt(args[0]); //no. of cores used in parallel runs
        //int nCores = 1;
        // want to do 1000 total reps where n_reps = nCores*nReps (increased to 2000)
        int nReps = 100; // changed 50 -> 100

        // log-norm distribution params
        final double scale_14pcres = 2.703747953786337, sigma_14pcres = 0.5690825284230452;
        final double scale_15pcres = 2.613325684685574, sigma_15pcres = 0.6260058161550592;
        final double scale_16pcres = 2.47772924764521,  sigma_16pcres = 0.7060073500033884;

        //time to failure params. these are varied to get ttf vs param plots
        String varied_param_key = "biofilm_threshold"; // the Map key corresponding to the parameter being varied todo - change the (Integer) flag back for non-K simulations
        double c_max = 5.; // default is 5
        double alpha = 0.01; //slope of biocide gradient. default value is 0.01
        double r_imm = 20; //immigration rate. default value is 20.  (do 18 -> 22 in steps of 1. -done)
        double g_max = 0.083; //max growth rate. default value is 0.083. (do 0.073 -> 0.093 in steps of 0.005)
        int K = 550; //carrying capacity. default value is 550 (do 500 -> 600 in steps of 25 -done)
        double biofilm_threshold = 0.65; //biofilm formation density (N = biofilm_threshold*K). default is 0.75 for phase 2 params
        double deterioration_ratio = 0.22; // r_det = deterioration_ratio*g_max. default is 0.22 for phase 2 params (do 0.2 -> 0.24 in steps of 0.01 -done)

        //we'll use an object array to store [directory_ID, N*, r_det_ratio]. default parameter values [0.75, 0.22]
        //Object[] phase2_params = new Object[]{"_phase2", biofilm_threshold, 0.22};

        // we'll use Maps to store the variables instead of this Object array nonsense
        // just do the 14 pc res params for now, can modify the keys etc later if needs be
        Map<String, Object> param_map = new HashMap<>();
        param_map.put("file_prefix", "timeToFailure-14_pc_res"); // file prefix for the ttf results
        param_map.put("phase_label", "_phase2"); // used in directory naming scheme, mainly obsolete
        param_map.put("varied_param_key", varied_param_key); // the key used to access the parameter that's being varied (nneds to match the param keys used in the map)
        param_map.put("n_cores", nCores); // no. of cores used in parallel processing
        param_map.put("n_reps", nReps); // no. of reps per core
        param_map.put("scale", scale_14pcres); // lognorm scale param (14 % resistant)
        param_map.put("sigma", sigma_14pcres); // lognorm sigma param (14 % resistant)
        param_map.put("c_max", c_max); // max biocide concn
        param_map.put("alpha", alpha);
        param_map.put("r_imm", r_imm); // immigration rate
        param_map.put("g_max", g_max); // max growth rate
        param_map.put("K", K); // carrying capacity
        param_map.put("biofilm_threshold", biofilm_threshold);
        param_map.put("deterioration_ratio", deterioration_ratio);

        BioSystem.timeToFailure_vs_x_param(param_map);


    }
}
