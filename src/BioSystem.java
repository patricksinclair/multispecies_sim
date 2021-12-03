import org.apache.commons.math3.distribution.PoissonDistribution;
import java.util.ArrayList;
import java.util.Map;
import java.util.Random;
import java.util.stream.IntStream;

class BioSystem {

    private Random rand = new Random();
    //need to initialise these in the performAction method in order to incorporate the tau halving
    private PoissonDistribution poiss_immigration;
    private PoissonDistribution poiss_deterioration;
    private PoissonDistribution poiss_migration;
    private PoissonDistribution poiss_migration_edge;

    private double alpha, c_max; //steepness and max val of antimicrobial concn (currently 0.01 and 10)
    private double scale, sigma; //mic distb shape parameters
    private ArrayList<Microhabitat> microhabitats;

    //exit time is the time it took for the biofilm to reach the thickness limit, if it did
    //failure time is the time it took for the system to "fail", in this case form a single biofilm section
    private double time_elapsed, exit_time, failure_time;
    private int immigration_index;

    private double biofilm_threshold; //= 0.75; these are now encapsulated in the phase_N parameter arrays which are passed to the constructor
    private double deterioration_rate; // = 0.0168; these are now encapsulated in the phase_N parameter arrays which are passed to the constructor
    // these "final" params are now initialised in the constructor
    private final int K; // = 550;
    private final double g_max;// = 0.083; //maximum value of the growth rate (2 per day)
    private final double immigration_rate;// = 20.;
    private final double migration_rate = 1.;
    private final double delta_x = 1.; //thickness of a microhabitat in microns
    private double tau = 0.2; //much larger value now that the bug is fixed
    //this is how big the system can get before we exit. should reduce overall simulation duration
    private int thickness_limit;
    //this is how thick the biofilm can get before the system is deemed to have "failed"
    private int failure_limit;
    private int detachment_counts = 0, death_counts = 0, replication_counts = 0, immigration_counts = 0, migration_counts = 0, n_tauHalves = 0;


    private BioSystem(double alpha, double c_max, double biofilm_threshold, double deterioration_ratio, double scale, double sigma){

        this.alpha = alpha;
        this.c_max = c_max;
        this.scale = scale;
        this.sigma = sigma;
        this.microhabitats = new ArrayList<>();
        this.time_elapsed = 0.;
        this.exit_time = 0.;
        this.failure_time = 0.;
        this.immigration_index = 0;

        // variable parameters (varied in the ttf sims, kept constant otherwise)
        this.immigration_rate = 20.;
        this.g_max = 0.083;
        this.deterioration_rate = deterioration_ratio*g_max; //this is now in terms of the g_max ratio, as seen in the biofilm threshold theory stuff
        this.K = 550;
        this.biofilm_threshold = biofilm_threshold;

        //added these initialisers here so we can change the thickness limit with another constructor for the time to failure runs
        //value of 40 for the geno distb ones.
        this.thickness_limit = 40;
        this.failure_limit = 1;

        microhabitats.add(new Microhabitat(K, calc_C_i(0, this.c_max, this.alpha, this.delta_x), scale, sigma, this.biofilm_threshold));
        microhabitats.get(0).setSurface();
        microhabitats.get(0).addARandomBacterium_x_N(5);
    }

    private BioSystem(double alpha, double c_max, double biofilm_threshold, double deterioration_ratio, double scale, double sigma, int thickness_limit){
        //this constructor is used for the time to failure vs c_max sims.  thickness limit is set to 1 as any amount of biofilm is classed as a failure.
        this.alpha = alpha;
        this.c_max = c_max;
        this.scale = scale;
        this.sigma = sigma;
        this.microhabitats = new ArrayList<>();
        this.time_elapsed = 0.;
        this.exit_time = 0.;
        this.failure_time = 0.;
        this.immigration_index = 0;

        // variable parameters (varied in the ttf sims, kept constant otherwise)
        this.immigration_rate = 20.;
        this.g_max = 0.083;
        this.deterioration_rate = deterioration_ratio*g_max; //this is now in terms of the g_max ratio, as seen in the biofilm threshold theory stuff
        this.K = 550;
        this.biofilm_threshold = biofilm_threshold;
        //added these initialisers here so we can change the thickness limit with another constructor for the time to failure runs
        this.thickness_limit = thickness_limit;
        this.failure_limit = thickness_limit;

        microhabitats.add(new Microhabitat(K, calc_C_i(0, this.c_max, this.alpha, this.delta_x), scale, sigma, this.biofilm_threshold));
        microhabitats.get(0).setSurface();
        microhabitats.get(0).addARandomBacterium_x_N(5);
    }

    private BioSystem(double alpha, double c_max, double biofilm_threshold, double deterioration_ratio, double scale, double sigma, int thickness_limit, double immigration_rate){
        //this constructor is used for the time to failure vs r_imm sims.  thickness limit is set to 1 as any amount of biofilm is classed as a failure.
        this.alpha = alpha;
        this.c_max = c_max;
        //this.biofilm_threshold = biofilm_threshold;
        //this.deterioration_rate = deterioration_ratio*g_max; //this is now in terms of the g_max ratio, as seen in the biofilm threshold theory stuff
        this.scale = scale;
        this.sigma = sigma;
        this.microhabitats = new ArrayList<>();
        this.time_elapsed = 0.;
        this.exit_time = 0.;
        this.failure_time = 0.;
        this.immigration_index = 0;

        // variable parameters (varied in the ttf sims, kept constant otherwise)
        this.immigration_rate = immigration_rate;
        this.g_max = 0.083;
        this.deterioration_rate = deterioration_ratio*g_max; //this is now in terms of the g_max ratio, as seen in the biofilm threshold theory stuff
        this.K = 550;
        this.biofilm_threshold = biofilm_threshold;

        //added these initialisers here so we can change the thickness limit with another constructor for the time to failure runs
        this.thickness_limit = thickness_limit;
        this.failure_limit = thickness_limit;

        microhabitats.add(new Microhabitat(K, calc_C_i(0, this.c_max, this.alpha, this.delta_x), scale, sigma, this.biofilm_threshold));

        microhabitats.get(0).setSurface();
        microhabitats.get(0).addARandomBacterium_x_N(5);
    }

    private BioSystem(double alpha, double c_max, double biofilm_threshold, double deterioration_ratio, double scale, double sigma, int thickness_limit,
                      double immigration_rate, double g_max, int K){
        //this constructor is used for the time to failure vs various param sims.
        // it takes all the possible variable parameters so they can be

        // thickness limit is set to 1 as any amount of biofilm is classed as a failure.
        this.alpha = alpha;
        this.c_max = c_max;
        this.scale = scale;
        this.sigma = sigma;
        this.microhabitats = new ArrayList<>();
        this.time_elapsed = 0.;
        this.exit_time = 0.;
        this.failure_time = 0.;
        this.immigration_index = 0;

        // variable parameters (varied in the ttf sims, kept constant otherwise)
        // their values will be passed via an array in TimeToFailureMain
        this.immigration_rate = immigration_rate; //default 20.;
        this.g_max = g_max; //default 0.083;
        this.K = K; // default 550;
        this.deterioration_rate = deterioration_ratio*0.083; //this is now (hard coded 0.083) in terms of the g_max ratio, as seen in the biofilm threshold theory stuff
        this.biofilm_threshold = biofilm_threshold;

        //added these initialisers here so we can change the thickness limit with another constructor for the time to failure runs
        this.thickness_limit = thickness_limit;
        this.failure_limit = thickness_limit;

        // use the constructor that ACTUALLY changes g_max
        microhabitats.add(new Microhabitat(this.K, calc_C_i(0, this.c_max, this.alpha, this.delta_x), scale, sigma, this.biofilm_threshold, this.g_max));

        microhabitats.get(0).setSurface();
        microhabitats.get(0).addARandomBacterium_x_N(5);
    }


    private double getDeterioration_rate(){ return deterioration_rate; }
    private double getBiofilm_threshold(){ return biofilm_threshold; }
    private int getDetachment_counts(){ return detachment_counts; }
    private int getDeath_counts(){ return death_counts; }
    private int getReplication_counts(){ return replication_counts; }
    private int getImmigration_counts(){ return immigration_counts; }
    private int getMigration_counts(){ return migration_counts; }
    private double getTimeElapsed(){return time_elapsed;}
    private double getExit_time(){return exit_time;}
    private int getSystemSize(){return microhabitats.size();}
    private double getFailure_time(){return failure_time;}


    private int getTotalN(){
        int runningTotal = 0;
        for(Microhabitat m : microhabitats) {
            runningTotal += m.getN();
        }
        return runningTotal;
    }

    private int getBiofilmEdge(){
        int edgeIndex = 0;
        for(int i = 0; i < microhabitats.size(); i++){
            if(microhabitats.get(i).isBiofilm_region()) edgeIndex = i;
        }
        return edgeIndex;
    }

    private int getBiofilmThickness(){
        int thickness = 0;
        for(int i = 0; i < microhabitats.size(); i++){
            if(microhabitats.get(i).isBiofilm_region()) thickness = i+1;
        }
        return thickness;
    }

    private ArrayList<ArrayList<Double>> getMicrohabPopulations(){

        ArrayList<ArrayList<Double>> mh_pops = new ArrayList<>();

        for(Microhabitat m : microhabitats) {
            ArrayList<Double> mh_pop = new ArrayList<>();
            for(Double geno : m.getPopulation()) {
                mh_pop.add(geno);
            }
            mh_pops.add(mh_pop);
        }

        return mh_pops;
    }



    private void immigrate(int mh_index, int n_immigrants){
        microhabitats.get(mh_index).addARandomBacterium_x_N(n_immigrants);
    }




    private static double calc_C_i(int i, double c_max, double alpha, double delta_x){
        //updated so it's now the concentration in the middle of the microhabitat
        double index = (i+1.)/2.;
        return c_max*Math.exp(-alpha*index*delta_x);
    }


    public void migrate(int mh_index, int bac_index){

        double migrating_bac = microhabitats.get(mh_index).getPopulation().get(bac_index);
        microhabitats.get(mh_index).removeABacterium(bac_index);

        if(microhabitats.get(mh_index).isSurface()){
            microhabitats.get(mh_index+1).addABacterium(migrating_bac);
        }else if(microhabitats.get(mh_index).isImmigration_zone()){
            microhabitats.get(mh_index-1).addABacterium(migrating_bac);
        }else{
            if(rand.nextBoolean()){
                microhabitats.get(mh_index+1).addABacterium(migrating_bac);
            }else{
                microhabitats.get(mh_index-1).addABacterium(migrating_bac);
            }
        }
    }


    private void updateBiofilmSize(){
        //once the edge microhabitat is sufficiently populated, this adds another microhabitat onto the system list
        //which is then used as the immigration zone
        //TODO- if we're doing time to failure stuff, then add failure limit things into the getSystemSize() check at the end

        if(microhabitats.get(immigration_index).atBiofilmThreshold()){

            microhabitats.get(immigration_index).setBiofilm_region();
            microhabitats.get(immigration_index).setImmigration_zone(false);

            int i = microhabitats.size();
            microhabitats.add(new Microhabitat(K, BioSystem.calc_C_i(i, c_max, alpha, delta_x), scale, sigma, biofilm_threshold));
            immigration_index = i;
            microhabitats.get(immigration_index).setImmigration_zone(true);
        }

        //todo - use this condition for the species composition simulations
        //this stops sims going onn unnecessarily too long. if the biofilm reaches the thickness limit then we record the
        //time this happened at and move on
        if(immigration_index==thickness_limit){
            exit_time = time_elapsed;
            time_elapsed = 9e9; //this way the time elapsed is now way above the duration value, so the simulation will stop

        }

        //todo - use this condition for the time to failure simulations
        //todo - THIS IF SECTION IS NOW OUTDATED (i think) ONLY USE THE ABOVE ONE INSTEAD.  i think i solved the issues with just using immigration_index and thickness_limit
        //if immigration index is the same as the failure limit, then we also move on
//        if(getSystemSize() == thickness_limit || immigration_index == failure_limit) {
//            exit_time = time_elapsed;
//            failure_time = time_elapsed;
//            System.out.println("Testy");
//            time_elapsed = 9e9; //this way the time elapsed is now way above the duration value, so the simulation will stop
//        }
    }


    public void performAction(){

        double tau_step = tau;

        int system_size = microhabitats.size();
        int[][] replication_allocations;
        int[][] death_allocations;
        int[][] migration_allocations;
        int[] detachment_allocations;
        int[] original_popsizes;
        int n_immigrants;

        whileloop:
        while(true) {
            poiss_immigration = new PoissonDistribution(immigration_rate*tau_step);
            poiss_deterioration = new PoissonDistribution(deterioration_rate*tau_step);
            poiss_migration = new PoissonDistribution(migration_rate*tau_step);
            poiss_migration_edge = new PoissonDistribution(0.5*migration_rate*tau_step);

            replication_allocations = new int[system_size][];
            death_allocations = new int[system_size][];
            migration_allocations = new int[system_size][];
            original_popsizes = new int[system_size];
            detachment_allocations = new int[microhabitats.get(immigration_index).getN()];

            for(int mh_index = 0; mh_index < system_size; mh_index++) {

                //we iterate through all the bacteria and calculate the events which they'll experience
                int mh_pop = microhabitats.get(mh_index).getN();
                int[] n_replications = new int[mh_pop];
                int[] n_deaths = new int[mh_pop];
                int[] n_migrations = new int[mh_pop];

                for(int bac_index = 0; bac_index < mh_pop; bac_index++) {
                    ///////// REPLICATIONS AND DEATHS ///////////////////
                    double[] g_and_d_rate = microhabitats.get(mh_index).replicationAndDeathRates(bac_index);
                    double g_rate = g_and_d_rate[0], d_rate = Math.abs(g_and_d_rate[1]);

                    if(g_rate > 0.) {
                        PoissonDistribution poiss_replication = new PoissonDistribution(g_rate*tau_step);
                        poiss_replication.reseedRandomGenerator(rand.nextLong());
                        n_replications[bac_index] = poiss_replication.sample();
                    }

                    //d_rate is always > 0 due to inclusion of uniform death rate, so no need for the if statements
                    //seen in earlier versions
                    PoissonDistribution poiss_death = new PoissonDistribution(d_rate*tau_step);
                    poiss_death.reseedRandomGenerator(rand.nextLong());
                    n_deaths[bac_index] = poiss_death.sample();

                    //Bacteria cannot replicate and die in the same timestep
                    if(n_deaths[bac_index] > 0 && n_replications[bac_index] > 0){
                        n_tauHalves++;
                        tau_step /= 2;
                        continue whileloop;
                    }

                    //bacteria can't die twice, so need to handle this
                    if(n_deaths[bac_index] > 1) {
                        n_tauHalves++;
                        tau_step /= 2;
                        continue whileloop;
                    }


                    ///////// MIGRATIONS AND DETACHMENTS //////////////////////
                    //only non-dead bacteria can migrate or detach
                    if(n_deaths[bac_index] == 0) {

                        //firstly work out the migrations
                        //do edge cases and bulk, then do detachments and set detaching migrations to 0
                        //only do migrations if there's multiple microhabs
                        if(immigration_index > 0) {
                            if(mh_index == 0 || mh_index == immigration_index) {
                                n_migrations[bac_index] = poiss_migration_edge.sample();
                            } else {
                                n_migrations[bac_index] = poiss_migration.sample();
                            }
                            //check for double events
                            if(n_migrations[bac_index] > 1) {
                                //tau_halves_counter++;
                                tau_step /= 2.;
                                continue whileloop;
                            }
                        }

                        //Now do detachments
                        //detaching bacteria can't migrate
                        if(mh_index == immigration_index){
                            detachment_allocations[bac_index] = poiss_deterioration.sample();
                            //check for double events
                            if(detachment_allocations[bac_index] > 1) {
                                //tau_halves_counter++;
                                tau_step /= 2.;
                                continue whileloop;
                            }
                            //bacteria can only migrate if it's not detaching
                            if(detachment_allocations[bac_index] > 0) {
                                n_migrations[bac_index] = 0;
                            }

                        }

                    }
                    //////////////////////////////////////////////////////
                }
                replication_allocations[mh_index] = n_replications;
                death_allocations[mh_index] = n_deaths;
                migration_allocations[mh_index] = n_migrations;
                original_popsizes[mh_index] = microhabitats.get(mh_index).getN();
            }
            n_immigrants = poiss_immigration.sample();
            break whileloop;
        }


        //now we carry out the actions
        for(int mh_index = 0; mh_index < system_size; mh_index++){
            //iterate backwards over the bacteria so we can remove them without getting index errors
            for(int bac_index = original_popsizes[mh_index]-1; bac_index >= 0; bac_index--){

                if(death_allocations[mh_index][bac_index]!= 0) {
                    microhabitats.get(mh_index).removeABacterium(bac_index);
                    death_counts++;
                }

                else{
                    microhabitats.get(mh_index).replicateABacterium_x_N(bac_index, replication_allocations[mh_index][bac_index]);
                    replication_counts += replication_allocations[mh_index][bac_index];

                    if(system_size > 1){
                        if(migration_allocations[mh_index][bac_index] != 0) migrate(mh_index, bac_index);
                    }

                    if(mh_index == immigration_index){
                        if(detachment_allocations[bac_index] != 0) {
                            microhabitats.get(mh_index).removeABacterium(bac_index);
                            detachment_counts++;
                        }
                    }
                }
            }
        }

        immigrate(immigration_index, n_immigrants);
        immigration_counts += n_immigrants;
        updateBiofilmSize();
        //update the time elapsed in the system by the value of tau used in the final events
        time_elapsed += tau_step;

    }



    static void getEventCountersAndRunPopulations(int nCores, int nBlocks, Object[] suscep_params, Object[] phase_params, int runID_offset){
        //this is the method for the big runs to get the population distributions over time
        //it returns a csv file that for each run contains the event counters
        //it also returns a folder full of the bacteria distributions for each run, sampled at regular intervals
        //these are used to make those pink and blue plots of the geno distbs over space and time
        long startTime = System.currentTimeMillis();
        int nRuns = nCores*nBlocks; //total number of simulations performed
        int nSamples = 100; //no of samples taken during the runs
        double duration = 26.*7.*24.; //26 week duration
        //double duration = 512.; //testing duration

        DataBox[] dataBoxes = new DataBox[nRuns]; //this will be used to store all the data

        //unpack the Object arrays here
        String phaseID = String.valueOf(phase_params[0]);
        String suscepID_andDate = String.valueOf(suscep_params[0]);
        //headers used for the counters file
        String[] headers = new String[]{"runID", "bf_thickness", "exit_time", "final_pop", "avg_pop", "n_deaths", "n_detachments", "n_immigrations", "n_migrations", "n_replications"};

        double scale = (double)suscep_params[1];
        double sigma = (double)suscep_params[2];
        double biofilm_threshold = (double)phase_params[1];
        double deterioration_ratio = (double)phase_params[2];

        //this is the master directory where all of the geno_distb runs are kept
        String results_directory = "/Disk/ds-sopa-personal/s1212500/multispecies-sims/geno_distbs_results"+phaseID;
        //this is the directory inside the master directory where all the data for a given run is stored
        String run_directory = results_directory+"/"+suscepID_andDate;
        //this is the filename for the counters file (this contains the counters for all the runs in this batch)
        String counters_filename = suscepID_andDate+"-event_counters-sigma="+String.format("%.5f", sigma)+"-t="+String.format("%.1f", duration);


        for(int j = 0; j < nBlocks; j++){
            System.out.println("section: "+j);

            IntStream.range(j*nCores, (j+1)*nCores).parallel().forEach(i ->
                    dataBoxes[i] = getEventCountersAndRunPops_Subroutine(duration, nSamples, i, biofilm_threshold, deterioration_ratio, scale, sigma, runID_offset));

        }

        Toolbox.writeDataboxEventCountersToFile(run_directory, counters_filename, headers, dataBoxes);

        for (DataBox dataBox : dataBoxes) {

            Toolbox.writeGenosOverTimeToCSV(run_directory, dataBox);

        }


        long finishTime = System.currentTimeMillis();
        String diff = Toolbox.millisToShortDHMS(finishTime - startTime);
        System.out.println("results written to file");
        System.out.println("Time taken: "+diff);

    }

    private static DataBox getEventCountersAndRunPops_Subroutine(double duration, int nSamples, int runID, double biofilm_threshold, double deterioration_ratio, double scale, double sigma, int runID_offset){

        int K = 550;
        //todo - in order to get some sort of growth ocurring, c_max has been lowered from 10 -> 6 -> 5.
        //todo - this is roughly equal to 2x the avg mic
        double c_max = 5.0, alpha = 0.01;
        double interval = duration/nSamples;
        boolean alreadyRecorded = false;
        //added in the runID_offset, this should make it somewhat easier to combine runs from successive dates
        int runID_adjusted = runID+runID_offset;

        BioSystem bs = new BioSystem(alpha, c_max, biofilm_threshold, deterioration_ratio, scale, sigma);
        ArrayList<ArrayList<ArrayList<Double>>> mh_pops_over_time = new ArrayList<>();
        ArrayList<Double> times = new ArrayList<>();
        ArrayList<Integer> totalN_overTime = new ArrayList<>(); //this is used to track the pop size over time, can then average it at the end and save as a counter

        while(bs.time_elapsed <= duration+0.02*interval){

            if((bs.getTimeElapsed()%interval <= 0.02*interval) && !alreadyRecorded){

                int max_poss_pop = bs.getSystemSize()*K;
                System.out.println("runID: "+runID_adjusted+"\tt: "+bs.getTimeElapsed()+"\tpop size: "+bs.getTotalN()+"/"+max_poss_pop+"\tbf_edge: "+bs.getBiofilmEdge()+"\tN*: "+bs.biofilm_threshold+"\tdet_r: "+bs.deterioration_rate);
                alreadyRecorded = true;

                times.add(bs.getTimeElapsed());
                mh_pops_over_time.add(bs.getMicrohabPopulations());
                totalN_overTime.add(bs.getTotalN());

            }
            if(bs.getTimeElapsed()%interval >= 0.1*interval) alreadyRecorded = false;

            bs.performAction();
        }
        if((int)bs.exit_time == 0) bs.exit_time = duration;


        int[] event_counters = new int[]{runID_adjusted, bs.getBiofilmThickness(), (int)bs.getExit_time(), bs.getTotalN(), Toolbox.averageArraylist(totalN_overTime), bs.getDeath_counts(), bs.getDetachment_counts(), bs.getImmigration_counts(), bs.getMigration_counts(), bs.getReplication_counts()};

        return new DataBox(runID_adjusted, event_counters, times, mh_pops_over_time);
    }

    static void timeToFailure_vs_r_imm(int nCores, int nReps, Object[] suscep_params, Object[] phase_params){
        //TODO - MAKE SURE THE SECTION IN THE UPDATE BIOFILM SIZE METHOD REGARDING FAILURE LIMIT IS UNCOMMENTED (THE SECOND IF STATEMENT).
        long startTime = System.currentTimeMillis();
        //this version of the time to failure routine is used to investigate how changing the immigration rate affects the time to failure.
        String results_directory = "/Disk/ds-sopa-personal/s1212500/multispecies-sims/time_to_failure_vs_r_imm";
        String file_ID = suscep_params[0]+String.format("-r_imm=%.2f", suscep_params[4])+"-session2";
        String[] headers = new String[]{"runID", "failure_time"};

        double duration = 365.*24.; //1 year duration.

        int nRuns = nCores*nReps;
        DataBox[] dataBoxes = new DataBox[nRuns];

        for(int j = 0; j < nReps; j++){
            IntStream.range(j*nCores, (j+1)*nCores).parallel().forEach(i ->
                    dataBoxes[i] = timeToFailure_vs_r_Imm_subroutine(duration, i, suscep_params, phase_params));
        }

        Toolbox.writeTimeToFailureDataToFile(results_directory, file_ID, headers, dataBoxes);

        long finishTime = System.currentTimeMillis();
        String diff = Toolbox.millisToShortDHMS(finishTime - startTime);
        System.out.println("results written to file");
        System.out.println("Time taken: "+diff);

    }

    private static DataBox timeToFailure_vs_r_Imm_subroutine(double duration, int i, Object[] suscep_params, Object[] phase_params){
        //extract parameters from the arrays
        double r_imm = (double)suscep_params[4];
        double c_max = (double)suscep_params[3];
        double alpha = 0.01; //this stays fixed.
        double scale = (double)suscep_params[1];
        double sigma = (double)suscep_params[2];
        double biofilm_threshold = (double)phase_params[1];
        double deterioration_ratio = (double)phase_params[2];
        int thickness_limit = 1; //immigration index limit, sim ends when it gets to this value

        BioSystem bs = new BioSystem(alpha, c_max, biofilm_threshold, deterioration_ratio, scale, sigma, thickness_limit, r_imm);
        System.out.println("run: "+i+"\t r_imm: "+bs.immigration_rate);
        //todo - took out the interval thing for now as it doesn't hugely matter for this method if we're not sampling over time
        while(bs.time_elapsed <= duration){
            //don't bother with the sampling things for now

            bs.performAction();
        }
        if((int)bs.exit_time == 0) bs.exit_time = duration;

        return new DataBox(i, bs.exit_time);
    }


    static void timeToFailure_vs_c_max(int nCores, int nReps, Object[] suscep_params, Object[] phase_params){
        //TODO - MAKE SURE THE SECTION IN THE UPDATE BIOFILM SIZE METHOD REGARDING FAILURE LIMIT IS UNCOMMENTED (THE SECOND IF STATEMENT).
        long startTime = System.currentTimeMillis();
        // Updated version of the old time to failure routine.  Here we will see how changing c_max affects the rate at which
        // the biofilm can form.  We'll do several values of c_max for each of the current percent resistant.
        // Will need to be careful with wording here.  As the x% resistant corresponds to a lognorm distb with c_max = 5,
        // varying c_max means that we'll need to be careful with terminology.  Maybe refer to % resistance as small, medium, large in the paper.

        // Method takes in no. of cores, no. of reps per core and an array containing the system parameters

        String results_directory = "/Disk/ds-sopa-personal/s1212500/multispecies-sims/time_to_failure_vs_cmax";
        String file_ID = suscep_params[0]+String.format("-c_max=%.2f", suscep_params[3])+"session-2";
        String[] headers = new String[]{"runID", "failure_time"};

        //double duration = 26.*7.*24.; //26 week duration
        double duration = 365.*24.; //1 year duration

        int nRuns = nCores*nReps;
        DataBox[] dataBoxes = new DataBox[nRuns];

        for(int j = 0; j < nReps; j++){

            IntStream.range(j*nCores, (j+1)*nCores).parallel().forEach(i ->
                    dataBoxes[i] = timeToFailure_vs_c_max_subroutine(duration, i, suscep_params, phase_params));
        }

        Toolbox.writeTimeToFailureDataToFile(results_directory, file_ID, headers, dataBoxes);

        long finishTime = System.currentTimeMillis();
        String diff = Toolbox.millisToShortDHMS(finishTime - startTime);
        System.out.println("results written to file");
        System.out.println("Time taken: "+diff);

    }

    private static DataBox timeToFailure_vs_c_max_subroutine(double duration, int i, Object[] suscep_params, Object[] phase_params){

        //extract parameters from the arrays
        double c_max = (double)suscep_params[3];
        double alpha = 0.01; //this stays fixed.
        double scale = (double)suscep_params[1];
        double sigma = (double)suscep_params[2];
        double biofilm_threshold = (double)phase_params[1];
        double deterioration_ratio = (double)phase_params[2];
        int thickness_limit = 1; //immigration index limit, sim ends when it gets to this value

        BioSystem bs = new BioSystem(alpha, c_max, biofilm_threshold, deterioration_ratio, scale, sigma, thickness_limit);
        System.out.println("run: "+i);
        //todo - took out the interval thing for now as it doesn't hugely matter for this method if we're not sampling over time
        while(bs.time_elapsed <= duration){
            //don't bother with the sampling things for now
//            if((bs.getTimeElapsed()%interval <= 0.02*interval) && !alreadyRecorded){
//
//                int max_poss_pop = bs.getSystemSize()*K;
//                System.out.println("runID: "+runID_adjusted+"\tt: "+bs.getTimeElapsed()+"\tpop size: "+bs.getTotalN()+"/"+max_poss_pop+"\tbf_edge: "+bs.getBiofilmEdge()+"\tN*: "+bs.biofilm_threshold+"\tdet_r: "+bs.deterioration_rate);
//                alreadyRecorded = true;
//
//                times.add(bs.getTimeElapsed());
//                mh_pops_over_time.add(bs.getMicrohabPopulations());
//                totalN_overTime.add(bs.getTotalN());
//
//            }
//            if(bs.getTimeElapsed()%interval >= 0.1*interval) alreadyRecorded = false;

            bs.performAction();
        }
        if((int)bs.exit_time == 0) bs.exit_time = duration;

        return new DataBox(i, bs.exit_time);
    }



    static void timeToFailure_vs_x_param(Map model_params){
        //TODO - MAKE SURE THE SECTION IN THE UPDATE BIOFILM SIZE METHOD REGARDING FAILURE LIMIT IS UNCOMMENTED (THE SECOND IF STATEMENT).
        long startTime = System.currentTimeMillis();
        // Updated updated version of the old time to failure routine.
        // Here we pass all the model params in a Map object, including a string which represents which param we want to vary.
        // Aiming for 1000 reps

        int nCores = (int)model_params.get("n_cores"), nReps = (int)model_params.get("n_reps");

        // Method takes in no. of cores, no. of reps per core and an array containing the system parameters
        String file_prefix = (String) model_params.get("file_prefix");
        String varied_param_string = (String) model_params.get("varied_param_key"); // Map key of the varied param, also used for file naming
        double varied_param_val = (double) model_params.get(varied_param_string); // numerical value of the model param being varied
        //double varied_param_val = (Integer) model_params.get(varied_param_string); // todo - change (Integer) back to (double) for non-K values

        String results_directory = "/Disk/ds-sopa-personal/s1212500/multispecies-sims/time_to_failure_vs_"+varied_param_string;
        //String results_directory = "testo";
        String file_ID = file_prefix+"-"+varied_param_string+String.format("=%.3f", varied_param_val);
        String[] headers = new String[]{"runID", "failure_time"};

        //double duration = 26.*7.*24.; //26 week duration
        double duration = 365.*24.; //1 year duration

        int nRuns = nCores*nReps;
        DataBox[] dataBoxes = new DataBox[nRuns];

        for(int j = 0; j < nReps; j++){

            IntStream.range(j*nCores, (j+1)*nCores).parallel().forEach(i ->
                    dataBoxes[i] = timeToFailure_vs_x_param_subroutine(duration, i, model_params));
        }

        Toolbox.writeTimeToFailureDataToFile(results_directory, file_ID, headers, dataBoxes);

        long finishTime = System.currentTimeMillis();
        String diff = Toolbox.millisToShortDHMS(finishTime - startTime);
        System.out.println("results written to file");
        System.out.println("Time taken: "+diff);

    }

    private static DataBox timeToFailure_vs_x_param_subroutine(double duration, int i, Map model_params){

        //extract parameters from the arrays
        double c_max = (double)model_params.get("c_max");
        double alpha = (double)model_params.get("alpha"); //0.01; //this stays fixed.
        double scale = (double)model_params.get("scale");
        double sigma = (double)model_params.get("sigma");
        double biofilm_threshold = (double)model_params.get("biofilm_threshold");
        double deterioration_ratio = (double)model_params.get("deterioration_ratio");
        double immigration_rate = (double)model_params.get("r_imm");
        double g_max = (double)model_params.get("g_max");
        int K = (int)model_params.get("K");

        int thickness_limit = 1; //immigration index limit, sim ends when it gets to this value

        //BioSystem bs = new BioSystem(alpha, c_max, biofilm_threshold, deterioration_ratio, scale, sigma, thickness_limit);
        BioSystem bs = new BioSystem(alpha, c_max, biofilm_threshold, deterioration_ratio, scale, sigma, thickness_limit, immigration_rate, g_max, K);
        System.out.println("run: "+i);
        while(bs.time_elapsed <= duration){
            //don't bother with the sampling things for now
            bs.performAction();
        }
        if((int)bs.exit_time == 0) bs.exit_time = duration;
        return new DataBox(i, bs.exit_time);
    }






}


