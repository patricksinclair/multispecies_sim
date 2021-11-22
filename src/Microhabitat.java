import org.apache.commons.math3.distribution.LogNormalDistribution;
import java.util.ArrayList;

class Microhabitat {

    private LogNormalDistribution MIC_distribution;

    private double c; //concn of antimicrobial
    private ArrayList<Double> population; //list of MICs of bacteria in microhab

    private int K; //karrying kapacity
    private boolean surface = false, biofilm_region, immigration_zone = false;
    private double g_max; //max growth rate =  2/day
    private double uniform_dRate = 0.018; //all bacteria have this death rate
    double biofilm_threshold; //fraction occupied needed to transition to biofilm


    Microhabitat(int K, double c, double scale, double sigma, double biofilm_threshold){
        this.K = K;
        this.c = c;
        double mu = Math.log(scale);
        this.population = new ArrayList<>(K);
        this.biofilm_threshold = biofilm_threshold;
        this.MIC_distribution = new LogNormalDistribution(mu, sigma);
        this.biofilm_region = false;

        this.g_max = 0.083;
    }

    Microhabitat(int K, double c, double scale, double sigma, double biofilm_threshold, double g_max){
        // this constructor is used to vary the value of g_max, for the ttf simulations
        this.K = K;
        this.c = c;
        double mu = Math.log(scale);
        this.population = new ArrayList<>(K);
        this.biofilm_threshold = biofilm_threshold;
        this.g_max = g_max;

        this.MIC_distribution = new LogNormalDistribution(mu, sigma);
        this.biofilm_region = false;
    }




    int getN(){
        return population.size();
    }

    boolean isSurface(){
        return surface;
    }

    boolean isBiofilm_region(){
        return biofilm_region;
    }

    boolean isImmigration_zone(){
        return immigration_zone;
    }

    ArrayList<Double> getPopulation(){
        return population;
    }

    void setSurface(){
        this.surface = true;
    }

    void setBiofilm_region(){
        this.biofilm_region = true;
    }

    void setImmigration_zone(boolean immigration_zone){
        this.immigration_zone = immigration_zone;
    }


    private double fractionFull(){
        return getN()/(double) K;
    }

    boolean atBiofilmThreshold(){
        return fractionFull() >= biofilm_threshold;
    }


    private double beta(int index){
        return population.get(index);
    }

    private double phi_c(int index){
        //pharmacodynamic function
        double cB = c/beta(index);
        return 1. - (6.*cB*cB)/(5. + cB*cB);
    }


    double[] replicationAndDeathRates(int index){
        //returns either the growth rate and the uniform death rate if the bacteria is resistant,
        //or the sums of the uniform and pharmacodyncamic death rates is the batceria is susceptible
        double phi_c_scaled = g_max *phi_c(index);
        double gRate = phi_c_scaled > 0. ? phi_c_scaled*(1. - (getN()/(double)K)) : 0.;
        double dRate = phi_c_scaled > 0. ? uniform_dRate : phi_c_scaled + uniform_dRate;

        return new double[]{gRate, dRate};
    }


    void addARandomBacterium_x_N(int n_bacteria){
        for(int i = 0; i < n_bacteria; i++) {
            population.add(MIC_distribution.sample());
        }
    }

    void replicateABacterium_x_N(int index, int nReps){
        for(int i = 0; i < nReps; i++) {
            population.add(population.get(index));
        }
    }

    void addABacterium(double MIC){
        population.add(MIC);
    }

    void removeABacterium(int index){
        population.remove(index);
    }


}
