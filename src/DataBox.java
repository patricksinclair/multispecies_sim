import java.util.ArrayList;

public class DataBox {

    private int runID;
    private int[] event_counters;
    private ArrayList<Double> times;
    private ArrayList<ArrayList<ArrayList<Double>>> all_microhab_pops;

    private double exit_time;

    public DataBox(int runID, int[] event_counters, ArrayList<Double> times, ArrayList<ArrayList<ArrayList<Double>>> all_microhab_pops){
        this.runID = runID;
        this.event_counters = event_counters;
        this.times = times;
        this.all_microhab_pops = all_microhab_pops;
    }

    public DataBox(int runID, double exit_time){
        // this constructor is used for the time to failure runs.
        this.runID = runID;
        this.exit_time = exit_time;
    }

    public int getRunID(){return runID;}
    public int[] getEvent_counters(){return event_counters;}
    public ArrayList<Double> getTimes(){return times;}
    public ArrayList<ArrayList<ArrayList<Double>>> getAll_microhab_pops(){return all_microhab_pops;}

    public double getExit_time(){return exit_time;}


}