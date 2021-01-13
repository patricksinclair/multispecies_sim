import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

public class Toolbox {


    public static void writeGenosOverTimeToCSV(String directoryName, DataBox dataBox){
        //new method of saving the data for the big geno distb runs.
        //instead of putting everything in one big file, each run will have a dedicated sub-directory
        //In each of these sub-directories there will be several csv files which contain the genos in each microhabitat
        //at a specific timestep.
        //For simplicity, we'll write these csv files such that each row represents a microhabitat,
        //then we can transpose the dataframes in python later


        //All the data in a run is stored in a databox.
        //Here, for each databox passed to this method, we iterate through the times list and the 3D arraylist containing all the geno data
        //use outer index of arraylist to iterate through geno(t) and the times list
        //use times list for filenames

        try{

            //each run has its own directory so we add the run ID here
            //subDirectoryName += dataBox.getRunID();
            directoryName += "/runID_"+dataBox.getRunID();
            File directory = new File(directoryName);
            if(!directory.exists()) directory.mkdirs();

            //iterate over all the timesteps stored
            for(int t = 0; t < dataBox.getAll_microhab_pops().size(); t++) {

                //create a new file for each timestep
                //
                String file_name = "geno_distb"+String.format("-t=%.2f", dataBox.getTimes().get(t));
                File file = new File(directoryName+"/"+file_name+".csv");
                //if(!file.exists()) file.createNewFile();

                FileWriter fw = new FileWriter(file.getAbsoluteFile());
                BufferedWriter bw = new BufferedWriter(fw);

                //iterate over all the microhabitats, new line for each one
                for(int mh = 0; mh < dataBox.getAll_microhab_pops().get(t).size(); mh++) {
                    String geno_distb = "mh_"+mh;

                    //now iterate over all the genos in each microhab
                    for(int g = 0; g < dataBox.getAll_microhab_pops().get(t).get(mh).size(); g++) {
                        geno_distb += String.format(",%.4f", dataBox.getAll_microhab_pops().get(t).get(mh).get(g));
                    }
                    bw.write(geno_distb);
                    bw.newLine();

                }
                bw.close();
            }

        }catch (IOException e){}


    }


    static void writeDataboxEventCountersToFile(String directoryName, String filename, String[] headers, DataBox[] dataBoxes){



        File directory = new File(directoryName);
        if(!directory.exists()) directory.mkdirs();

        File file = new File(directoryName+"/"+filename+".csv");
        //if(!file.exists()) file.createNewFile();


        try{

            FileWriter fw = new FileWriter(file.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);

            int ncols = headers.length;
            int string_length = Math.max(12, Toolbox.largestHeaderLength(headers)+3);
            String file_header = headers[0];
            for(int h = 1; h < headers.length; h++){
                file_header += ","+headers[h];
            }
//
//            String head_start = "#"+headers[0]+",";
//            String file_header = String.format("%-"+string_length+"s", head_start);
//            for(int i = 1; i < headers.length-1; i++){
//                String heado = headers[i]+",";
//                file_header += String.format("%-"+string_length+"s", heado);
//            }
//            String heado = headers[headers.length-1];
//            file_header += String.format("%-"+string_length+"s", heado);
            bw.write(file_header);
            bw.newLine();


            for(int i = 0; i < dataBoxes.length; i++){

                //runID is also included in the event_counters array, so it's printed to file here too
                int[] event_counters = dataBoxes[i].getEvent_counters();
                String output = String.format("%d", event_counters[0]);

                for(int nc = 1; nc < ncols; nc++){
                    output += String.format(",%d", event_counters[nc]);
                }

                bw.write(output);
                bw.newLine();
            }
            bw.close();

        }catch (IOException e){}

    }

    static void writeTimeToFailureDataToFile(String directoryID, String fileID, String[] headers, DataBox[] dataBoxes){
        //Method to write the many time to failure times.  Two columns with runID and exit time.

        File directory = new File(directoryID);
        if(!directory.exists()) directory.mkdirs();

        File file = new File(directoryID+"/"+fileID+".csv");

        try{
            FileWriter fw = new FileWriter(file.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);

            //write the headers to the file
            String file_header = "";
            for(int i = 0; i < headers.length-1; i++){
                file_header += headers[i]+",";
            }
            file_header += headers[headers.length-1];
            bw.write(file_header);

            //write all data from seperate runs to the file
            //only 2 entries so do it manually
            for(DataBox db : dataBoxes){
                bw.newLine();
                String output = String.format("%d,%.3f", db.getRunID(), db.getExit_time());
                bw.write(output);
            }

            bw.close();

        }catch (IOException e){}

    }

    public static int averageArraylist(ArrayList<Integer> pops){
        //calculates the average value of the population size
        //can add this to the event counters
        double runningTotal = 0;
        for(int p : pops){
            runningTotal += p;
        }
        return (int)(runningTotal/(double)pops.size());
    }




    private static int largestHeaderLength(String[] headers){
        int biggun = 0;
        for(String s : headers){
            if(s.length() > biggun) biggun = s.length();
        }
        return biggun;
    }

    private static double averageOfArrayList(ArrayList<Double> listo){

        if(listo.size() > 0) {
            double sum = 0.;

            for(Double d : listo) {
                sum += d;
            }

            return sum/(double) listo.size();
        }else{
            return 0.;
        }
    }

    private static double[] averageAndStDevOfArray(double[] results){
        double sum = 0.;
        for(double d : results){
            sum += d;
        }
        double mean = sum/results.length;

        double sumSq = 0.;
        for(double d : results){
            sumSq += (d-mean)*(d-mean);
        }
        double stDev = Math.sqrt(sumSq/(results.length - 1.));

        return new double[]{mean, stDev};
    }






    public static String millisToShortDHMS(long duration) {
        String res = "";
        long days  = TimeUnit.MILLISECONDS.toDays(duration);
        long hours = TimeUnit.MILLISECONDS.toHours(duration)
                - TimeUnit.DAYS.toHours(TimeUnit.MILLISECONDS.toDays(duration));
        long minutes = TimeUnit.MILLISECONDS.toMinutes(duration)
                - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(duration));
        long seconds = TimeUnit.MILLISECONDS.toSeconds(duration)
                - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(duration));
        if (days == 0) {
            res = String.format("%02d:%02d:%02d", hours, minutes, seconds);
        }
        else {
            res = String.format("%dd%02d:%02d:%02d", days, hours, minutes, seconds);
        }
        return res;
    }


}