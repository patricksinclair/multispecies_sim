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

        //todo - check how the data is stored in the simulation here
        //All the data in a run is stored in a databox.
        //Here, for each databox passed to this method, we iterate through the times list and the 3D arraylist containing all the geno data
        //use outer index of arraylist to iterate through geno(t) and the times list
        //use times list for filenames

        try{

            //each run has its own directory so we add the run ID here
            //subDirectoryName += dataBox.getRunID();
            directoryName += dataBox.getRunID();
            File directory = new File(directoryName);
            if(!directory.exists()) directory.mkdirs();

            //iterate over all the timesteps stored
            for(int t = 0; t < dataBox.getAll_microhab_pops().size(); t++) {

                //create a new file for each timestep
                //
                String file_name = "geno_distb"+String.format("-t=%.5E", dataBox.getTimes().get(t));
                File file = new File(directoryName+"/"+file_name+".csv");
                //if(!file.exists()) file.createNewFile();

                FileWriter fw = new FileWriter(file.getAbsoluteFile());
                BufferedWriter bw = new BufferedWriter(fw);

                //iterate over all the microhabitats, new line for each one
                for(int mh = 0; mh < dataBox.getAll_microhab_pops().get(t).size(); mh++) {
                    String geno_distb = "mh_"+mh;

                    //now iterate over all the genos in each microhab
                    for(int g = 0; g < dataBox.getAll_microhab_pops().get(t).get(mh).size(); g++) {
                        geno_distb += String.format(",%.5E", dataBox.getAll_microhab_pops().get(t).get(mh).get(g));
                    }
                    bw.write(geno_distb);
                    bw.newLine();

                }
                bw.close();
            }

        }catch (IOException e){}


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
            sumSq +=(d-mean)*(d-mean);
        }
        double stDev = Math.sqrt(sumSq/(results.length - 1.));

        return new double[]{mean, stDev};
    }






    public static void writeHistoArrayToFile(String filename, int[] inputData){
        try {
            File file = new File(filename+".txt");
            if(!file.exists()) file.createNewFile();

            FileWriter fw = new FileWriter(file.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);

            int nReadings = inputData.length;

            for(int i = 0; i < nReadings-1; i++){

                String output = String.format("%d", inputData[i]);

                bw.write(output);
                bw.newLine();
            }
            String output = String.format("%d", inputData[nReadings-1]);
            bw.write(output);
            bw.close();
        }catch (IOException e){}
    }










    public static void writeMultipleColumnsToFile(String filename, String[] headers, double[][] results){

        try{
            File file = new File(filename+".txt");
            if(!file.exists()) file.createNewFile();

            FileWriter fw = new FileWriter(file.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);

            int ncols = headers.length;
            int string_length = Math.max(12, Toolbox.largestHeaderLength(headers)+3);
            String head_start = "#"+headers[0]+",";
            String file_header = String.format("%-"+string_length+"s", head_start);
            for(int i = 1; i < headers.length-1; i++){
                String heado = headers[i]+",";
                file_header += String.format("%-"+string_length+"s", heado);
            }
            String heado = headers[headers.length-1];
            file_header += String.format("%-"+string_length+"s", heado);
            bw.write(file_header);
            bw.newLine();


            for(int i = 0; i < results[0].length; i++){

                String output = "";

                for(int nc = 0; nc < ncols-1; nc++){
                    String num_val = String.format("%.4E", results[nc][i])+",";
                    output += String.format("%-"+string_length+"s", num_val);
                }
                String num_val = String.format("%.4E", results[ncols-1][i]);
                output += String.format("%-"+string_length+"s", num_val);

                bw.write(output);
                bw.newLine();
            }
            bw.close();

        }catch (IOException e){}
    }




    static void writeCountersToFile(String filename, String[] headers, int[][] counters){
        try{
            File file = new File(filename+".txt");
            if(!file.exists()) file.createNewFile();

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


            for(int i = 0; i < counters.length; i++){

                String output = String.format("%d", counters[i][0]);

                for(int nc = 0; nc < ncols; nc++){
                    output += String.format(",%d", counters[i][nc]);
                }


//                for(int nc = 0; nc < ncols-1; nc++){
//                    String num_val = String.format("%d", counters[i][nc])+",";
//                    output += String.format("%-"+string_length+"s", num_val);
//                }
//                String num_val = String.format("%d", counters[i][ncols-1]);
//                output += String.format("%-"+string_length+"s", num_val);

                bw.write(output);
                bw.newLine();
            }
            bw.close();

        }catch (IOException e){}

    }


    static void writeDataboxEventCountersToFile(String directoryName, String filename, String[] headers, DataBox[] dataBoxes){



        File directory = new File(directoryName);
        if(!directory.exists()) directory.mkdirs();

        File file = new File(directoryName+"/"+filename+".txt");
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

                for(int nc = 0; nc < ncols; nc++){
                    output += String.format(",%d", event_counters[nc]);
                }


//                for(int nc = 0; nc < ncols-1; nc++){
//                    String num_val = String.format("%d", counters[i][nc])+",";
//                    output += String.format("%-"+string_length+"s", num_val);
//                }
//                String num_val = String.format("%d", counters[i][ncols-1]);
//                output += String.format("%-"+string_length+"s", num_val);

                bw.write(output);
                bw.newLine();
            }
            bw.close();

        }catch (IOException e){}

    }


    static void writeDataboxMicrohabPopsToFile(String directoryName, String filename, DataBox dataBox){

        File directory = new File(directoryName);
        if(!directory.exists()) directory.mkdirs();

        File file = new File(directoryName+"/"+filename+".txt");

        try{

            FileWriter fw = new FileWriter(file.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);

            ArrayList<Double> times = dataBox.getTimes();
            ArrayList<ArrayList<ArrayList<Double>>> mh_pops_over_time = dataBox.getAll_microhab_pops();


            int string_length = 12;

            for(int t = 0; t < mh_pops_over_time.size(); t++){

                bw.write("#t = "+String.format("%.3E", times.get(t)));
                bw.newLine();

                for(int mh = 0; mh < mh_pops_over_time.get(t).size(); mh++){

                    String output = "";
                    int nbac = mh_pops_over_time.get(t).get(mh).size();

                    for(int b = 0; b < nbac-1; b++){

                        String geno_val = String.format("%.5E", mh_pops_over_time.get(t).get(mh).get(b))+",";
                        output += String.format("%-"+string_length+"s", geno_val);
                    }
                    //need to handle if the array is of length 0 - think this does that
                    if(nbac > 0){
                        String geno_val = String.format("%.5E", mh_pops_over_time.get(t).get(mh).get(nbac-1));
                        output += String.format("%-"+string_length+"s", geno_val);
                    }

                    if(nbac==0){
                        //adds a negative value to avoid empty strings
                        String empty_val = String.format("%.5E", -9999.0);
                        output += String.format("%-"+string_length+"s", empty_val);
                    }

                    bw.write(output);
                    bw.newLine();
                }
            }

            bw.close();

        }catch (IOException e){}
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