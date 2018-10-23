package statalign;

import java.net.URL;
import java.net.URLClassLoader;
import java.util.Arrays;

import mpi.*;
import statalign.base.MainManager;
import statalign.base.Utils;

public class StatAlignParallel {

    public static void main(String[] args) {

    	// The MPJ framework appends (at least) three parameters in front of the normal 
    	// parameters when it is initializing the processes. We want to separate those
    	// from the normal parameters (mpiArgs = [0, 1, 'mpidev'], realArgs = [....]).    	
    	if (args.length < 3) {
            boolean parallel = true;
            CommandLine cl = new CommandLine(parallel);
            MainManager manager = new MainManager(null,parallel);
            System.out.println(cl.getUsageString(null, manager));
    	}
        String[] mpiArgs = new String[3];
        String[] realArguments = new String[args.length - 3];
        for (int i = 0; i < 3; i++) {
            mpiArgs[i] = args[i];
        }
        for (int i = 3; i < args.length; i++) {
            realArguments[i - 3] = args[i];
        }
        
        // Initializes the MPJ framework.
        MPI.Init(mpiArgs);
       
        int rank = MPI.COMM_WORLD.Rank();
        int noOfProcesses = MPI.COMM_WORLD.Size();
        
        if (rank==0) {
            System.out.println("\n    StatAlign " + StatAlign.version + " - parallel version\n");
            System.out.println("    Arguments: " + Arrays.toString(realArguments));
            System.out.println("    No of processes: " + noOfProcesses + "");
        }
        if (rank==1) {
        	//Utils.DEBUG = true;
        }
       
        // in case we're running on a headless server
        System.setProperty("java.awt.headless", "true");

        /*
         * Set java.class.path so that we can locate plugins.
         * This is necessary when StatAlign is invoked via MPJ's starter.jar.
         */
//        System.out.println(System.getProperty("java.class.path"));
//        URL[] urls = ((URLClassLoader) (Thread.currentThread().getContextClassLoader())).getURLs();
//        for (URL url : urls) {
//        	System.out.println(url.toString());
//        }
        String statalign_bin = ((URLClassLoader) (Thread.currentThread().getContextClassLoader())).getURLs()[1].toString();
        statalign_bin = statalign_bin.replaceFirst("file:",""); 
        System.setProperty("java.class.path",statalign_bin);

        boolean parallel = true;
        MainManager manager = new MainManager(null,parallel);
        CommandLine cl = new CommandLine(parallel);
        
        // Only get INFO messages (errors and warnings excluded) from the master.
        if (rank==0) {
        	cl.setVerbose(true);
        }
        
        // Configure the program.
        if (cl.fillParams(realArguments, manager) > 0) {
            System.exit(1);
        }
        
        // TODO cf. comment in StatAlign.java
		manager.init(cl.pluginParameters);
		//manager.init(Postprocess.pluginParameters);

        // Sets up a barrier to synchronize the processes here.
        MPI.COMM_WORLD.Barrier();
        
        manager.start(noOfProcesses, rank);

        MPI.Finalize();
    }

}
