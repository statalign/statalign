package statalign;

import java.util.Arrays;

import mpi.MPI;
import statalign.base.MainManager;
import statalign.postprocess.Postprocess;

public class StatAlignParallel {

    public static void main(String[] args) {

    	// The MPJ framework appends (at least) three parameters in front of the normal 
    	// parameters when it is initializing the processes. We want to separate those
    	// from the normal parameters (mpiArgs = [0, 1, 'mpidev'], realArgs = [....]).
    	
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
        
        if (MPIUtils.isMaster(rank)) {
            System.out.println("\n    StatAlign " + StatAlign.version + " - parallel version\n");
            System.out.println("    Arguments: " + Arrays.toString(realArguments));
            System.out.println("    No of processes: " + noOfProcesses + "");
        }

        MainManager manager = new MainManager(null);
        CommandLine cl = new CommandLine(true);
        
        // Only get INFO messages (errors and warnings excluded) from the master.
        if (MPIUtils.isMaster(rank)) {
        	cl.setVerbose(true);
        }
        
        // Configure the program.
        if (cl.fillParams(realArguments, manager) > 0) {
            System.exit(1);
        }
        
        // TODO cf. comment in StatAlign
		manager.init(cl.pluginParameters);
		//manager.init(Postprocess.pluginParameters);

        // Sets up a barrier to synchronize the processes here.
        MPI.COMM_WORLD.Barrier();
        
        manager.startParallel(noOfProcesses, rank);

        MPI.Finalize();
    }

}
