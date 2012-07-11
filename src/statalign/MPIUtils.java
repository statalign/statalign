package statalign;

public class MPIUtils {
	
	public static void println(int rank, String message) {
		String rankString = (rank == 0 ? "Master" : "Worker " + rank);
		System.out.printf("%8s: %s\n", rankString, message);
	}
	
	/**
	 * Determines if a process with a certain <tt>rank</tt> is the master. By default only processes 
	 * with the rank of 0 are masters. 
	 * 
	 * @param rank the rank of the process
	 * @return whether the process is a master.
	 */
	public static boolean isMaster(int rank) {
		return rank == 0;
	}
	
}
