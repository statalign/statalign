package com.ppfold.algo;

/**
 * Abstract class for asynchronous job execution.
 * 
 * @author M.Vaerum
 */

public abstract class AsynchronousJobExecutor {
	// This is the interface
	public abstract void startExecution(CYKJob cYKJob, JobListener listener);

	public abstract void startExecution(PhyloJob job, JobListener listener);
	
	public abstract void startExecution(PhyloJobFuzzy job, JobListener listener);

	public abstract String getDescription();

	public abstract String getId();
	
	//This will tell the Executor that all threads should stop 
	public abstract void shutDown();
	
	//Returns true when all the jobs are terminated.
	public abstract boolean isTerminated();
	
	
	// Use this to get the actual instance based on a param instance
	// public abstract AsynchronousJobExecutor createJobExecutor(AlgoParameters
	// param);

	@Override
	public String toString() {
		return getDescription();
	}

	private static AsynchronousJobExecutor localExecutor = null;
	private static AsynchronousJobExecutor remoteExecutor = null;
	private static int nrcores = Runtime.getRuntime().availableProcessors();

	// Static methods

	static {
		AsynchronousJobExecutorThreadPool executor = new AsynchronousJobExecutorThreadPool(
				Runtime.getRuntime().availableProcessors());
		setLocalExecutor(executor);
	}

	public static void setLocalExecutor(AsynchronousJobExecutor local) {
		localExecutor = local;
	}

	public static void setRemoteExecutor(AsynchronousJobExecutor remote) {
		remoteExecutor = remote;
	}

	
	public static AsynchronousJobExecutor getLocalExecutor() {
		return localExecutor;
	}

	//ZS method to enable user selection of nr of cores to use
	public static AsynchronousJobExecutor getLocalExecutor(int inputnr) {
		nrcores = inputnr; 
		AsynchronousJobExecutorThreadPool executor = new AsynchronousJobExecutorThreadPool(nrcores);
		setLocalExecutor(executor);
		return localExecutor;
	}

	
	public static AsynchronousJobExecutor getRemoteExecutor() {
		return remoteExecutor;
	}
	
	
	

}
