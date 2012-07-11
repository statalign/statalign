package com.ppfold.algo;

/**
 * Interface for a progress bar
 * 
 * @author M.Vaerum
 */

public interface Progress {

	/**
	 * Get the current progress (from 0.0 to 1.0).
	 * 
	 * @return the current progress
	 */
	public abstract double getProgress();

	/**
	 * Set the progress. Use Double.NaN to indicate that the progress is
	 * unknown.
	 * 
	 * @param progress
	 *            the progress (between 0.0 and 1.0)
	 */
	public abstract void setProgress(double progress);

	/**
	 * Get a suitable Activity for a child process.
	 * 
	 * @param fraction
	 *            the fraction of the overall progress that the child process
	 *            represents
	 * @return a child Activity
	 */
	public abstract Progress getChildProgress(double fraction);

	/**
	 * Check whether the process is suspended.
	 * 
	 * @return Whether the process is suspended.
	 */
	public abstract boolean isSuspended();

	/**
	 * Set whether the process should suspend.
	 * 
	 * @param isSuspended
	 *            whether the process should suspend
	 */
	public abstract void setSuspended(boolean isSuspended);

	/**
	 * Check whether we should stop the ongoing process. The thread will sleep
	 * if the process is suspended.
	 * 
	 * @return whether we should stop
	 */
	public abstract boolean shouldStop();

	/**
	 * Check whether we should stop the ongoing process. The thread will sleep
	 * if the process is suspended.
	 * 
	 * @throws InterruptedException
	 *             if the process should be stopped.
	 */
	public abstract void checkStop() throws InterruptedException;

	/**
	 * Sets a short description of what is going on right now.
	 * 
	 * @param activity
	 */
	public abstract void setCurrentActivity(String activity);

	/**
	 * Returns a short description of what is going on right now.
	 */
	public abstract String getCurrentActivity();

}
