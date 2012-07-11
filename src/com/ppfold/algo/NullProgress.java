package com.ppfold.algo;

/**
 * Class is used for standalone version, where no actual progress bar is
 * implemented.
 * 
 * @author M.Vaerum
 */

public class NullProgress implements Progress {
	public static final Progress INSTANCE = new NullProgress();

	private NullProgress() {

	}

	public double getProgress() {
		return 0.0;
	}

	public void setProgress(double progress) {
	}

	public Progress getChildProgress(double fraction) {
		return this;
	}

	public boolean isSuspended() {
		return false;
	}

	public void setSuspended(boolean isSuspended) {

	}

	public boolean shouldStop() {
		return false;
	}

	public void checkStop() throws InterruptedException {

	}

	public void setCurrentActivity(String activity) {

	}

	public String getCurrentActivity() {
		return null;
	}

}
