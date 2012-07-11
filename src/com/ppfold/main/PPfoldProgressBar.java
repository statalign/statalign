package com.ppfold.main;

import java.util.ArrayList;

import javax.swing.JProgressBar;

import com.ppfold.algo.Progress;

public class PPfoldProgressBar implements Progress {

	private JProgressBar barGUI;
	private PPfoldProgressBar parent;
	private double contribution; 
	private double progress; 
	private String message; 
	
	public PPfoldProgressBar(JProgressBar bar, PPfoldProgressBar parin, double contrib){
		barGUI = bar;
		parent = parin;
		contribution = contrib;
	}

	public double getProgress() {
		return progress;
	}
	
	public void setProgress(double newprogress) {
		double oldprogress = this.progress;
		this.progress = newprogress;
		if(barGUI!=null){ //if it has a bar then it has no parent
			barGUI.setValue((int)(progress*100));
			barGUI.setString(barGUI.getValue() + "% (" + message + ")");
		}
		else{ //in this case it has a parent
			this.parent.childReportActivity(oldprogress, newprogress, contribution);
		}
	}

	private void childReportActivity(double oldprogress, double newprogress, double contribution){
		this.setProgress(this.progress - oldprogress*contribution + newprogress*contribution);
	}
	
	public Progress getChildProgress(double fraction) {
		PPfoldProgressBar child = new PPfoldProgressBar(null, this, fraction);
		return child;
	}

	public boolean isSuspended() {
		//Not suspending
		return false;
	}

	public void setSuspended(boolean isSuspended) {
		//Not suspending
	}

	public boolean shouldStop() {
		// Returns true if folding should stop 
		return PPfoldMain.shouldstop;
	}

	public void checkStop() throws InterruptedException {
		if(shouldStop()){
			throw new InterruptedException("Process has been interrupted!");
		}
	}

	public void setCurrentActivity(String activity) {
		message = activity;
		if(barGUI!=null){
			barGUI.setString(barGUI.getValue() + "% (" + message + ")");
		}
		else{
			parent.setCurrentActivity(activity);
		}

	}

	public String getCurrentActivity() {
		return message;
	}
	
}
