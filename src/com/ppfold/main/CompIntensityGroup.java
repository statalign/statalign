package com.ppfold.main;

import java.awt.Dimension;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

public class CompIntensityGroup extends JPanel {

	private static final long serialVersionUID = -6531182945482227985L;
	JSpinner compIntensitySpinner; 
	JSpinner iterationSpinner; 
	
	
	public CompIntensityGroup(){
		
        SpinnerNumberModel numberSpinnerModel = new SpinnerNumberModel(Runtime.getRuntime().availableProcessors(), 1, Runtime.getRuntime().availableProcessors(), 1); 
		compIntensitySpinner = new JSpinner(numberSpinnerModel); 
		compIntensitySpinner.setPreferredSize(new Dimension(50,20));
		add(new JLabel("Computational intensity (cores):"));
		compIntensitySpinner.setToolTipText("Number of processor cores to use for the calculation. Higher numbers will give faster results but your computer will be less responsive meanwhile.");
		add(compIntensitySpinner);
		
		compIntensitySpinner.addChangeListener(new ChangeListener(){
			public void stateChanged(ChangeEvent arg0) {
				PPfoldMain.nrprocessors = (Integer) compIntensitySpinner.getValue();
			}
		});
		
        SpinnerNumberModel numberSpinnerModel2 = new SpinnerNumberModel(PPfoldMain.iterlimit, 1, 200, 1); 
		iterationSpinner = new JSpinner(numberSpinnerModel2); 
		iterationSpinner.setPreferredSize(new Dimension(50,20));
		iterationSpinner.setToolTipText("The maximum number of times branch lengths are iterated for optimization");
		add(new JLabel("Tree optimization iteration limit:"));
		add(iterationSpinner);
		
		iterationSpinner.addChangeListener(new ChangeListener(){
			public void stateChanged(ChangeEvent arg0) {
				PPfoldMain.iterlimit = (Integer) iterationSpinner.getValue();
			}
		});
		
	}
	
	@Override
	public void setEnabled(boolean value){
		compIntensitySpinner.setEnabled(value);
		iterationSpinner.setEnabled(value);
	}
	
}
