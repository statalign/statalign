package com.ppfold.main;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

public class SequenceExportsGroup extends JPanel {

	private static final long serialVersionUID = -6196506765119020438L;
	private JCheckBox sequenceExportBox;
	private JTextField sequenceNameField;
	
	public SequenceExportsGroup(){
		this.setLayout(new BoxLayout(this,BoxLayout.LINE_AXIS));
		
		sequenceExportBox = new JCheckBox("Also export data for this sequence: ");
		sequenceExportBox.setToolTipText("If a sequence is of special interest in the alignment, this will spare you the work of converting" +
				" alignment coordinates to sequence coordinates.");
		
		
		sequenceNameField = new JTextField("");
		sequenceNameField.setEnabled(false);
		sequenceNameField.setToolTipText("If a sequence is of special interest in the alignment, this will spare you the work of converting" +
				" alignment coordinates to sequence coordinates.");
		sequenceExportBox.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent arg0) {
				if(sequenceExportBox.isSelected()){
					sequenceNameField.setEnabled(true);
				}
				else{
					sequenceNameField.setEnabled(false);
				}
				updateFileHandle();
			}
			
		});
		
		sequenceNameField.getDocument().addDocumentListener(new DocumentListener(){
			public void changedUpdate(DocumentEvent arg0) {updateFileHandle();}
			public void insertUpdate(DocumentEvent arg0) {updateFileHandle();}
			public void removeUpdate(DocumentEvent arg0) {updateFileHandle();}
		});
		
		this.add(sequenceExportBox);
		this.add(sequenceNameField);
		
	}


	@Override 
	public void setEnabled(boolean value){
		sequenceExportBox.setEnabled(value);
		if(value==true && sequenceExportBox.isSelected()){
			sequenceNameField.setEnabled(true);
		}
		else{
			sequenceNameField.setEnabled(false);
		}
	}
	
	public void updateFileHandle(){
		if(sequenceExportBox.isSelected() && !sequenceNameField.getText().trim().isEmpty()){
			PPfoldMain.seqexportname = sequenceNameField.getText().trim();
		}
		else{		
			PPfoldMain.seqexportname = null;
		}
	}

}
