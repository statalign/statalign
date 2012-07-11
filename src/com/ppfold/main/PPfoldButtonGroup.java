package com.ppfold.main;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.border.Border;

public abstract class PPfoldButtonGroup extends JPanel implements Updatable {
	  // members:
	  private JLabel title; 
	  private JButton dialogButton;
	  private JButton resetButton;
	  private JTextField dataname;
	  
	  public static int FILECHOOSER = 0;
	  public static int FOLDERCHOOSER = 1;
	  
	  // constructors:
	  public PPfoldButtonGroup(final String titleString, final int choosertype) {

		this.setLayout(new BoxLayout(this, BoxLayout.LINE_AXIS));
		this.setSize(100,20);
		this.setPreferredSize(new Dimension(100,20));
		
	    // create components
		title = new JLabel(titleString);
		dataname = new JTextField("<No file selected>");
		dataname.setEditable(false);
	    dialogButton = new JButton("Browse...");
	    dialogButton.setToolTipText("Browse");
	    title.setVisible(true);
	    dataname.setVisible(true);
	    //dialogButton.setPreferredSize(new Dimension(90,30));
	    
	    //dataname.setPreferredSize(new Dimension(230,5));
	   // title.setPreferredSize(new Dimension(130,10));
	    
	    dialogButton.setVisible(true);
	    
	    resetButton = new JButton("Remove");
	    resetButton.setToolTipText("Remove");
	    //resetButton.setPreferredSize(new Dimension(80,30));
	    resetButton.setVisible(true); 
	    resetButton.setEnabled(false);
	    
	    // add buttons to current panel
	    dataname.setAlignmentX(Component.LEFT_ALIGNMENT);
	    title.setAlignmentX(Component.LEFT_ALIGNMENT);
	    title.setHorizontalAlignment(JLabel.CENTER);
	    dialogButton.setAlignmentX(Component.RIGHT_ALIGNMENT);
	    resetButton.setAlignmentX(Component.RIGHT_ALIGNMENT);
	    dataname.setAlignmentY(Component.CENTER_ALIGNMENT);
	    title.setAlignmentY(Component.CENTER_ALIGNMENT);
	    dialogButton.setAlignmentY(Component.CENTER_ALIGNMENT);
	    resetButton.setAlignmentY(Component.CENTER_ALIGNMENT);
	    
	    add(title);
	    add(dataname);
		this.add(Box.createRigidArea(new Dimension(5, 0)));
	    add(dialogButton);  // add button to current panel
		this.add(Box.createRigidArea(new Dimension(5, 0)));
	    add(resetButton);	    
	    
	    // register the current panel as listener for the buttons
	    dialogButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
            			final JFileChooser fc = new JFileChooser();
            			if(PPfoldGUIMainWindow.directory!=null){
            				File theDirectory = new File(PPfoldGUIMainWindow.directory);
            			fc.setCurrentDirectory(theDirectory);
            			}
            			if(choosertype == FOLDERCHOOSER){
            				fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
            			}
            			int returnVal = fc.showOpenDialog(null);
            	        if (returnVal == JFileChooser.APPROVE_OPTION) {
            	            File file = fc.getSelectedFile();
            	            dataname.setText(file.getAbsolutePath());
                    		resetButton.setEnabled(true);
                    		PPfoldGUIMainWindow.directory = file.getParent(); ;
                    		//System.out.println(PPfoldGUIMainWindow.directory);
                    		updateModel();
            	        }  //otherwise user cancelled, do nothing }     			
            		}
        	});
	    
	    resetButton.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent arg0) {
				dataname.setText("<No file selected>");
				resetButton.setEnabled(false);
				updateModel();
	    }
	  });
	    
	  } 
	  public String getFileName(){
		  return dataname.getText();
	  }
	  public void setText(String text){
		  if(!text.startsWith("<No file selected>")){
			  resetButton.setEnabled(true);
		  }
		  dataname.setText(text);
	  }
	  
	  @Override 
	  public void setEnabled(boolean val){
		  dialogButton.setEnabled(val);
		  resetButton.setEnabled(val);
		  dataname.setEnabled(val);
		  //Double-check that Reset is not enabled when not wanted
		  if(dataname.getText().startsWith("<No file selected>")){
			resetButton.setEnabled(false);}
	  }
	  
	  
	} 


