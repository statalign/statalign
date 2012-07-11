package com.ppfold.main;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.PrintStream;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import com.ppfold.algo.Progress;

public class MainButtonPanel extends JPanel {
	  // members:
	  private JButton startButton;
	  private JButton exitButton;
	  private JButton checkInputButton; 
	  private PPfoldGUIMainWindow parent;
	  // constructors:
	  public MainButtonPanel(final PPfoldGUIMainWindow inparent) {
		this.parent = inparent;
		this.setLayout(new GridLayout());

	    // create buttons
		checkInputButton = new JButton("Check input");
	    startButton = new JButton("Start");
	    exitButton = new JButton("Exit");

	    // add buttons to current panel
	    add(checkInputButton,BorderLayout.WEST);
	    add(startButton,BorderLayout.WEST);  // add button to current panel
	    add(exitButton,BorderLayout.WEST); // add button to current panel
	    
	    // register the current panel as listener for the buttons
	    startButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
            		if(PPfoldMain.alignmentfilename==null){
            			JOptionPane.showMessageDialog(parent,
            				    "Alignment must be given!",
            				    "Missing input",
            				    JOptionPane.ERROR_MESSAGE);
            		}
            		else if(PPfoldMain.outputdir==null||PPfoldMain.outputdir==""||PPfoldMain.outputdir.startsWith("<No file selected>")){
            			JOptionPane.showMessageDialog(parent,
            				    "Not a valid directory!",
            				    "Missing input",
            				    JOptionPane.ERROR_MESSAGE);
            		}
            		else{            			
            			PPfoldMain.userfinished = true;
            			parent.disableAll();
            			parent.checkAllData(false);
            			if(!CheckAllData.nothingfailed){
            				Object[] options = {"Yes",
                            "No, show me the problems"};
            				int n = JOptionPane.showOptionDialog(parent,
            						"Potential problems were detected in your input, which are likely to \n" +
            						"cause errors or incorrect prediction results. \n" +
            						"Do you want to run PPfold anyway?",
            						"Potential problems",
            						JOptionPane.YES_NO_OPTION,
            						JOptionPane.QUESTION_MESSAGE,
            						null,     //do not use a custom Icon
            						options,  //the titles of buttons
            						options[1]); //default button title
            				if(n==JOptionPane.NO_OPTION){
            					parent.enableAll();
            					parent.checkAllData(true);
            					return;
            				}
            			}
            			//we fold in on thread
            			parent.resetOutput();
            			final PPfoldProgressBar activity = new PPfoldProgressBar(parent.progressbar, null, 1.0);
            			PPfoldMain.setProgressBar(activity);
            			//TODO: instead of clearing extradata only, set all parameters after what is in the GUI
            			//Just to double-check that everything is being set correctly!!!
            			final PPfoldMain thread = new PPfoldMain();
                    	Thread foldingthread = new Thread(thread);
                    	foldingthread.start();
                    	//and listen for finishing in another thread
                    	//so the GUI doesn't block.
                    	final ExecutorService listener;
                    	listener = Executors.newSingleThreadExecutor();
                    	listener.execute(new Runnable(){
							public void run() {
					        	while(!thread.foldingfinished && !PPfoldMain.isstopping){
				        		//Wait until folding is finished
					        		try {
					        			Thread.sleep(100);
									} catch (InterruptedException e) {
										System.out.println("User has quit");
									}
					        	}
					        	parent.enableAll();
					        	if(PPfoldMain.shouldstop){					        		
					        		PPfoldMain.shouldstop = false;
					        		PPfoldMain.isstopping = false;
					        		activity.setCurrentActivity("Folding stopped by user.");
					        	}
					        	else if(thread.success()){
					        		JOptionPane.showMessageDialog(parent, "PPfold successfully finished!",
					        				"Success",
					        			    JOptionPane.INFORMATION_MESSAGE);
					        	}
					        	else{
					        		String message = "An error occured! Try checking your data with the 'Check input' button\n"+
					        				         "and follow any instructions given below. If the problem persists or you\n" +
					        				         "have a question, please email the following text to zs@mb.au.dk and \n" +
					        				         "describe what you did.";
					            	String errormessage = "PPfold version " + PPfoldMain.versionnumber + "\n";
					        		try{
					        			errormessage += "Platform: " + System.getProperty("os.name") + ", version " + System.getProperty("os.version") + "\n";
					        			errormessage += "JVM: " + System.getProperty("java.vm.vendor") + ", JRE version "+ System.getProperty("java.version") + "\n";
					        		}
					        		catch(SecurityException e1){
					        			errormessage += "System check permission denied \n";    			
					        		}
					        		errormessage += "\n";
					        		errormessage += thread.errormessage;
									new CustomTextDialog("Error",message,errormessage);
					        	}
					        	listener.shutdown();
					        	thread.cleanUp(); //To remove lingering data to avoid memory leaks.
							}
                    	});
            		}
            }
        });  
	    exitButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
        		System.setOut(PPfoldGUIMainWindow.ps_orig);
        		System.setErr(PPfoldGUIMainWindow.ps_orig_err);
            	System.out.println("Quitting...");
                System.exit(0);
            }
        });
	    
	    checkInputButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
            	parent.checkAllData(true);
            }
        }); 
	    
	    checkInputButton.setToolTipText("Checks input data before you run PPfold");

	  }
	  
	  @Override 
	  public void setEnabled(boolean value){
		  checkInputButton.setEnabled(value);
		  startButton.setEnabled(value);
	  }
	  
	} 


