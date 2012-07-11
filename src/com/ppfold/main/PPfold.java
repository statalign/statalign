package com.ppfold.main;

import java.awt.HeadlessException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.swing.JOptionPane;

public class PPfold {
	public static void main(String[] args) {
		System.out.println("**************");
        System.out.println("Running PPfold version " + PPfoldMain.versionnumber);       
        System.out.println("Checking input arguments...");

        PPfoldGUIMainWindow mainframe = null;
		try {
			mainframe = ArgumentParser.parseArgs(args);
		} 
		catch(HeadlessException e){
    		System.out.println("GUI could not be initialized. Error stack trace: ");
    		e.printStackTrace();
    		System.exit(-1);
		}
		catch (Exception e1) {
			System.err.println("There was a problem opening some of the input files.");
			e1.printStackTrace();
			if(mainframe==null){
				System.exit(-1);
			}
		}        
                
        if(mainframe!=null){
        	PPfoldMain.gui=true;
            System.out.println("Welcome to PPfold version " + PPfoldMain.versionnumber + "!");   
            System.out.println("**********************************");
        	System.out.println("Please select an alignment and any additional data.");
        	System.out.println("Optional: You can check your data before running PPfold using the 'Check input' button.");

        	
        	while(true){
        		//Wait until user has finished choosing inputs
        		try {
					Thread.sleep(100);
				} catch (InterruptedException e) {
					// This happens when the program quits because the user pressed quit in the GUI
				}
        	}
        }
        else{
        	try{
        		final PPfoldMain thread = new PPfoldMain();
        		new Thread(thread).start();
        		final ExecutorService listener;
            	listener = Executors.newSingleThreadExecutor();
            	listener.execute(new Runnable(){
					public void run() {
			        	while(!thread.foldingfinished && !PPfoldMain.isstopping){
		        		//Wait until folding is finished or the user stopped it
			        		try {
			        			Thread.sleep(100);
							} catch (InterruptedException e) {
								System.out.println("User has quit");
							}
			        	}
			        	if(PPfoldMain.shouldstop){					        		
			        		PPfoldMain.shouldstop = false;
			        		PPfoldMain.isstopping = false;
			        	}
			        	else if(thread.success()){
			        		System.out.println("PPfold finished successfully.");
			        	}
			        	else{
			        		String message = "An error occured! Information is given below. If the problem persists or you\n" +
			        				         "have a question, please email the following text to " +
			        				         "zs@mb.au.dk and describe \nwhat you did.";
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
							System.out.println(message);
							System.out.println(errormessage);
							System.exit(-1);
			        	}
			        	listener.shutdown();
					}
            	});
        	}
        	catch(Exception e){
        		//Catch all exceptions not caught elsewhere
    			System.err.println("There was a problem while running PPfold.");
        		e.printStackTrace();
        	}
        }
        
	}
}
