package com.ppfold.main;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.border.Border;
import javax.swing.text.DefaultCaret;

public class PPfoldGUIMainWindow extends JFrame {

	MainButtonPanel exitStartButtons;
	PPfoldButtonGroup alignmentGroup;
	PPfoldButtonGroup treeGroup;
	PPfoldButtonGroup outputGroup;
	DataButtonGroup databuttons; 
	ExportsGroup exportbuttons; 
	SequenceExportsGroup seqexportbuttons; 
	JTextArea textArea;
	CompIntensityGroup intensityArea;
	JProgressBar progressbar; 
	JButton suspendbutton;
	JButton stopbutton;
	
	//public static ArrayList<String> dataIDs = new ArrayList<String>(); 
	
	public static PrintStream ps_orig; //output stream
	public static PrintStream ps_orig_err; //Error stream
	
	
	public static String directory = null;
	
	private static final long serialVersionUID = 1L;

	public PPfoldGUIMainWindow(){
		createAndShowGUI();
	}
	
	private void createAndShowGUI(){
		setTitle("PPfold version " + PPfoldMain.versionnumber);
		setSize(400,400);
		setLocation(200,10);
		setResizable( false );
		
		Container contentPane = this.getContentPane();
		contentPane.setLayout(new BoxLayout(contentPane, BoxLayout.PAGE_AXIS));
		
		JPanel mainpanel = new JPanel();
		Border border = BorderFactory.createEmptyBorder(10,10,10,10);
		mainpanel.setBorder(border);
		
		outputGroup = new PPfoldButtonGroup("Output folder: ", PPfoldButtonGroup.FOLDERCHOOSER){
			private static final long serialVersionUID = 1L;
			public void updateModel(){
				if(outputGroup.getFileName().startsWith("<No file selected>")){
					PPfoldMain.outputdir = null;}
				else{
					PPfoldMain.outputdir = outputGroup.getFileName();
				}
			}
		};
		
		alignmentGroup = new PPfoldButtonGroup("Alignment file: ", PPfoldButtonGroup.FILECHOOSER){
			private static final long serialVersionUID = 1L;
			public void updateModel(){
				if(alignmentGroup.getFileName().startsWith("<No file selected>")){
					PPfoldMain.alignmentfilename = null;}
				else{
					PPfoldMain.alignmentfilename = alignmentGroup.getFileName();
					if(PPfoldMain.outputdir==null){
						PPfoldMain.outputdir = directory;
						outputGroup.setText(directory);
					}
					File alignfile = new File (PPfoldMain.alignmentfilename);
					String exportfilehandle = alignfile.getName();
		        	int lastdotpos = exportfilehandle.lastIndexOf('.');
		        	exportfilehandle = exportfilehandle.substring(0,lastdotpos);
					exportbuttons.setExportPrefix(exportfilehandle);
				}				
			}
		};
		alignmentGroup.setAlignmentX(Component.CENTER_ALIGNMENT);
		
		treeGroup = new PPfoldButtonGroup("Tree file: ", PPfoldButtonGroup.FILECHOOSER){
			private static final long serialVersionUID = 1L;
			public void updateModel(){
				if(treeGroup.getFileName().startsWith("<No file selected>")){
					PPfoldMain.treefilename = null;
					PPfoldMain.createtree=true;
					PPfoldMain.optimizetree=true;
				}
				else{
					PPfoldMain.treefilename = treeGroup.getFileName();
					PPfoldMain.createtree=false;
					PPfoldMain.optimizetree=false;
				}
			}
		};
		treeGroup.setAlignmentX(Component.CENTER_ALIGNMENT);
	
		mainpanel.setLayout(new BoxLayout(mainpanel, BoxLayout.PAGE_AXIS)); 
		mainpanel.add(alignmentGroup);
		mainpanel.add(Box.createRigidArea(new Dimension(0, 5)));
		mainpanel.add(treeGroup);
		mainpanel.add(Box.createRigidArea(new Dimension(0, 5)));
		mainpanel.add(outputGroup);		
		mainpanel.add(Box.createRigidArea(new Dimension(0, 5)));

		contentPane.add( mainpanel); 

		
		addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e){
					PPfoldMain.userfinished = true; 
					System.exit(0);
			} // Terminate the program
		});
	
		databuttons = new DataButtonGroup();
		databuttons.setAlignmentX(Component.CENTER_ALIGNMENT);
		mainpanel.add(databuttons);
		
		exportbuttons = new ExportsGroup();
		mainpanel.add(exportbuttons);
		mainpanel.add(Box.createRigidArea(new Dimension(0,5)));
		seqexportbuttons = new SequenceExportsGroup();
		mainpanel.add(seqexportbuttons);
		
		intensityArea = new CompIntensityGroup();
		mainpanel.add(intensityArea);
		
		exitStartButtons = new MainButtonPanel(this);
		exitStartButtons.setBorder(border);
		contentPane.add( exitStartButtons); 
		
		JPanel outputpanel = new JPanel();
		outputpanel.setBorder(border);
		
		textArea = new AutoScrollingTextArea();
		float[] somecolor = new float[3];
		somecolor = Color.RGBtoHSB(175,238,238, null);
		textArea.setBackground(Color.getHSBColor(somecolor[0], somecolor[1], somecolor[2]));
		textArea.setLineWrap(true);
		textArea.setWrapStyleWord(true);
		textArea.setToolTipText("Stdout and stderr are redirected here.");
		System.out.println("Redirecting output and errors to GUI window.");
		
		ps_orig = System.out;//This is for backup
		ps_orig_err = System.err;//This is for backup

		try{
			System.setOut(new PrintStream(new TextAreaOutputStream(textArea)));
			//System.setErr(new PrintStream(new TextAreaOutputStream(textArea)));
		}
		catch(SecurityException e){
			System.err.println("Permission denied, output and errors remain on command-line.");;
		}

		outputpanel.setLayout(new BoxLayout(outputpanel, BoxLayout.PAGE_AXIS));


		JLabel outputlabel = new JLabel("PPfold messages:");
		outputlabel.setAlignmentX(Component.CENTER_ALIGNMENT);
		outputpanel.add(outputlabel); 
		outputpanel.add(textArea);
		
		JScrollPane scrollpane = new JScrollPane(textArea);
		scrollpane.setVerticalScrollBarPolicy(
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
		scrollpane.setPreferredSize(new Dimension(400, 200));
		
		outputpanel.add(scrollpane);		
		
		JButton copyButton = new JButton("Copy to clipboard");
		 copyButton.addActionListener(new ActionListener() {
	            public void actionPerformed(ActionEvent e) {
	        		String selection = textArea.getText();
	       		    StringSelection data = new StringSelection(selection);
	       		   Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
	       		   clipboard.setContents(data, data);
	            }
	        }); 
		 
		copyButton.setAlignmentX(Component.CENTER_ALIGNMENT);
		outputpanel.add(copyButton);
		
		contentPane.add(outputpanel);

		JPanel progresspanel = new JPanel();
		progresspanel.setLayout(new BoxLayout(progresspanel, BoxLayout.LINE_AXIS));
		progresspanel.setBorder(border);
		
		progressbar = new JProgressBar(0, 100);
		//JLabel progresslabel = new JLabel("Folding progress: ");
		//progresslabel.setAlignmentX(Component.CENTER_ALIGNMENT);
		progressbar.setAlignmentX(Component.CENTER_ALIGNMENT);
		progressbar.setStringPainted(true);
		progressbar.setString("(Folding not started)");
		progressbar.setToolTipText("Folding progress is shown here");
		//progresspanel.add(progresslabel);
		progresspanel.add(progressbar);
		stopbutton = new JButton("Stop");
		stopbutton.setEnabled(false);
		stopbutton.setAlignmentX(Component.CENTER_ALIGNMENT);
		stopbutton.setToolTipText("Stops the current folding.");
		
		progresspanel.add(Box.createRigidArea(new Dimension(5, 0)));
		progresspanel.add(stopbutton);
		stopbutton.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent arg0) {
				PPfoldMain.shouldstop = true;
				stopbutton.setEnabled(false); //Isn't allowed to press twice!
			}
		});
		
		contentPane.add(progresspanel);	
		
		mainpanel.setVisible(true);
		outputpanel.setVisible(true);
		contentPane.setVisible(true); 

		pack();
		setVisible(true);
		
		}
	
	public void checkAllData(boolean messageDialog){
		try{
			String testresults = CheckAllData.checkData();
			String message = "The results of the test are given below. If any of the tests failed,\n" +
					"PPfold is unlikely to run correctly.";
			if(messageDialog){
				@SuppressWarnings("unused")
				JDialog testResults = new CustomTextDialog("Test results",message,testresults);
			}
		}
		catch(Exception e){
			String message = "An error occured while checking the data! This definitely should not have happened.\n" +
					"Please email the following debug text to zs@mb.au.dk and describe what you did.";
        	String text = new String("");
        	for(StackTraceElement s:e.getStackTrace()){
        		text = text.concat(s.toString()+"\n");
        	}
        	String errormessage = "PPfold version " + PPfoldMain.versionnumber + "\n";
    		try{
    			errormessage += "Platform: " + System.getProperty("os.name") + ", version " + System.getProperty("os.version") + "\n";
    			errormessage += "JVM: " + System.getProperty("java.vm.vendor") + ", JRE version "+ System.getProperty("java.version") + "\n";
    		}
    		catch(SecurityException e1){
    			errormessage += "System check permission denied \n";    			
    		}
    		errormessage += "\n";
    		errormessage += e.toString()+":\n"+text;
    		@SuppressWarnings("unused")
			JDialog failure = new CustomTextDialog("Error",message,errormessage);
    		e.printStackTrace();
		}
		
	}
	
	public void resetOutput(){
		textArea.setText("");
	}
	
	public void enableAll(){
		alignmentGroup.setEnabled(true);
		treeGroup.setEnabled(true);
		outputGroup.setEnabled(true);
		databuttons.setEnabled(true);
		exitStartButtons.setEnabled(true); //exit will always be enabled
		exportbuttons.setEnabled(true);
		seqexportbuttons.setEnabled(true);
		intensityArea.setEnabled(true);
		stopbutton.setEnabled(false);

	}
	
	public void disableAll(){
		alignmentGroup.setEnabled(false);
		treeGroup.setEnabled(false);
		outputGroup.setEnabled(false);
		databuttons.setEnabled(false);
		exitStartButtons.setEnabled(false); //exit will always be enabled 
		exportbuttons.setEnabled(false);
		seqexportbuttons.setEnabled(false);
		intensityArea.setEnabled(false);
		stopbutton.setEnabled(true);
	}

}
