package com.ppfold.main;

import java.awt.BorderLayout;
import java.awt.FileDialog;
import java.awt.Frame;
import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

public class ArgumentParser {
	
		
	public static PPfoldGUIMainWindow parseArgs(String[] args) throws Exception{
		
		if(args.length == 0){
        	System.out.println("None found. Attempting to trigger GUI...");
        	try{
        		PPfoldGUIMainWindow PPfoldGUI = new PPfoldGUIMainWindow();
        		PPfoldGUI.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        		return PPfoldGUI;
        	}
        	catch(Exception e){
        		System.out.println("GUI could not be initialized. Please run in non-graphical mode!");
        		e.printStackTrace();
        		System.exit(-1);
        	}
        }
        else if(args.length == 1 && (args[0].startsWith("--help") || args[0].startsWith("-h"))){
        	quitWithUsageMessage();
        }
        else if(args.length >= 1){
        	PPfoldMain.alignmentfilename = args[0];
        	
        	for(int i = 1; i<args.length; i++){
        		String argument = args[i];
        		//parse the rest of the arguments
        		if(argument.startsWith("--usetree") || argument.startsWith("-t")){
        			if(args.length < i+2){
        				System.out.println("Not enough arguments: tree file missing");
        				quitWithUsageMessage();
        			}
        			else{
        				String newargument = args[i+1];
        				i++;
        				PPfoldMain.treefilename = newargument;
        				PPfoldMain.createtree=false;
        				PPfoldMain.optimizetree=false;
        			}
        		}
        		else if(argument.startsWith("--entropy")){
        			if(args.length < i+2){
        				System.out.println("Not enough arguments: sequence name (entropy) missing");
        				quitWithUsageMessage();
        			}
        			else{
        				String newargument = args[i+1];
        				i++;
        				PPfoldMain.entropycalc = true;
        				PPfoldMain.entropyString = newargument; 
        			}
        		}
        		else if(argument.startsWith("--usedata") || argument.startsWith("-d")){
        			System.out.println(argument);
        			if(args.length < i+1){
        				System.out.println("Not enough arguments: data file and/or sequence name missing");
        				quitWithUsageMessage();
        			}
        			else{
        				String newargument = args[i+1];
        				DataInfo newdatainfo = new DataInfo();
        				newdatainfo.setFileName(newargument);
        				System.out.println(newargument);
        				PPfoldMain.auxdata=true;
        				newargument = args[i+2];
        				newdatainfo.setSequenceName(newargument);
        				i = i+2;
        				System.out.println(newargument);
        				if(args.length>i+1 && args[i+1].startsWith("--dist")){
        					newargument = args[i+2];
        					if(newargument.toUpperCase().startsWith("DEFAULT")){
            					//Default dist file
            					newdatainfo.setDistFileName(PPfoldMain.defaultDataDistfile);
                				newdatainfo.setType(0); //Bar dist type.
                				newdatainfo.setContactDistance(-1);
            				}
        					else{
        						//Not default dist file
        						newdatainfo.setDistFileName(newargument);
        						System.out.println(newargument);
                				newdatainfo.setContactDistance(0); //Bar dist type.
                				newdatainfo.setContactDistance(-1);
        					}
        					i = i+2; 
        				}
        				else if(args.length>i+1 && args[i+1].startsWith("--force")){
        					newdatainfo.setType(2); //Force constraints
            				//add empty string so lengths will fit
            				newdatainfo.setDistFileName(new String(""));
            				newdatainfo.setContactDistance(-1);
        					i = i+1; 
        				}
        				else if(args.length>i+1 && args[i+1].startsWith("--direct")){
            				//add empty string so lengths will fit
        					newdatainfo.setDistFileName(new String(""));
        					newdatainfo.setType(1); //Force constraint
                			newdatainfo.setContactDistance(-1);
            				i = i+1;
        				}
        				else{
        					System.out.println("Problem: confused when reading arguments");
        					quitWithUsageMessage();
        				}
        				PPfoldMain.datainfo.add(newdatainfo);
        			}
        		}
        		else if(argument.startsWith("--contactdistance") || argument.startsWith("-cd")){
        			DataInfo newdatainfo = new DataInfo();
        			if(args.length < i+2){
        				System.out.println("Not enough arguments: contact distance missing");
        				quitWithUsageMessage();
        			}
        			else{
        				String newargument = args[i+1];
        				i++;
        				PPfoldMain.auxdata=true;
        				newdatainfo.setFileName(null); //file name is null. will just generate an empty forced constraints
        				newdatainfo.setDistFileName(new String(""));
            			newdatainfo.setType(2); //Force constraint
        				newdatainfo.setContactDistance(Integer.parseInt(newargument));
        				newdatainfo.setSequenceName(null); //Sequence name is null. will just generate an empty forced constraints
        			}
        		}
        		
        		else if(argument.startsWith("--mle") || argument.startsWith("-m")){
        			PPfoldMain.optimizetree=true;
        		}
        		else if(argument.startsWith("--iterlim")||argument.startsWith("-l")){
        			if(args.length < i+2){
        				System.out.println("Not enough arguments: number of iterations missing");
        				quitWithUsageMessage();
        			}
        			else{
        				String newargument = args[i+1];
        				i++;
        				PPfoldMain.iterlimit = Integer.parseInt(newargument);
        			}
        		}
        		else if(argument.startsWith("--onlyCT")){
        			PPfoldMain.onlyCT=true;
        		}
        		else if(argument.startsWith("-f")){
        			String newargument = args[i+1];
    				i++;
    				PPfoldMain.exportfilehandle= newargument;
        			PPfoldMain.specialname=true;
        		}
        		
        		else if(argument.startsWith("--help")){
        			quitWithUsageMessage();
        		}
        		else if(argument.startsWith("--verbose")||argument.startsWith("-v")){
        			PPfoldMain.verbose = true; 
        		}
        		else if(argument.startsWith("--exports")||argument.startsWith("-e")){
        			PPfoldMain.exportson = true;
        		}
        		else if(argument.startsWith("--paramfl")){
        			if(args.length < i+2){
        				System.out.println("Not enough arguments: parameter file missing");
        				quitWithUsageMessage();
        			}
        			else{
        				String newargument = args[i+1];
        				i++;
        				PPfoldMain.paramfilename = newargument;
        			}
        		}
        		else if(argument.startsWith("--scfgjnr")||argument.startsWith("-s")){
        			if(args.length < i+2){
        				System.out.println("Not enough arguments: number of SCFG divisions missing");
        				quitWithUsageMessage();
        			}
        			else{
        				String newargument = args[i+1];
        				i++;
        				PPfoldMain.scfgdivisions = Integer.parseInt(newargument);
        			}
        		}
        		else if(argument.startsWith("--phyljnr")||argument.startsWith("-p")){
        			if(args.length < i+2){
        				System.out.println("Not enough arguments: number of phylogenetic divisions missing");
        				quitWithUsageMessage();
        			}
        			else{
        				String newargument = args[i+1];
        				i++;
            			PPfoldMain.phylodivisions = Integer.parseInt(newargument);	
        			}
        		}
        		else if(argument.startsWith("--proccnt")||argument.startsWith("-c")){
        			if(args.length < i+2){
        				System.out.println("Not enough arguments: number of processors to use missing");
        				quitWithUsageMessage();
        			}
        			else{
        				String newargument = args[i+1];
        				i++;
        				PPfoldMain.nrprocessors = Integer.parseInt(newargument);	
            			if(PPfoldMain.nrprocessors > Runtime.getRuntime().availableProcessors()){
            				PPfoldMain.nrprocessors = Runtime.getRuntime().availableProcessors();
            				System.out.println("WARNING: too large processor count requested. " +
            						"Using maximum available number of processors (" + PPfoldMain.nrprocessors + ").");
            			}	
        			}
        		}
        		else if(argument.startsWith("--outputd")||argument.startsWith("-o")){
        			if(args.length < i+2){
        				System.out.println("Not enough arguments: output directory missing");
        				quitWithUsageMessage();
        			}
        			else{
        				String newargument = args[i+1];
        				i++;
        				PPfoldMain.outputdir = newargument;
            			try {
            				new File(PPfoldMain.outputdir).mkdir();
            			}
            			catch (Exception e1) {
            				System.out.println("Couldn't open output directory! Quitting...");
            				throw new Exception(e1);
            			}
        			}
        			
        		}
        		else{
        			System.out.println("Incorrect input argument: " + argument);
        			quitWithUsageMessage(); 
        		}
        	}
        }
        else{
        	System.out.println("Incorrect input arguments.");
        	quitWithUsageMessage(); 
        }
		
		return null;
		
	}
	
	
	public static void quitWithUsageMessage(){
		System.out.println();
    	System.out.println("Correct usage: ");
    	System.out.println();
    	System.out.println("java -jar PPfold.jar [ALIGNMENTFILE] [OPTIONS]");
    	System.out.println("");
    	System.out.println("ALIGNMENTFILE: optional, must be in FASTA format.");
    	
    	System.out.println();
    	
    	System.out.println("OPTIONS: ");
    	System.out.println("--help or -h:      triggers this usage message");
    	
    	System.out.println("--paramfl or -p:   optional argument specifying the input parameter file.");
    	System.out.println("                   There must be a space between the argument and the filename.");
    	System.out.println("                   Please refer to PPfold documentation for the correct format.");
    	System.out.println("                   If none specified, default values (pfold) will be used.");
    	
    	System.out.println("--tree or -t:      optional argument specifying the input tree file.");
    	System.out.println("                   There must be a space between the argument and the filename.");
    	System.out.println("                   The nodes of the tree must exactly match the sequences in the alignment.");
    	System.out.println("                   The tree must be in Newick format.");
    	System.out.println("                   If none specified, the MLE tree will be calculated from the alignment.");
    	
    	System.out.println("--mle or -m:       Trigger maximum likelihood estimation of the branch lengths");
    	System.out.println("                   Only needed if there is a non-MLE input tree.");
    	System.out.println("                   (The default options will calculate the MLE tree from the alignment,");
    	System.out.println("                   but if a tree is given, it will by default not be optimized.)");
    	
    	System.out.println("--iterlim or -l:   Maximum number of iterations in maximum likelihood estimation of the tree.");
    	
    	System.out.println("--verbose or -v:   optional argument, produces lots of extra output. ");   	
    	
    	System.out.println("--exports or -e:   optional argument, turns all exports on.");
    	System.out.println("                   If none specified, the following files will be exported: ");
    	System.out.println("                   .st, .seq, .lseq, .ct and .newick");

    	System.out.println("--outputd or -o:   Specifies output folder ");
    	System.out.println("                   There must be a space between the argument and the folder name.");
    	System.out.println("                   Do not add / or \\ at the end!");
    	
    	System.out.println("-f:                Specifies an export file handle different from the default.");
    	
		System.out.println("--scfgjnr or -s:   Number of SCFG divisions (default=available cores*8)");
		
		System.out.println("--phyljnr or -p:   Number of phylogenetic divisions (default=available cores*2)" );
		
		System.out.println("--proccnt or -c:   Number of processors to use (default=max)");
    	
    	System.out.println("--entropy:         calculates the information entropy for the named sequence");
    	
    	System.out.println("--usedata or -d:   Uses the chemical probing data file specified in the predictions.");
    	System.out.println("                   Usage: --usedata [shapefile] [sequencename]");
    	System.out.println("                   Must also specify data type, as the next argument, as: ");
    	System.out.println("                   --dist DEFAULT   : interprets the data as data from SHAPE experiments.");
    	System.out.println("                   --dist filename  : interprets the data as from a generic probing experiment," );
    	System.out.println("                                      uses filename as the structure distribution data.");
    	System.out.println("                   --force          : interpets the data as hard constraints (mfold-style)");
    	System.out.println("                   --direct         : interprets the data as posterior probabilities (advanced)");
    	System.out.println("                   For more details, see the PPfold website.");
    	
    	System.out.println("--contactdistance or -cd: limits the contact distance to the specified number");
    	System.out.println("                          (Might not be exactly fulfilled in the final structure due to ");
    	System.out.println("                          column removal)");
    	
    	System.out.println("--onlyCT:          Only exports a CT file.");
    	
    	System.out.println();
    	System.out.println("If no arguments are given at all, PPfold will attempt to trigger the GUI.");

    	System.out.println();
    	System.out.println("Quitting...");
    	System.exit(0);
	}
	
	
}
