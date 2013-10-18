package statalign;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;

import javax.swing.JOptionPane;

import statalign.base.MainManager;
import statalign.postprocess.Postprocess;
import statalign.ui.ErrorMessage;
import statalign.ui.MainFrame;

/**
 * <p>The entry point of the application. If the program is called from command line and
 * arguments are given, it runs in terminal mode.
 * 
 * <p>If no arguments are given or launched from a jar file, it opens a main frame.
 *
 * @author miklos, novak
 * 
 */
public class StatAlign{

	/**
	 * StatAlign version data.
	 */

	public static final int majorVersion = 3;
	public static final int minorVersion = 0;
	public static final String version = "v3.0";
	
	public static boolean allowVersionCheck = true;
	
	public static final String webPageURL = "http://statalign.github.io/";

	/**
	 * Only method of the class.
	 * If a Fasta input file is given as argument, it runs in terminal mode
	 * (without graphical interface), otherwise it launches the main GUI
	 * of the program.
	 * 
	 * @param args [0]: the input file name containing the sequences in Fasta format
	 * @throws IOException
	 */
	public static void main(String args[]) {
			
		System.out.println("StatAlign "+version);
		
		if(args.length != 0) {
			// console mode
			MainManager manager = new MainManager(null);
			CommandLine cl = new CommandLine(false);
			cl.setVerbose(true);
			if(cl.fillParams(args, manager) > 0)
				System.exit(1);
			
			manager.init(cl.pluginParameters);			
			manager.start();

		} else {
			// GUI mode			
			MainFrame mf = null;
			try {
				mf = new MainFrame();				
			} catch(Exception e) {
				ErrorMessage.showPane(mf, e, true);
			}
		}
	}

}