package statalign;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.nio.file.Files;

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
	public static final int minorVersion = 2;
	public static final String version = "v3.2";
	
	public static boolean allowVersionCheck = true;
	
	public static final String webPageURL = "http://statalign.github.io/";

	/** 
	 * If command line arguments are provided, terminal mode is launched
	 * (without graphical interface). Running with no arguments launches 
	 * the GUI version of the program.
	 * 
	 * Information on command line options can be obtained by running
	 * 
	 * java -jar StatAlign.jar -help
	 * 
	 * @param args (optional)
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
