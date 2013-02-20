package statalign.postprocess.plugins;

//import java.awt.*;

import java.awt.BorderLayout;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JPanel;

import statalign.base.InputData;
import statalign.base.McmcStep;
import statalign.base.State;
import statalign.postprocess.gui.LogLikelihoodTraceGUI;
import statalign.postprocess.utils.LogLikelihoodTraceContainer;

/**
 * This is the implementation of the postprocess that handles loglikelihoods.
 *
 * @author miklos, novak
 *
 */
public class LogLikelihoodTrace extends statalign.postprocess.Postprocess{

	JPanel pan;
	public ArrayList<LogLikelihoodTraceContainer> list;
	int current;
	int step;
	int count;
	List<Double> loglikelihoods;
	private LogLikelihoodTraceGUI gui;

	/**
	 * It constructs a postprocess handler that is  screenable (=can be shown on the GUI),
	 * outputable (= can be written into the log file) and postprocessable
	 * (= can generate its own output file).
	 *
	 * The default setting is that it writes loglikelihoods both into the logfile and its own
	 * output file.
	 */
	public LogLikelihoodTrace(){
		screenable = true;
		outputable = true;
		postprocessable = true;
		sampling = true;
		postprocessWrite = true;
		rnaAssociated = false;
	}

	/**
	 * Returns the name of the postprocess, 'Loglikelihood'.
	 */
	@Override
	public String getTabName(){
		//return null;
		return "LogLikelihood";
	}

	/**
	 * Returns the icon generated from the picture in file icons/loglikelihood1.gif.
	 */
	@Override
	public Icon getIcon(){
		return new ImageIcon(ClassLoader.getSystemResource("icons/loglikelihood1.gif"));
//		return new ImageIcon("icons/loglikelihood1.gif");
	}

	/**
	 * It generates a panel and returns it.
	 */
	@Override
	public JPanel getJPanel(){
		pan = new JPanel(new BorderLayout());
		return pan;
	}

	/**
	 * Returns the tip of the tabulated panel 'Log-likelihood trace'
	 * (shown when mouse is moved over the panel)
	 */
	@Override
	public String getTip(){
		return "Log-likelihood values";
	}

	@Override
	public double getTabOrder() {
		return 4.0d;
	}
	
	@Override
	public String getFileExtension() {
		return "ll";
	}

	/**
	 * Writes the loglikelihoods into a logfile, if sampling mode is switched on.
	 */
	@Override
	public void newSample(State state, int no, int total) {
		double logLike = state.logLike;
		if(postprocessWrite)
			loglikelihoods.add(logLike);
		if(sampling){
			try {
				file.write("Sample "+no+"\tLoglikelihood:\t"+logLike+"\n");
			} catch (IOException e) {
			}
		}
		else{
			//System.out.println("Not sampling loglikelihood!!!");
		}
	}

	/**
	 * It makes a dynamically updated list containing <300 samples evenly distributed from the steps
	 * already made. The graphical interface draws the loglikelihood chain based on this list.
	 */
	@Override
	public void newStep(McmcStep mcmcStep) {
		LogLikelihoodTraceContainer container = new LogLikelihoodTraceContainer(
				mcmcStep.newLogLike, mcmcStep.burnIn);
		if(list.size() < 300){
			list.add(container);
		} else {
			count++;
			if(count == step){
				count = 0;
				list.remove(current);
				list.add(container);
				current += 2;
				if(current > 150){
					current = 0;
					step *= 2;
				}
			}
		}
		if(show)
			gui.repaint();
	}

	/**
	 * Initialises the graphical interface.
	 */
	@Override
	public void beforeFirstSample(InputData input) {
		if(show) {
			pan.removeAll();
			gui = new LogLikelihoodTraceGUI(pan, this);
			pan.add(gui);
			pan.getParent().getParent().getParent().validate();
		}
		
		list = new ArrayList<LogLikelihoodTraceContainer>();
		current = 0;
		step = 2;
		count = 0;
		loglikelihoods = new ArrayList<Double>();
	}

	/**
	 * It writes the sampled loglikelihood values into a separate file.
	 */
	@Override
	public void afterLastSample() {
		if (postprocessWrite) {
			try {
				for (double ll : loglikelihoods) {
					outputFile.write(ll + "\n");
				}
			} catch (IOException e) {
			}
		}
	}
	
	@Override
	public List<String> getCreatedFileDescriptions() {
		return Arrays.asList("Log-likelihood trace of the samples");
	}

	/**
	 * It switches on and off the sampling mode.
	 */
	@Override
	public void setSampling(boolean enabled) {
		sampling = enabled;

	}

}