package statalign.postprocess.plugins;

import java.awt.BorderLayout;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import statalign.base.InputData;
import statalign.base.State;
import statalign.postprocess.gui.CNetworkView;
import statalign.postprocess.plugins.contree.CNetwork;
import statalign.postprocess.plugins.contree.CTMain;
/**
 * This is the postprocessmanager for showing the consensus network.
 * 
 * @author miklos, novak, wood
 *
 */
public class ConsensusNetwork extends statalign.postprocess.Postprocess{

	private int updateFrequency = 100;
	private int updateCnt;
    private CTMain main;
	public JPanel pan;
	CNetworkView gui;
	
	/**
	 * It construct a postprocess that is screenable (=can be shown on the GUI),
	 * outputable (= can be written into th elog file) but not postprocessable
	 * (= cannot generate its own output file).
	 */
	public ConsensusNetwork(){
		screenable = true;
		outputable = true;
		postprocessable = false;
		rnaAssociated = false;
	}
	
	/**
	 * It constructs a new JPanel, and returns with it
	 */
	@Override
	public JPanel getJPanel() {
		pan = new JPanel(new BorderLayout());
		return pan;
	}

	/**
	 * It generates a new icon based on the figure in file icons/calignment.gif
	 */
	@Override
	public Icon getIcon() {		
//		return new ImageIcon("icons/calignment.gif");
		return new ImageIcon(ClassLoader.getSystemResource("icons/calignment.gif"));
	}

	/**
	 * It returns with its tab name, 'Alignment'
	 */
	@Override
	public String getTabName() {
		return "Network";
	}

	/**
	 * It returns with the tip of the tabulated panel 'Current alignment in the Markov chain'
	 * (shown when mouse is moved over the panel)
	 */
	@Override
	public String getTip() {
		return "Consensus network of tree samples";
	}
	
	@Override
    public double getTabOrder() {
        return 7.0d;
    }


	@Override
	public String getFileExtension() {
		return "cnw";
	}

	/**
	 * Initializes the graphical interface.
	 */
	@Override
	public void beforeFirstSample(InputData inputData) {

		if(show) {
			pan.removeAll();
			JScrollPane scroll = new JScrollPane();
			scroll.setViewportView(gui = new CNetworkView(scroll));//, mcmc.tree.printedAlignment()));
			pan.add(scroll, BorderLayout.CENTER);
			//System.out.println("Consensus network parent: " + pan.getParent());
			pan.getParent().validate();
		}
		updateCnt = 0;

		
	}

	/**
	 * At each MCMC sampling point, the current alignment is shown on the GUI.
	 * Also, it writes the current alignment into the log file when in
	 * sampling mode.
	 */
	@Override
	public void newSample(State state, int no, int total) {
		if(show) {
			CNetwork outputNetwork = mcmc.tree.network;
			gui.network = outputNetwork;
			// gui.redraw();
			gui.repaint();
		}
		//if (show) refreshGUI();
	}
	
	/**
	 * Toggles sampling mode.
	 */
	@Override
	public void setSampling(boolean enabled){
		sampling = enabled;
	}
	
	/**
	 * Empty function, since there is nothing to do after the last sample in this postprocess thread.
	 */
	@Override
	public void afterLastSample() {
	}

}