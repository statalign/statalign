package statalign.postprocess.plugins;

import java.awt.BorderLayout;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import fr.orsay.lri.varna.exceptions.ExceptionNonEqualLength;

import statalign.base.InputData;
import statalign.base.State;
import statalign.postprocess.gui.StructureGUI;
import statalign.postprocess.gui.TestGUI;
import statalign.postprocess.utils.RNAFoldingTools;

public class Structure extends statalign.postprocess.Postprocess {

	private JPanel pan;
	private StructureGUI gui;
	
	private static float[][] probMatrix;
	private static String currentSequence;
	private static String currentDotBracketStructure;
	
	public Structure() {
		sampling = true;
		screenable = true;
		outputable = true;
		postprocessable = true;
		postprocessWrite = true;
	}
	
	
	@Override
	public String getTabName() {
		return "Consensus Structure";
	}

	@Override
	public Icon getIcon() {
		// TODO Auto-generated method stub
		return new ImageIcon("icons/test1.png");
	}

	@Override
	public JPanel getJPanel() {
		// TODO Auto-generated method stub
		pan = new JPanel(new BorderLayout());
		return pan;
	}

	@Override
	public String getTip() {
		// TODO Auto-generated method stub
		return "Consensus structure";
	}

	@Override
	public void setSampling(boolean enabled) {
		// TODO Auto-generated method stub
		sampling = enabled;
		
	}
	
	public static void updateMatrix(float[][] newMatrix) {
		probMatrix = newMatrix;
	}
	
	public static void updateStructure() {
		RNAFoldingTools rnaTools = new RNAFoldingTools();
		int[] pairedSites = rnaTools.getPosteriorDecodingConsensusStructureMultiThreaded(probMatrix);
		currentDotBracketStructure = RNAFoldingTools.getDotBracketStringFromPairedSites(pairedSites);
	}
	
	public static void updateSequence() {
		currentSequence = PPFold.getSequenceByName(PPFold.getSequences(), PPFold.getRefName()).replaceAll("-", "");;
	}
	
	@Override
	public void beforeFirstSample(InputData input) {
		if(show) {
			pan.removeAll();
			JScrollPane scroll = new JScrollPane();
			try {
				scroll.setViewportView(gui = new StructureGUI(currentSequence, currentDotBracketStructure));
			} catch (ExceptionNonEqualLength e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}// ,
																				// mcmc.tree.printedAlignment()));
			pan.add(scroll, BorderLayout.CENTER);
			pan.getParent().validate();
		}
	}
	
	@Override
	public void newSample(State state, int no, int total) {
		
		if(sampling) {
			//updateSequence();
			//updateStructure();
			//gui.updateAndDraw(currentSequence, currentDotBracketStructure);
		}
	}
	
	@Override
	public void afterLastSample() {
		
	}

}
