package statalign.postprocess.plugins;

import java.io.IOException;
import java.util.ArrayList;

import javax.swing.Icon;
import javax.swing.ImageIcon;

import statalign.base.InputData;
import statalign.base.State;
import statalign.model.ext.ModelExtManager;
import statalign.postprocess.Postprocess;
import statalign.postprocess.TreeVisualizer;
import statalign.postprocess.gui.treeviews.HorizontalCladogramTreeView;
import statalign.postprocess.gui.treeviews.HorizontalPhylogramTreeView;
import statalign.postprocess.gui.treeviews.TreeView;
import statalign.postprocess.plugins.contree.CNetwork;
import statalign.postprocess.plugins.contree.CTMain;
import statalign.postprocess.plugins.contree.CTree;
import statalign.postprocess.utils.NewickParser;

public class ConsensusTree extends TreeVisualizer {

    // Variables

    private CTMain main;
	public CNetwork network;

    //private ArrayList<String> consensusTrees;

    // Functions

    public ConsensusTree() {
        // Configuration
        outputable = true;
        postprocessable = true;
        postprocessWrite = false;
        rnaAssociated = false;
        noOfSamples = 1;
		useInParallelMode = false;
    }
    
    @Override
    public void init(ModelExtManager manager) {
    	if(show) {
			super.init(new TreeView[]{
    				new HorizontalCladogramTreeView(),
    				new HorizontalPhylogramTreeView()
    		});	
		}
    }

    @Override
    public Icon getIcon() {
        return new ImageIcon(ClassLoader.getSystemResource("icons/tree.gif"));
//        return new ImageIcon("icons/tree.gif");
    }

    @Override
    public String getTabName() {
        return "Consensus tree";
    }

    @Override
    public String getTip() {
        return "The consensus tree based on all of the samples";
    }

    @Override
	public String getFileExtension() {
		return "ctree";
	}

	@Override
    public double getTabOrder() {
        return 6.0d;
    }

    @Override
    public void newSample(State state, int no, int total) {
		if (!active) return;
    	if (state.isBurnin) return;
    	
    	++noOfSamples;
    	String treeString = state.getNewickStringWithLabels();
        NewickParser parser = new NewickParser(treeString);
        TreeNode root = parser.parse();
        
        // Initialization (only happens in the beginning)
        if (!main.initialized) { // The first sample has arrived - initialize!
            main.initialize(root, total);
            main.addNewTree(root);
            return;
        }
        main.addNewTree(root);
        CTree output = main.constructMajorityTree();
        TreeNode outputRoot = output.getRoot();
        String outputRootString = outputRoot.toStringWithProbs(noOfSamples);       
        
        network = main.constructNetwork(output); 
        
        // Logging.
        if (sampling) {
            try {
                file.write("Sample " + no + "\tConsensus tree:\t" + outputRootString + "\n");
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        if (postprocessWrite) {
            try {
           		outputFile.write(outputRootString + "\n");
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        if (show) {
	        for (TreeView view : treeViews) {
	            view.newSample(outputRoot);
	        }
	
	        refreshGUI();
        }
    }

    @Override
    public void beforeFirstSample(InputData input) {
        super.beforeFirstSample(input); // Mandatory.

        main = new CTMain();
        //consensusTrees = new ArrayList<String>();
    }

    @Override
    public void afterLastSample() {
        super.afterLastSample();
        
		if (postprocessWrite) {
			try {			
				outputFile.flush();
				outputFile.close();
			} catch (IOException ioex) {
				ioex.printStackTrace();
			}
		}
    }

}
