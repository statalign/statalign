package statalign.postprocess.plugins;

import java.io.IOException;
import java.util.ArrayList;

import javax.swing.Icon;
import javax.swing.ImageIcon;

import statalign.base.InputData;
import statalign.base.State;
import statalign.postprocess.TreeVisualizer;
import statalign.postprocess.gui.treeviews.HorizontalCladogramTreeView;
import statalign.postprocess.gui.treeviews.HorizontalPhylogramTreeView;
import statalign.postprocess.gui.treeviews.TreeView;
import statalign.postprocess.plugins.contree.CTMain;
import statalign.postprocess.plugins.contree.CTree;
import statalign.postprocess.utils.NewickParser;

public class ConsensusTreeVisualizer extends TreeVisualizer {

    // Variables

    private CTMain main;
    private ArrayList<String> consensusTrees;

    // Functions

    public ConsensusTreeVisualizer() {
        // Configuration
        outputable = true;
        postprocessable = true;
        postprocessWrite = true;
    }
    
    @Override
    public void init() {
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
        String treeString = state.getNewickString();
        NewickParser parser = new NewickParser(treeString);
        TreeNode root = parser.parse();
        
        // Initialization (only happens in the beginning)
        if (no == 0) { // The first sample has arrived - initialize!
            main.initialize(root, total);
            main.addNewTree(root);
            return;
        }
        main.addNewTree(root);
        CTree output = main.constructMajorityTree();
        mcmc.tree.network = main.constructNetwork(output); 
        TreeNode outputRoot = output.getRoot();
        String outputRootString = outputRoot.toString();
        
        // Logging.
        consensusTrees.add(outputRootString);
        if (sampling) {
            try {
                file.write("Sample " + no + "\tConsensus tree:\t" + outputRootString + "\n");
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
        consensusTrees = new ArrayList<String>();
    }

    @Override
    public void afterLastSample() {
        super.afterLastSample();
        
		if (postprocessWrite) {
			try {
				for (String tree : consensusTrees) {
					outputFile.write(tree + "\n");
				}
				
				// TODO: should be handled in the PostprocessManager 
				//       in my opinion.
				outputFile.close();
			} catch (IOException ioex) {
				ioex.printStackTrace();
			}
		}
    }

}
