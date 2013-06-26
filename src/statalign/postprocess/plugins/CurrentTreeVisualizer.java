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
import statalign.postprocess.gui.treeviews.NetworkTreeView;
import statalign.postprocess.gui.treeviews.TreeView;
import statalign.postprocess.gui.treeviews.VerticalPhylogramTreeView;
import statalign.postprocess.utils.NewickParser;

public class CurrentTreeVisualizer extends TreeVisualizer {

    // Functions

    public CurrentTreeVisualizer() {
        // Configuration
        outputable = true;
        postprocessable = true;
        postprocessWrite = true;
        sampling = true;
        rnaAssociated = false;
    }
    
    @Override
    public void init() {
    	if(show) {
			super.init(new TreeView[]{
	                new VerticalPhylogramTreeView(),
	                new HorizontalCladogramTreeView(),
	                new HorizontalPhylogramTreeView(),
	                new NetworkTreeView()
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
        return "Tree";
    }

    @Override
    public String getTip() {
        return "Current tree in the Markov chain";
    }

    @Override
    public double getTabOrder() {
        return 3.0d;
    }

    @Override
    public void newPeek(State state) {
//    	System.out.println("CurrentTreeVisualizer: New peek.");
        String treeString = state.getNewickString();
        NewickParser parser = new NewickParser(treeString);
        TreeNode node = parser.parse();

        if(show) {
        	for (TreeView view : treeViews) {
        		view.newSample(node);
        	}

        	refreshGUI();
        }
    }

    @Override
    public String getFileExtension() {
        return "tree";
    }

    @Override
    public void beforeFirstSample(InputData input) {
        super.beforeFirstSample(input); // Mandatory.

        String[] names = new String[input.seqs.seqNames.size()];
		
		for (int i = 0; i < names.length; i++) {
			names[i] = input.seqs.seqNames.get(i).replaceAll(" ", "");
		}
		if (postprocessWrite) {
			try {
		        outputFile.write("#nexus\n\nBEGIN Taxa;\nDIMENSIONS ntax="
						+ names.length + ";\nTAXLABELS\n");
				for (int i = 0; i < names.length; i++) {
					outputFile.write("[" + (i + 1) + "] '" + names[i] + "'\n");
				}
				outputFile.write(";\nEND; [Taxa]\nBEGIN Trees;\n");
			} catch (IOException e) {
	            e.printStackTrace();
	        }
		}
    }

    @Override
    public void newSample(State state, int no, int total) {
//    	System.out.println("CurrentTreeVisualizer: New sample.");
        String sampledTree = state.getNewickString();

        // TODO: START OF FIX
        NewickParser parser = new NewickParser(sampledTree);
        TreeNode node = parser.parse();

        if(show) {
        	for (TreeView view : treeViews) {
        		view.newSample(node);
        	}

        	refreshGUI();
        }
        // TODO: END OF FIX 
        
        if (sampling) {
            try {
                file.write("Sample " + no + "\tTree:\t" + sampledTree + "\n");
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        if (postprocessWrite) {
            try {
                outputFile.write("[" + no + "] tree 'tree_" + no + "'= "
        				+ sampledTree + "\n");
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    @Override
    public void afterLastSample() {
		if (postprocessWrite) {
			try {								
				outputFile.write("END; [Trees]\n" + "BEGIN st_Assumptions;\n"
						+ "\ttreestransform=TreeSelector;\n"
						+ "\tsplitstransform=EqualAngle;\n"
						+ "\tSplitsPostProcess filter=dimension value=4;\n"
						+ "\tautolayoutnodelabels;\nEND; [st_Assumptions]");
				outputFile.flush();
				outputFile.close();
			} catch (IOException ioex) {
				ioex.printStackTrace();
			}
		}
    }

}
