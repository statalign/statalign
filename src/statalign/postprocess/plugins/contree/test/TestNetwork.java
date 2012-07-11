package statalign.postprocess.plugins.contree.test;
import java.awt.BorderLayout;
import java.util.ArrayList;
import java.util.BitSet;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import statalign.postprocess.gui.CNetworkView;
import statalign.postprocess.plugins.contree.CNetwork;
import statalign.postprocess.plugins.contree.CTMain;
import statalign.postprocess.plugins.contree.CTree;
/**
 * Creates trees that create a network that exhibits many incompatible splits in order to test the network drawing feature.
 * NOTE: make sure trees only contain compatible splits!
 * Run from the ConsensusNetworkTester file...
 * 
 * @author wood
 *
 */
public class TestNetwork {
	// Consensus Tree
	private CTree conTree;
	// Consensus Network (derived from Consensus Tree!)
	private CNetwork conNetwork;
	// Main Consensus Thread for calculating consensus tree and network
	private CTMain main;
	// Array of trees
	private ArrayList<ArrayList<BitSet>> trees;
	// Panel to store GUI in
	private JPanel pan;
	// GUI of network output
	private CNetworkView gui;
	// Frame to store panel in
	private JFrame frame;
	/**
	 * Creates the test network by creating trees and storing them in the hash table using test functions in CTMain, running the consensus tree and network creating algorithms and then outputting the network in a window.
	 * 
	 *  */
	public void TestNetwork(){
			// specify the trees by splits...
			trees = new ArrayList<ArrayList<BitSet>>();
			// first tree
				ArrayList<BitSet> tree1 = new ArrayList<BitSet>();
				
					//first split
					BitSet split1a = new BitSet(10);
					split1a.set(2,true);
					split1a.set(4,true);
					tree1.add(split1a);
					
					//second split
					BitSet split1b = new BitSet(10);
					split1b.set(4,true);
					split1b.set(5,true);
					split1b.set(2,true);
					tree1.add(split1b);
					
					//third split
					BitSet split1c = new BitSet(10);
					split1c.set(6,true);
					split1c.set(7,true);
					split1c.set(8,true);
					tree1.add(split1c);
					
					//fourth split
					BitSet split1d = new BitSet(10);
					split1d.set(6,true);
					split1d.set(7,true);
					tree1.add(split1d);
					
				trees.add(tree1);
				// second tree
				ArrayList<BitSet> tree2 = new ArrayList<BitSet>();
					
					//first split
					BitSet split2a = new BitSet(10);
					split2a.set(1,true);
					split2a.set(3,true);
					split2a.set(4,true);
					tree2.add(split2a);
					//second split
					BitSet split2b = new BitSet(10);
					split2b.set(4,true);
					split2b.set(3,true);
					tree2.add(split2b);
					// third split
					BitSet split2c = new BitSet(10);
					split2c.set(4,true);
					split2c.set(3,true);
					split2c.set(1,true);
					split2c.set(2,true);
					tree2.add(split2c);
					// fourth split
					BitSet split2d = new BitSet(10);
					split2d.set(4,true);
					split2d.set(3,true);
					split2d.set(1,true);
					split2d.set(2,true);
					split2d.set(5,true);
					tree2.add(split2d);
					// fifth split
					BitSet split2e = new BitSet(10);
					split2e.set(6,true);
					split2e.set(9,true);
					tree2.add(split2e);
					
				trees.add(tree2);
				
				/*
				// third tree
				ArrayList<BitSet> tree3 = new ArrayList<BitSet>();
					//first split
					BitSet split3a = new BitSet(10);
					split3a.set(1,true);
					split3a.set(2,true);
					split3a.set(3,true);
					split3a.set(4,true);
					tree3.add(split3a);
					//second split
					BitSet split3b = new BitSet(10);
					split3b.set(4,true);
					split3b.set(3,true);
					tree3.add(split3b);
				trees.add(tree3);
				*/
			// Load the CTMain up...
			main = new CTMain();
			main.setSeed(1);
			main.InitialiseNetworkTester(10,trees.size());
			// add the trees using the test methods in CTMain...
			for(ArrayList<BitSet> tree : trees){
				main.addTestTree(tree);
			}
			// perform the calculations...
			conTree = main.constructMajorityTree();
			conNetwork = main.constructNetwork(conTree);
			//visualise!	
			pan = new JPanel(new BorderLayout());
			JScrollPane scroll = new JScrollPane();
			scroll.setViewportView(gui = new CNetworkView(scroll));
			gui.network = conNetwork;
			pan.add(scroll, BorderLayout.CENTER);
			pan.setVisible(true);
			frame = new JFrame("Network test results");
			frame.setSize(1000, 500);
			frame.setLocation(300,200);
			frame.setVisible(true);
			frame.add(pan);
		}
}
