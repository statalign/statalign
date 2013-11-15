package statalign.postprocess.gui;

import javax.swing.*;

import org.apache.commons.math3.util.Pair;

import java.awt.*;

import statalign.model.subst.SubstitutionModel;
import statalign.postprocess.Postprocess;
import statalign.postprocess.Track;

/**
 * This is the graphical interface for showing alignments.
 * 
 * @author miklos,novak
 *
 */
public class AlignmentGUI extends JPanel{

/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//	JPanel pan;
//	CurrentAlignment owner;
	SubstitutionModel subst;
	Postprocess owner = null;
	/**
	 *  The alignment that is written onto the screen is stored in this array
	 */
	public String[] alignment = null;
	public String[] sequenceNames = null;
	/**
	 * The title of the current analysis
	 */
	public String title;
	/**
	 * These are the posterior decoding values
	 */
	public double[] decoding = null;
	private int aligNum = 15;
	private int aligLen = 200;

	/** Values greater than SCALE_FACTOR * mean will be cropped */
	static double SCALE_FACTOR = 2.5;
	static final int COLUMN_WIDTH = 9;
	static final int FONT_HEIGHT = 15;
	static final int OFFSET_X = 10;
	static final int OFFSET_Y = 10;
	static final int TITLE_Y = 20;

	//static final String magentaCharacters = "RHK";
	//static final String     redCharacters = "AVFPMILW";
	//static final String    blueCharacters = "DE";
	//static final String   greenCharacters = "STYCNGQ";

/**
 * It initializes a new panel for alignments.
 * 
 * @param title This is the title of the current analysis
 * @param subst The substitution model defines the background colors of characters
 */
	public AlignmentGUI(String title, SubstitutionModel subst) {//CurrentAlignment inp, String[] s){
//		String a = "";
//		for(int i = 0; i < s.length; i++){
//			a += s[i]+"\n";
//		}
//		this.pan = pan;
//		this.owner = inp;
		this.title = title;
		this.subst = subst;
	}

	/**
	 * It initializes a new panel for alignments.
	 * 
	 * @param title This is the title of the current analysis
	 * @param subst The substitution model defines the background colors of characters
	 * @param _owner The Postprocess plugin that owns this GUI.
	 */
	public AlignmentGUI(String title, SubstitutionModel subst, Postprocess _owner) {//CurrentAlignment inp, String[] s){
//			String a = "";
//			for(int i = 0; i < s.length; i++){
//				a += s[i]+"\n";
//			}
//			this.pan = pan;
//			this.owner = inp;
		this.title = title;
		this.subst = subst;
		owner = _owner;
	}
//	private static boolean allTab(String[] s, int p){
//		boolean b = true;
//		for(int i = 0; i < s.length && b; i++){
//			b = s[i].charAt(p) == '\t';
//		}
//		return b;
//	}

	/**
	 * This function updates the graphics
	 */
	public void paintComponent(Graphics gr){
		// System.out.println("Updating current alignment");
//		super.paintComponent(gr);
		//System.out.println("Updating current alignment");
		double maxScore = 1.0;
		//setSize((int)(owner.pan.getSize().width*0.9),(int)(owner.pan.getSize().height*0.9));
		// text.setRows(pan.getHeight()/13);
		// text.setColumns((int)(pan.getWidth()/6.6));
		if(alignment != null){
			//  text.replaceRange(owner.alignment, 0, length);
			//length = owner.alignment.length();
		}
		Graphics2D g = (Graphics2D)gr;

		g.setBackground(Color.WHITE);
		g.clearRect(0, 0, this.getWidth(), this.getHeight());
//		alignment = owner.allAlignment;
//		title = owner.title;
		if(alignment != null && alignment[0] != null) {

			//int colHeight = (decoding == null ? 0 :  Math.min(200, g.getClipBounds().height - (2 * OFFSET_Y + TITLE_Y + alignment.length * FONT_HEIGHT)));
			
			// Changed this to fixed height:
			boolean extraTracks = owner != null && owner.getTracks() != null && owner.getTracks().size() > 0;
			int colHeight = ((decoding == null && !extraTracks) ? 0 : 130);
			int maxNameLength = 0;
			aligNum = alignment.length;			
			for (int i=0; i<aligNum; i++) {
				//System.out.println(sequenceNames[i]+"...");
				int l = sequenceNames[i].length();
				if (l > 0) {
					maxNameLength = (l > maxNameLength) ? l : maxNameLength;
				}
			}
			aligLen = alignment[0].length() + maxNameLength + 5;
			setSize(OFFSET_X + COLUMN_WIDTH * aligLen + 30, OFFSET_Y + TITLE_Y + FONT_HEIGHT * aligNum+30);
			//Find the maximum length of the sequence names
//            int[] nameLength = new int[alignment.length];
//            int maxSeqNameLength = 0;
//            for (int i=0; i<alignment.length; i++) {
//            	int j=0;
//            	while(alignment[i].charAt(j++) != '\t') { nameLength[i] = j + 1; }
//            	maxSeqNameLength = j > maxSeqNameLength ? j : maxSeqNameLength; 
//            }
//            int maxSNL = maxSeqNameLength; // Shorthand 
            
//			int tab =alignment[0].length()-2;
//			while(!allTab(alignment, tab)){
//				tab--;
//			}
			//System.out.println("Tab: "+tab);

			g.setColor(Color.BLACK);
			g.setFont(new Font("SANS_SERIF", Font.BOLD, 16));
			//g.drawString(title, OFFSET_X, TITLE_Y);
			g.setFont(new Font("MONOSPACED", Font.PLAIN, 12));
//			g.setClip(0,0, COLUMN_WIDTH * alignment[0].length() + 2 * OFFSET_X, g.getClipBounds().height);			

			for (int i = 0; i < alignment.length; i++) {
				//System.out.println(alignment[i]);
				g.drawString(sequenceNames[i],OFFSET_X,colHeight + OFFSET_Y + TITLE_Y + FONT_HEIGHT * (i+1));
//				if (sequenceNames[i].length() < maxNameLength) {
//					for (int j=0; j<maxNameLength-sequenceNames[i].length(); j++) {
//						g.drawString(' ',OFFSET_X + COLUMN_WIDTH * j);
//					}
//				}
			//	System.out.println("i = "+i+", maxNameLength = "+maxNameLength);

				for (int j = 0; j < alignment[i].length(); j++) {
					char ch = alignment[i].charAt(j);
//					int jOffset = j;
//					if (j==nameLength[i]) { 
//						jOffset = maxSNL;
//					}
//					if (j>=maxNameLength){
						g.setColor(subst.getColor(ch));
						g.fillRect(OFFSET_X + COLUMN_WIDTH * (j+maxNameLength),
								colHeight + OFFSET_Y + TITLE_Y + FONT_HEIGHT * i + 3, 
								COLUMN_WIDTH,
								FONT_HEIGHT);
						//System.out.println((OFFSET_X + COLUMN_WIDTH * j)+" "+(colHeight + OFFSET_Y + TITLE_Y + FONT_HEIGHT * i + 3));
					//}
					g.setColor(Color.BLACK);
					g.drawString(ch + "",
							OFFSET_X + COLUMN_WIDTH * (j+maxNameLength),
							colHeight + OFFSET_Y + TITLE_Y + FONT_HEIGHT * (i+1));

				}
			}
			if(decoding != null){
//				System.out.println(decoding[1]);
				//System.out.println("Printing decoding.");
				Track decodingTrack = new Track(Color.BLUE,decoding,1,0,1);				
				plotScore(g,decodingTrack,maxNameLength,colHeight);								
			}
			if (extraTracks) {
				//System.out.println("Printing "+owner.getTracks().size()+" extra track(s).");
				for (Track track : owner.getTracks()) {
					plotScore(g,track,maxNameLength,colHeight);
				}
			}

		}
		else{
			g.setColor(Color.BLACK);
			g.setFont(new Font("SANS_SERIF", Font.BOLD, 16));
			g.drawString("Waiting for data...", OFFSET_X, TITLE_Y);
		}
//		this.updateUI();
	}

	private void plotScore(Graphics2D g, Track track, int maxNameLength, int colHeight) {
		//int colHeight = Math.min(100, g.getClipBounds().height - (2 * OFFSET_Y + TITLE_Y + alignment.length * FONT_HEIGHT));
		//System.out.println("Score length = "+score.length);
		Color oldColor = g.getColor();
		g.setColor(track.color);
		double maxScore = track.max;
		if (maxScore > SCALE_FACTOR * track.mean) maxScore = SCALE_FACTOR * track.mean;
		for(int i = 1; i < track.scores.length; i++){
			if (Double.isNaN(track.scores[i-1])) { continue; }
			if (Double.isNaN(track.scores[i])) { ++i; continue; }
			g.drawLine( OFFSET_X + (maxNameLength + i - 1) * COLUMN_WIDTH + COLUMN_WIDTH / 2,
					(int)(OFFSET_Y + TITLE_Y + colHeight - Math.round(track.scores[i-1] * (colHeight) / maxScore)),
//					(int)(OFFSET_Y + TITLE_Y - Math.round(decoding[i-1] * (TITLE_Y) / max)),
					OFFSET_X + (maxNameLength + i) * COLUMN_WIDTH + COLUMN_WIDTH / 2,
					(int)(OFFSET_Y + TITLE_Y + colHeight - Math.round(track.scores[i] * (colHeight) / maxScore))
//					(int)(OFFSET_Y + TITLE_Y - Math.round(decoding[i] * (TITLE_Y) / max))
					);
		}
		g.setColor(oldColor);
	}
	/**
	 * This function tells the minimum size of the panel
	 */
	public Dimension getMinimumSize(){
		return getPreferredSize();
	}

	/**
	 * This function tells the preferred size of the panel
	 */
	public Dimension getPreferredSize() {
		return new Dimension(OFFSET_X + COLUMN_WIDTH * aligLen + 30, OFFSET_Y + TITLE_Y + FONT_HEIGHT * aligNum+30);
//		if (alignment != null && alignment[0] != null) {
//			return new Dimension((alignment[0].length() + 3) * COLUMN_WIDTH + 6 * OFFSET_X, 100);
//		} else {
//			return new Dimension(0,100);
//		}
	}


}

/*


package statalign.postprocess.gui;

import java.text.*;
import java.awt.*;
import java.awt.font.*;
import java.awt.geom.*;
import javax.swing.*;


import statalign.postprocess.plugins.CurrentAlignment;

public class CurrentAlignmentGUI extends JPanel{

	CurrentAlignment currentAlignment;
	JPanel panel;
	public JTextArea text;

	public CurrentAlignmentGUI(JPanel panel, CurrentAlignment currentAlignment) {
		this.currentAlignment = currentAlignment;
		this.panel = panel;
		text = new JTextArea(20,20);
	}

	private String[] alignment = null;
	private String title;
	private int numCols;
	private int numRows;

	static final int COLUMN_WIDTH = 9;
	static final int FONT_HEIGHT = 15;
	static final int OFFSET_X = 10;
	static final int OFFSET_Y = 10;
	static final int TITLE_Y = 20;

	static final String magentaCharacters = "RHK";
	static final String     redCharacters = "AVFPMILW";
	static final String    blueCharacters = "DE";
	static final String   greenCharacters = "STYCNGQ";

	public CurrentAlignmentGUI() {
	}

	public void updateData(String[] alignment, String title) {
		this.alignment = alignment;
		this.title = title;
		numCols = alignment[0].length() + 1;
		numRows = alignment.length;
	 }

	  public void paintComponent(Graphics gr){
		  super.paintComponent(gr);
		  text.setColumns((int)(panel.getWidth()/6.6));
		  text.setRows(panel.getHeight()/13);
	  }
	/*
	public void paintComponent(Graphics gr) {
		System.out.println("Updating current alignment");
		super.paintComponent(gr);
		Graphics2D g = (Graphics2D)gr.create();

		g.setBackground(Color.WHITE);
		g.clearRect(0, 0, this.getWidth(), this.getHeight());
		alignment = currentAlignment.alignment;
		if(alignment != null) {
		    int colHeight = Math.min(100, g.getClipBounds().height - (2 * OFFSET_Y + TITLE_Y + numRows * FONT_HEIGHT));
		   // double max = 1;
		    //for (int i = 0; i < numCols; i++) {
			//max = Math.max(max, rs.getValue(i));
		   // }
		    //System.out.println("Maximum "+max);
	//	    g.setColor(Color.LIGHT_GRAY);
		//    for (int i = 0; i < 5; i++) {
			//g.drawLine( OFFSET_X+24*COLUMN_WIDTH,
				//    OFFSET_Y + TITLE_Y + Math.round(i * (float)colHeight / 4),
				  //  OFFSET_X + numCols * COLUMN_WIDTH,
				    //OFFSET_Y + TITLE_Y + Math.round(i * (float)colHeight / 4));
		    //}

		    g.setColor(Color.BLACK);
		    g.setFont(new Font("SANS_SERIF", Font.BOLD, 16));
		    g.drawString(title, OFFSET_X, TITLE_Y);
		    g.setFont(new Font("MONOSPACED", Font.PLAIN, 12));
		    g.setClip(0,0, COLUMN_WIDTH * numCols + 2 * OFFSET_X, g.getClipBounds().height);			

		    for (int i = 0; i < numRows; i++) {
		    	for (int j = 0; j < alignment[i].length(); j++) {
		    		if(j>23){
		    			g.setColor(Color.LIGHT_GRAY);
		    			if(magentaCharacters.indexOf(alignment[i].charAt(j)) != -1){
		    				g.setColor(Color.MAGENTA);
		    			}
		    			else if(redCharacters.indexOf(alignment[i].charAt(j)) != -1){
		    				g.setColor(Color.RED);
		    			}
		    			else if(blueCharacters.indexOf(alignment[i].charAt(j)) != -1){
		    				g.setColor(Color.BLUE);
		    			}
		    			else if(greenCharacters.indexOf(alignment[i].charAt(j)) != -1){
		    				g.setColor(Color.GREEN);
		    			}
		    			g.fillRect(OFFSET_X + COLUMN_WIDTH * i,colHeight + OFFSET_Y + TITLE_Y + FONT_HEIGHT * j + 3,COLUMN_WIDTH,FONT_HEIGHT);
		    		}
		    		g.setColor(Color.BLACK);
		    		g.drawString(alignment[i].charAt(j) + "",
		    				OFFSET_X + COLUMN_WIDTH * j,
		    				colHeight + OFFSET_Y + TITLE_Y + FONT_HEIGHT * (i+1));

		    	}
		    }
		}

		//	if(i>24){
			//    g.setColor(Color.BLUE);
			  //  g.drawLine( OFFSET_X + (j - 1) * COLUMN_WIDTH + COLUMN_WIDTH / 2,
				//	(int)(OFFSET_Y + TITLE_Y + colHeight - Math.round(rs.getValue(j-1) * (colHeight) / max)),
					//OFFSET_X + j * COLUMN_WIDTH + COLUMN_WIDTH / 2,
		//			(int)(OFFSET_Y + TITLE_Y + colHeight - Math.round(rs.getValue(j) * (colHeight) / max))
			//		);
		//	}
		  //  }
		    //}

		this.updateUI();
	}

	public Dimension getMinimumSize(){
		return getPreferredSize();
	}

	public Dimension getPreferredSize() {
		if (alignment != null) {
			return new Dimension((alignment[0].length() + 3) * COLUMN_WIDTH, 100);
		} else {
			return new Dimension(0,100);
		}
	}

}
 */