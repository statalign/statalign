/*
 VARNA is a tool for the automated drawing, visualization and annotation of the secondary structure of RNA, designed as a companion software for web servers and databases.
 Copyright (C) 2008  Kevin Darty, Alain Denise and Yann Ponty.
 electronic mail : Yann.Ponty@lri.fr
 paper mail : LRI, bat 490 University Paris-Sud 91405 Orsay Cedex France

 This file is part of VARNA version 3.1.
 VARNA version 3.1 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 VARNA version 3.1 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with VARNA version 3.1.
 If not, see http://www.gnu.org/licenses.
 */
package fr.orsay.lri.varna.applications;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import fr.orsay.lri.varna.VARNAPanel;
import fr.orsay.lri.varna.controlers.ControleurInterpolator;
import fr.orsay.lri.varna.exceptions.ExceptionDrawingAlgorithm;
import fr.orsay.lri.varna.exceptions.ExceptionFileFormatOrSyntax;
import fr.orsay.lri.varna.exceptions.ExceptionModeleStyleBaseSyntaxError;
import fr.orsay.lri.varna.exceptions.ExceptionNAViewAlgorithm;
import fr.orsay.lri.varna.exceptions.ExceptionNonEqualLength;
import fr.orsay.lri.varna.exceptions.ExceptionParameterError;
import fr.orsay.lri.varna.exceptions.ExceptionUnmatchedClosingParentheses;
import fr.orsay.lri.varna.exceptions.MappingException;
import fr.orsay.lri.varna.models.VARNAConfig;
import fr.orsay.lri.varna.models.rna.Mapping;
import fr.orsay.lri.varna.models.rna.ModeleBase;
import fr.orsay.lri.varna.models.rna.ModeleStyleBase;
import fr.orsay.lri.varna.models.rna.RNA;

import fr.orsay.lri.varna.interfaces.InterfaceVARNAListener;;

public class NussinovDemo extends JFrame implements InterfaceVARNAListener {

	/**
	 * 
	 */
	private static final long serialVersionUID = -790155708306987257L;

	private static final String SEQUENCE_A =        "AGGCACGUCU";
	private static final String SEQUENCE_B =  "GAGUAGCCUC";
	private static final String SEQUENCE_C = "GCAUAGCUGC";
	
	private static final String SEQUENCE_BIG = "AAAACAAAAACACCAUGGUGUUUUCACCCAAUUGGGUGAAAACAGAGAUCUCGAGAUCUCUGUUUUUGUUUU"; 

	private static final String DEFAULT_STRUCTURE = "..........";
	// private static final String DEFAULT_STRUCTURE1 = "((((....))))";
	// private static final String DEFAULT_STRUCTURE2 =
	// "((((..(((....)))..))))";

	private VARNAPanel _vpMaster;

	private JPanel _tools = new JPanel();
	private JPanel _input = new JPanel();

	private JPanel _seqPanel = new JPanel();
	private JPanel _structPanel = new JPanel();
	private JLabel _info = new JLabel();
	private JLabel _struct = new JLabel(DEFAULT_STRUCTURE);
	private JComboBox _seq1 = new JComboBox();
	private JLabel _structLabel = new JLabel("Predicted Secondary Structure");
	private JLabel _seqLabel = new JLabel("RNA sequence");
	private JButton _goButton = new JButton("Fold");
	private JButton _switchButton = new JButton("Reset");

	private static String errorOpt = "error";
	@SuppressWarnings("unused")
	private boolean _error;

	private Color _backgroundColor = Color.white;

	@SuppressWarnings("unused")
	private int _algoCode;


	public static ModeleStyleBase createStyle(String txt) 
	{
		ModeleStyleBase result = new ModeleStyleBase();
		try {
			result.assignParameters(txt);
		} catch (ExceptionModeleStyleBaseSyntaxError e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ExceptionParameterError e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return result;
	}
	
	public void applyTo(VARNAPanel vp, ModeleStyleBase mb, int[] indices)
	{
		for(int i=0;i<indices.length;i++)
		{ 
			ModeleBase m = vp.getRNA().getBaseAt(indices[i]);
			m.setStyleBase(mb);
			if (m.getElementStructure()!=-1)
			{
				vp.getRNA().getBaseAt(m.getElementStructure()).setStyleBase(mb);
			}
		}
		vp.repaint();
	}
	
	
	public NussinovDemo() {
		super();
		try {
			_vpMaster = new VARNAPanel(getSeq(), "");
		} catch (ExceptionNonEqualLength e) {
			_vpMaster.errorDialog(e);
		}
		_vpMaster.setPreferredSize(new Dimension(600, 600));
		RNAPanelDemoInit();
	}

	private void RNAPanelDemoInit() {
		int marginTools = 250;
		Font textFieldsFont = Font.decode("MonoSpaced-BOLD-16");
		Font labelsFont = _seqLabel.getFont().deriveFont(16f);

		_seq1.setFont(textFieldsFont);
		String[] seqs = {SEQUENCE_A,SEQUENCE_B,SEQUENCE_C,SEQUENCE_BIG};
		_seq1.setModel(new DefaultComboBoxModel(seqs));
		_seq1.setEditable(true);

		setBackground(_backgroundColor);
		_vpMaster.setBackground(_backgroundColor);
		_vpMaster.addVARNAListener(this);
		//_vpSlave.setModifiable(false);


		_seqLabel.setHorizontalTextPosition(JLabel.LEFT);
		_seqLabel.setPreferredSize(new Dimension(marginTools, 15));
		_seqLabel.setFont(labelsFont);
		_structLabel.setFont(labelsFont);

		_goButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				showSolution();
				onStructureRedrawn();
			}
		});

		_switchButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
							RNA r = new RNA();
							r.setRNA("", "");
							_struct.setText("");
							_vpMaster.setTitle("");
							_vpMaster.showRNA(r);
							onStructureRedrawn();
					} 
				catch (ExceptionFileFormatOrSyntax e2) {
					e2.printStackTrace();
				} catch (ExceptionUnmatchedClosingParentheses e2) {
					// TODO Auto-generated catch block
					e2.printStackTrace();
				}
				_vpMaster.repaint();
			}
		});

		_seqPanel.setLayout(new BorderLayout());
		_seqPanel.add(_seqLabel, BorderLayout.WEST);
		_seqPanel.add(_seq1, BorderLayout.CENTER);

		_structLabel.setPreferredSize(new Dimension(marginTools, 15));
		_structLabel.setHorizontalTextPosition(JLabel.LEFT);
		_struct.setFont(textFieldsFont);
		_structPanel.setLayout(new BorderLayout());
		_structPanel.add(_structLabel, BorderLayout.WEST);
		_structPanel.add(_struct, BorderLayout.CENTER);

		_input.setLayout(new GridLayout(2, 0));
		_input.add(_seqPanel);
		_input.add(_structPanel);

		JPanel goPanel = new JPanel();
		goPanel.setLayout(new BorderLayout());

		_tools.setLayout(new BorderLayout());
		_tools.add(_input, BorderLayout.CENTER);
		_tools.add(_info, BorderLayout.SOUTH);
		_tools.add(goPanel, BorderLayout.EAST);

		goPanel.add(_goButton, BorderLayout.CENTER);
		goPanel.add(_switchButton, BorderLayout.SOUTH);

		getContentPane().setLayout(new BorderLayout());
		JPanel VARNAs = new JPanel();
		VARNAs.setLayout(new GridLayout(1,1));
		VARNAs.add(_vpMaster);
		getContentPane().add(VARNAs, BorderLayout.CENTER);
		getContentPane().add(_tools, BorderLayout.SOUTH);

		setVisible(true);

		_vpMaster.getVARNAUI().UIRadiate();
		_vpMaster.setTitleFontSize(26f);
		_vpMaster.setTitleFontStyle(Font.PLAIN);
		
		this.setTitle("RNA Folding Demo - Simple Matching Algorithm");
		
		onStructureRedrawn();
	}

	private void showSolution()
	{
		RNA rflat = getRNA();
		rflat.drawRNALine();
		_vpMaster.setTitle("#Secondary Structures: "+count(getSeq()));
		_struct.setText(getStruct());

		_vpMaster.drawRNA(rflat);
		_vpMaster.showRNAInterpolated(rflat);

		RNA rfolded = getRNA();
		rfolded.drawRNARadiate(_vpMaster.getConfig());
        _vpMaster.showRNAInterpolated(rfolded);
		
        RNA rLinear = getRNA();
		rLinear.drawRNALine();
		_vpMaster.showRNAInterpolated(rLinear);
		
	}
	

	
	public RNA getRNA() {
		RNA r = new RNA();
		try {
			r.setRNA(getSeq(), getStruct());
		} catch (ExceptionUnmatchedClosingParentheses e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ExceptionFileFormatOrSyntax e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return r;
	}


 	public String getSeq()
	{
		return ""+_seq1.getSelectedItem();
	}

 	private boolean canBasePair(char a, char b)
 	{
 		if ((a=='G')&&(b=='C'))
 			return true;
 		if ((a=='C')&&(b=='G'))
 			return true;
 		if ((a=='U')&&(b=='A'))
 			return true;
 		if ((a=='A')&&(b=='U'))
 			return true;
 		return false;
 	}
 	
 	public int[][] fillMatrix(String seq)
 	{
		int n = seq.length();
		int[][] tab = new int[n][n];
		for(int m=1;m<=n;m++)
		{
			for(int i=0;i<n-m+1;i++)
			{
				int j = i+m-1;
				tab[i][j] = 0;
				if (i<j)
				{ 
					tab[i][j] = Math.max(tab[i][j], tab[i+1][j]); 
					for (int k=i+1;k<=j;k++)
					{
						if (canBasePair(seq.charAt(i),seq.charAt(k)))
						{
							int fact1 = 0;
							if (k>i+1)
							{
								fact1 = tab[i+1][k-1];
							}
							int fact2 = 0;
							if (k<j-1)
							{
								fact2 = tab[k+1][j];
							}
							tab[i][j] = Math.max(tab[i][j],1+fact1+fact2);
						} 
					}
				}
			}			
		}
 		return tab;
 	}

 	private String backtrack(int[][] tab, String seq)
 	{
 		return backtrack(tab,seq, 0, seq.length()-1);
 	}

 	private String backtrack(int[][] tab, String seq, int i, int j)
 	{
 		if (i<j)
		{ 
			if (tab[i][j] == tab[i+1][j])
			{
				return "."+backtrack(tab, seq, i+1,j);
			}

			for (int k=i+1;k<=j;k++)
			{
				if (canBasePair(seq.charAt(i),seq.charAt(k)))
				{
					int fact1 = 0;
					if (k>i+1)
					{
						fact1 = tab[i+1][k-1];
					}
					int fact2 = 0;
					if (k<j-1)
					{
						fact2 = tab[k+1][j];
					}
					if (tab[i][j]==1+fact1+fact2)
					{ 
						return "("+backtrack(tab, seq, i+1,k-1)+")"+backtrack(tab, seq, k+1,j);
					}
				} 
			}
		}
		else if  (i==j)
		{
			return ".";
		}
		return "";
 	}
 	
 	public long count(String seq)
 	{
		int n = seq.length();
		long[][] tab = new long[n][n];
		for(int m=1;m<=n;m++)
		{
			for(int i=0;i<n-m+1;i++)
			{
				int j = i+m-1;
				tab[i][j] = 0;
				if (i<j)
				{ 
					tab[i][j] += tab[i+1][j]; 
					for (int k=i+1;k<=j;k++)
					{
						if (canBasePair(seq.charAt(i),seq.charAt(k)))
						{
							long fact1 = 1;
							if (k>i+1)
							{
								fact1 = tab[i+1][k-1];
							}
							long fact2 = 1;
							if (k<j-1)
							{
								fact2 = tab[k+1][j];
							}
							tab[i][j] += fact1*fact2;
						} 
					}
				}
				else
				{
					tab[i][j] = 1;
				}
			}			
		}
 		return tab[0][n-1];
 	}
 	
 	
	public String getStruct() {
		String seq = getSeq();
		seq = seq.toUpperCase();
		int n = seq.length();
		int[][] mfe = fillMatrix(seq);
		String back = backtrack(mfe,seq);
		return back;
	}

	private String cleanStruct(String struct) {
		struct = struct.replaceAll("[:-]", "");
		return struct;
	}



	public void init() {
		_vpMaster.setBackground(_backgroundColor);
		_error = true;
	}

	@SuppressWarnings("unused")
	private Color getSafeColor(String col, Color def) {
		Color result;
		try {
			result = Color.decode(col);
		} catch (Exception e) {
			try {
				result = Color.getColor(col, def);
			} catch (Exception e2) {
				return def;
			}
		}
		return result;
	}

	public VARNAPanel get_varnaPanel() {
		return _vpMaster;
	}

	public void set_varnaPanel(VARNAPanel surface) {
		_vpMaster = surface;
	}


	public JLabel get_info() {
		return _info;
	}

	public void set_info(JLabel _info) {
		this._info = _info;
	}

	public static void main(String[] args) {
		NussinovDemo d = new NussinovDemo();
		d.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		d.pack();
		d.setVisible(true);
	}

	public void onStructureRedrawn() {
		_vpMaster.repaint();
	}

	public void onWarningEmitted(String s) {
		// TODO Auto-generated method stub
		
	}

	public void onLoad(String path) {
		// TODO Auto-generated method stub
		
	}

	public void onLoaded() {
		// TODO Auto-generated method stub
		
	}

	public void onUINewStructure(VARNAConfig v, RNA r) {
		// TODO Auto-generated method stub
		
	}
}
