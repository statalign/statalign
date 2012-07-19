package statalign.postprocess.gui;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import javax.swing.undo.UndoManager;

import fr.orsay.lri.varna.VARNAPanel;
import fr.orsay.lri.varna.controlers.ControleurBlinkingThread;
import fr.orsay.lri.varna.controlers.ControleurClicMovement;
import fr.orsay.lri.varna.controlers.ControleurDraggedMolette;
import fr.orsay.lri.varna.controlers.ControleurInterpolator;
import fr.orsay.lri.varna.controlers.ControleurMolette;
import fr.orsay.lri.varna.controlers.ControleurVARNAPanelKeys;
import fr.orsay.lri.varna.exceptions.ExceptionNonEqualLength;
import fr.orsay.lri.varna.models.VARNAConfig;
import fr.orsay.lri.varna.models.rna.RNA;
import fr.orsay.lri.varna.views.VueUI;

public class StructureGUI extends VARNAPanel {

	/**
	 * This class implements a graphical interface for displaying the current consensus structure. RNA
	 * graphics are provided by VARNA.
	 * 
	 * @author Preeti Arunapuram
	 */
	private static final long serialVersionUID = 1L;
	private static String sequence;
	private static String structure;
	
	private Graphics2D g2;
	public String title;
	

	private VARNAConfig _conf = new VARNAConfig();
	
	UndoManager _manager;
	
	private ControleurBlinkingThread _blink;
	


	// private Point _positionRelativeSouris;
	private Point _translation;
	private boolean _horsCadre;
	private boolean _premierAffichage;


	private ControleurInterpolator _interpolator;

	private VueUI _UI = new VueUI(this);

	
	public StructureGUI(String title) throws ExceptionNonEqualLength {
		sequence = null;
		structure = null;
		this.title = title;
	}
	
	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		g2 = (Graphics2D)g;
		
		if(structure == null || sequence == null) {
			g2.drawString("Waiting for data..", 30, 30);
            return;
		}
		
		g2.setPaint(Color.BLACK);
		g2.setFont(new Font("SANS_SERIF", Font.BOLD, 16));
		g2.drawString(title, 10, 20);
		g2.setFont(new Font("MONOSPACED", Font.PLAIN, 12));
		
	}
	
	public void updateAndDraw(String seq, String str) {
		sequence = seq;
		structure = str;
		
		repaint();
		
		super.drawRNAInterpolated(sequence, structure, RNA.DRAW_MODE_RADIATE);
		
		setBackground(VARNAConfig.DEFAULT_BACKGROUND_COLOR);
		_manager = new UndoManager();
		_manager.setLimit(10000);
		_UI.addUndoableEditListener(_manager);
		

		_blink = new ControleurBlinkingThread(this,
				ControleurBlinkingThread.DEFAULT_FREQUENCY, 0, 1.0, 0.0, 0.2);
		_blink.start();

		_premierAffichage = true;
		_translation = new Point(0, 0);

		_horsCadre = false;
		this.setFont(_conf._fontBasesGeneral);

		// ajout des controleurs au VARNAPanel
		ControleurClicMovement controleurClicMovement = new ControleurClicMovement(	this);
		this.addMouseListener(controleurClicMovement);
		this.addMouseMotionListener(controleurClicMovement);
		this.addMouseWheelListener(new ControleurMolette(this));

		ControleurDraggedMolette ctrlDraggedMolette = new ControleurDraggedMolette(
				this);
		this.addMouseMotionListener(ctrlDraggedMolette);
		this.addMouseListener(ctrlDraggedMolette);

		ControleurVARNAPanelKeys ctrlKey = new ControleurVARNAPanelKeys(this);
		this.addKeyListener(ctrlKey);
		this.addFocusListener(ctrlKey);
		
		_interpolator = new ControleurInterpolator(this);
		_interpolator.start();
	}
	}
	
	
	
	
	
	
	

