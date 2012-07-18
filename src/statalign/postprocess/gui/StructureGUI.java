package statalign.postprocess.gui;

import java.awt.Dimension;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.geom.Point2D;
import java.util.ArrayList;

import javax.swing.JPanel;
import javax.swing.undo.UndoManager;

import fr.orsay.lri.varna.VARNAPanel;
import fr.orsay.lri.varna.controlers.ControleurBlinkingThread;
import fr.orsay.lri.varna.controlers.ControleurClicMovement;
import fr.orsay.lri.varna.controlers.ControleurDraggedMolette;
import fr.orsay.lri.varna.controlers.ControleurInterpolator;
import fr.orsay.lri.varna.controlers.ControleurMolette;
import fr.orsay.lri.varna.controlers.ControleurVARNAPanelKeys;
import fr.orsay.lri.varna.exceptions.ExceptionNonEqualLength;
import fr.orsay.lri.varna.interfaces.InterfaceVARNAListener;
import fr.orsay.lri.varna.interfaces.InterfaceVARNARNAListener;
import fr.orsay.lri.varna.interfaces.InterfaceVARNASelectionListener;
import fr.orsay.lri.varna.models.BaseList;
import fr.orsay.lri.varna.models.VARNAConfig;
import fr.orsay.lri.varna.models.annotations.TextAnnotation;
import fr.orsay.lri.varna.models.rna.ModeleBase;
import fr.orsay.lri.varna.models.rna.RNA;
import fr.orsay.lri.varna.views.VueMenu;
import fr.orsay.lri.varna.views.VueUI;

public class StructureGUI extends VARNAPanel {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private static String sequence;
	private static String structure;
	
	private boolean _debug = false;

	private VARNAConfig _conf = new VARNAConfig();
	
	private ArrayList<InterfaceVARNAListener> _VARNAListeners = new ArrayList<InterfaceVARNAListener>();
	private ArrayList<InterfaceVARNASelectionListener> _selectionListeners = new ArrayList<InterfaceVARNASelectionListener>(); 
	private ArrayList<InterfaceVARNARNAListener> _RNAListeners = new ArrayList<InterfaceVARNARNAListener>();
	
	UndoManager _manager;

	
	private Point2D.Double[] _realCoords = new Point2D.Double[0];
	private Point2D.Double[] _realCenters = new Point2D.Double[0];
	private double _scaleFactor = 1.0;
	private Point2D.Double _offsetPanel = new Point2D.Double();
	private Point2D.Double _offsetRNA = new Point2D.Double();

	private double _offX;
	private double _offY;
	
	
	private ControleurBlinkingThread _blink;
	private BaseList _selectedBases = new BaseList("selection");
	private ArrayList<ModeleBase> _backupSelection = new ArrayList<ModeleBase>();
	private Integer _nearestBase = null;
	private Point2D.Double _lastSelectedCoord = new Point2D.Double(0.0,0.0);
	
	private Point2D.Double _linkOrigin = null;
	private Point2D.Double _linkDestination = null;
	
	private Rectangle _selectionRectangle = null;

	private boolean _highlightAnnotation = false;

	private int _titleHeight;
	private Dimension _border = new Dimension(0,0);

	private boolean _drawBBox = false;
	private boolean _drawBorder = false;
	


	// private Point _positionRelativeSouris;
	private Point _translation;
	private boolean _horsCadre;
	private boolean _premierAffichage;


	private ControleurInterpolator _interpolator;
	
	private VueMenu _popup = new VueMenu(this);

	private VueUI _UI = new VueUI(this);


	private TextAnnotation _selectedAnnotation;
	
	public StructureGUI(String seq, String str) throws ExceptionNonEqualLength {
		sequence = seq;
		structure = str;
	}
	
	public void updateAndDraw(String seq, String str) {
		sequence = seq;
		structure = str;
		
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
	
	
	
	
	
	
	

